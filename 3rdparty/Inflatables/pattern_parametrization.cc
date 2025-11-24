#include "pattern_parametrization.hh"

#include <MeshFEM/SparseMatrices.hh>
#include <MeshFEM/Laplacian.hh>
#include <MeshFEM/GlobalBenchmark.hh>
#include <MeshFEM/MeshIO.hh>
#include <complex>
#include <set>

#include <MeshFEM/ParallelAssembly.hh>
#include "circular_mean.hh"
#include "subdivide_triangle.hh"
#include "curvature.hh"

namespace parametrization {

////////////////////////////////////////////////////////////////////////////////
// RegularizedGenericParametrizer: Global nonlinear energy with auxiliary variables with bounds on both singular values
////////////////////////////////////////////////////////////////////////////////
// Initialize from the local-global parametrizer
RegularizedPatternParametrizer::RegularizedPatternParametrizer(LocalGlobalGenericParametrizer &lggparam, MatInfoFunction matInfo, Eigen::VectorXd default_pattern_params, size_t input_num_pattern_vars)
    : Base(lggparam),
      materialInfo(matInfo)
{
    Base::setStretchRegW(0.0);
    Base::setDiffRegW(0.0);

    if (size_t(default_pattern_params.size()) != size_t(mesh().numTris() * input_num_pattern_vars)) {
        throw std::runtime_error("Invalid number of default pattern parameters!");
    }

    m_num_pattern_vars = input_num_pattern_vars;
    m_pattern_params = default_pattern_params;

    m_set_matInfoArgs();

    for (size_t i = 0; i < materialInfo.size() / 3; ++i) {
        if (size_t(materialInfo[i * 3](getMatInfoArgs()).size()) != size_t(mesh().numTris())) {
            throw std::runtime_error("The information from splines does not match the number of pattern parameters! Index: " + std::to_string(i));
        }
    }

    setVars(getVars());
    m_bendingSensitivityCache.clear();

    m_pp_bounds = MX2d::Zero(num_pattern_vars(), 2); // Upper and lower bounds on the pattern parameters
    //Initialize the bounds eigen vectors of negative and positive infinity
    m_pp_bounds.col(0) = Eigen::VectorXd::Constant(num_pattern_vars(), -std::numeric_limits<double>::infinity());
    m_pp_bounds.col(1) = Eigen::VectorXd::Constant(num_pattern_vars(), std::numeric_limits<double>::infinity());

    // Default normalization factors
    m_param_normalization_factors = Eigen::VectorXd::Ones(num_pattern_vars());
     
}

Real RegularizedPatternParametrizer::energy(PatternEnergyType etype) const {
    BENCHMARK_SCOPED_TIMER_SECTION timer("RegularizedPatternParametrizer.energy");
    // rgp fitting energy + phi regularization
    // with new bending energy and pattern parameters smoothness regularization 
    Real result = 0.0;
    if ((etype == PatternEnergyType::Full) || (etype == PatternEnergyType::RGP)) {
        result += Base::energy(RegularizedGenericParametrizer::EnergyType::Full);
    }
    if (etype == PatternEnergyType::DEBUG_Fitting) result += Base::energy(RegularizedGenericParametrizer::EnergyType::Fitting);
    if (etype == PatternEnergyType::DEBUG_PhiRegularization) result += Base::energy(RegularizedGenericParametrizer::EnergyType::PhiRegularization);
    if (etype == PatternEnergyType::DEBUG_StretchRegularization) result += Base::energy(RegularizedGenericParametrizer::EnergyType::StretchRegularization);
    const size_t nt = mesh().numTris();
    if (m_bending_w != 0.0 && ((etype == PatternEnergyType::Full) || (etype == PatternEnergyType::Bending))) {
        Real bending_energy = 0.0;
        m_bendingSensitivityCache.update(*this);
        for (const auto &tri : mesh().elements()) {
            const size_t ti = tri.index();
            bending_energy += m_bendingSensitivityCache.sensitivityForBending.at(ti).objective * tri->volume() *  m_bending_w;
        }
        if (Base::scaleInvariantFittingEnergy) bending_energy /= mesh().volume();
        result += bending_energy;
    }
    // Edge-based regularization terms
    // Dual Laplacian-based regularization terms
    if (((etype == PatternEnergyType::Full) || (etype == PatternEnergyType::PatternRegularization)) && m_pattern_reg_w != 0.0) {
        dualLaplacianStencil.visit_edges([this, &result, nt](size_t i, size_t j, Real w_ij) {
            const size_t ti = i,
                         tj = j;
            for (size_t pi = 0; pi < this->num_pattern_vars(); ++pi) {
                result += w_ij * std::pow(std::abs(this->getPatternParams()[pi * nt + ti] - this->getPatternParams()[pi * nt + tj]) / m_param_normalization_factors[pi], m_pattern_reg_p) * m_pattern_reg_w;
            }
        }
        );
    }

    if (etype == PatternEnergyType::Full || etype == PatternEnergyType::PatternBoundConstraint) {
        result += energyBoundConstraint();
    }
    return result;
}

Eigen::VectorXd RegularizedPatternParametrizer::gradient(PatternEnergyType etype) const {
    BENCHMARK_SCOPED_TIMER_SECTION timer("RegularizedPatternParametrizer.gradient");
    Eigen::VectorXd result(numVars());
    result.setZero();
    const size_t nt = mesh().numTris();

    Eigen::VectorXd rgp_grad = Base::gradient(RegularizedGenericParametrizer::EnergyType::Full);
    if ((etype == PatternEnergyType::Full) || (etype == PatternEnergyType::RGP)) {
        result.head(patternOffset()) += rgp_grad.head(patternOffset());
        result.segment(patternOffset(), num_pattern_vars() * nt) += get_pattern_params_gradient(rgp_grad.tail(Base::numStretchVars()));
    }
    // Add gradients of the bending energy over pattern parameters and psi, and smoothness energy over pattern parameters.
    if (((etype == PatternEnergyType::Full) || (etype == PatternEnergyType::Bending)) && m_bending_w != 0.0) {
        m_bendingSensitivityCache.update(*this);
        // Gradient of bending energy
        for (const auto &tri : mesh().elements()) {
            const size_t ti = tri.index();
            Eigen::VectorXd psi_p =  m_bendingSensitivityCache.sensitivityForBending.at(ti).psi_p_gradient * m_bending_w;
            result[Base::psiOffset() + ti] += psi_p[0] * tri->volume() / mesh().volume();
            for (size_t pi = 0; pi < num_pattern_vars(); ++pi) {
                result[patternOffset() + pi * nt + ti] += psi_p[pi + 1] * tri->volume() / mesh().volume();
            }
        }
    }
    if (((etype == PatternEnergyType::Full) || (etype == PatternEnergyType::PatternRegularization)) && m_pattern_reg_w != 0.0) {
        // Gradient of edge-based regularization terms
        dualLaplacianStencil.visit_edges([this, &result, nt](size_t ti, size_t tj, Real w_ij) {    
        // Smoothness regularization
            // Using std::copysign(1.0, stretch_diff) doesn't work since it gives bad derivatives around stretch_diff = 0.
            // We get better results explicitly setting the derivative equal to zero in this case.
            Real sign = 0.0;
            for (size_t pi = 0; pi < num_pattern_vars(); ++pi) {
                Real p_diff = this->getPatternParams()[pi * nt + ti] - this->getPatternParams()[pi * nt + tj];
                if (p_diff > 0) sign =  1.0;
                if (p_diff < 0) sign = -1.0;

                Real val;
                if (m_pattern_reg_p == 1.0) { val = m_pattern_reg_w * sign; }
                else                        { val = m_pattern_reg_w * std::pow(std::abs(p_diff), m_pattern_reg_p - 1.0) * sign; }
                val *= w_ij;

                result[patternOffset() + pi * nt + ti] += val;
                result[patternOffset() + pi * nt + tj] -= val;
            }
        }
        );
    }

    if (etype == PatternEnergyType::Full || etype == PatternEnergyType::PatternBoundConstraint) {
        addBoundConstraintGradient(result);
    }
    
    // Debug
    if (etype == PatternEnergyType::DEBUG_Fitting) {
        Eigen::VectorXd fitting_grad = Base::gradient(RegularizedGenericParametrizer::EnergyType::Fitting);
        result.head(patternOffset()) += fitting_grad.head(patternOffset());
        result.segment(patternOffset(), num_pattern_vars() * nt) += get_pattern_params_gradient(fitting_grad.tail(Base::numStretchVars()));
    }
    if (etype == PatternEnergyType::DEBUG_PhiRegularization) {
        Eigen::VectorXd phi_reg_grad = Base::gradient(RegularizedGenericParametrizer::EnergyType::PhiRegularization);
        result.head(patternOffset()) += phi_reg_grad.head(patternOffset());
        result.segment(patternOffset(), num_pattern_vars() * nt) += get_pattern_params_gradient(phi_reg_grad.tail(Base::numStretchVars()));
    }
    if (etype == PatternEnergyType::DEBUG_StretchRegularization) {
        Eigen::VectorXd stretch_reg_grad = Base::gradient(RegularizedGenericParametrizer::EnergyType::StretchRegularization);
        result.head(patternOffset()) += stretch_reg_grad.head(patternOffset());
        result.segment(patternOffset(), num_pattern_vars() * nt) += get_pattern_params_gradient(stretch_reg_grad.tail(Base::numStretchVars()));
    }
    return result;
}

SuiteSparseMatrix RegularizedPatternParametrizer::baseHessianSparsityPattern() const {
    BENCHMARK_SCOPED_TIMER_SECTION("RegularizedPatternParametrizer::baseHessianSparsityPattern");
    if (m_cachedBaseHessianSparsity) return *m_cachedBaseHessianSparsity;
    auto baseHsp = Base::hessianSparsityPattern();
    m_cachedBaseHessianSparsity = std::make_unique<SuiteSparseMatrix>(baseHsp);

    return baseHsp;
}


SuiteSparseMatrix RegularizedPatternParametrizer::hessianSparsityPattern(Real val) const {
    BENCHMARK_SCOPED_TIMER_SECTION timer("RegularizedPatternParametrizer.hessianSparsityPattern");
    if (m_cachedHessianSparsity) {
        if (m_cachedHessianSparsity->Ax[0] != val) m_cachedHessianSparsity->fill(val);
        return *m_cachedHessianSparsity;
    }
    SuiteSparseMatrix result(numVars(), numVars());
    result.symmetry_mode = SuiteSparseMatrix::SymmetryMode::UPPER_TRIANGLE;
    result.Ap.reserve(numVars() + 1);
    auto &Ap = result.Ap;
    auto &Ai = result.Ai;
    SuiteSparseMatrix rgp_hsp = baseHessianSparsityPattern();

    Ap.resize(patternOffset() + 1);
    std::copy(rgp_hsp.Ap.begin(), rgp_hsp.Ap.begin() + patternOffset() + 1, Ap.begin());
    size_t rgp_nnz = rgp_hsp.Ap[patternOffset() + 1];
    Ai.resize(rgp_nnz);
    std::copy(rgp_hsp.Ai.begin(), rgp_hsp.Ai.begin() + rgp_nnz, Ai.begin());


    auto addIdx = [&](const size_t idx) { Ai.push_back(idx); };

    auto finalizeCol = [&]() {
        const size_t colStart = Ap.back();
        const size_t colEnd = Ai.size();
        Ap.push_back(colEnd);
        std::sort(Ai.begin() + colStart, Ai.begin() + colEnd);
    };

    const auto &m = mesh();
    const size_t nv = m.numVertices();
    const size_t nt = mesh().numTris();

    // Tri field columns: interact with corner vertices, neighbors, and selves
    const size_t phio = Base::phiOffset(), psio = Base::psiOffset();
    const size_t numTriFields = num_pattern_vars();
    // In this case fieldOffset only iterates over pattern parameter indices.
    for (size_t fieldOffset = 0; fieldOffset < numTriFields; ++fieldOffset) {
        for (const auto &tri : m.elements()) {
            const size_t tj = tri.index();
            for (const auto &v : tri.vertices()) {
                addIdx(v.index());      // u variable
                addIdx(v.index() + nv); // v variable
            }
            addIdx(phio + tj); // phi-phi/phi-psi/phi-alpha/phi-beta interaction
            addIdx(psio   + tj); // psi-psi/psi-alpha/psi-beta interaction

            for (size_t pi = 0; pi < fieldOffset + 1; ++pi) {
                addIdx(patternOffset() + nt * pi + tj);
            }

            // Laplacian-style regularization (upper triangle)
            dualLaplacianStencil.visit(tri.index(), [this, &m, &addIdx, fieldOffset, phio, tj, nt](size_t /* j */, size_t i, Real /* w_ij */) {
                auto tri_i = m.tri(i);
                const size_t ti = tri_i.index();
                if (!tri_i) return;
                if (Base::useBarrier()) {
                    //Barrier interactions
                    addIdx(phio + ti); // phi-p_{fieldOffset} interaction between neighboring triangles
                    for (size_t pi = 0; pi < fieldOffset; ++pi) {
                        addIdx(patternOffset() + nt * pi + ti); // p_{fieldOffset}-p_{pi} interaction between neighboring triangles
                    }
                }
                if (ti > tj) return;
                // Each pattern parameter interacts only with the neighboring corresponding pattern parameter. 
                addIdx(patternOffset() + nt * fieldOffset + ti); // p_{fieldOffset}_-p_{fieldOffset} interaction between neighboring triangles
            });

            finalizeCol();
        }
    }

    result.nz = result.Ai.size();
    result.Ax.assign(result.nz, val);

    m_cachedHessianSparsity = std::make_unique<SuiteSparseMatrix>(result);
    return result;
}

SuiteSparseMatrix RegularizedPatternParametrizer::hessian(PatternEnergyType etype) const {
    SuiteSparseMatrix H = hessianSparsityPattern();
    hessian(H, etype);
    return H;
}

template<typename Real_>
struct d_stretch_dp_entry {
    typename CSCMatrix<SuiteSparse_long, Real_>::index_type first;
    typename CSCMatrix<SuiteSparse_long, Real_>::value_type second;
};
template<typename Real_>
using d_stretch_dp_type = std::vector<std::vector<d_stretch_dp_entry<Real_>>>;

template<typename Real_>
void get_d_stretch_dp(d_stretch_dp_type<Real_> &d_stretch_dp, size_t nv, size_t nt, size_t stretchOffset, const Eigen::VectorXd  &d_alpha_dp, const Eigen::VectorXd &d_beta_dp, size_t num_pattern_vars) {
    using index_type = SuiteSparse_long;

    d_stretch_dp.resize(nv);
    for (auto &row: d_stretch_dp) { row.clear(); row.reserve(num_pattern_vars); } // only the pattern parameters for a triangle affect the stretch variables for that triangle.
    // For non stretch variables, the jacobian is identity.
    for (size_t i = 0; i < stretchOffset; ++i) d_stretch_dp[i].push_back({index_type(i), 1.0});
    // For stretch variables, use the spline jacobians.
    for (size_t k = 0; k < nt; ++k) {
        for (size_t pi = 0; pi < num_pattern_vars; ++pi) {
            size_t linear_index = pi * nt + k;
            d_stretch_dp[stretchOffset + k].push_back({index_type(stretchOffset + linear_index), d_alpha_dp[linear_index]});
            d_stretch_dp[stretchOffset + k + nt].push_back({index_type(stretchOffset + linear_index), d_beta_dp[linear_index]});
        }
    }

}

void RegularizedPatternParametrizer::hessian(SuiteSparseMatrix &H, PatternEnergyType etype) const {
    BENCHMARK_SCOPED_TIMER_SECTION timer("RegularizedPatternParametrizer.hessian");
    if (H.symmetry_mode != SuiteSparseMatrix::SymmetryMode::UPPER_TRIANGLE) throw std::runtime_error("Hessian must be upper triangular");

    // Copy the beginning rgp block over. 
    auto baseH = Base::hessianSparsityPattern();
    Base::hessian(baseH, RegularizedGenericParametrizer::EnergyType::Full);

    const size_t psio = Base::psiOffset(), po = patternOffset();
    const size_t nt = mesh().numTris();

    if ((etype == PatternEnergyType::Full) || (etype == PatternEnergyType::RGP)) {
        // Convert hessian over stretching variables to hessian over pattern parameters.
        Eigen::VectorXd d_alpha_dp = materialInfo[1](getMatInfoArgs());
        Eigen::VectorXd d_beta_dp = materialInfo[1 + num_info_per_mat_var](getMatInfoArgs());

        d_stretch_dp_type<Real> d_stretch_dp;
        get_d_stretch_dp(d_stretch_dp, Base::numVars(), nt, Base::stretchOffset(), d_alpha_dp, d_beta_dp, num_pattern_vars());
        // Accumulate contribution of each (upper triangle) entry in baseH to the
        // full Hessian term:
        //      d_stretch_dpi baseH_kl d_stretch_dj
        using Idx = typename SuiteSparseMatrix::index_type;
        Idx idx = 0, idx2 = 0;
        Idx ncol = baseH.n, colbegin = baseH.Ap[0];
        // Sparse matrix multiplication: A^T baseH A where A is the block diagonal and contains the sensitivity of d_x_hat_d_x.
        for (Idx l = 0; l < ncol; ++l) {
            const Idx colend = baseH.Ap[l + 1];
            for (auto entry = colbegin; entry < colend; ++entry) {
                const Idx k = baseH.Ai[entry];
                const auto v = baseH.Ax[entry];
                if (v == 0.0) continue;
                if (!(k <= l)) throw std::runtime_error("Base hessian must be upper triangular");
                const auto &dvk_dr = d_stretch_dp[k];
                const auto &dvl_dr = d_stretch_dp[l];
                for (const auto &dvl_drj : dvl_dr) {
                    const Idx j = dvl_drj.first;
                    if (dvk_dr.size() == 0) continue;
                    const auto val = dvl_drj.second * v;
                    {
                        const Idx i = dvk_dr[0].first;
                        if (i > j) continue;
                        idx = H.template addNZAtLoc</* _knownGood = */ false>(i, j, val * dvk_dr[0].second, idx);
                    }
                    for (size_t ii = 1; ii < dvk_dr.size(); ++ii) {
                        const Idx i = dvk_dr[ii].first;
                        if (i > j) break;
                        idx = H.template addNZAtLoc</* _knownGood = */ true>(i, j, val * dvk_dr[ii].second, idx);
                    }
                }
                if (k != l) {
                    // Contribution from (l, k), if it falls in the upper triangle of H; capture all the missed entries from the previous loop due to the (i>j) check.
                    for (const auto &dvk_drj : dvk_dr) {
                        const Idx j = dvk_drj.first;
                        if (dvl_dr.size() == 0) continue;
                        const auto val = dvk_drj.second * v;
                        {
                            const Idx i = dvl_dr[0].first;
                            if (i > j) continue;
                            idx2 = H.template addNZAtLoc</* _knownGood = */ false>(i, j, val * dvl_dr[0].second, idx2);
                        }
                        for (size_t ii = 1; ii < dvl_dr.size(); ++ii) {
                            const Idx i = dvl_dr[ii].first;
                            if (i > j) break;
                            idx2 = H.template addNZAtLoc</* _knownGood = */ true>(i, j, val * dvl_dr[ii].second, idx2);
                        }
                    }
                }
            }
            colbegin = colend;
        }

        // Accumulate contribution of the Hessian of x_hat wrt the bending transformation variables.
        //  dE/stretch^j (d^2 stretch^j / dp_k dp_l)
        const Eigen::VectorXd rgp_grad = Base::gradient(RegularizedGenericParametrizer::EnergyType::Full);
        Real alphaOffset = patternOffset();

        Eigen::VectorXd d_alpha_d_p_d_p = materialInfo[2](getMatInfoArgs());

        std::vector<Eigen::VectorXd> alpha_hessians(num_pattern_vars() * num_pattern_vars());
        for (size_t pi = 0; pi < num_pattern_vars(); ++pi) {
            for (size_t pj = 0; pj < num_pattern_vars(); ++pj) {
                size_t linear_index = pi * num_pattern_vars() + pj;
                alpha_hessians[pi * num_pattern_vars() + pj] = d_alpha_d_p_d_p.segment(linear_index * nt, nt);
            }
        }
        for (size_t ti = 0; ti < nt; ++ti) {
            for (size_t pi = 0; pi < num_pattern_vars(); ++pi) {
                for (size_t pj = pi; pj < num_pattern_vars(); ++pj) {
                    Real val = rgp_grad[alphaOffset + ti] * alpha_hessians[pi * num_pattern_vars() + pj][ti];
                    if (val != 0.0) H.addNZ(patternOffset() + pi * nt + ti, patternOffset() + pj * nt + ti, val);
                }
            }
        }

        size_t betaOffset = patternOffset() + nt;
        Eigen::VectorXd d_beta_d_p_d_p = materialInfo[2 + num_info_per_mat_var](getMatInfoArgs());

        std::vector<Eigen::VectorXd> beta_hessians(num_pattern_vars() * num_pattern_vars());
        for (size_t pi = 0; pi < num_pattern_vars(); ++pi) {
            for (size_t pj = pi; pj < num_pattern_vars(); ++pj) {
                size_t linear_index = pi * num_pattern_vars() + pj;
                beta_hessians[pi * num_pattern_vars() + pj] = d_beta_d_p_d_p.segment(linear_index * nt, nt);
            }
        }
        for (size_t ti = 0; ti < nt; ++ti) {
            for (size_t pi = 0; pi < num_pattern_vars(); ++pi) {
                for (size_t pj = pi; pj < num_pattern_vars(); ++pj) {
                    Real val = rgp_grad[betaOffset + ti] * beta_hessians[pi * num_pattern_vars() + pj][ti];
                    if (val != 0.0) H.addNZ(patternOffset() + pi * nt + ti, patternOffset() + pj * nt + ti, val);
                }
            }
        }
    }

    if (((etype == PatternEnergyType::Full) || (etype == PatternEnergyType::Bending)) && m_bending_w != 0.0) {
        m_bendingSensitivityCache.update(*this);
        for (const auto &tri : mesh().elements()) {
            const size_t ti = tri.index();
            Real A = tri->volume() / mesh().volume();
            auto psi_p = (m_bendingSensitivityCache.sensitivityForBending.at(ti).psi_p_hessian * m_bending_w).eval();
            H.addNZ(psio + ti, psio + ti, psi_p(0, 0) * A);
            for (size_t pi = 0; pi < num_pattern_vars(); ++pi) {
                H.addNZ(psio + ti, patternOffset() + pi * nt + ti, psi_p(0, pi + 1) * A);
                for (size_t pj = pi; pj < num_pattern_vars(); ++pj) {
                    H.addNZ(patternOffset() + pi * nt + ti, patternOffset() + pj * nt + ti, psi_p(pi + 1, pj + 1) * A);
                }
            }
        }
    }

    if (((etype == PatternEnergyType::Full) || (etype == PatternEnergyType::PatternRegularization)) && m_pattern_reg_w != 0.0) {
        dualLaplacianStencil.visit_edges([this, &H, nt, po](size_t ti, size_t tj, Real w_ij) {
            if ((m_pattern_reg_p == 1.0)) return;
            if (ti > tj) std::swap(ti, tj); // Visit each stencil edge (unordered triangle pair) exactly once

            for (size_t pi = 0; pi < num_pattern_vars(); ++pi) {
                const Real p_diff = m_pattern_params[pi * nt + ti] - m_pattern_params[pi * nt + tj];
                if (m_pattern_reg_p < 2.0 && std::abs(p_diff) < 1e-14) continue;
                Real val = w_ij * m_pattern_reg_w * (m_pattern_reg_p - 1.0) * std::pow(std::abs(p_diff), m_pattern_reg_p - 2.0);
                if (val != 0.0) {
                    H.addNZ(po + pi * nt + ti, po + pi * nt + tj, -val);
                    H.addNZ(po + pi * nt + tj, po + pi * nt + tj,  val);
                    H.addNZ(po + pi * nt + ti, po + pi * nt + ti,  val);
                }
            }
        }
        );
    }
    if (etype == PatternEnergyType::Full || etype == PatternEnergyType::PatternBoundConstraint) {
        addBoundConstraintHessian(H);
    }
}


std::vector<Eigen::VectorXd> RegularizedPatternParametrizer::perVertexPatternParams() const {
    const auto &m = mesh();
    std::vector<Eigen::VectorXd> results;
    results.resize(num_pattern_vars());
    const std::vector<Eigen::VectorXd> &params = getMatInfoArgs();

    for (size_t pi = 0; pi < num_pattern_vars(); ++pi) {
        results[pi] = Eigen::VectorXd::Zero(m.numVertices());
        for (auto v : m.vertices()) {
            double &curr_param = results[pi](v.index());
            curr_param = 0;
            size_t tri_valence = 0;
            for (auto he : v.incidentHalfEdges()) {
                if (!he.tri()) continue;
                curr_param += params[pi][he.tri().index()];
                ++tri_valence;
            }
            curr_param /= tri_valence;
        }
    }

    return results;
}

std::tuple<std::shared_ptr<Mesh>, Eigen::VectorXd, std::vector<Eigen::VectorXd>>
RegularizedPatternParametrizer::upsampledVertexLeftStretchAnglesAndPatternParameters(size_t nsubdiv, double agreementThreshold) const {
    std::tuple<std::shared_ptr<Mesh>, Eigen::VectorXd, std::vector<Eigen::VectorXd>> result;

    std::vector<MeshIO::IOVertex > subVertices;
    std::vector<MeshIO::IOElement> subElements;

    auto coarseVertexAngles = Base::perVertexLeftStretchAngles(agreementThreshold);
    std::vector<Eigen::VectorXd> coarseVertexPatternParams = perVertexPatternParams();

    std::vector<double> subAngles;
    std::vector<std::vector<double>> subPatternParams;
    subPatternParams.resize(num_pattern_vars());

    // Until we have implemented a weighted angle averaging algorithm,
    // implement the rational barycentric coordinate weights by duplicating the
    // corresponding angles (inefficient).
    const size_t barycentricDenominator = nsubdiv + 1;
    std::vector<double> cornerAngleVec(barycentricDenominator);

    const auto &m = mesh();
    PointGluingMap indexForPoint;
    // size_t triIdx = 0;
    for (const auto &tri : m.elements()) {
        // bool verbose = (triIdx++ == 6467);

        auto newPt = [&](const Point3D &p, double lambda_0, double lambda_1, double lambda_2) {
            cornerAngleVec.clear();
            auto replicateAngle = [&](size_t corner, double lambda) {
                const size_t numerator = std::round(lambda * barycentricDenominator);
                double angle = 2.0 * coarseVertexAngles[tri.vertex(corner).index()]; // average 2x the angle to account for 2-RoSy
                for (size_t i = 0; i < numerator; ++i) cornerAngleVec.push_back(angle);
            };

            replicateAngle(0, lambda_0);
            replicateAngle(1, lambda_1);
            replicateAngle(2, lambda_2);

            assert(cornerAngleVec.size() == barycentricDenominator);

            subVertices.emplace_back(p);
            subAngles.push_back(0.5 * circularMean(cornerAngleVec));
            for (size_t pi = 0; pi < num_pattern_vars(); ++pi) {
                subPatternParams[pi].push_back(lambda_0 * coarseVertexPatternParams[pi][tri.vertex(0).index()] +
                                                lambda_1 * coarseVertexPatternParams[pi][tri.vertex(1).index()] +
                                                lambda_2 * coarseVertexPatternParams[pi][tri.vertex(2).index()]);
            }
            
            // if (verbose) {
            //     std::cout << "pt " << p.transpose() << " mean " << subAngles.back() << " from";
            //     for (double v : cornerAngleVec)
            //         std::cout << "\t" << v;
            //     std::cout << std::endl;
            // }

            return subVertices.size() - 1;
        };

        subdivide_triangle(nsubdiv,
                padTo3D(m_uv.row(tri.vertex(0).index()).transpose().eval()),
                padTo3D(m_uv.row(tri.vertex(1).index()).transpose().eval()),
                padTo3D(m_uv.row(tri.vertex(2).index()).transpose().eval()),
                indexForPoint,
                newPt, [&](size_t i0, size_t i1, size_t i2) { subElements.emplace_back(i0, i1, i2); });
    }

    std::get<0>(result) = std::make_shared<Mesh>(subElements, subVertices);
    std::get<1>(result) = Eigen::Map<Eigen::VectorXd>(subAngles.data(), subAngles.size());
    std::vector<Eigen::VectorXd> subPatternParamsResults;
    subPatternParamsResults.resize(num_pattern_vars());
    for (size_t pi = 0; pi < num_pattern_vars(); ++pi) {
        subPatternParamsResults[pi] = Eigen::Map<Eigen::VectorXd>(subPatternParams[pi].data(), subPatternParams[pi].size());
    }
    std::get<2>(result) = subPatternParamsResults;
    return result;
}

}