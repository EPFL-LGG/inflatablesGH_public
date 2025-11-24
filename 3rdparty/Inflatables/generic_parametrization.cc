#include "generic_parametrization.hh"

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
RegularizedGenericParametrizer::RegularizedGenericParametrizer(LocalGlobalGenericParametrizer &lggparam)
    : Parametrizer(lggparam.meshPtr()),
      dualLaplacianStencil(*lggparam.meshPtr()),
      m_alphaMin(lggparam.alphaMin()),
      m_alphaMax(lggparam.alphaMax()),
      m_betaMin(lggparam.betaMin()),
      m_betaMax(lggparam.betaMax())
{
    const size_t nt = mesh().numTris();
    m_phi.resize(nt);
    m_psi.resize(nt);
    m_stretch.resize(2*nt);

    // Initialize the variable fields from the local-global parametrizer
    for (size_t i = 0; i < nt; ++i) {
        m_phi[i]   = lggparam. leftStretchAngle(i);
        m_psi[i]   = lggparam.rightStretchAngle(i);
        m_stretch[i] = lggparam.getAlpha(i);
        m_stretch[nt + i] = lggparam.getBeta(i);
    }

    setUV(lggparam.uv());

    // Cache the (constant) Laplacian block of the Hessian.
    m_laplacian = SuiteSparseMatrix(Laplacian::construct(mesh()));

    dualLaplacianStencil.type = DualLaplacianStencil<Mesh>::Type::DualMeshIDT;

    // Construct per-triangle averaged shape operator
    {
        const auto &m = mesh();
        CurvatureInfo cinfo(mesh());

        m_shapeOperators.reserve(m.numTris());
        m_shapeOperators.clear();

        for (auto tri : m.tris()) {
            const size_t ti = tri.index();
            M2d S(M2d::Zero());
            for (auto corner : tri.vertices()) {
                M2d d;
                V2d k(cinfo.kappa_1[corner.index()],
                      cinfo.kappa_2[corner.index()]);
                for (size_t i = 0; i < 2; ++i) {
                    // Project curvature direction onto triangle tangent plane and re-normalize.
                    d.col(i) = (cinfo.d(i).row(corner.index()) * m_B[ti]).normalized().eval();
                }
                S += d * k.asDiagonal() * d.transpose();
            }
            S /= 3.0;
            m_shapeOperators.push_back(S);
        }
    }
}


void RegularizedGenericParametrizer::m_evalIterate() {
    const size_t nt = mesh().numTris();
    m_M.resize(nt);
    m_U.resize(nt);
    m_V.resize(nt);
    m_dU_dphi.resize(nt);
    m_dV_dpsi.resize(nt);

    for (const auto &tri : mesh().elements()) {
        const size_t ti = tri.index();
        auto &M = m_M[ti];
        auto &U = m_U[ti];
        auto &V = m_V[ti];
        auto &dU_dphi = m_dU_dphi[ti];
        auto &dV_dpsi = m_dV_dpsi[ti];

        U = Eigen::Rotation2D<Real>(m_phi[ti]).matrix();
        V = Eigen::Rotation2D<Real>(m_psi[ti]).matrix();

        dU_dphi = Eigen::Rotation2D<Real>(m_phi[ti] + M_PI / 2).matrix();
        dV_dpsi = Eigen::Rotation2D<Real>(m_psi[ti] + M_PI / 2).matrix();

        M = U * Vector2D(m_stretch[ti], m_stretch[nt + ti]).asDiagonal() * V.transpose();
    }
}

Real RegularizedGenericParametrizer::energy(EnergyType etype) const {
    // Accumulate energy contributions in temporaries so that other threads don't read intermediate values.
    Real fittingEnergy = 0,
         stretchRegEnergy = 0,
         phiRegEnergy = 0,
         diffRegEnergy = 0,
         phiRegEnergyWeighted = 0;
    
    const size_t nt = mesh().numTris();
    
    for (const auto &tri : mesh().elements()) {
        const size_t ti = tri.index();
        const M2d JB = m_J[ti] * m_B[ti];
        fittingEnergy += (JB - m_M[ti]).squaredNorm() * tri->volume();
    }
    // Edge-based regularization terms
    // Dual Laplacian-based regularization terms
    dualLaplacianStencil.visit_edges([this, &phiRegEnergy, &stretchRegEnergy, &diffRegEnergy, nt, etype](size_t i, size_t j, Real w_ij) {
        const size_t ti = i,
                     tj = j;

        
        const Real diff_i = m_stretch[ti] - m_stretch[nt + ti];
        const Real diff_j = m_stretch[tj] - m_stretch[nt + tj];
        std::vector<Real> barrier_data;
        m_barrier(diff_i, barrier_data, BarrierEvalLevel::Val);
        const Real barrier_i = barrier_data[0];
        m_barrier(diff_j, barrier_data, BarrierEvalLevel::Val);
        const Real barrier_j = barrier_data[0];
        if (((etype == EnergyType::Full) || (etype == EnergyType::PhiRegularization)) && m_phi_reg_w != 0.0) {
            if (!m_useBarrier)
                phiRegEnergy += w_ij * std::pow(std::abs(sin(m_phi[ti] - m_phi[tj])), m_phi_reg_p);
            else {
                // Barrier for umbilical points
                phiRegEnergy += w_ij * std::pow(std::abs(sin(m_phi[ti] - m_phi[tj])), m_phi_reg_p) * barrier_i * barrier_j;
            }
        }
        if (((etype == EnergyType::Full) || (etype == EnergyType::DiffRegularization)) && m_diff_reg_w != 0.0) {
            diffRegEnergy += w_ij * (std::exp( -std::pow(diff_i, 2.0)) + std::exp( -std::pow(diff_j, 2.0)));
            // diffRegEnergy += - w_ij * std::pow(diff_i, 2.0) - std::pow(diff_j, 2.0);
        }
        if (m_variableStretch && ((etype == EnergyType::Full) || (etype == EnergyType::StretchRegularization)) && m_stretch_reg_w != 0.0) {
            stretchRegEnergy += w_ij * std::pow(std::abs(m_stretch[ti] - m_stretch[tj]), m_stretch_reg_p);
            stretchRegEnergy += w_ij * std::pow(std::abs(m_stretch[nt + ti] - m_stretch[nt + tj]), m_stretch_reg_p);
        }
    }
    );

    if (etype != EnergyType::Full) {
        if (etype != EnergyType::Fitting              ) fittingEnergy = 0.0;
        if (etype != EnergyType::StretchRegularization) stretchRegEnergy = 0.0;
        if (etype != EnergyType::PhiRegularization    ) phiRegEnergy = 0.0;
        if (etype != EnergyType::DiffRegularization   ) diffRegEnergy = 0.0;
    }
    if (scaleInvariantFittingEnergy) fittingEnergy /= mesh().volume();
    return (0.5 * fittingEnergy) + (m_stretch_reg_w / m_stretch_reg_p) * stretchRegEnergy + (m_phi_reg_w / m_phi_reg_p) * (phiRegEnergy) + diffRegEnergy * m_diff_reg_w;
}



Eigen::VectorXd RegularizedGenericParametrizer::gradient(EnergyType etype) const {
    Eigen::VectorXd result(numVars());
    result.setZero();
    Eigen::Map<UVMap> grad_uv(result.data(), m_uv.rows(), m_uv.cols());

    const size_t nt = mesh().numTris();
        

    // Gradient of fitting energy
    for (const auto &tri : mesh().elements()) {
        const size_t ti = tri.index();

        if ((etype == EnergyType::Full) || (etype == EnergyType::Fitting)) {
            const M2d JB = m_J[ti] * m_B[ti];
            M2d scaled_dist = tri->volume() * (JB - m_M[ti]);
            if (scaleInvariantFittingEnergy) scaled_dist /= mesh().volume();

            // Gradient wrt parametrization:
            for (const auto &v : tri.vertices())
                grad_uv.row(v.index()) += (scaled_dist * m_B[ti].transpose()) * tri->gradBarycentric().col(v.localIndex());

            Vector2D tgt_sigma(m_stretch[ti], m_stretch[nt + ti]);
            // Gradient wrt rotations:
            result[phiOffset() + ti] -= (scaled_dist.transpose() * m_dU_dphi[ti] * tgt_sigma.asDiagonal() *       m_V[ti].transpose()).trace();
            result[psiOffset() + ti] -= (scaled_dist.transpose() *       m_U[ti] * tgt_sigma.asDiagonal() * m_dV_dpsi[ti].transpose()).trace();

            if (!m_variableStretch) continue;

            result[alphaOffset() + ti] -= ((scaled_dist.transpose() * m_U[ti].col(0)) * m_V[ti].col(0).transpose()).trace();
            result[betaOffset() + ti] -= ((scaled_dist.transpose() * m_U[ti].col(1)) * m_V[ti].col(1).transpose()).trace();
        }
    }

    // Gradient of edge-based regularization terms
    dualLaplacianStencil.visit_edges([this, etype, &result, nt](size_t ti, size_t tj, Real w_ij) {    
        // Phi regularization
        if (((etype == EnergyType::Full) || (etype == EnergyType::PhiRegularization)) && m_phi_reg_w != 0.0) {
            Real phi_diff = m_phi[ti] - m_phi[tj];
            Real s = sin(phi_diff),
                c = cos(phi_diff);
            // Using std::copysign(1.0, s) doesn't work since it gives bad derivatives around phi_diff = 0.
            // We get better results explicitly setting the derivative equal to zero in this case.
            Real sign = 0.0;
            if (s > 0) sign =  1.0;
            if (s < 0) sign = -1.0;

            // Real val = m_phi_reg_w * std::pow(1 - cos(2 * phi_diff), m_phi_reg_p * 0.5 - 1.0) * sin(2 * phi_diff);

            Real val;
            if (m_phi_reg_p == 1.0) { val = m_phi_reg_w * c * sign; }
            else                    { val = m_phi_reg_w * std::pow(std::abs(s), m_phi_reg_p - 1.0) * c * sign; } // This is well-behaved for p > 0 (finite, non-nan value)
            val *= w_ij;

            
            
            const Real diff_i = m_stretch[ti] - m_stretch[nt + ti];
            const Real diff_j = m_stretch[tj] - m_stretch[nt + tj];
            std::vector<Real> barrier_data;
            m_barrier(diff_i, barrier_data, BarrierEvalLevel::Grad);
            const Real barrier_i = barrier_data[0], barrier_grad_i = barrier_data[1];
            m_barrier(diff_j, barrier_data, BarrierEvalLevel::Grad);
            const Real barrier_j = barrier_data[0], barrier_grad_j = barrier_data[1];
            if (!m_useBarrier) {
                result[phiOffset() + ti] += val;
                result[phiOffset() + tj] -= val;
            }
            else {

                result[phiOffset() + ti] += val * barrier_i * barrier_j;
                result[phiOffset() + tj] -= val * barrier_i * barrier_j;

                const Real phiRegEnergy = w_ij * (m_phi_reg_w / m_phi_reg_p) * std::pow(std::abs(s), m_phi_reg_p);

                result[stretchOffset() + ti] += phiRegEnergy * barrier_grad_i * barrier_j;
                result[stretchOffset() + nt + ti] -= phiRegEnergy * barrier_grad_i * barrier_j;
                result[stretchOffset() + tj] += phiRegEnergy * barrier_grad_j * barrier_i;
                result[stretchOffset() + nt + tj] -= phiRegEnergy * barrier_grad_j * barrier_i;
            }                
        }

        if (((etype == EnergyType::Full) || (etype == EnergyType::DiffRegularization)) && m_diff_reg_w != 0.0) {
            const Real diff_i = m_stretch[ti] - m_stretch[nt + ti];
            const Real diff_j = m_stretch[tj] - m_stretch[nt + tj];
            const Real diff_grad_i = - w_ij * 2.0 * diff_i * std::exp( -std::pow(diff_i, 2.0)) * m_diff_reg_w;
            const Real diff_grad_j = - w_ij * 2.0 * diff_j * std::exp( -std::pow(diff_j, 2.0)) * m_diff_reg_w;
            
            // const Real diff_grad_i = - 2.0 * diff_i;
            // const Real diff_grad_j = - 2.0 * diff_j;
            result[stretchOffset() + ti] += diff_grad_i;
            result[stretchOffset() + nt + ti] -= diff_grad_i;
            result[stretchOffset() + tj] += diff_grad_j;
            result[stretchOffset() + nt + tj] -= diff_grad_j;
        }

        // Stretch regularization
        if (m_variableStretch && ((etype == EnergyType::Full) || (etype == EnergyType::StretchRegularization)) && m_stretch_reg_w != 0.0) {
            Real stretch_diff = m_stretch[ti] - m_stretch[tj];

            // Using std::copysign(1.0, stretch_diff) doesn't work since it gives bad derivatives around stretch_diff = 0.
            // We get better results explicitly setting the derivative equal to zero in this case.
            Real sign = 0.0;
            if (stretch_diff > 0) sign =  1.0;
            if (stretch_diff < 0) sign = -1.0;

            Real val;
            if (m_stretch_reg_p == 1.0) { val = m_stretch_reg_w * sign; }
            else                      { val = m_stretch_reg_w * std::pow(std::abs(stretch_diff), m_stretch_reg_p - 1.0) * sign; }
            val *= w_ij;

            result[stretchOffset() + ti] += val;
            result[stretchOffset() + tj] -= val;



            stretch_diff = m_stretch[nt + ti] - m_stretch[nt + tj];
            if (stretch_diff > 0) sign =  1.0;
            if (stretch_diff < 0) sign = -1.0;
            if (m_stretch_reg_p == 1.0) { val = m_stretch_reg_w * sign; }
            else                      { val = m_stretch_reg_w * std::pow(std::abs(stretch_diff), m_stretch_reg_p - 1.0) * sign; }
            val *= w_ij;

            result[stretchOffset() + nt + ti] += val;
            result[stretchOffset() + nt + tj] -= val;
        }
    }
    );
    return result;
}

SuiteSparseMatrix RegularizedGenericParametrizer::hessianSparsityPattern(Real val) const {
    SuiteSparseMatrix result(numVars(), numVars());
    result.symmetry_mode = SuiteSparseMatrix::SymmetryMode::UPPER_TRIANGLE;
    result.Ap.reserve(numVars() + 1);

    auto &Ap = result.Ap;
    auto &Ai = result.Ai;

    auto addIdx = [&](const size_t idx) { Ai.push_back(idx); };

    auto finalizeCol = [&]() {
        const size_t colStart = Ap.back();
        const size_t colEnd = Ai.size();
        Ap.push_back(colEnd);
        std::sort(Ai.begin() + colStart, Ai.begin() + colEnd);
    };

    // Build the sparsity pattern in compressed form one column (variable) at a time.
    result.Ap.push_back(0);

    // Laplacian blocks: each vertex value interacts with itself and its neighbors
    const auto &m = mesh();
    const size_t nv = m.numVertices();
    const size_t nt = mesh().numTris();
    for (size_t uvo = 0; uvo < 2; ++uvo) { // 0: u variables, 1: v variables
        for (const auto &v : m.vertices()) {
            size_t vi = v.index() + uvo * nv;
            addIdx(vi);
            for (const auto &he : v.incidentHalfEdges()) {
                size_t ui = he.tail().index() + uvo * nv;
                if (ui < vi) addIdx(ui);
            }
            finalizeCol();
        }
    }

    
    // Tri field columns: interact with corner vertices, neighbors, and selves
    const size_t phio = phiOffset(), psio = psiOffset(), alphao = alphaOffset(), betao = betaOffset();
    const size_t numTriFields = variableStretch() ? 4 : 2;
    for (size_t fieldOffset = 0; fieldOffset < numTriFields; ++fieldOffset) { // 0: phi variables, 1: psi variables, 2: alpha variables, 3: beta variables
        for (const auto &tri : m.elements()) {
            const size_t tj = tri.index();
            for (const auto &v : tri.vertices()) {
                addIdx(v.index());      // u variable
                addIdx(v.index() + nv); // v variable
            }
            addIdx(phio + tj); // phi-phi/phi-psi/phi-alpha/phi-beta interaction
            if (fieldOffset > 0) addIdx(psio   + tj); // psi-psi/psi-alpha/psi-beta interaction
            if (fieldOffset > 1) addIdx(alphao + tj); // alpha-alpha/alpha-beta interaction
            if (fieldOffset > 2) addIdx(betao + tj); // beta-beta interaction

            // Laplacian-style regularization (upper triangle)
            dualLaplacianStencil.visit(tri.index(), [this, &m, &addIdx, fieldOffset, phio, alphao, betao, tj](size_t /* j */, size_t i, Real /* w_ij */) {
                auto tri_i = m.tri(i);
                const size_t ti = tri_i.index();
                if (!tri_i) return;
                if (m_useBarrier) {
                    //Barrier interactions
                    if (fieldOffset > 1) { addIdx(phio + ti);}
                    if (fieldOffset > 2) { addIdx(alphao + ti);}    
                }
                if (ti > tj) return;
                
                if (fieldOffset == 0) { addIdx(phio   + ti); } //   phi regularization interaction
                
                if (fieldOffset == 2) { addIdx(alphao + ti); } // alpha regularization interaction
                if (fieldOffset == 3) { addIdx(betao + ti); } // beta regularization interaction
            });

            finalizeCol();
        }
    }

    result.nz = result.Ai.size();
    result.Ax.assign(result.nz, val);
    return result;
}

SuiteSparseMatrix RegularizedGenericParametrizer::hessian(EnergyType etype) const {
    SuiteSparseMatrix H = hessianSparsityPattern();
    hessian(H, etype);
    return H;
}

void RegularizedGenericParametrizer::hessian(SuiteSparseMatrix &H, EnergyType etype) const {
    const size_t uo = uOffset(),
                 vo = vOffset(),
                 phio   = phiOffset(),
                 psio   = psiOffset(),
                 alphao = alphaOffset(),
                 betao = betaOffset();

    const size_t nt = mesh().numTris();
    if ((etype == EnergyType::Full) || (etype == EnergyType::Fitting)) {
        // u-u, v-v (Laplacian)
        for (const auto &entry : m_laplacian) {
            auto val = entry.v;
            if (scaleInvariantFittingEnergy) val /= mesh().volume();
            H.addNZ(uo + entry.i, uo + entry.j, val);
            H.addNZ(vo + entry.i, vo + entry.j, val);
        }
    }
    for (const auto &tri : mesh().elements()) {
        const size_t ti = tri.index();
        if ((etype == EnergyType::Full) || (etype == EnergyType::Fitting)) {
            Real A = tri->volume();
            if (scaleInvariantFittingEnergy) A /= mesh().volume();
            Vector2D tgt_sigma(m_stretch[ti], m_stretch[nt + ti]);

            // target_fields-u, target_fields-v
            M2d dM_dphi   = m_dU_dphi[ti] * tgt_sigma.asDiagonal() *       m_V[ti].transpose(),
                dM_dpsi   =       m_U[ti] * tgt_sigma.asDiagonal() * m_dV_dpsi[ti].transpose(),
                dM_dalpha = m_U[ti].col(0) * m_V[ti].col(0).transpose(),
                dM_dbeta = m_U[ti].col(1) * m_V[ti].col(1).transpose();
            for (const auto &v : tri.vertices()) {
                Vector2D dE_duv_dphi = -A * ((dM_dphi * m_B[ti].transpose()) * tri->gradBarycentric().col(v.localIndex()));
                H.addNZ(uo + v.index(), phio + tri.index(), dE_duv_dphi[0]);
                H.addNZ(vo + v.index(), phio + tri.index(), dE_duv_dphi[1]);

                Vector2D dE_duv_dpsi = -A * ((dM_dpsi * m_B[ti].transpose()) * tri->gradBarycentric().col(v.localIndex()));
                H.addNZ(uo + v.index(), psio + ti, dE_duv_dpsi[0]);
                H.addNZ(vo + v.index(), psio + ti, dE_duv_dpsi[1]);

                Vector2D dE_duv_dalpha = -A * ((dM_dalpha * m_B[ti].transpose()) * tri->gradBarycentric().col(v.localIndex()));
                H.addNZ(uo + v.index(), alphao + ti, dE_duv_dalpha[0]);
                H.addNZ(vo + v.index(), alphao + ti, dE_duv_dalpha[1]);

                Vector2D dE_duv_dbeta = -A * ((dM_dbeta * m_B[ti].transpose()) * tri->gradBarycentric().col(v.localIndex()));
                H.addNZ(uo + v.index(), betao + ti, dE_duv_dbeta[0]);
                H.addNZ(vo + v.index(), betao + ti, dE_duv_dbeta[1]);
            }

            M2d dist = m_J[ti] * m_B[ti] - m_M[ti];
            M2d d2M_dphi_dpsi = m_dU_dphi[ti] * tgt_sigma.asDiagonal() * m_dV_dpsi[ti].transpose();

            // psi-phi
            H.addNZ(phio + ti, psio + ti, A * ((dM_dphi.transpose() * dM_dpsi).trace() - (dist.transpose() * d2M_dphi_dpsi).trace()));

            // phi-phi, psi-psi
            // Note: d^2U/dphi^2 = -U, so d^2M/dphi^2 = -M = d^2M/dpsi^2
            Real dist_contract_neg_d2M_dangle2 = (dist.transpose() * m_M[ti]).trace();
            H.addNZ(phio + ti, phio + ti, A * (dM_dphi.squaredNorm() + dist_contract_neg_d2M_dangle2));
            H.addNZ(psio + ti, psio + ti, A * (dM_dpsi.squaredNorm() + dist_contract_neg_d2M_dangle2));

            if (!variableStretch()) continue;

            M2d d2M_dphi_dalpha = m_dU_dphi[ti].col(0) * m_V[ti].col(0).transpose(),
                d2M_dpsi_dalpha = m_U[ti].col(0) * m_dV_dpsi[ti].col(0).transpose();
            M2d d2M_dphi_dbeta = m_dU_dphi[ti].col(1) * m_V[ti].col(1).transpose(),
                d2M_dpsi_dbeta = m_U[ti].col(1) * m_dV_dpsi[ti].col(1).transpose();

            // phi-alpha, psi-alpha
            H.addNZ(phio + ti, alphao + ti, A * ((dM_dphi.transpose() * dM_dalpha).trace() - (dist.transpose() * d2M_dphi_dalpha).trace()));
            H.addNZ(psio + ti, alphao + ti, A * ((dM_dpsi.transpose() * dM_dalpha).trace() - (dist.transpose() * d2M_dpsi_dalpha).trace()));

            // phi-beta, psi-beta
            H.addNZ(phio + ti, betao + ti, A * ((dM_dphi.transpose() * dM_dbeta).trace() - (dist.transpose() * d2M_dphi_dbeta).trace()));
            H.addNZ(psio + ti, betao + ti, A * ((dM_dpsi.transpose() * dM_dbeta).trace() - (dist.transpose() * d2M_dpsi_dbeta).trace()));

            // alpha-alpha (note d2M_dalpha_dalpha = 0)
            H.addNZ(alphao + ti, alphao + ti, A * dM_dalpha.squaredNorm());

            // beta-beta (note d2M_dbeta_dbeta = 0)
            H.addNZ(betao + ti, betao + ti, A * dM_dbeta.squaredNorm());

            // alpha-beta (note d2M_dalpha_dbeta = 0)
            H.addNZ(alphao + ti, betao + ti, A * (dM_dalpha.transpose() * dM_dbeta).trace());
        }
    }

    dualLaplacianStencil.visit_edges([this, etype, &H, nt, phio, alphao, betao](size_t ti, size_t tj, Real w_ij) {
        if (ti > tj) std::swap(ti, tj); // Visit each stencil edge (unordered triangle pair) exactly once
        
        if (((etype == EnergyType::Full) || (etype == EnergyType::PhiRegularization)) && m_phi_reg_w != 0.0) {
            Real phi_diff = m_phi[ti] - m_phi[tj];
            // Real c = cos(2 * phi_diff),
            //      s = sin(2 * phi_diff);
            // Real val = std::pow(1 - c, m_phi_reg_p * 0.5 - 2.0) * s * s * (m_phi_reg_p * 0.5 - 1.0) +
            //            std::pow(1 - c, m_phi_reg_p * 0.5 - 1.0) * c;
            // val *= 2 * m_phi_reg_w;

            Real s = sin(phi_diff),
                c = cos(phi_diff);
            // Using std::copysign(1.0, s) doesn't work since it gives bad derivatives around phi_diff = 0.
            // We get better results explicitly setting the derivative equal to zero in this case.
            Real sign = 0.0;
            if (s > 0) sign =  1.0;
            if (s < 0) sign = -1.0;

            // Note: second derivatives in both the "p = 1" and "p = 2" cases are well behaved,
            // but they blow up around phi_diff = 0 when "1 < p < 2". We discard Hessian
            // contributions near this blowup.
            Real val = 0.0;
            const Real p = m_phi_reg_p;
            if      (m_phi_reg_p == 1.0) { val = -std::abs(s); }
            else if (m_phi_reg_p == 2.0) { val = c * c - s * s; } // equivalently: cos(2 * phi_diff)
            else if (std::abs(s) > 1e-4) { val = (p - 1) * std::pow(std::abs(s), p - 2.0) * c * c
                                                        - std::pow(std::abs(s), p - 1.0) * s * sign; }
            else                         { val = 0.0; } // discard Hessian contributions in cases that blow up

            val *= m_phi_reg_w * w_ij;

            const Real diff_i = m_stretch[ti] - m_stretch[nt + ti];
            const Real diff_j = m_stretch[tj] - m_stretch[nt + tj];
            std::vector<Real> barrier_data;
            m_barrier(diff_i, barrier_data, BarrierEvalLevel::Hess);
            const Real barrier_i = barrier_data[0], barrier_grad_i = barrier_data[1], barrier_hess_i = barrier_data[2];
            m_barrier(diff_j, barrier_data, BarrierEvalLevel::Hess);
            const Real barrier_j = barrier_data[0], barrier_grad_j = barrier_data[1], barrier_hess_j = barrier_data[2];

            if (!m_useBarrier) {
                if (val != 0.0) {
                    H.addNZ(phio + ti, phio + tj, -val);
                    H.addNZ(phio + tj, phio + tj,  val);
                    H.addNZ(phio + ti, phio + ti,  val);
                    // H.addNZ(phio + tj, phio + ti, -val); (lower triangle)
                }
            }
            else {
                    // Barrier for umbilical points
                    const Real phiRegEnergy = w_ij * (m_phi_reg_w / m_phi_reg_p) * std::pow(std::abs(s), m_phi_reg_p);

                    if (val * barrier_i * barrier_j != 0.0) {
                        H.addNZ(phio + ti, phio + tj, -val * barrier_i * barrier_j);
                        H.addNZ(phio + tj, phio + tj,  val * barrier_i * barrier_j);
                        H.addNZ(phio + ti, phio + ti,  val * barrier_i * barrier_j);
                        // H.addNZ(phio + tj, phio + ti, -val * barrier_i * barrier_j); (lower triangle)
                    }
                    
                    val = phiRegEnergy;
                    if (val * barrier_hess_i * barrier_j != 0.0) {
                        H.addNZ(alphao + ti, alphao + ti, val * barrier_hess_i * barrier_j);
                        H.addNZ(betao + ti, betao + ti, val * barrier_hess_i * barrier_j);
                        H.addNZ(alphao + ti, betao + ti, -val * barrier_hess_i * barrier_j);
                        // H.addNZ(betao + ti, alphao + ti, -val * barrier_hess_i * barrier_j); (lower triangle)
                    }
                    if (val * barrier_hess_j * barrier_i != 0.0) {
                        H.addNZ(alphao + tj, alphao + tj, val * barrier_hess_j * barrier_i);
                        H.addNZ(betao + tj, betao + tj, val * barrier_hess_j * barrier_i);
                        H.addNZ(alphao + tj, betao + tj, -val * barrier_hess_j * barrier_i);
                        // H.addNZ(betao + tj, alphao + tj, -val * barrier_hess_j * barrier_i); (lower triangle)
                    }
                    val *= barrier_grad_i;
                    val *= barrier_grad_j;
                    if (val != 0.0) {
                        H.addNZ(alphao + ti, alphao + tj, val);
                        H.addNZ(betao + ti, betao + tj, val);
                        H.addNZ(alphao + ti, betao + tj, -val);
                        H.addNZ(alphao + tj, betao + ti, -val);  
                    }
                    if (m_phi_reg_p == 1.0) { val = m_phi_reg_w * c * sign; }
                    else                    { val = m_phi_reg_w * std::pow(std::abs(s), m_phi_reg_p - 1.0) * c * sign; }
                    val *= w_ij;
                    
                    if (val * barrier_grad_i * barrier_j != 0.0) {
                        H.addNZ(phio + ti, alphao + ti, val * barrier_grad_i * barrier_j);
                        H.addNZ(phio + tj, alphao + ti, -val * barrier_grad_i * barrier_j);
                        H.addNZ(phio + ti, betao + ti, -val * barrier_grad_i * barrier_j);
                        H.addNZ(phio + tj, betao + ti, val * barrier_grad_i * barrier_j);
                    }
                    if (val * barrier_grad_j * barrier_i != 0.0) {
                        H.addNZ(phio + ti, alphao + tj, val * barrier_grad_j * barrier_i);
                        H.addNZ(phio + tj, alphao + tj, -val * barrier_grad_j * barrier_i);
                        H.addNZ(phio + ti, betao + tj, -val * barrier_grad_j * barrier_i);
                        H.addNZ(phio + tj, betao + tj, val * barrier_grad_j * barrier_i);
                    }
            }
        }
        if (((etype == EnergyType::Full) || (etype == EnergyType::DiffRegularization)) && (m_diff_reg_w != 0.0)) {
            const Real diff_i = m_stretch[ti] - m_stretch[nt + ti];
            const Real diff_j = m_stretch[tj] - m_stretch[nt + tj];
            const Real diff_hess_i = w_ij * (4.0 * diff_i * diff_i - 2.0 )* std::exp( -std::pow(diff_i, 2.0)) * m_diff_reg_w;
            const Real diff_hess_j = w_ij * (4.0 * diff_j * diff_j - 2.0 )* std::exp( -std::pow(diff_j, 2.0)) * m_diff_reg_w;
            // const Real diff_hess_i = -2.0;
            // const Real diff_hess_j = -2.0;
            H.addNZ(alphao + ti, alphao + ti, diff_hess_i);
            H.addNZ(betao + ti, betao + ti, diff_hess_i);
            H.addNZ(alphao + ti, betao + ti, - diff_hess_i);

            H.addNZ(alphao + tj, alphao + tj, diff_hess_j);
            H.addNZ(betao + tj, betao + tj, diff_hess_j);
            H.addNZ(alphao + tj, betao + tj, - diff_hess_j);
        }
        
        if (((etype == EnergyType::Full) || (etype == EnergyType::StretchRegularization)) && (m_stretch_reg_w != 0.0)) {
            const Real alpha_diff = m_stretch[ti] - m_stretch[tj], beta_diff = m_stretch[nt + ti] - m_stretch[nt + tj];
            if (!m_variableStretch || (m_stretch_reg_p == 1.0)) return;
            
            if (!((m_stretch_reg_p < 2.0) && (std::abs(alpha_diff) < 1e-14))) {
                Real val = w_ij * m_stretch_reg_w * (m_stretch_reg_p - 1.0) * std::pow(std::abs(alpha_diff), m_stretch_reg_p - 2.0);
                H.addNZ(alphao + ti, alphao + tj, -val);
                H.addNZ(alphao + tj, alphao + tj,  val);

                // H.addNZ(alphao + tj, alphao + ti, -val); (lower triangle)
                H.addNZ(alphao + ti, alphao + ti,  val);
            }
            if (!((m_stretch_reg_p < 2.0) && (std::abs(beta_diff) < 1e-14))) {
                Real val = w_ij * m_stretch_reg_w * (m_stretch_reg_p - 1.0) * std::pow(std::abs(beta_diff), m_stretch_reg_p - 2.0);
                H.addNZ(betao + ti, betao + tj, -val);
                H.addNZ(betao + tj, betao + tj,  val);

                // H.addNZ(betao + tj, betao + ti, -val); (lower triangle)
                H.addNZ(betao + ti, betao + ti,  val);
            }
        }
    }
    );

}

}