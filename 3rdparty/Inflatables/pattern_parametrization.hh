#ifndef PATTERNPARAMETRIZATION_HH
#define PATTERNPARAMETRIZATION_HH


#include "generic_parametrization.hh"
#include "BendingStiffnessIntegralSensitivityPattern.hh"
#include "ConstraintBarrier.hh"

namespace parametrization {

////////////////////////////////////////////////////////////////////////////////
// RegularizedPatternParametrizer: Global nonlinear energy with auxiliary variables and pattern parameters
////////////////////////////////////////////////////////////////////////////////
// Perform the nonlinear minimization:
//      min 0.5 * ||grad f - M Bt||^2 + w_phi ||grad phi||_{p_phi} + w_alpha ||grad alpha||_{p_alpha}
// with a Newton-type method, where alpha, beta and phi are the target metric
// stretching factors/orientation.
// This relies on a good initialization (e.g., computed from a local-global method).
// The target metric M is of the form:
//      U(phi) diag(alpha, beta) V(psi)^T
// Where, e.g., U(phi) := [cos(phi) -sin(phi); sin(phi) cos(phi)]
// Use the pattern parameters as auxiliary variables to compute the distortion at each point.
struct RegularizedPatternParametrizer : public RegularizedGenericParametrizer {
    using Base = RegularizedGenericParametrizer;
    using VXd = Eigen::VectorXd;
    using MatInfoFunction = std::vector<std::function<VXd(std::vector<VXd>)>>;
    using MX2d = Eigen::Matrix<Real, Eigen::Dynamic, 2>;

    static constexpr size_t num_material_vars = 9; // two max / min stretching factors, two x, y stretching factors (F_bar inv; but since we only consider patches with reflectional symmetry, we only have two parameters),  five stiffness coefficients

    enum class PatternEnergyType { Full, RGP, Bending, PatternRegularization, PatternBoundConstraint, DEBUG_Fitting, DEBUG_PhiRegularization, DEBUG_StretchRegularization};

    size_t num_info_per_mat_var = 3; // objective, gradient, hessian: three lists of functions (depending on the number of pattern parameters, the number of partial gradient functions and hessian functions are different)
    size_t num_total_info = num_material_vars * num_info_per_mat_var;

    size_t num_stiffness_var = 5; 
    size_t num_non_stiffness_var = num_material_vars - num_stiffness_var;

    // Initialize from the local-global parametrizer
    RegularizedPatternParametrizer(LocalGlobalGenericParametrizer &lggparam, MatInfoFunction newMaterialInfo, VXd default_pattern_params, size_t input_num_pattern_vars);

    size_t patternOffset() const { return Base::stretchOffset(); }
    size_t numVars() const { return patternOffset() + m_pattern_params.rows(); }
    
    size_t num_pattern_vars() const { return m_num_pattern_vars; }

    VXd getVars() const {
        BENCHMARK_SCOPED_TIMER_SECTION timer("RegularizedPatternParametrizer.getVars");
        VXd result(numVars());
        VXd rgp_vars = Base::getVars();
        // Copy over everything except the stretch variables at the end of rgp_vars.
        result.head(patternOffset()) = rgp_vars.head(patternOffset());
        result.segment(patternOffset(), m_pattern_params.rows()) = m_pattern_params;
        return result;
    }

    void setVars(const VXd &vars) {
        BENCHMARK_SCOPED_TIMER_SECTION timer("RegularizedPatternParametrizer.setVars");
        if (size_t(vars.rows()) != numVars()) throw std::runtime_error("Invalid variable count");
        // Set rgp vars.
        VXd rgp_vars(Base::numVars());
        rgp_vars.head(patternOffset()) = vars.head(patternOffset());
        m_pattern_params = vars.segment(patternOffset(), num_pattern_vars() * mesh().numTris());
        m_set_matInfoArgs();
        rgp_vars.tail(Base::numStretchVars()) = get_stretch_from_pattern_params();
        Base::setVars(rgp_vars);
        // Set pattern params.
        m_bendingSensitivityCache.clear();
    }

    VXd get_stretch_from_pattern_params() const {
        const size_t nt = mesh().numTris();
        VXd result = VXd::Zero(Base::numStretchVars());
        if (nt * 2 != Base::numStretchVars()) throw std::runtime_error("Invalid stretch variable count");
        result.head(nt) = materialInfo[0](getMatInfoArgs());
        result.tail(nt) = materialInfo[num_info_per_mat_var](getMatInfoArgs());
        return result;
    }

    VXd get_stretch_angle_offset_from_pattern_params(std::vector<VXd> query_pattern_params) const {
        const size_t np = query_pattern_params[0].size();
        // std::cout<<"materialInfo.size() = "<<materialInfo.size()<<std::endl;
        VXd x_scale_factors = materialInfo[num_info_per_mat_var * 2](query_pattern_params);
        VXd y_scale_factors = materialInfo[num_info_per_mat_var * 3](query_pattern_params);

        VXd result = VXd::Zero(np);
        for (size_t i = 0; i < np; ++i) {
            if (x_scale_factors[i] > y_scale_factors[i]) result[i] = 0.0;
            else result[i] = M_PI / 2.0;
        }
        return result;
    }

    VXd get_pattern_params_gradient(const VXd &stretch_gradient) const {
        const size_t nt = mesh().numTris();
        if (size_t(stretch_gradient.rows()) != nt * 2) throw std::runtime_error("Invalid stretch gradient size");
        VXd alpha_grad = stretch_gradient.head(nt);
        VXd beta_grad = stretch_gradient.tail(nt);
        VXd result = VXd::Zero(num_pattern_vars() * nt);

        VXd d_alpha_dp = materialInfo[1](getMatInfoArgs());
        VXd d_beta_dp = materialInfo[1 + num_info_per_mat_var](getMatInfoArgs());

        for (size_t i = 0; i < num_pattern_vars(); ++i) {
            result.segment(i * nt, nt) = d_alpha_dp.segment(i * nt, nt).cwiseProduct(alpha_grad) + d_beta_dp.segment(i * nt, nt).cwiseProduct(beta_grad);
        }
        return result;
    }

    virtual Real energy() const override { return energy(PatternEnergyType::Full); }
    Real energy(PatternEnergyType etype) const;
    

    VXd gradient(PatternEnergyType etype = PatternEnergyType::Full) const;

    size_t hessianNNZ() const { return hessianSparsityPattern().nz; } // TODO: predict without constructing

    SuiteSparseMatrix baseHessianSparsityPattern() const;

    SuiteSparseMatrix hessianSparsityPattern(Real val = 0.0) const;

    void              hessian(SuiteSparseMatrix &H, PatternEnergyType etype = PatternEnergyType::Full) const; // accumulate Hessian to H
    SuiteSparseMatrix hessian(                      PatternEnergyType etype = PatternEnergyType::Full) const; // construct and return Hessian
    void              hessian(SuiteSparseMatrix &H, bool /* projectionMask */) const { hessian(H, PatternEnergyType::Full); }


    // Regularization parameters
    void setPatternRegW(Real val) { m_pattern_reg_w = val; }
    void setPatternRegP(Real val) { m_pattern_reg_p = val; }
    void    setBendRegW(Real val) {     m_bending_w = val; }

    Real patternRegW() const { return m_pattern_reg_w; }
    Real patternRegP() const { return m_pattern_reg_p; }
    Real    bendRegW() const { return m_bending_w; }

    void setPhiRegW(Real val) { Base::setPhiRegW(val); }
    void setPhiRegP(Real val) { Base::setPhiRegP(val); }
    Real phiRegW() const { return Base::phiRegW(); }

    void setPatternParamBounds(const MX2d &val) { 
        if (size_t(val.rows()) != num_pattern_vars()) throw std::runtime_error("Invalid pattern parameter count");
        for (size_t i = 0; i < num_pattern_vars(); ++i) {
            if (val(i, 0) > val(i, 1)) throw std::runtime_error("Invalid pattern parameter bounds. Lower bound must be <= upper bound.");
        }
        m_pp_bounds = val; 
    }
    MX2d patternParamBounds() const { return m_pp_bounds; }

    void setPatternParamNormalizationFactors(const VXd &val) { 
        if (size_t(val.rows()) != num_pattern_vars()) throw std::runtime_error("Invalid pattern parameter count");
        m_param_normalization_factors = val; 
    }
    VXd patternParamNormalizationFactors() const { return m_param_normalization_factors; }




    const VXd &getPatternParams() const { return m_pattern_params; }

    VXd getPatternParams(size_t pi) const {
        if (pi >= num_pattern_vars()) throw std::runtime_error("Invalid pattern parameter index");
        return m_pattern_params.segment(pi * mesh().numTris(), mesh().numTris());
    }

    // Get mat info args
    const std::vector<VXd> &getMatInfoArgs() const { 
        if (m_matInfoArgs.size() != num_pattern_vars()) throw std::runtime_error("Invalid mat info args");
        return m_matInfoArgs; 
    }

    // We have two stretching factors and five coefficients in the bending stiffness polynomial. We need objective, gradient, hessian for each of them.
    // The order of materialInfo is the objective, gradient, hessian for the first material, then the second material, etc.
    MatInfoFunction materialInfo;

    template<class F>
    void visitPatternParams(const F &f) const {
        size_t nt = mesh().numTris();
        for (size_t pi = 0; pi < num_pattern_vars(); ++pi) {
            for (size_t ti = 0; ti < mesh().numTris(); ++ti) {
                // Index of the pattern var, lower bound, upper bound.
                f(pi * nt + ti, m_pp_bounds(pi, 0), m_pp_bounds(pi, 1));
            }
        }
    }

    // Pattern variable bound constraint barrier.
    Real energyBoundConstraint() const {
        Real result = 0.0;
        if (not useConstraintBarrier) return result;
        visitPatternParams([&](size_t i, Real lower, Real upper) {
                result += m_constraintBarrier.eval(m_pattern_params(i), lower, upper);
            });
        return result;
    }

    void addBoundConstraintGradient(VXd &g) const {
        if (not useConstraintBarrier) return;
        visitPatternParams([&](size_t i, Real lower, Real upper) {
                g[i + patternOffset()] += m_constraintBarrier.deval(m_pattern_params(i), lower, upper);
            });
    }

    void addBoundConstraintHessian(SuiteSparseMatrix &H) const {
        if (not useConstraintBarrier) return;
        visitPatternParams([&](size_t i, Real lower, Real upper) {
                size_t var = i + patternOffset();
                H.addNZ(var, var, m_constraintBarrier.d2eval(m_pattern_params(i), lower, upper));
            });
    }

    virtual ~RegularizedPatternParametrizer() { }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    std::vector<VXd> perVertexPatternParams() const;

    std::tuple<std::shared_ptr<Mesh>, VXd, std::vector<VXd>>
    upsampledVertexLeftStretchAnglesAndPatternParameters(size_t nsubdiv, double agreementThreshold) const;

    void setConstraintBarrier(bool use) { useConstraintBarrier = use; }
    bool getConstraintBarrier() const { return useConstraintBarrier; }
    
protected:
    struct BendingSensitivityCache {
        // Cache of x_hat' Jacobians and Hessians
        // (to accelerate repeated calls to elastic energy Hessian/gradient).
        std::vector<BendingStiffnessIntegralSensitivityPattern> sensitivityForBending;
        void update(const RegularizedPatternParametrizer &rgpwp) {
            BENCHMARK_SCOPED_TIMER_SECTION timer("RegularizedPatternParametrizer.BendingSensitivityCache.update");

            if (evaluated) return;
            if (rgpwp.materialInfo.size() != rgpwp.num_total_info) throw std::runtime_error("Invalid material info size");
            
            const size_t nt = rgpwp.Base::mesh().numTris();
            sensitivityForBending.resize(nt);

            std::vector<VXd> evaluated_bsiInfo;
            evaluated_bsiInfo.resize(rgpwp.num_info_per_mat_var * rgpwp.num_stiffness_var);
            // Extract information about the bending stiffness coefficients.

            // The first two variables are stretching factors, so they are not included in the bending information.
            for (size_t i = 0; i < rgpwp.num_stiffness_var * rgpwp.num_info_per_mat_var; ++i) {
                evaluated_bsiInfo[i] = rgpwp.materialInfo[i + rgpwp.num_info_per_mat_var * rgpwp.num_non_stiffness_var](rgpwp.getMatInfoArgs());
            }

            auto processTri = [this, &rgpwp, &evaluated_bsiInfo, nt](size_t ti) {
                // Add print statement for everyline of the following code so we can see how far we are.
                std::vector<VXd> curr_bsi_info;
                // stiffness_coefficients, coefficient_gradient_alpha, coefficient_gradient_beta, coefficient_hessian_alpha, coefficient_hessian_beta, coefficient_hessian_alpha_beta
                size_t num_derivative_info = 1 + rgpwp.num_pattern_vars() + rgpwp.num_pattern_vars() * rgpwp.num_pattern_vars();
                curr_bsi_info.resize(num_derivative_info);
                size_t num_bending_var = rgpwp.num_stiffness_var;
                // Objective
                curr_bsi_info[0] = VXd::Zero(num_bending_var);
                for (size_t i = 0; i < num_bending_var; ++i) {
                    curr_bsi_info[0][i] = evaluated_bsiInfo[i * rgpwp.num_info_per_mat_var](ti);
                }
                // Gradient
                for (size_t i = 0; i < rgpwp.num_pattern_vars(); ++i) {
                    curr_bsi_info[i + 1] = VXd::Zero(num_bending_var);
                    for (size_t j = 0; j < num_bending_var; ++j) {
                        curr_bsi_info[i + 1][j] = evaluated_bsiInfo[j * rgpwp.num_info_per_mat_var + 1](i * nt + ti);
                    }
                }
                // Hessian
                for (size_t i = 0; i < rgpwp.num_pattern_vars(); ++i) {
                    for (size_t j = 0; j < rgpwp.num_pattern_vars(); ++j) {
                        size_t linear_index = i * rgpwp.num_pattern_vars() + j;
                        curr_bsi_info[linear_index + rgpwp.num_pattern_vars() + 1] = VXd::Zero(num_bending_var);
                        for (size_t k = 0; k < num_bending_var; ++k) {
                            curr_bsi_info[linear_index + rgpwp.num_pattern_vars() + 1][k] = evaluated_bsiInfo[k * rgpwp.num_info_per_mat_var + 2](linear_index * nt + ti);
                        }
                    }
                }

                Real delta = rgpwp.Base::kappaAngle(ti, 0) - rgpwp.Base::getPsis()(ti);

                sensitivityForBending[ti].num_psi_p = rgpwp.num_pattern_vars() + 1;
                sensitivityForBending[ti].update(delta, rgpwp.Base::kappa(ti, 0), rgpwp.Base::kappa(ti, 1), curr_bsi_info);
            }; 

            parallel_for_range(nt, processTri);
            // for (size_t ti = 0; ti < nt; ++ti) {
            //     processTri(ti);
            // }

            evaluated = true;
        }

        bool filled() const { return !sensitivityForBending.empty(); }
        const BendingStiffnessIntegralSensitivityPattern &lookup(size_t vi) const { return sensitivityForBending.at(vi); }
        bool evaluated = false;
        void clear() { evaluated = false;  sensitivityForBending.clear();}
    };

    mutable BendingSensitivityCache m_bendingSensitivityCache;
    ConstraintBarrier m_constraintBarrier;
    bool useConstraintBarrier = true;

private:
    size_t m_num_pattern_vars;
    VXd m_pattern_params; // Current per-triangle variables controlling the fusing pattern parameters
    std::vector<VXd> m_matInfoArgs;
    void m_set_matInfoArgs() {
        size_t nt = mesh().numTris();
        m_matInfoArgs.clear();
        m_matInfoArgs.resize(m_num_pattern_vars);
        for (size_t pi = 0; pi < m_num_pattern_vars; ++pi) {
            m_matInfoArgs[pi] = m_pattern_params.segment(pi * nt, nt);
        }
    }

    Real m_pattern_reg_w = 1.0;
    Real m_pattern_reg_p = 2.0; // must be >= 1.0!

    Real m_bending_w = 1.0;

    MX2d m_pp_bounds; // Upper and lower bounds on the pattern parameters
    VXd m_param_normalization_factors;
    mutable std::unique_ptr<SuiteSparseMatrix> m_cachedHessianSparsity, m_cachedBaseHessianSparsity;

};

}

#endif