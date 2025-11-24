#ifndef GENERICPARAMETRIZATION_HH
#define GENERICPARAMETRIZATION_HH


#include "local_global_parametrization.hh"
#include "BendingStiffnessIntegralSensitivityPattern.hh"

namespace parametrization {

////////////////////////////////////////////////////////////////////////////////
// RegularizedGenericParametrizer: Global nonlinear energy with auxiliary variables
////////////////////////////////////////////////////////////////////////////////
// Perform the nonlinear minimization:
//      min 0.5 * ||grad f - M Bt||^2 + w_phi ||grad phi||_{p_phi} + w_alpha ||grad alpha||_{p_alpha}
// with a Newton-type method, where alpha, beta and phi are the target metric
// stretching factors/orientation.
// This relies on a good initialization (e.g., computed from a local-global method).
// The target metric M is of the form:
//      U(phi) diag(alpha, beta) V(psi)^T
// Where, e.g., U(phi) := [cos(phi) -sin(phi); sin(phi) cos(phi)]
struct RegularizedGenericParametrizer : public Parametrizer {
    enum class EnergyType { Full, Fitting, PhiRegularization, StretchRegularization, DiffRegularization, BendingRegularization};

    // Initialize from the local-global parametrizer
    RegularizedGenericParametrizer(LocalGlobalGenericParametrizer &lggparam);

    Eigen::VectorXd getVars() const {
        Eigen::VectorXd result(numVars());
        Eigen::Map<UVMap>(result.data(), m_uv.rows(), m_uv.cols()) = m_uv;
        result.segment(phiOffset(), m_phi.rows()) = m_phi;
        result.segment(psiOffset(), m_psi.rows()) = m_psi;
        if (variableStretch()) result.segment(stretchOffset(), m_stretch.rows()) = m_stretch;
        return result;
    }

    void setVars(const Eigen::VectorXd &vars) {
        if (size_t(vars.rows()) != numVars()) throw std::runtime_error("Invalid variable count");
        m_uv = Eigen::Map<const UVMap>(vars.data(), m_uv.rows(), m_uv.cols());
        m_phi = vars.segment(phiOffset(), m_phi.rows());
        m_psi = vars.segment(psiOffset(), m_psi.rows());
        if (variableStretch()) m_stretch = vars.segment(stretchOffset(), m_stretch.rows());
        setUV(m_uv);

    }

    // Offsets of field variables within the full vector of parametrization variables.
    size_t    uvOffset() const { return 0; }
    size_t     uOffset() const { return uvOffset(); }
    size_t     vOffset() const { return uvOffset() + m_uv.rows(); }
    size_t   phiOffset() const { return numUVVars(); }
    size_t   psiOffset() const { return phiOffset() + m_phi.size(); }
    size_t alphaOffset() const { return psiOffset() + m_psi.size(); }
    size_t betaOffset() const { return alphaOffset() + m_stretch.size()/2; }
    size_t stretchOffset() const { return psiOffset() + m_psi.size(); }
    

    Real alphaMin() const { return m_alphaMin; }
    Real alphaMax() const { return m_alphaMax; }
    Real betaMin() const { return m_betaMin; }
    Real betaMax() const { return m_betaMax; }

    void setAlphaMin(Real alphaMin) { m_alphaMin = alphaMin; m_stretch.segment(0, m_stretch.rows()/2) = m_stretch.segment(0, m_stretch.rows()/2).array().max(alphaMin); parametrizationUpdated(); }
    void setAlphaMax(Real alphaMax) { m_alphaMax = alphaMax; m_stretch.segment(0, m_stretch.rows()/2) = m_stretch.segment(0, m_stretch.rows()/2).array().min(alphaMax); parametrizationUpdated(); }

    void setBetaMin(Real betaMin) { m_betaMin = betaMin; m_stretch.segment(m_stretch.rows()/2, m_stretch.rows()/2) = m_stretch.segment(m_stretch.rows()/2, m_stretch.rows()/2).array().max(betaMin); parametrizationUpdated(); }
    void setBetaMax(Real betaMax) { m_betaMax = betaMax; m_stretch.segment(m_stretch.rows()/2, m_stretch.rows()/2) = m_stretch.segment(m_stretch.rows()/2, m_stretch.rows()/2).array().min(betaMax); parametrizationUpdated(); }
    
    void setVariableStretch(bool varStretch) { m_variableStretch = varStretch; }
    bool    variableStretch() const { return m_variableStretch; }

    void setUseBarrier(bool useBarrier) { m_useBarrier = useBarrier; }
    bool    useBarrier() const { return m_useBarrier; }

    const Eigen::VectorXd getAlphas() const override { return m_stretch.segment(0, m_stretch.rows()/2); }
    const Eigen::VectorXd getBetas() const override { return m_stretch.segment(m_stretch.rows()/2, m_stretch.rows()/2); }
    const Eigen::VectorXd getBarriers() const { 
        Eigen::VectorXd result(m_stretch.rows()/2);
        
        const size_t nt = m_stretch.rows()/2;
        for (size_t ti = 0; ti < nt; ++ti) {
            const Real e_minus_bx_i = std::exp(-m_b_b * (m_stretch[ti] - m_stretch[nt + ti]));
            const Real e_bx_i = std::exp(m_b_b * (m_stretch[ti] - m_stretch[nt + ti]));
            
            result[ti] = m_b_a / (m_b_a + e_minus_bx_i) + m_b_a / (m_b_a + e_bx_i);
        }
        return result; 
    }
    const Eigen::VectorXd &getStretch() const { return m_stretch; }
    const Eigen::VectorXd &getPhis()   const { return m_phi; }
    const Eigen::VectorXd &getPsis()   const { return m_psi; }

    const M2d &jacobian_M(size_t i) const { return m_M[i]; }

    // Angle between the local frame's x axis and left/right singular vectors
    virtual Real  leftStretchAngle(size_t i) const override { return m_phi[i]; }
    virtual Real rightStretchAngle(size_t i) const override { return m_psi[i]; }

    virtual Real energy() const override { return energy(EnergyType::Full); }
    Real energy(EnergyType etype) const;
    

    Eigen::VectorXd gradient(EnergyType etype = EnergyType::Full) const;

    size_t hessianNNZ() const { return hessianSparsityPattern().nz; } // TODO: predict without constructing

    SuiteSparseMatrix hessianSparsityPattern(Real val = 0.0) const;

    void              hessian(SuiteSparseMatrix &H, EnergyType etype = EnergyType::Full) const; // accumulate Hessian to H
    SuiteSparseMatrix hessian(                      EnergyType etype = EnergyType::Full) const; // construct and return Hessian
    void              hessian(SuiteSparseMatrix &H, bool /* projectionMask */) const { hessian(H, EnergyType::Full); }

    size_t numVars()      const { return 2 * m_uv.rows() + m_phi.rows() + m_psi.rows() + numStretchVars(); }
    size_t numStretchVars() const { return m_variableStretch ? m_stretch.rows() : 0; }

    // Cache energy
    virtual void parametrizationUpdated() override { m_evalIterate(); }

    // Regularization parameters
    void setStretchRegW(Real val) { m_stretch_reg_w = val; parametrizationUpdated(); }
    void setStretchRegP(Real val) { m_stretch_reg_p = val; parametrizationUpdated(); }
    void   setPhiRegW(Real val) {   m_phi_reg_w = val; parametrizationUpdated(); }
    void   setPhiRegP(Real val) {   m_phi_reg_p = val; parametrizationUpdated(); }
    void   setDiffRegW(Real val) {   m_diff_reg_w = val; parametrizationUpdated(); }
    void setBarrierA(Real val) { m_b_a = val; m_barrier = SigmoidBarrier(m_b_a, m_b_b); parametrizationUpdated(); }
    void setBarrierB(Real val) { m_b_b = val; m_barrier = SigmoidBarrier(m_b_a, m_b_b); parametrizationUpdated(); }

    Real stretchRegW() const { return m_stretch_reg_w; }
    Real stretchRegP() const { return m_stretch_reg_p; }
    Real   phiRegW() const { return   m_phi_reg_w; }
    Real   phiRegP() const { return   m_phi_reg_p; }
    Real   diffRegW() const { return   m_diff_reg_w; }
    Real barrierA() const { return m_b_a; }
    Real barrierB() const { return m_b_b; }

    void setLines(const Eigen::MatrixXd &lines) { 
        m_lines = lines;
        // Normalize lines
        auto lines_normalized = m_lines;
        for (int i = 0; i < lines_normalized.rows(); ++i) {
            lines_normalized.row(i) /= lines_normalized.row(i).segment(0, 2).norm();
        }
        m_lines = lines_normalized;
    }

    bool isPointInLines(const V2d &pt) {
        // Extend pt to V3d
        V3d pt3d(pt[0], pt[1], 1.0);
        auto result = m_lines * pt3d;
        return result.maxCoeff() <= 0.0;
    }
    void projectPointInLines(V2d &pt) {
        // Normalize lines
        auto lines_normalized = m_lines;
        // for (int i = 0; i < lines_normalized.rows(); ++i) {
        //     lines_normalized.row(i) /= lines_normalized.row(i).segment(0, 2).norm();
        // }
        Eigen::VectorXd distances = lines_normalized * V3d(pt[0], pt[1], 1.0);
        // clamp distances to positive
        distances.cwiseMax(0.0);
        pt += ((-distances).asDiagonal() * lines_normalized).colwise().sum();
    }

    size_t numLinearInequalityConstraints() const { return m_lines.rows() * numStretchVars() / 2; }

    std::vector<LinearInequality> getLinearInequalityConstraints() const {
        const size_t nc = numLinearInequalityConstraints();
        const int nt = mesh().numTris();
        std::vector<LinearInequality> result(nc);

        for (int li = 0; li < m_lines.rows(); ++li) {
            for (int vi = 0; vi < nt; ++vi) {
                result[li * nt + vi].vars = std::vector<int> {(int)alphaOffset() +  vi, (int)betaOffset() + vi};
                result[li * nt + vi].coeffs = std::vector<Real> {m_lines(li, 0), m_lines(li, 1)};
                result[li * nt + vi].constPart = m_lines(li, 2);
            }     
        }
        return result;
    }

    // get principal curvature "kappa_i" and 2D direction "d_i" on triangle "tri"
    // (2D vector of components in tri's tangent space basis).
    std::pair<Real, V2d> curvature(size_t tri, size_t i) const {
        if (i > 2) throw std::runtime_error("Kappa subscript i out of bounds");
        Eigen::SelfAdjointEigenSolver<M2d> solver;
        solver.compute(m_shapeOperators.at(tri));

        // Note: Eigen does not guarantee a sorting order; we want k_0 to be
        // the largest (algebraic) eigenvalue.
        if (solver.eigenvalues()[1] > solver.eigenvalues()[0])
            i = 1 - i;

        return { solver.eigenvalues()[i], solver.eigenvectors().col(i) };
    }

    // get principal curvature "kappa_i" and 3D direction "d_i" on triangle "tri"
    // (3D vector lying in the same plane as "tri")
    std::pair<Real, V3d> curvature3d(size_t tri, size_t i) const {
        auto c2d = curvature(tri, i);
        return { c2d.first, m_B.at(tri) * c2d.second };
    }


    // Helper functions for querying the principal curvatures.
    Real kappa(size_t tri, size_t i) const { return curvature(tri, i).first; }
    Real kappaAngle(size_t tri, size_t i) const { 
        V2d vec =  curvature(tri, i).second; 
        if (vec.y() == 0 && vec.x() == 0) throw std::runtime_error("kappaAngle: zero vector");
        return std::atan2(vec.y(), vec.x()); 
    }

    const Eigen::VectorXd getKappaAngle() const {
        Eigen::VectorXd result(mesh().numTris());
        for (size_t ti = 0; ti < mesh().numTris(); ++ti) {
            result[ti] = kappaAngle(ti, 0);
        }
        return result;
    }

    const Eigen::VectorXd getKappa(size_t kappaIdx) const {
        Eigen::VectorXd result(mesh().numTris());
        for (size_t ti = 0; ti < mesh().numTris(); ++ti) {
            result[ti] = kappa(ti, kappaIdx);
        }
        return result;
    }

    virtual ~RegularizedGenericParametrizer() { }

    DualLaplacianStencil<Mesh> dualLaplacianStencil;
    bool scaleInvariantFittingEnergy = true;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
    Real m_alphaMin = 1.0, m_alphaMax = 1.0;
    Real m_betaMin = 1.0, m_betaMax = 1.0;
    Eigen::VectorXd m_phi, m_psi, m_stretch; // Current per-triangle variables controlling the target metric stretch and orientation.
    Eigen::MatrixXd m_lines; // Linear inequality constraints on the stretch variables
    aligned_std_vector<M2d> m_M, m_U, m_V, m_dU_dphi, m_dV_dpsi;
    SuiteSparseMatrix m_laplacian;

    // Averaged shape operators on each triangle
    aligned_std_vector<M2d> m_shapeOperators;


    // Real m_b_a = 1.0, m_b_b = 0.1, m_b_c = (1.0 - std::exp(-10.0)); // Barrier function parameters
    Real m_b_a = 0.1, m_b_b = 100; // Barrier function parameters
    SigmoidBarrier m_barrier = SigmoidBarrier(m_b_a, m_b_b);

    // Real m_energy = 0.0;

    Real m_stretch_reg_w = 1.0;
    Real m_stretch_reg_p = 2.0; // must be >= 1.0!
    Real   m_phi_reg_w = 1.0;
    Real   m_phi_reg_p = 2.0; // must be >= 1.0!
    Real m_diff_reg_w = 1.0;

    bool m_variableStretch = true, m_useBarrier = false;

    void m_evalIterate();

    // We have five coefficients in the bending stiffness polynomial and we also need gradient and hessian information for them. So there are 30 functions that we get from the spline fitting result from python, order as coefficients (5), gradient(10: partial alpha, partial beta for each), hessian(15: partial partial alpha, partial partial beta, partial alpha partial beta for each).
    std::vector<std::function<Eigen::VectorXd(Eigen::VectorXd, Eigen::VectorXd)>> bendingStiffnessIntegralInfo;
};

}

#endif