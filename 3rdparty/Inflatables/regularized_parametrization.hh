#ifndef REGULARIZED_PARAMETRIZATION_HH
#define REGULARIZED_PARAMETRIZATION_HH

#include "parametrization.hh"

#include "local_global_parametrization.hh"

#include "DualLaplacianStencil.hh"

#include <MeshFEM/FEMMesh.hh>
#include <MeshFEM/SparseMatrices.hh>
#include <memory>
#include <utility>

#include "SVDSensitivity.hh"

#include "BendingStiffnessIntegralSensitivity.hh"

namespace parametrization {

////////////////////////////////////////////////////////////////////////////////
// RegularizedParametrizer: Global nonlinear energy with auxiliary variables
////////////////////////////////////////////////////////////////////////////////
// Perform the nonlinear minimization:
//      min 0.5 * ||grad f - M Bt||^2 + w_phi ||grad phi||_{p_phi} + w_alpha ||grad alpha||_{p_alpha}
// with a Newton-type method, where alpha and phi are the target metric
// stretching factors/orientation.
// This relies on a good initialization (e.g., computed from a local-global method).
// The target metric M is of the form:
//      U(phi) diag(alpha, 1) V(psi)^T
// Where, e.g., U(phi) := [cos(phi) -sin(phi); sin(phi) cos(phi)]
struct RegularizedParametrizer : public Parametrizer {
    enum class EnergyType { Full, Fitting, PhiRegularization, AlphaRegularization };

    // Initialize from the local-global parametrizer
    RegularizedParametrizer(LocalGlobalParametrizer &lgparam);

    Eigen::VectorXd getVars() const {
        Eigen::VectorXd result(numVars());
        Eigen::Map<UVMap>(result.data(), m_uv.rows(), m_uv.cols()) = m_uv;
        result.segment(phiOffset(), m_phi.rows()) = m_phi;
        result.segment(psiOffset(), m_psi.rows()) = m_psi;
        if (variableAlpha()) result.segment(alphaOffset(), m_alpha.rows()) = m_alpha;
        return result;
    }

    void setVars(const Eigen::VectorXd &vars) {
        if (size_t(vars.rows()) != numVars()) throw std::runtime_error("Invalid variable count");
        m_uv = Eigen::Map<const UVMap>(vars.data(), m_uv.rows(), m_uv.cols());
        m_phi = vars.segment(phiOffset(), m_phi.rows());
        m_psi = vars.segment(psiOffset(), m_psi.rows());
        if (variableAlpha()) m_alpha = vars.segment(alphaOffset(), m_alpha.rows());
        setUV(m_uv);
    }

    // Offsets of field variables within the full vector of parametrization variables.
    size_t    uvOffset() const { return 0; }
    size_t     uOffset() const { return uvOffset(); }
    size_t     vOffset() const { return uvOffset() + m_uv.rows(); }
    size_t   phiOffset() const { return numUVVars(); }
    size_t   psiOffset() const { return phiOffset() + m_phi.size(); }
    size_t alphaOffset() const { return psiOffset() + m_psi.size(); }

    Real alphaMin() const { return m_alphaMin; }
    Real alphaMax() const { return m_alphaMax; }

    // Adjust bounds on alpha, clamping existing variables.
    void setAlphaMin(Real alphaMin) { m_alphaMin = alphaMin; m_alpha = m_alpha.array().max(alphaMin); parametrizationUpdated(); }
    void setAlphaMax(Real alphaMax) { m_alphaMax = alphaMax; m_alpha = m_alpha.array().min(alphaMax); parametrizationUpdated(); }

    void setVariableAlpha(bool varAlpha) { m_variableAlpha = varAlpha; }
    bool    variableAlpha() const { return m_variableAlpha; }
    

    const Eigen::VectorXd getAlphas() const override { return m_alpha; }
    const Eigen::VectorXd getPhis()   const { return m_phi; }

    // Angle between the local frame's x axis and left/right singular vectors
    virtual Real  leftStretchAngle(size_t i) const override { return m_phi[i]; }
    virtual Real rightStretchAngle(size_t i) const override { return m_psi[i]; }

    virtual Real energy() const override { return m_energy; }

    Eigen::VectorXd gradient(EnergyType etype = EnergyType::Full) const;

    size_t hessianNNZ() const { return hessianSparsityPattern().nz; } // TODO: predict without constructing

    SuiteSparseMatrix hessianSparsityPattern(Real val = 0.0) const;

    void              hessian(SuiteSparseMatrix &H, EnergyType etype = EnergyType::Full) const; // accumulate Hessian to H
    SuiteSparseMatrix hessian(                      EnergyType etype = EnergyType::Full) const; // construct and return Hessian
    void              hessian(SuiteSparseMatrix &H, bool /* projectionMask */) const { hessian(H, EnergyType::Full); }

    size_t numVars()      const { return 2 * m_uv.rows() + m_phi.rows() + m_psi.rows() + numAlphaVars(); }
    size_t numAlphaVars() const { return m_variableAlpha ? m_alpha.rows() : 0; }

    // Cache energy
    virtual void parametrizationUpdated() override { m_evalIterate(); }

    // Regularization parameters
    void setAlphaRegW(Real val) { m_alpha_reg_w = val; parametrizationUpdated(); }
    void setAlphaRegP(Real val) { m_alpha_reg_p = val; parametrizationUpdated(); }
    void   setPhiRegW(Real val) {   m_phi_reg_w = val; parametrizationUpdated(); }
    void   setPhiRegP(Real val) {   m_phi_reg_p = val; parametrizationUpdated(); }

    Real alphaRegW() const { return m_alpha_reg_w; }
    Real alphaRegP() const { return m_alpha_reg_p; }
    Real   phiRegW() const { return   m_phi_reg_w; }
    Real   phiRegP() const { return   m_phi_reg_p; }

    virtual ~RegularizedParametrizer() { }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
    Real m_alphaMin = 1.0, m_alphaMax = 1.0;
    Eigen::VectorXd m_phi, m_psi, m_alpha; // Current per-triangle variables controlling the target metric stretch and orientation.
    aligned_std_vector<M2d> m_M, m_U, m_V, m_dU_dphi, m_dV_dpsi;
    SuiteSparseMatrix m_laplacian;

    Real m_energy = 0.0;

    Real m_alpha_reg_w = 1.0;
    Real m_alpha_reg_p = 2.0; // must be >= 1.0!
    Real   m_phi_reg_w = 1.0;
    Real   m_phi_reg_p = 2.0; // must be >= 1.0!

    bool m_variableAlpha = true;

    void m_evalIterate();
};


////////////////////////////////////////////////////////////////////////////////
// RegularizedParametrizerSVD: Global nonlinear energy with uv variables only
////////////////////////////////////////////////////////////////////////////////
// Perform the nonlinear minimization to find a mapping f:
// min    0.5 * ||sigma_1 - 1.0||^2
//  f   + 0.5 * ||(sigma_0 - alpha_min)_{-}||^2
//      + 0.5 * ||(sigma_0 - alpha_max)_{+}||^2
//      + w_phi ||grad phi(u_0)||_{p_phi}
//      + w_alpha ||grad sigma_0||_{p_alpha}
// where sigma_0 and sigma_1 are the largest and second largest singular values of f,
// u_0 is the left singular vector corresponding to sigma_0, and phi(u_0) is the
// angle this singular vector makes with the "x" axis.
// This relies on a good initialization (e.g., computed from a local-global method).
// The target metric M is of the form:
//      U(phi) diag(alpha, 1) V(psi)^T
// Where, e.g., U(phi) := [cos(phi) -sin(phi); sin(phi) cos(phi)]
struct RegularizedParametrizerSVD : public Parametrizer {
    enum class EnergyType { Full, Fitting, PhiRegularization, AlphaRegularization, BendingRegularization };

    RegularizedParametrizerSVD(LocalGlobalParametrizer &lgparam)
        : RegularizedParametrizerSVD(lgparam.meshPtr(), lgparam.uv(), lgparam.alphaMin(), lgparam.alphaMax(), true) { }
    RegularizedParametrizerSVD(const std::shared_ptr<Mesh> &inMesh, const UVMap &uvInit, Real amin = 1.0, Real amax = 1.0, bool transformForRigidMotionConstraint = false);

    Eigen::VectorXd getVars() const {
        Eigen::VectorXd result(numVars());
        Eigen::Map<UVMap>(result.data(), m_uv.rows(), m_uv.cols()) = m_uv;
        return result;
    }

    void setVars(const Eigen::VectorXd &vars) { setUV(Eigen::Map<const UVMap>(vars.data(), m_uv.rows(), m_uv.cols())); }

    // Offsets of field variables within the full vector of parametrization variables.
    size_t    uvOffset() const { return 0; }
    size_t     uOffset() const { return uvOffset(); }
    size_t     vOffset() const { return uvOffset() + m_uv.rows(); }

    const std::array<size_t, 3> &rigidMotionPinVars() const { return m_rigidMotionPinVars; }

    Real alphaMin() const { return m_alphaMin; }
    Real alphaMax() const { return m_alphaMax; }

    // Adjust bounds on alpha
    void setAlphaMin(Real alphaMin) { m_alphaMin = alphaMin; } // energy must be re-evaluated when alphamin/max are changed...
    void setAlphaMax(Real alphaMax) { m_alphaMax = alphaMax; }

    Real getAlpha(size_t i) const { return m_alpha[i]; }
    Real   getPhi(size_t i) const { return m_phi[i]; }

    const Eigen::VectorXd getAlphas() const override { return m_alpha; }
    const Eigen::VectorXd getPhis()   const { return m_phi; }

    // For debugging/analyzing alpha regularization only...
    void setAlphas(const Eigen::VectorXd &val) { m_alpha = val; }

    Eigen::VectorXd getMinSingularValues() const {
        const size_t nt = mesh().numTris();
        Eigen::VectorXd result(nt);
        for (size_t ti = 0; ti < nt; ++ti) result[ti] = m_svds[ti].sigma(1);
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

    // Angle between the local frame's x axis and left/right singular vectors
    virtual Real  leftStretchAngle(size_t i) const override { return m_phi[i]; }
    virtual Real rightStretchAngle(size_t i) const override { const auto &v0 = m_svds.at(i).v(0); return std::atan2(v0[1], v0[0]); }

    // Get the 3D orientations of the implied air channel tube for each triangle
    // (i.e. the right stretch vector corresponding sigma_1 ~= 1.0).
    Eigen::MatrixX3d tubeDirections() const {
        const size_t nt = mesh().numTris();
        Eigen::MatrixX3d result(nt, 3);
        for (size_t i = 0; i < nt; ++i)
            result.row(i) = m_B[i] * m_svds.at(i).v(1);
        return result;
    }

    virtual Real energy() const override { return energy(EnergyType::Full); }
    Real energy(EnergyType etype) const;

    Eigen::VectorXd gradient(EnergyType etype = EnergyType::Full) const;

    size_t hessianNNZ() const { return hessianSparsityPattern().nz; } // TODO: predict without constructing

    SuiteSparseMatrix hessianSparsityPattern(Real val = 0.0) const;

    void              hessian(SuiteSparseMatrix &H, EnergyType etype = EnergyType::Full, bool projectionMask = false) const; // accumulate Hessian to H
    SuiteSparseMatrix hessian(                      EnergyType etype = EnergyType::Full, bool projectionMask = false) const; // construct and return Hessian
    void              hessian(SuiteSparseMatrix &H, bool projectionMask) const { hessian(H, EnergyType::Full, projectionMask); }


    size_t numVars()   const { return 2 * m_uv.rows(); }

    // Cache energy
    virtual void parametrizationUpdated() override { m_evalIterate(); }

    // Fitting energy parameters
    void setStretchDeviationP(Real val) { m_stretch_deviation_p = val; }
    Real stretchDeviationP() const { return m_stretch_deviation_p; }

    // Regularization parameters
    void setAlphaRegW(Real val) { m_alpha_reg_w = val; }
    void setAlphaRegP(Real val) { m_alpha_reg_p = val; }
    void   setPhiRegW(Real val) {   m_phi_reg_w = val; }
    void   setPhiRegP(Real val) {   m_phi_reg_p = val; }
    void  setBendRegW(Real val) {  m_bend_reg_w = val; }

    Real alphaRegW() const { return m_alpha_reg_w; }
    Real alphaRegP() const { return m_alpha_reg_p; }
    Real   phiRegW() const { return   m_phi_reg_w; }
    Real   phiRegP() const { return   m_phi_reg_p; }
    Real  bendRegW() const { return  m_bend_reg_w; }

    virtual ~RegularizedParametrizerSVD() { }

    DualLaplacianStencil<Mesh> dualLaplacianStencil;
    bool scaleInvariantFittingEnergy = false; // `true` is more desirable, but we keep `false` by default for backwards compatibility

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
    Real m_alphaMin = 1.0, m_alphaMax = 1.0;

    Eigen::VectorXd m_alpha, m_phi; // stretch factor and orientation from the current parametrization

    // Averaged shape operators on each triangle
    aligned_std_vector<M2d> m_shapeOperators;

    aligned_std_vector<SVDSensitivity> m_svds;

    std::array<size_t, 3> m_rigidMotionPinVars;

    Real m_alpha_reg_w = 1.0;
    Real m_alpha_reg_p = 2.0; // must be >= 1.0!
    Real   m_phi_reg_w = 1.0;
    Real   m_phi_reg_p = 2.0; // must be >= 1.0!
    Real  m_bend_reg_w = 0.0;

    Real m_stretch_deviation_p = 2.0; // Note: set to > 2.0 for C^2 fitting term (at 1.0 it's only C^1).

    void m_evalIterate();
};

////////////////////////////////////////////////////////////////////////////////
// RegularizedGenericParametrizerSVD: Global nonlinear energy with uv variables only
////////////////////////////////////////////////////////////////////////////////
// Perform the nonlinear minimization to find a mapping f:
// min    0.5 * ||sigma_1 - beta_min)_{-}||^2
//  f   + 0.5 * ||sigma_1 - beta_max)_{+}||^2
//      + 0.5 * ||(sigma_0 - alpha_min)_{-}||^2
//      + 0.5 * ||(sigma_0 - alpha_max)_{+}||^2
//      + w_phi ||grad phi(u_0)||_{p_phi}
//      + w_alpha ||grad sigma_0||_{p_alpha}
//      + w_beta ||grad sigma_1||_{p_beta}
// where sigma_0 and sigma_1 are the largest and second largest singular values of f,
// u_0 is the left singular vector corresponding to sigma_0, and phi(u_0) is the
// angle this singular vector makes with the "x" axis.
// This relies on a good initialization (e.g., computed from a local-global method).
// The target metric M is of the form:
//      U(phi) diag(alpha, beta) V(psi)^T
// Where, e.g., U(phi) := [cos(phi) -sin(phi); sin(phi) cos(phi)]
struct RegularizedGenericParametrizerSVD : public Parametrizer {
    enum class EnergyType { Full, Fitting, PhiRegularization, AlphaRegularization, BetaRegularization, BendingRegularization };

    RegularizedGenericParametrizerSVD(LocalGlobalGenericParametrizer &lgparam)
        : RegularizedGenericParametrizerSVD(lgparam.meshPtr(), lgparam.uv(), lgparam.alphaMin(), lgparam.alphaMax(), lgparam.betaMin(), lgparam.betaMax(), true) { }
    RegularizedGenericParametrizerSVD(const std::shared_ptr<Mesh> &inMesh, const UVMap &uvInit, Real amin = 1.0, Real amax = 1.0, Real bmin = 1.0, Real bmax = 1.0,  bool transformForRigidMotionConstraint = false);

    Eigen::VectorXd getVars() const {
        Eigen::VectorXd result(numVars());
        Eigen::Map<UVMap>(result.data(), m_uv.rows(), m_uv.cols()) = m_uv;
        return result;
    }

    void setVars(const Eigen::VectorXd &vars) { setUV(Eigen::Map<const UVMap>(vars.data(), m_uv.rows(), m_uv.cols())); }

    // Offsets of field variables within the full vector of parametrization variables.
    size_t    uvOffset() const { return 0; }
    size_t     uOffset() const { return uvOffset(); }
    size_t     vOffset() const { return uvOffset() + m_uv.rows(); }

    const std::array<size_t, 3> &rigidMotionPinVars() const { return m_rigidMotionPinVars; }

    Real alphaMin() const { return m_alphaMin; }
    Real alphaMax() const { return m_alphaMax; }

    Real betaMin() const { return m_betaMin; }
    Real betaMax() const { return m_betaMax; }

    // Adjust bounds on alpha
    void setAlphaMin(Real alphaMin) { m_alphaMin = alphaMin; } // energy must be re-evaluated when alphamin/max are changed...
    void setAlphaMax(Real alphaMax) { m_alphaMax = alphaMax; }

    // Adjust bounds on beta
    void setBetaMin(Real betaMin) { m_betaMin = betaMin; } // energy must be re-evaluated when betamin/max are changed...
    void setBetaMax(Real betaMax) { m_betaMax = betaMax; }

    Real getAlpha(size_t i) const { return m_alpha[i]; }
    Real getBeta(size_t i) const { return m_beta[i]; }
    Real   getPhi(size_t i) const { return m_phi[i]; }

    const Eigen::VectorXd getAlphas() const override { return m_alpha; }
    const Eigen::VectorXd getBetas() const override { return m_beta; }
    const Eigen::VectorXd getPhis()   const { return m_phi; }

    // For debugging/analyzing alpha/beta regularization only...
    void setAlphas(const Eigen::VectorXd &val) { m_alpha = val; }
    void setBetas(const Eigen::VectorXd &val) { m_beta = val; }

    Eigen::VectorXd getMinSingularValues() const {
        const size_t nt = mesh().numTris();
        Eigen::VectorXd result(nt);
        for (size_t ti = 0; ti < nt; ++ti) result[ti] = m_svds[ti].sigma(1);
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

    // Angle between the local frame's x axis and left/right singular vectors
    virtual Real  leftStretchAngle(size_t i) const override { return m_phi[i]; }
    virtual Real rightStretchAngle(size_t i) const override { const auto &v0 = m_svds.at(i).v(0); return std::atan2(v0[1], v0[0]); }

    // Get the 3D orientations of the implied air channel tube for each triangle
    // (i.e. the right stretch vector corresponding sigma_1 = beta).
    Eigen::MatrixX3d tubeDirections() const {
        const size_t nt = mesh().numTris();
        Eigen::MatrixX3d result(nt, 3);
        for (size_t i = 0; i < nt; ++i)
            result.row(i) = m_B[i] * m_svds.at(i).v(1);
        return result;
    }

    virtual Real energy() const override { return energy(EnergyType::Full); }
    Real energy(EnergyType etype) const;

    Eigen::VectorXd gradient(EnergyType etype = EnergyType::Full) const;

    size_t hessianNNZ() const { return hessianSparsityPattern().nz; } // TODO: predict without constructing

    SuiteSparseMatrix hessianSparsityPattern(Real val = 0.0) const;

    void              hessian(SuiteSparseMatrix &H, EnergyType etype = EnergyType::Full, bool projectionMask = false) const; // accumulate Hessian to H
    SuiteSparseMatrix hessian(                      EnergyType etype = EnergyType::Full, bool projectionMask = false) const; // construct and return Hessian
    void              hessian(SuiteSparseMatrix &H, bool projectionMask) const { hessian(H, EnergyType::Full, projectionMask); }


    size_t numVars()   const { return 2 * m_uv.rows(); }

    // Cache energy
    virtual void parametrizationUpdated() override { m_evalIterate(); }

    // Fitting energy parameters
    void setStretchDeviationP(Real val) { m_stretch_deviation_p = val; }
    Real stretchDeviationP() const { return m_stretch_deviation_p; }

    // Regularization parameters
    void setStretchRegW(Real val) { m_stretch_reg_w = val; }
    void setStretchRegP(Real val) { m_stretch_reg_p = val; }
    void   setPhiRegW(Real val) {   m_phi_reg_w = val; }
    void   setPhiRegP(Real val) {   m_phi_reg_p = val; }
    void  setBendRegW(Real val) { m_bend_reg_w = val; }

    Real stretchRegW() const { return m_stretch_reg_w; }
    Real stretchRegP() const { return m_stretch_reg_p; }
    Real   phiRegW() const { return   m_phi_reg_w; }
    Real   phiRegP() const { return   m_phi_reg_p; }
    Real  bendRegW() const { return  m_bend_reg_w; }

    virtual ~RegularizedGenericParametrizerSVD() { }

    DualLaplacianStencil<Mesh> dualLaplacianStencil;
    bool scaleInvariantFittingEnergy = false; // `true` is more desirable, but we keep `false` by default for backwards compatibility

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
    Real m_alphaMin = 1.0, m_alphaMax = 1.0;
    Real m_betaMin = 1.0, m_betaMax = 1.0;

    Eigen::VectorXd m_alpha, m_beta, m_phi; // stretch factors and orientation from the current parametrization

    // Averaged shape operators on each triangle
    aligned_std_vector<M2d> m_shapeOperators;

    aligned_std_vector<SVDSensitivity> m_svds;

    std::array<size_t, 3> m_rigidMotionPinVars;

    Real m_stretch_reg_w = 1.0;
    Real m_stretch_reg_p = 2.0; // must be >= 1.0!
    Real   m_phi_reg_w = 1.0;
    Real   m_phi_reg_p = 2.0; // must be >= 1.0!
    Real  m_bend_reg_w = 0.0;

    Real m_stretch_deviation_p = 2.0; // Note: set to > 2.0 for C^2 fitting term (at 1.0 it's only C^1).

    void m_evalIterate();
};

}

#endif