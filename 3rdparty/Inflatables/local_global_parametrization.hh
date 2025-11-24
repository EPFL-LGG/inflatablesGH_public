#ifndef LOCAL_GLOBAL_PARAMETRIZATION_HH
#define LOCAL_GLOBAL_PARAMETRIZATION_HH

#include "parametrization.hh"
#include "DualLaplacianStencil.hh"

#include <MeshFEM/FEMMesh.hh>
#include <MeshFEM/SparseMatrices.hh>
#include <memory>
#include <utility>

#include "SVDSensitivity.hh"

#include "BendingStiffnessIntegralSensitivity.hh"

namespace parametrization {

struct LocalGlobalParametrizer : public Parametrizer {
    LocalGlobalParametrizer(const std::shared_ptr<Mesh> &inMesh, const UVMap &uvInit);

    Real alphaMin() const { return m_alphaMin; }
    Real alphaMax() const { return m_alphaMax; }

    void setAlphaMin(Real alphaMin) { m_alphaMin = alphaMin;           parametrizationUpdated(); }
    void setAlphaMax(Real alphaMax) { m_alphaMax = alphaMax;           parametrizationUpdated(); }
    void setAlpha(Real alpha)       { m_alphaMin = m_alphaMax = alpha; parametrizationUpdated(); } // use constant target stretch value of "alpha"

    void runIteration();

    // Principal directions scaled by the current stretch factors (singular values)
    std::pair<MX3d, MX3d> scaledPrincipalDirections(Domain d) const {
        const size_t nt = mesh().numTris();
        std::pair<MX3d, MX3d> result;
        result.first.resize(nt, 3);
        result.second.resize(nt, 3);

        for (size_t ti = 0; ti < nt; ++ti) {
            M23d frame_transpose;
            if (d == Domain::UV) { frame_transpose << m_lambda.row(ti).asDiagonal().inverse() * (m_U[ti] * m_R[ti]).transpose(), Eigen::Vector2d::Zero(); }
            else                 { frame_transpose =  m_lambda.row(ti).asDiagonal()           * (m_B[ti] * m_R[ti]).transpose(); }

            result.first .row(ti) = frame_transpose.row(0);
            result.second.row(ti) = frame_transpose.row(1);
        }

        return result;
    }

    const aligned_std_vector<M2d> &getR()      const { return m_R; }
    const aligned_std_vector<M2d> &getU()      const { return m_U; }
    const Eigen::VectorXd  getAlphas() const override { return m_alpha; }

    const M2d &getR(size_t i) const { return m_R.at(i); }
    const M2d &getU(size_t i) const { return m_U.at(i); }
    Real   getAlpha(size_t i) const { return m_alpha[i]; }

    // Angle between the local frame's x axis and left/right singular vectors
    virtual Real  leftStretchAngle(size_t i) const override;
    virtual Real rightStretchAngle(size_t i) const override;

    virtual Real energy() const override { return m_energy; }

    virtual void parametrizationUpdated() override { m_localStep(); }

    virtual ~LocalGlobalParametrizer(); // Out-of-line destructor needed since SPSDSystemSolver is incomplete type

    std::unique_ptr<SPSDSystemSolver> L;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
    void m_localStep();

    Real m_energy = 0.0;

    // Range of admissible stretch factors
    Real m_alphaMin = 1.0, m_alphaMax = 1.0;

    // Cached quantities computed from the local step
    aligned_std_vector<M23d> m_M_Bt;
    aligned_std_vector<M2d> m_R, m_U;
    Eigen::VectorXd m_alpha; // Controls anisotropy of the target singular values for the Jacobian ([alpha, 1]). We assume alpha > 1.
    Eigen::Matrix<Real, Eigen::Dynamic, 2> m_lambda;
};


struct LocalGlobalGenericParametrizer : public Parametrizer {
    LocalGlobalGenericParametrizer(const std::shared_ptr<Mesh> &inMesh, const UVMap &uvInit);

    Real alphaMin() const { return m_alphaMin; }
    Real alphaMax() const { return m_alphaMax; }
    Real betaMin() const { return m_betaMin; }
    Real betaMax() const { return m_betaMax; }

    void setStretchMin(Real alphaMin, Real betaMin) { m_alphaMin = alphaMin; m_betaMin = betaMin;           parametrizationUpdated(); }
    void setStretchMax(Real alphaMax, Real betaMax) { m_alphaMax = alphaMax; m_betaMax = betaMax;           parametrizationUpdated(); }

    void setAlphaMin(Real alphaMin) { m_alphaMin = alphaMin; parametrizationUpdated(); }
    void setAlphaMax(Real alphaMax) { m_alphaMax = alphaMax; parametrizationUpdated(); }
    void setBetaMin(Real betaMin) { m_betaMin = betaMin; parametrizationUpdated(); }
    void setBetaMax(Real betaMax) { m_betaMax = betaMax; parametrizationUpdated(); }
    
    void setStretch(Real alpha, Real beta)       { m_alphaMin = m_alphaMax = alpha; m_betaMin = m_betaMax = beta; parametrizationUpdated(); } // use constant target stretch values

    void runIteration();

    // Principal directions scaled by the current stretch factors (singular values)
    std::pair<MX3d, MX3d> scaledPrincipalDirections(Domain d) const {
        const size_t nt = mesh().numTris();
        std::pair<MX3d, MX3d> result;
        result.first.resize(nt, 3);
        result.second.resize(nt, 3);

        for (size_t ti = 0; ti < nt; ++ti) {
            M23d frame_transpose;
            if (d == Domain::UV) { frame_transpose << m_lambda.row(ti).asDiagonal().inverse() * (m_U[ti] * m_R[ti]).transpose(), Eigen::Vector2d::Zero(); }
            else                 { frame_transpose =  m_lambda.row(ti).asDiagonal()           * (m_B[ti] * m_R[ti]).transpose(); }

            result.first .row(ti) = frame_transpose.row(0);
            result.second.row(ti) = frame_transpose.row(1);
        }

        return result;
    }

    const aligned_std_vector<M2d> &getR()      const { return m_R; }
    const aligned_std_vector<M2d> &getU()      const { return m_U; }
    const Eigen::VectorXd getAlphas() const override { return m_stretch.segment(0, m_stretch.size() / 2); }
    const Eigen::VectorXd getBetas() const override { return m_stretch.segment(m_stretch.size() / 2, m_stretch.size() / 2); }
    
    const M2d &getR(size_t i) const { return m_R.at(i); }
    const M2d &getU(size_t i) const { return m_U.at(i); }
    Real   getAlpha(size_t i) const { return getAlphas()[i]; }
    Real   getBeta(size_t i) const { return getBetas()[i]; }

    // Angle between the local frame's x axis and left/right singular vectors
    virtual Real  leftStretchAngle(size_t i) const override;
    virtual Real rightStretchAngle(size_t i) const override;

    virtual Real energy() const override { return m_energy; }

    virtual void parametrizationUpdated() override { m_localStep(); }

    Eigen::MatrixXd getLines() const { return m_lines; }
    
    void setLines(const Eigen::MatrixXd &lines) { 
        m_lines = lines;
        // Normalize lines
        auto lines_normalized = m_lines;
        for (int i = 0; i < lines_normalized.rows(); ++i) {
            lines_normalized.row(i) /= lines_normalized.row(i).segment(0, 2).norm();
        }
        m_lines = lines_normalized;
        computeCorners();
    }

    void computeCorners() {
        // go through pairs of consecutive lines in m_lines
        // find intersection of each pair
        // store intersection points in m_corners
        m_corners.clear();
        for(auto i = 0; i < m_lines.rows(); ++i)
        {
            auto j = (i + 1) % m_lines.rows();
            V3d a = m_lines.row(i).transpose(), b = m_lines.row(j).transpose();
            V3d cross_product = a.cross(b);
            if (cross_product[2] == 0.0)
                throw std::runtime_error("Invalid cross product of two lines in plane");
            m_corners.push_back(V2d(cross_product[0] / cross_product[2], cross_product[1] / cross_product[2]));
        }

    }

    bool isPointInLines(const V2d &pt, const Real eps = 0.0) {
        // Extend pt to V3d
        V3d pt3d(pt[0], pt[1], 1.0);
        auto result = m_lines * pt3d;
        return result.maxCoeff() <= eps;
    }
    V2d projectPointInLines(const V2d &pt) {
        V2d result = pt;
        if (isPointInLines(pt)) return result;
        //assumed that m_lines is normalized when setLines is called
        auto lines_normalized = m_lines;
        Eigen::VectorXd distances = lines_normalized * V3d(pt[0], pt[1], 1.0);
        // loop through projections and see if there is a valid one
        for (int i = 0; i < lines_normalized.rows(); ++i) {
            if (distances[i] <= 0.0) continue;
            V2d projection = pt + (-distances[i] * lines_normalized.row(i).segment(0, 2).transpose());
            if (isPointInLines(projection, 1e-6)) {
                result = projection;
                return result;
            }
        }
        // if no valid projection found, project onto closest corner
        // find closest corner
        size_t closest_corner = 0;
        double closest_distance = (m_corners[0] - pt).norm();
        for (size_t i = 1; i < m_corners.size(); ++i) {
            double distance = (m_corners[i] - pt).norm();
            if (distance < closest_distance) {
                closest_corner = i;
                closest_distance = distance;
            }
        }
        return m_corners[closest_corner];
    }

    virtual ~LocalGlobalGenericParametrizer(); // Out-of-line destructor needed since SPSDSystemSolver is incomplete type

    std::unique_ptr<SPSDSystemSolver> L;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
    void m_localStep();

    Real m_energy = 0.0;

    // Range of admissible stretch factors
    Real m_alphaMin = 1.0, m_alphaMax = 1.0;
    Real m_betaMin = 1.0, m_betaMax = 1.0;

    // Cached quantities computed from the local step
    aligned_std_vector<M23d> m_M_Bt;
    aligned_std_vector<M2d> m_R, m_U;
    Eigen::VectorXd m_stretch; // Controls anisotropy of the target singular values for the Jacobian ([alpha, beta]).
    Eigen::Matrix<Real, Eigen::Dynamic, 2> m_lambda;
    Eigen::MatrixXd m_lines; // Linear inequality constraints on the stretch variables
    std::vector<V2d> m_corners; // Corners of the polygon defined by m_lines
};

} // namespace parametrization
#endif /* end of include guard: LOCAL_GLOBAL_PARAMETRIZATION_HH */
