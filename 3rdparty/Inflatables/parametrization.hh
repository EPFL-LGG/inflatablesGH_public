#ifndef PARAMETRIZATION_HH
#define PARAMETRIZATION_HH

#include "DualLaplacianStencil.hh"

#include <MeshFEM/FEMMesh.hh>
#include <MeshFEM/SparseMatrices.hh>
#include <memory>
#include <utility>

#include "SVDSensitivity.hh"

#include "BendingStiffnessIntegralSensitivity.hh"
#include "circular_mean.hh"

namespace parametrization {


struct LinearInequality {
    std::vector<int>  vars;
    std::vector<Real> coeffs;
    Real              constPart;
};

using Mesh = FEMMesh<2, 1, Vector3D>; // Piecewise linear triangle mesh embedded in R^3
using UVMap = Eigen::Matrix<Real, Eigen::Dynamic, 2, Eigen::ColMajor>;
using NDMap = Eigen::MatrixXd;

struct SPSDSystemSolver : public SPSDSystem<Real> {
    using Base = SPSDSystem<Real>;
    using Base::Base;
};

// Compute a least-squares conformal parametrization with the global scale factor
// chosen to minimize the L2 norm of the pointwise area distortion.
UVMap lscm(const Mesh &mesh);

// Compute a harmonic map with prescribed boundary positions (in 2D or 3D)
NDMap harmonic(const Mesh &mesh, NDMap &boundaryData);

enum class Domain { UV, XYZ };

enum class BarrierEvalLevel { Val, Grad, Hess};
struct Barrier {
    virtual void operator()(const Real x, std::vector<Real> & result, const BarrierEvalLevel evalLevel) const = 0;
    virtual ~Barrier() { }
};
struct SigmoidBarrier : public Barrier {
    SigmoidBarrier(Real a, Real b) : m_a(a), m_b(b) { 
        m_offset = 2 * m_a / (m_a + 1.0);
        m_scale = (1.0 + m_a)/(1.0 - m_a);
    }
    virtual void operator()(const Real x, std::vector<Real> & result, const BarrierEvalLevel evalLevel = BarrierEvalLevel::Val) const override {
        result.clear();
        switch (evalLevel)
        {
        case BarrierEvalLevel::Val: {
            const Real e_minus_bx = std::exp(-m_b * x);
            const Real e_bx = std::exp(m_b * x);
            const Real val = (m_a / (m_a + e_minus_bx) + m_a / (m_a + e_bx) - m_offset)*m_scale;
            result.push_back(val);
            break;
        }
            
        case BarrierEvalLevel::Grad: {
            const Real e_minus_bx = std::exp(-m_b * x);
            const Real e_bx = std::exp(m_b * x);
            
            const Real val = (m_a / (m_a + e_minus_bx) + m_a / (m_a + e_bx) - m_offset)*m_scale;
            const Real grad = m_a * m_b * (e_minus_bx/std::pow(m_a + e_minus_bx, 2.0) - e_bx/std::pow(m_a + e_bx, 2.0)) * m_scale;
            result.push_back(val);
            result.push_back(grad);
            break;
        }
            
        case BarrierEvalLevel::Hess:{
            const Real e_minus_bx = std::exp(-m_b * x);
            const Real e_bx = std::exp(m_b * x);
            const Real val = (m_a / (m_a + e_minus_bx) + m_a / (m_a + e_bx) - m_offset)*m_scale;
            const Real grad = m_a * m_b * (e_minus_bx/std::pow(m_a + e_minus_bx, 2.0) - e_bx/std::pow(m_a + e_bx, 2.0)) * m_scale;
            const Real hess = (- m_a * m_b * m_b * (e_minus_bx * (m_a - e_minus_bx)/std::pow(m_a + e_minus_bx, 3.0) + e_bx * (m_a - e_bx)/std::pow(m_a + e_bx, 3.0))) * m_scale;
            result.push_back(val);
            result.push_back(grad);
            result.push_back(hess);
            break;
        }
        default:
            break;
        }
    }
private:
    Real m_a, m_b, m_offset, m_scale;
};

struct SigmoidBarrier2 : public Barrier {
    SigmoidBarrier2(Real a, Real b) : m_a(a), m_b(b) { 
        m_offset = 2 * m_a / (m_a + 1.0);
        m_scale = (1.0 + m_a)/(1.0 - m_a);
    }
    virtual void operator()(const Real x, std::vector<Real> & result, const BarrierEvalLevel evalLevel = BarrierEvalLevel::Val) const override {
        result.clear();
        switch (evalLevel)
        {
        case BarrierEvalLevel::Val: {
            const Real e_minus_bx2 = std::exp(-m_b * x * x);
            const Real e_bx2 = std::exp(m_b * x * x);
            const Real val = (m_a / (m_a + e_minus_bx2) + m_a / (m_a + e_bx2) - m_offset)*m_scale;
            result.push_back(val);
            break;
        }
            
        case BarrierEvalLevel::Grad: {
            const Real e_minus_bx2 = std::exp(-m_b * x * x);
            const Real e_bx2 = std::exp(m_b * x * x);
            
            const Real val = (m_a / (m_a + e_minus_bx2) + m_a / (m_a + e_bx2) - m_offset)*m_scale;
            const Real grad = 2.0 * m_a * m_b * x * (e_minus_bx2/std::pow(m_a + e_minus_bx2, 2.0) - e_bx2/std::pow(m_a + e_bx2, 2.0)) * m_scale;
            result.push_back(val);
            result.push_back(grad);
            break;
        }
            
        case BarrierEvalLevel::Hess:{
            const Real e_minus_bx2 = std::exp(-m_b * x * x);
            const Real e_bx2 = std::exp(m_b * x * x);
            
            const Real val = (m_a / (m_a + e_minus_bx2) + m_a / (m_a + e_bx2) - m_offset)*m_scale;
            const Real grad = 2.0 * m_a * m_b * x * (e_minus_bx2/std::pow(m_a + e_minus_bx2, 2.0) - e_bx2/std::pow(m_a + e_bx2, 2.0)) * m_scale;
            const Real hess = (2.0 * m_a * m_b * (e_minus_bx2/std::pow(m_a + e_minus_bx2, 2.0) - e_bx2/std::pow(m_a + e_bx2, 2.0)) - 4.0 * m_a * m_b * m_b * x * x * (e_minus_bx2 * (m_a - e_minus_bx2)/std::pow(m_a + e_minus_bx2, 3.0) + e_bx2 * (m_a - e_bx2)/std::pow(m_a + e_bx2, 3.0))) * m_scale;
            result.push_back(val);
            result.push_back(grad);
            result.push_back(hess);
            break;
        }
        default:
            break;
        }
    }
private:
    Real m_a, m_b, m_offset, m_scale;
};
struct AbsLogBarrier : public Barrier {
    AbsLogBarrier(Real a, Real b) : m_a(a), m_b(b) { 
    // const Real exp_abs_diff_i = std::exp(std::abs(m_stretch[ti] - m_stretch[nt + ti]));
    // const Real exp_abs_diff_j = std::exp(std::abs(m_stretch[tj] - m_stretch[nt + tj]));
    // const Real barrier_i = m_b_a + m_b_b * std::log(exp_abs_diff_i - m_b_c);
    // const Real barrier_j = m_b_a + m_b_b * std::log(exp_abs_diff_j - m_b_c);
    // Real b_sign_i = 0.0, b_sign_j = 0.0;
    // // b_sign_i = 1.0; b_sign_j = 1.0;
    // if (m_stretch[ti] - m_stretch[nt + ti] > 0) b_sign_i =  1.0;
    // if (m_stretch[ti] - m_stretch[nt + ti] < 0) b_sign_i = -1.0;
    // if (m_stretch[tj] - m_stretch[nt + tj] > 0) b_sign_j =  1.0;
    // if (m_stretch[tj] - m_stretch[nt + tj] < 0) b_sign_j = -1.0;
    // const Real barrier_grad_i = m_b_b * exp_abs_diff_i/(exp_abs_diff_i - m_b_c) * b_sign_i;
    // const Real barrier_grad_j = m_b_b * exp_abs_diff_j/(exp_abs_diff_j - m_b_c) * b_sign_j;
    // const Real barrier_hess_i = - barrier_grad_i * barrier_j * m_b_c /(exp_abs_diff_i - m_b_c) * b_sign_i;
    // const Real barrier_hess_j = - barrier_grad_j * barrier_i * m_b_c /(exp_abs_diff_j - m_b_c) * b_sign_j;
                    
    }
    virtual void operator()(const Real x, std::vector<Real> & result, const BarrierEvalLevel evalLevel = BarrierEvalLevel::Val) const override {
        result.clear();
        switch (evalLevel)
        {
        case BarrierEvalLevel::Val: {
            const Real exp_abs = std::exp(std::abs(x));
            const Real val = m_a + m_b * std::log(exp_abs - m_c);
            result.push_back(val);
            break;
        }
            
        case BarrierEvalLevel::Grad: {
            const Real exp_abs = std::exp(std::abs(x));
            Real sign = (x > 0) ? 1.0 : -1.0;
            sign = (x == 0) ? 0.0 : sign;
            
            const Real val = m_a + m_b * std::log(exp_abs - m_c);
            const Real grad = m_b * exp_abs/(exp_abs - m_c) * sign;
            result.push_back(val);
            result.push_back(grad);
            break;
        }
            
        case BarrierEvalLevel::Hess:{


            const Real exp_abs = std::exp(std::abs(x));
            Real sign = (x > 0) ? 1.0 : -1.0;
            sign = (x == 0) ? 0.0 : sign;
            
            const Real val = m_a + m_b * std::log(exp_abs - m_c);
            const Real grad = m_b * exp_abs/(exp_abs - m_c) * sign;
            const Real hess = - grad * m_c /(exp_abs - m_c) * sign;
            result.push_back(val);
            result.push_back(grad);
            result.push_back(hess);
            break;
        }
        default:
            break;
        }
    }
private:
    Real m_a, m_b, m_c;
};



// Abstract base class for parametrization algorithms requiring a per-triangle tangent space basis.
struct Parametrizer {
    using  V2d = Eigen::Vector2d;
    using  V3d = Eigen::Vector3d;
    using  V4d = Eigen::Vector4d;
    using  M2d = Eigen::Matrix2d;
    using  M3d = Eigen::Matrix3d;
    using M23d = Eigen::Matrix<Real, 2, 3>;
    using MX3d = Eigen::Matrix<Real, Eigen::Dynamic, 3>;
    using  M4d = Eigen::Matrix4d;

    Parametrizer(std::shared_ptr<Mesh> inMesh) : m_mesh(inMesh) {
        const auto &m = mesh();
        const size_t nt = m.numTris();
        m_B.resize(nt);
        for (auto tri : m.elements()) {
            Vector3D b0 = (tri.node(1)->p - tri.node(0)->p).normalized();
            Vector3D b1 = tri->normal().cross(b0);
            m_B[tri.index()].col(0) = b0;
            m_B[tri.index()].col(1) = b1;
        }
    }

    // Replace the parametrization, also trigging an update to derived class' cached data
    // (e.g., re-running local step of LocalGlobalParametrizer).
    void setUV(Eigen::Ref<const UVMap> uv, bool debug = false);
    void setUVDebug(Eigen::Ref<const UVMap> uv);

    const UVMap &uv() { return m_uv; }

    size_t numUVVars() const { return m_uv.size(); }

          Mesh &mesh()       { return *m_mesh; }
    const Mesh &mesh() const { return *m_mesh; }

    virtual Real energy() const = 0;
    size_t numFlips() const { return m_flipCount; }

    const aligned_std_vector<Eigen::Matrix<Real, 3, 2>> &B() const { return m_B; }

    // Access the mesh shared pointer from this instance
    std::shared_ptr<Mesh> meshPtr() { return m_mesh; }

    // Jacobian of the uv map on triangle i
    const M23d &jacobian(size_t i) const { return m_J[i]; }

    // Angle between the local frame's x axis and left/right singular vectors
    // The basis for the left stretch vector is the "u" axis of the UV domain, and
    // the basis for the right stretch vector is m_B[i].col(0).
    virtual Real  leftStretchAngle(size_t i) const = 0;
    virtual Real rightStretchAngle(size_t i) const = 0;
    virtual const Eigen::VectorXd getAlphas() const = 0;
    virtual const Eigen::VectorXd getBetas() const { return Eigen::VectorXd::Ones(mesh().numElements());}

    Eigen::VectorXd leftStretchAngles() const {
        const size_t nt = mesh().numTris();
        Eigen::VectorXd result(nt);
        for (size_t i = 0; i < nt; ++i) result[i] = leftStretchAngle(i);
        return result;
    }

    // Average the stretching directions from the triangles onto the vertices
    // in a smoothness-aware way (only average if the directions are reasonably
    // coherent; otherwise pick the single incident direction that "agrees" with
    // the majority of the rest)
    Eigen::VectorXd perVertexLeftStretchAngles(double agreementThreshold = M_PI / 4) const;
    Eigen::VectorXd perVertexAlphas() const;
    Eigen::VectorXd perVertexBetas() const;

    // Return a higher resolution flattened mesh with smootly interpolated uvs at its vertices.
    std::tuple<std::shared_ptr<Mesh>, UVMap> upsampledUV(size_t nsubdiv = 2) const;

    // Return a higher resolution flattened mesh with smootly interpolated stretching angles (phis) and magnitudes (alphas) at its vertices.
    std::tuple<std::shared_ptr<Mesh>, Eigen::VectorXd, Eigen::VectorXd> upsampledVertexLeftStretchAnglesAndMagnitudes(size_t nsubdiv = 2, double agreementThreshold = M_PI / 4) const;

    Eigen::VectorXd rightStretchAngles() const {
        const size_t nt = mesh().numTris();
        Eigen::VectorXd result(nt);
        for (size_t i = 0; i < nt; ++i) result[i] = rightStretchAngle(i);
        return result;
    }

    virtual void parametrizationUpdated() = 0;

    size_t numLinearInequalityConstraints() const { return 0; }
    std::vector<LinearInequality> getLinearInequalityConstraints() const { return std::vector<LinearInequality>(); }


    virtual ~Parametrizer() { }

protected:
    UVMap m_uv;         // Current parametrization (|V| x 2 matrix)

    // Orthonormal basis for each triangle's tangent space.
    aligned_std_vector<Eigen::Matrix<Real, 3, 2>> m_B;
    aligned_std_vector<M23d> m_J;

    std::shared_ptr<Mesh> m_mesh;

    // Cached quantities computed from current iterate
    size_t m_flipCount = 0;
};

} // namespace parametrization

#endif /* end of include guard: PARAMETRIZATION_HH */