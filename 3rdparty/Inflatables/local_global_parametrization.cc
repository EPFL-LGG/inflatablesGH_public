#include "local_global_parametrization.hh"

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

LocalGlobalParametrizer::LocalGlobalParametrizer(const std::shared_ptr<Mesh> &inMesh, const UVMap &uvInit)
    : Parametrizer(inMesh)
{
    setUV(uvInit);
    // Constant Laplacian matrix used throughout local/global iterations
    L = std::make_unique<SPSDSystemSolver>(Laplacian::construct(mesh()));
    L->fixVariables(std::vector<size_t>{0}, std::vector<Real>{0.0}); // fix first vertex's coordinate in u or v axis.
}

// Replace the parametrization, updating the local-global energy energy (running the local step)
void LocalGlobalParametrizer::m_localStep() {
    BENCHMARK_START_TIMER_SECTION("Local step");
    const auto &m = mesh();
    const size_t ne = m.numElements();
    m_J.resize(ne);
    m_M_Bt.resize(ne);
    m_R.resize(ne);
    m_U.resize(ne);
    m_alpha.resize(ne);
    m_lambda.resize(ne, 2);

    // Local step: compute closest admissible Jacobian U R(theta) [alpha 0; 0 1] R(theta)^T
    // and construct the RHS for the global step.
    auto process_tri = [&](const size_t ti) {
        M2d JB = m_J[ti] * m_B[ti];

        // Decompose JB = U R Lambda R^T where "U" is a post-stretch rotation
        // in the parametric domain, and R Lambda R^T describes how material is
        // stretched in the (b0, b1) tangent plane.
        // Column j of R is the principal stretch vector stretched by lambda[j]
        auto &R = m_R[ti];
        auto &U = m_U[ti];
        auto lambda = m_lambda.row(ti);
        {
            Eigen::JacobiSVD<M2d> svd(JB, Eigen::ComputeFullU | Eigen::ComputeFullV);
            M2d tmp = svd.matrixU();
            lambda = svd.singularValues();
            // Note: we want to make sure both U *and* R are true rotations, not reflections.
            // If det(JB) < 0, a singular value needs to be flipped negative
            // (along with its column in tmp), which will guarantee a positive
            // determinant of U = tmp * V^T.
            // But R could still be a reflection; we negate its last column in
            // this case (which leaves the mapping U R Lambda R^T unchanged).
            if (JB.determinant() < 0) {
                tmp.col(1) *= -1;
                lambda[1]  *= -1;
            }
            U = tmp * svd.matrixV().transpose(); // positive determinant
            R = svd.matrixV();
            if (R.determinant() < 0) { R.col(1) *= -1; }
        }
        m_alpha[ti] = std::min(m_alphaMax, std::max(m_alphaMin, lambda[0]));
        Vector2D lambda_target(m_alpha[ti], 1.0);

        M2d M = U * (R * (lambda_target.asDiagonal() * R.transpose()));
        m_M_Bt[ti] = M * m_B[ti].transpose();
    };

    const size_t nt = m.numTris();
#if MESHFEM_WITH_TBB
    tbb::parallel_for(tbb::blocked_range<size_t>(0, nt), [&](const tbb::blocked_range<size_t> &b) { for (size_t ti = b.begin(); ti < b.end(); ++ti) process_tri(ti); });
#else
    for (size_t ti = 0; ti < nt; ++ti) process_tri(ti);
#endif

    // Update the energy
    {
        // Accumulate in temporary so other threads don't read intermediate values.
        Real energy = 0;
        for (size_t ti = 0; ti < nt; ++ti)
            energy += 0.5 * (m_J[ti] - m_M_Bt[ti]).squaredNorm() * m.element(ti)->volume();
        m_energy = energy;
    }

    BENCHMARK_STOP_TIMER_SECTION("Local step");
}

void LocalGlobalParametrizer::runIteration() {
    const auto &m = mesh();
    const size_t nv = m.numVertices();

    // Global step
    BENCHMARK_START_TIMER_SECTION("Global step");
    // Compute RHS vectors
    UVMap rhs_uv = UVMap::Zero(nv, 2);
    for (auto tri : m.elements()) {
        const size_t ti = tri.index();
        for (auto v : tri.vertices())
            rhs_uv.row(v.index()) += (m_M_Bt[ti] * tri->gradBarycentric().col(v.localIndex())) * tri->volume();
    }

    // Solve two Poisson equations
    UVMap uv_new(m_uv.rows(), 2);
    Eigen::VectorXd soln;
    L->solve(rhs_uv.col(0), soln);
    uv_new.col(0) = soln;
    L->solve(rhs_uv.col(1), soln);
    uv_new.col(1) = soln;
    BENCHMARK_STOP_TIMER_SECTION("Global step");

    setUV(uv_new); // Update Jacobians and run the next local step, allowing us to evaluate energy.
}

Real LocalGlobalParametrizer::leftStretchAngle(size_t i) const {
    Eigen::Rotation2D<Real> UR(getU(i) * getR(i));
    return UR.angle();
}

Real LocalGlobalParametrizer::rightStretchAngle(size_t i) const {
    return Eigen::Rotation2D<Real>(getR(i)).angle();
}

LocalGlobalParametrizer::~LocalGlobalParametrizer() { }


LocalGlobalGenericParametrizer::LocalGlobalGenericParametrizer(const std::shared_ptr<Mesh> &inMesh, const UVMap &uvInit)
    : Parametrizer(inMesh)
{
    setUV(uvInit);
    // Constant Laplacian matrix used throughout local/global iterations
    L = std::make_unique<SPSDSystemSolver>(Laplacian::construct(mesh()));
    L->fixVariables(std::vector<size_t>{0}, std::vector<Real>{0.0}); // fix first vertex's coordinate in u or v axis.
}

// Replace the parametrization, updating the local-global energy energy (running the local step)
void LocalGlobalGenericParametrizer::m_localStep() {
    BENCHMARK_START_TIMER_SECTION("Local step");
    const auto &m = mesh();
    const size_t ne = m.numElements();
    m_J.resize(ne);
    m_M_Bt.resize(ne);
    m_R.resize(ne);
    m_U.resize(ne);
    m_stretch.resize(2*ne);
    m_lambda.resize(ne, 2);

    // Local step: compute closest admissible Jacobian U R(theta) [alpha 0; 0 beta] R(theta)^T
    // and construct the RHS for the global step.
    auto process_tri = [&](const size_t ti) {
        M2d JB = m_J[ti] * m_B[ti];

        // Decompose JB = U R Lambda R^T where "U" is a post-stretch rotation
        // in the parametric domain, and R Lambda R^T describes how material is
        // stretched in the (b0, b1) tangent plane.
        // Column j of R is the principal stretch vector stretched by lambda[j]
        auto &R = m_R[ti];
        auto &U = m_U[ti];
        auto lambda = m_lambda.row(ti);
        {
            Eigen::JacobiSVD<M2d> svd(JB, Eigen::ComputeFullU | Eigen::ComputeFullV);
            M2d tmp = svd.matrixU();
            lambda = svd.singularValues();
            // Note: we want to make sure both U *and* R are true rotations, not reflections.
            // If det(JB) < 0, a singular value needs to be flipped negative
            // (along with its column in tmp), which will guarantee a positive
            // determinant of U = tmp * V^T.
            // But R could still be a reflection; we negate its last column in
            // this case (which leaves the mapping U R Lambda R^T unchanged).
            if (JB.determinant() < 0) {
                tmp.col(1) *= -1;
                lambda[1]  *= -1;
            }
            U = tmp * svd.matrixV().transpose(); // positive determinant
            R = svd.matrixV();
            if (R.determinant() < 0) { R.col(1) *= -1; }
        }
        
        m_stretch[ti] = std::min(m_alphaMax, std::max(m_alphaMin, lambda[0]));
        m_stretch[ne + ti] = std::min(m_betaMax, std::max(m_betaMin, lambda[1]));
        
        if ((m_lines.size() > 0) && !isPointInLines(V2d(m_stretch[ti], m_stretch[ne + ti]))) {

            auto feasible_stretch = projectPointInLines(V2d(m_stretch[ti], m_stretch[ne + ti]));
            m_stretch[ti] = feasible_stretch[0];
            m_stretch[ne + ti] = feasible_stretch[1];    
        }
        
        Vector2D lambda_target(m_stretch[ti], m_stretch[ne + ti]);

        M2d M = U * (R * (lambda_target.asDiagonal() * R.transpose()));
        m_M_Bt[ti] = M * m_B[ti].transpose();
    };

    const size_t nt = m.numTris();
#if MESHFEM_WITH_TBB
    tbb::parallel_for(tbb::blocked_range<size_t>(0, nt), [&](const tbb::blocked_range<size_t> &b) { for (size_t ti = b.begin(); ti < b.end(); ++ti) process_tri(ti); });
#else
    for (size_t ti = 0; ti < nt; ++ti) process_tri(ti);
#endif

    // Update the energy
    {
        // Accumulate in temporary so other threads don't read intermediate values.
        Real energy = 0;
        for (size_t ti = 0; ti < nt; ++ti)
            energy += 0.5 * (m_J[ti] - m_M_Bt[ti]).squaredNorm() * m.element(ti)->volume();
        m_energy = energy;
    }

    BENCHMARK_STOP_TIMER_SECTION("Local step");
}

void LocalGlobalGenericParametrizer::runIteration() {
    const auto &m = mesh();
    const size_t nv = m.numVertices();

    // Global step
    BENCHMARK_START_TIMER_SECTION("Global step");
    // Compute RHS vectors
    UVMap rhs_uv = UVMap::Zero(nv, 2);
    for (auto tri : m.elements()) {
        const size_t ti = tri.index();
        for (auto v : tri.vertices())
            rhs_uv.row(v.index()) += (m_M_Bt[ti] * tri->gradBarycentric().col(v.localIndex())) * tri->volume();
    }

    // Solve two Poisson equations
    UVMap uv_new(m_uv.rows(), 2);
    Eigen::VectorXd soln;
    L->solve(rhs_uv.col(0), soln);
    uv_new.col(0) = soln;
    L->solve(rhs_uv.col(1), soln);
    uv_new.col(1) = soln;
    BENCHMARK_STOP_TIMER_SECTION("Global step");

    setUV(uv_new); // Update Jacobians and run the next local step, allowing us to evaluate energy.
}

Real LocalGlobalGenericParametrizer::leftStretchAngle(size_t i) const {
    Eigen::Rotation2D<Real> UR(getU(i) * getR(i));
    return UR.angle();
}

Real LocalGlobalGenericParametrizer::rightStretchAngle(size_t i) const {
    return Eigen::Rotation2D<Real>(getR(i)).angle();
}

LocalGlobalGenericParametrizer::~LocalGlobalGenericParametrizer() { }

}