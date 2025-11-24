#include "parametrization.hh"

#include <MeshFEM/SparseMatrices.hh>
#include <MeshFEM/Laplacian.hh>
#include <MeshFEM/GlobalBenchmark.hh>
#include <MeshFEM/MeshIO.hh>
#include <complex>
#include <set>

#include <MeshFEM/ParallelAssembly.hh>
#include "subdivide_triangle.hh"
#include "curvature.hh"

namespace parametrization {

// Compute a least-squares conformal parametrization with the global scale factor
// chosen to minimize the L2 norm of the pointwise area distortion.
UVMap lscm(const Mesh &mesh) {
    const size_t nv = mesh.numVertices();
    UVMap uv(nv, 2);

    TripletMatrix<> K(2 * nv, 2 * nv);
    K.symmetry_mode = TripletMatrix<>::SymmetryMode::UPPER_TRIANGLE;

    // Assemble (upper triangle of) LSCM matrix K =
    // [L   A] = [L   A]
    // [A^T L]   [-A  L]
    // where L_ij = int grad phi_i . grad phi_j dA          is the P1 FEM Laplacian and
    //       A_ij = int n . (grad phi_j x grad phi_i) dA    is the skew symmetric "parametric area calculator"
    //            = int s_ij (1 / 2A) dA = sum_T s_ij|_T / 2
    //       s_ij|_T = 1 if local(i) == local(j) + 1, -1 if local(i) == local(j) - 1, 0 otherwise   (this is as evaluated on a particular triangle T)
    // This is the quadratic form for [u v] giving the LSCM energy.
    // Note: the interior edge contributions to the "area calculator" matrix cancel out, and it can be written as an integral over the boundary.
    // However, if we want to support varying triangle weights as recommended in Spectral Conformal Parametrization,
    // we need to compute the per-triangle contribution. (This seems to actually be a bad idea though--probably they just need to incorporate a mass matrix in their generalized eigenvalue problem.)
    for (auto tri : mesh.elements()) {
        const auto &gradLambda = tri->gradBarycentric();
        for (auto ni : tri.nodes()) {
            for (auto nj : tri.nodes()) {
                if (ni.index() > nj.index()) continue; // lower triangle
                // Symmetric Laplacian blocks
                const Real val = gradLambda.col(ni.localIndex()).dot(gradLambda.col(nj.localIndex())) * tri->volume();
                K.addNZ(     ni.index(),      nj.index(), val); // (u, u) block
                K.addNZ(nv + ni.index(), nv + nj.index(), val); // (v, v) block

                // Skew symmetric A block (u, v)
                if (ni.localIndex() == nj.localIndex()) continue;
                int s = (ni.localIndex() == (nj.localIndex() + 1) % 3) ? 1.0 : -1.0;
                K.addNZ(ni.index(), nv + nj.index(),  0.5 * s);
                K.addNZ(nj.index(), nv + ni.index(), -0.5 * s);
            }
        }
    }

    SPSDSystemSolver Ksys(K);

    // Pin down the null-space (scale, rotation) by fixing two vertices' UVs: vertex 0 and the vertex furthest from it.
    {
        Point3D p0 = mesh.node(0)->p;
        Real furthestDist = 0;
        size_t furthestIdx = 0;

        for (auto n : mesh.nodes()) {
            Real dist = (n->p - p0).norm();
            if (dist > furthestDist) {
                furthestDist = dist;
                furthestIdx = n.index();
            }
        }

        std::vector<size_t>    fixedVars = {0, furthestIdx, nv, nv + furthestIdx};
        std::vector<Real> fixedVarValues = {0.0, furthestDist, 0.0, 0.0};

        Ksys.fixVariables(fixedVars, fixedVarValues);
    }

    Eigen::VectorXd soln;
    Ksys.solve(Eigen::VectorXd::Zero(2 * nv), soln);
    Eigen::Map<Eigen::VectorXd>(uv.data(), 2 * nv) = soln;

    // Compute per-triangle areas before and after parametrization
    Eigen::VectorXd origArea(mesh.numTris()), paramArea(mesh.numTris());
    for (const auto t : mesh.elements()) {
        origArea[t.index()] = t->volume();
        std::array<Point2D, 3> poly;
        for (auto v : t.vertices())
            poly[v.localIndex()] = uv.row(v.index()).transpose();
        paramArea[t.index()] = area(poly);
    }

    // Scale the full parametrization to minimize the squared difference in areas
    // min_s 1/2 ||s paramArea - origArea||^2 ==> (s paramArea - origArea) . paramArea = 0 ==> s = (origArea . paramArea) / ||paramArea||^2
    uv *= std::sqrt(origArea.dot(paramArea) / paramArea.squaredNorm());

    return uv;
}

NDMap harmonic(const Mesh &mesh, NDMap &boundaryData) {
    const size_t nbn = mesh.numBoundaryNodes(),
                 nn  = mesh.numNodes();
    if (size_t(boundaryData.rows()) != nbn) throw std::runtime_error("Invalid boundary data size");
    size_t numComponents = boundaryData.cols();

    NDMap result(nn, numComponents);

    auto L = Laplacian::construct(mesh);
    L.sumRepeated();
    L.needs_sum_repeated = false;
    SPSDSystemSolver Lsys(L);

    // Avoid resetting the SPSDSystemSolver and fixing variables anew for each component solve
    // by always fixing the boundary variables to "0" and directly computing the "load"
    // contributed by these constraints
    std::vector<size_t> fixedVars(nbn);
    for (auto bn : mesh.boundaryNodes())
        fixedVars[bn.index()] = bn.volumeNode().index();
    Lsys.fixVariables(fixedVars, std::vector<double>(nbn, 0.0));
    std::vector<double> negDirichletValues(nn, 0.0);
    std::vector<double> soln;

    for (size_t c = 0; c < numComponents; ++c) {
        for (auto bn : mesh.boundaryNodes())
            negDirichletValues[bn.volumeNode().index()] = -boundaryData(bn.index(), c);
        auto rhs = L.apply(negDirichletValues);
        Lsys.solve(rhs, soln);

        for (auto n : mesh.nodes()) {
            auto bn = n.boundaryNode();
            result(n.index(), c) = bn ? boundaryData(bn.index(), c) : soln[n.index()];
        }
    }

    return result;
}

void Parametrizer::setUV(Eigen::Ref<const UVMap> uv, bool debug /*= false*/) {
    const auto &m = mesh();
    if (size_t(uv.rows()) != m.numVertices()) throw std::runtime_error("Invalid parametrization size");
    m_uv = uv;

    // Update the cached Jacobians and count flips
    const size_t nt = m.numTris();
    m_J.resize(nt);
    M23d f_restrict_T;
    m_flipCount = 0;
    for (const auto &tri : m.elements()) {
        f_restrict_T.col(0) = m_uv.row(tri.vertex(0).index());
        f_restrict_T.col(1) = m_uv.row(tri.vertex(1).index());
        f_restrict_T.col(2) = m_uv.row(tri.vertex(2).index());
        auto &J = m_J[tri.index()];
        
        if (debug) {
            std::cout<<"f_restrict_T"<<std::endl;
            std::cout<<f_restrict_T<<std::endl;
            std::cout<<"tri->gradBarycentric().transpose()"<<std::endl;
            std::cout<<tri->gradBarycentric().transpose()<<std::endl;
        }
        

        J = f_restrict_T * tri->gradBarycentric().transpose();
        if ((J * m_B[tri.index()]).determinant() < 0) ++m_flipCount;
    }

    parametrizationUpdated(); // Notify derived class that the parametrization has been updated (invalidate cache)
}

void Parametrizer::setUVDebug(Eigen::Ref<const UVMap> uv) {
    const auto &m = mesh();
    if (size_t(uv.rows()) != m.numVertices()) throw std::runtime_error("Invalid parametrization size");
    m_uv = uv;

    // Update the cached Jacobians and count flips
    const size_t nt = m.numTris();
    m_J.resize(nt);
    M23d f_restrict_T;
    m_flipCount = 0;
    for (const auto &tri : m.elements()) {
        f_restrict_T.col(0) = m_uv.row(tri.vertex(0).index());
        f_restrict_T.col(1) = m_uv.row(tri.vertex(1).index());
        f_restrict_T.col(2) = m_uv.row(tri.vertex(2).index());
        auto &J = m_J[tri.index()];
        J = f_restrict_T * tri->gradBarycentric().transpose();
        if ((J * m_B[tri.index()]).determinant() < 0) ++m_flipCount;
    }

    // parametrizationUpdated(); // Notify derived class that the parametrization has been updated (invalidate cache)
}

Eigen::VectorXd Parametrizer::perVertexLeftStretchAngles(double /* agreementThreshold */) const {
    const auto &m = mesh();
    Eigen::VectorXd result(m.numVertices());

    std::vector<double> twiceIncidentAngles;
    for (const auto &v : m.vertices()) {
        twiceIncidentAngles.clear();
        for (const auto &he : v.incidentHalfEdges())
            if (he.tri()) twiceIncidentAngles.push_back(2 * leftStretchAngle(he.tri().index()));
        result[v.index()] = 0.5 * circularMean(twiceIncidentAngles);
    }

    return result;
}

Eigen::VectorXd Parametrizer::perVertexAlphas() const {
    const auto &m = mesh();
    Eigen::VectorXd result(m.numVertices());
    const Eigen::VectorXd &alphas = getAlphas();

    for (auto v : mesh().vertices()) {
        double &alpha = result[v.index()];
        alpha = 0;
        size_t tri_valence = 0;
        for (auto he : v.incidentHalfEdges()) {
            if (!he.tri()) continue;
            alpha += alphas[he.tri().index()];
            ++tri_valence;
        }
        alpha /= tri_valence;
    }

    return result;
}

Eigen::VectorXd Parametrizer::perVertexBetas() const {
    const auto &m = mesh();
    Eigen::VectorXd result(m.numVertices());
    const Eigen::VectorXd &betas = getBetas();

    for (auto v : mesh().vertices()) {
        double &beta = result[v.index()];
        beta = 0;
        size_t tri_valence = 0;
        for (auto he : v.incidentHalfEdges()) {
            if (!he.tri()) continue;
            beta += betas[he.tri().index()];
            ++tri_valence;
        }
        beta /= tri_valence;
    }

    return result;
}

std::tuple<std::shared_ptr<Mesh>, UVMap>
Parametrizer::upsampledUV(size_t nsubdiv) const {
    std::tuple<std::shared_ptr<Mesh>, UVMap> result;

    std::vector<MeshIO::IOVertex > subVertices;
    std::vector<MeshIO::IOElement> subElements;
    aligned_std_vector<V2d> subUV;

    const auto &m = mesh();
    PointGluingMap indexForPoint;
    for (const auto &tri : m.elements()) {
        auto newPt = [&](const Point3D &p, double lambda_0, double lambda_1, double lambda_2) {
            subVertices.emplace_back(p);
            subUV.push_back(lambda_0 * m_uv.row(tri.vertex(0).index()) +
                            lambda_1 * m_uv.row(tri.vertex(1).index()) +
                            lambda_2 * m_uv.row(tri.vertex(2).index()));
            return subVertices.size() - 1;
        };

        subdivide_triangle(nsubdiv,
                tri.vertex(0).node()->p,
                tri.vertex(1).node()->p,
                tri.vertex(2).node()->p,
                indexForPoint,
                newPt, [&](size_t i0, size_t i1, size_t i2) { subElements.emplace_back(i0, i1, i2); });
    }

    std::get<0>(result) = std::make_shared<Mesh>(subElements, subVertices);
    auto &fineUV = std::get<1>(result);
    fineUV.resize(subUV.size(), 2);

    for (size_t i = 0; i < subUV.size(); ++i)
        fineUV.row(i) = subUV[i];

    return result;
}

// Return Stretches [*alphas , *betas] if beta is not constant, otherwise just return Alphas
std::tuple<std::shared_ptr<Mesh>, Eigen::VectorXd, Eigen::VectorXd>
Parametrizer::upsampledVertexLeftStretchAnglesAndMagnitudes(size_t nsubdiv, double agreementThreshold) const {
    std::tuple<std::shared_ptr<Mesh>, Eigen::VectorXd, Eigen::VectorXd> result;

    std::vector<MeshIO::IOVertex > subVertices;
    std::vector<MeshIO::IOElement> subElements;

    std::cout<<"upsampledVertexLeftStretchAnglesAndMagnitudes"<<std::endl;
    auto coarseVertexAngles = perVertexLeftStretchAngles(agreementThreshold);
    Eigen::VectorXd coarseVertexStretches = perVertexAlphas();
    if (getBetas() != Eigen::VectorXd::Ones(mesh().numElements())) {
        const auto nv = coarseVertexStretches.size();
        coarseVertexStretches.resize(2*nv);
        coarseVertexStretches.head(nv) = perVertexAlphas();
        coarseVertexStretches.tail(nv) = perVertexBetas();
    }

    std::vector<double> subAngles, subAlphas, subBetas;

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

            if (cornerAngleVec.size() != barycentricDenominator) throw std::runtime_error("Invalid barycentric coordinate weights");

            subVertices.emplace_back(p);
            subAngles.push_back(0.5 * circularMean(cornerAngleVec));
            subAlphas.push_back(lambda_0 * coarseVertexStretches[tri.vertex(0).index()] +
                                lambda_1 * coarseVertexStretches[tri.vertex(1).index()] +
                                lambda_2 * coarseVertexStretches[tri.vertex(2).index()]);
            const auto nv = mesh().numVertices();
            if (coarseVertexStretches.size() > nv)
                subBetas.push_back(lambda_0 * coarseVertexStretches[nv + tri.vertex(0).index()] +
                                    lambda_1 * coarseVertexStretches[nv + tri.vertex(1).index()] +
                                    lambda_2 * coarseVertexStretches[nv + tri.vertex(2).index()]);

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
    std::cout<<"subAlphas.size() = "<<subAlphas.size()<<std::endl;
    std::cout<<"subBetas.size() = "<<subBetas.size()<<std::endl;
    if (subBetas.size() > 0)
        subAlphas.insert( subAlphas.end(), subBetas.begin(), subBetas.end() );
    std::get<2>(result) = Eigen::Map<Eigen::VectorXd>(subAlphas.data(), subAlphas.size());
    return result;
}
} // namespace parametrization