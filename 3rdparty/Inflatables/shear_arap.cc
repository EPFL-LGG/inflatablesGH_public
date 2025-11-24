#include "parametrization.hh"
#include "local_global_parametrization.hh"
#include "regularized_parametrization.hh"
#include <MeshFEM/MeshIO.hh>
#include <MeshFEM/MSHFieldWriter.hh>
#include <MeshFEM/GlobalBenchmark.hh>

#include "fd_tests.hh"

using Mesh = parametrization::Mesh;

int main(int argc, const char *argv[]) {
    if (argc != 2) {
        std::cerr << "usage: shear_arap in_mesh.obj" << std::endl;
        exit(-1);
    }

    std::string inPath = argv[1];

    std::vector<MeshIO::IOVertex > vertices;
    std::vector<MeshIO::IOElement> elements;
    MeshIO::load(inPath, vertices, elements);
    auto mesh = std::make_shared<Mesh>(elements, vertices);

    auto writeFlattenedMesh = [&](const parametrization::UVMap &uv, const std::string &path) {
        std::vector<MeshIO::IOVertex> flatVertices(vertices.size());
        for (size_t i = 0; i < flatVertices.size(); ++i)
            flatVertices[i].point << uv(i, 0), uv(i, 1), 0;
        MeshIO::save(path, flatVertices, elements);
    };

    // Initialize the parametrization with a conformal map.
    auto f = parametrization::lscm(*mesh);
    writeFlattenedMesh(f, "lscm_flattened.msh");

    parametrization::LocalGlobalParametrizer param(mesh, f);
    param.setAlpha(M_PI / 2);
    // param.setAlpha(1.0);

    std::cout.precision(19);

    // Local-global iterations
    const size_t nit = 2;
    for (size_t it = 0; it < nit; ++it) {
        std::cout << "ARAP energy:\t" << param.energy() << std::endl;
        param.runIteration();
        // writeFlattenedMesh(param.uv(), "it_" + std::to_string(it) + ".msh");
    }
    std::cout << "ARAP energy:\t" << param.energy() << std::endl;

    // Regularized parametrizer iterations
    parametrization::RegularizedParametrizer rparam(param);

    fd_tests(rparam);

    // BENCHMARK_REPORT_NO_MESSAGES();

#if 0
    MSHFieldWriter writer("debug.msh", vertices, elements);
    VectorField<Real, 3> b0(nt), b1(nt);
    for (size_t i = 0; i < nt; ++i) {
        b0(i) = B[i].col(0);
        b1(i) = B[i].col(1);
    }

    writer.addField("b0", b0, DomainType::PER_ELEMENT);
    writer.addField("b1", b1, DomainType::PER_ELEMENT);

    // Compute the finite strain measure
    using M3D = Eigen::Matrix3d;
    SymmetricMatrixField<Real, 3> greenStrain(nt), corotatedStrain(nt), stretch(nt);
    for (const auto triRest : param.mesh().elements()) {
        const size_t ti = triRest.index();
        const auto triDefo = defoMesh.element(ti);
        // pt(lambda(x)) = [p0 | p1 | p2](l0, l1, l2)^T
        M3D P;
        P.col(0) = triDefo.node(0)->p;
        P.col(1) = triDefo.node(1)->p;
        P.col(2) = triDefo.node(2)->p;
        M3D F = P * triRest->gradBarycentric().transpose() // each column of gradBarycentric is the gradient of a barycentric coordinate
              + triDefo->normal() * triRest->normal().transpose();
#endif

    return 0;
}
