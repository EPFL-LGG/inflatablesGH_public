#include <iostream>
#include <iomanip>
#include <sstream>
#include <utility>
#include <memory>
#include <functional>

#include <MeshFEM/GlobalBenchmark.hh>
#include <MeshFEM/MeshIO.hh>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/iostream.h>
#include <pybind11/functional.h>
namespace py = pybind11;

#include "../DualLaplacianStencil.hh" // Must be included before MeshFEM's Triplet is defined...

#include "../parametrization.hh"
#include "../local_global_parametrization.hh"
#include "../regularized_parametrization.hh"
#include "../parametrization_knitro.hh"
#include "../parametrization_newton.hh"
#include "../generic_parametrization.hh"
#include "../pattern_parametrization.hh"

#include "../circular_mean.hh"
#include "../EllipsePointParameters.hh"

template <typename T>
std::string to_string_with_precision(const T &val, const int n = 6) {
    std::ostringstream ss;
    ss << std::setprecision(n) << val;
    return ss.str();

}
// Conversion of std::tuple to and from a py::tuple, since pybind11 doesn't seem to provide this...
template<typename... Args, size_t... Idxs>
py::tuple to_pytuple_helper(const std::tuple<Args...> &args, std::index_sequence<Idxs...>) {
    return py::make_tuple(std::get<Idxs>(args)...);
}

template<typename... Args>
py::tuple to_pytuple(const std::tuple<Args...> &args) {
    return to_pytuple_helper(args, std::make_index_sequence<sizeof...(Args)>());
}

template<class OutType>
struct FromPytupleImpl;

template<typename... Args>
struct FromPytupleImpl<std::tuple<Args...>> {
    template<size_t... Idxs>
    static auto run_helper(const py::tuple &t, std::index_sequence<Idxs...>) {
        return std::make_tuple((t[Idxs].cast<Args>())...);
    }
    static auto run(const py::tuple &t) {
        if (t.size() != sizeof...(Args)) throw std::runtime_error("Mismatched tuple size for py::tuple to std::tuple conversion.");
        return run_helper(t, std::make_index_sequence<sizeof...(Args)>());
    }
};

template<class OutType>
OutType from_pytuple(const py::tuple &t) {
    return FromPytupleImpl<OutType>::run(t);
}

using namespace parametrization;

PYBIND11_MODULE(inflatables_parametrization, m) {
    m.doc() = "Shear metric parametrization";

    py::module::import("py_newton_optimizer");
    py::module::import("mesh_utilities");

    ////////////////////////////////////////////////////////////////////////////////
    // Mesh construction (for mesh type used by parametrization routines)
    ////////////////////////////////////////////////////////////////////////////////
    // WARNING: Mesh's holder type is a shared_ptr; returning a unique_ptr will lead to a dangling pointer in the current version of Pybind11
    m.def("Mesh", [](const std::string &path) { return std::shared_ptr<Mesh>(Mesh::load(path)); }, py::arg("path"));
    m.def("Mesh", [](const Eigen::MatrixX3d &V, const Eigen::MatrixX3i &F) { return std::make_shared<Mesh>(F, V); }, py::arg("V"), py::arg("F"));

    ////////////////////////////////////////////////////////////////////////////////
    // Free-standing functions
    ////////////////////////////////////////////////////////////////////////////////
    m.def("lscm",     &lscm,     py::arg("mesh"),                          "Compute least-squares conformal parametrization");
    m.def("harmonic", &harmonic, py::arg("mesh"), py::arg("boundaryData"), "Compute harmonic map with given dirichlet data");
    m.def("regularized_parametrization_knitro", &regularized_parametrization_knitro<RegularizedParametrizer   >, py::arg("rparam"), py::arg("maxIter"), py::arg("fixedVars"), py::arg("gradTol") = 1e-8);
    m.def("pattern_parametrization_knitro", &regularized_parametrization_knitro<RegularizedPatternParametrizer   >, py::arg("rparam"), py::arg("maxIter"), py::arg("fixedVars"), py::arg("gradTol") = 1e-8);
    m.def("regularized_parametrization_knitro", &regularized_parametrization_knitro<RegularizedGenericParametrizer   >, py::arg("rparam"), py::arg("maxIter"), py::arg("fixedVars"), py::arg("gradTol") = 1e-8);
    m.def("regularized_parametrization_knitro", &regularized_parametrization_knitro<RegularizedParametrizerSVD>, py::arg("rparam"), py::arg("maxIter"), py::arg("fixedVars"), py::arg("gradTol") = 1e-8);

    m.def("regularized_parametrization_newton", &regularized_parametrization_newton<RegularizedParametrizer   >, py::arg("rparam"), py::arg("fixedVars"), py::arg("options") = NewtonOptimizerOptions());
    m.def("regularized_parametrization_newton", &regularized_parametrization_newton<RegularizedGenericParametrizer   >, py::arg("rparam"), py::arg("fixedVars"), py::arg("options") = NewtonOptimizerOptions());
    m.def("regularized_parametrization_newton", &regularized_parametrization_newton<RegularizedParametrizerSVD>, py::arg("rparam"), py::arg("fixedVars"), py::arg("options") = NewtonOptimizerOptions());
    m.def("regularized_parametrization_newton", &regularized_parametrization_newton<RegularizedGenericParametrizerSVD>, py::arg("rparam"), py::arg("fixedVars"), py::arg("options") = NewtonOptimizerOptions());
    m.def("pattern_parametrization_newton", regularized_parametrization_newton<RegularizedPatternParametrizer>, py::arg("rparam"), py::arg("fixedVars"), py::arg("options") = NewtonOptimizerOptions());

    using AngleVec = std::vector<double>;
    m.def("circularDistance",       &circularDistance      <double          >, py::arg("alpha"), py::arg("beta"));
    m.def("sumSquaredCircularDist", &sumSquaredCircularDist<double, AngleVec>, py::arg("alpha"), py::arg("angles"));
    m.def("circularMean",           &circularMean          <        AngleVec>, py::arg("angles"));

    m.def("ellipsePointParameters", [](double s, double a, double b) {
        std::vector<double> tvalues;
        double pointAreas;
        ellipsePointParameters(s, a, b, tvalues, pointAreas);
        return std::make_tuple(tvalues, pointAreas);
    }, py::arg("s"), py::arg("a"), py::arg("b"));
    
    ////////////////////////////////////////////////////////////////////////////////
    // Parametrization base class
    ////////////////////////////////////////////////////////////////////////////////
    py::class_<Parametrizer>(m, "Parametrizer")
        .def("setUV",    &Parametrizer::setUV, py::arg("uv"), py::arg("debug") = false)
        .def("setUVDebug",    &Parametrizer::setUVDebug, py::arg("uv"))
        .def("uv",       &Parametrizer::uv)
        .def("mesh",     py::overload_cast<>(&Parametrizer::mesh), py::return_value_policy::reference)
        .def("jacobian", &Parametrizer::jacobian, py::arg("tri_idx"))
        .def("B", &Parametrizer::B)
        .def("leftStretchAngles",                             &Parametrizer::leftStretchAngles)
        .def("rightStretchAngle",                            &Parametrizer::rightStretchAngle)
        .def("perVertexLeftStretchAngles",                    &Parametrizer::perVertexLeftStretchAngles,                                        py::arg("agreementThreshold") = M_PI / 8)
        .def("perVertexAlphas",                               &Parametrizer::perVertexAlphas)
        .def("perVertexBetas",                               &Parametrizer::perVertexBetas)
        // .def("perVertexQuantity",                               &Parametrizer::perVertexQuantity, py::arg("quantity"))
        .def("upsampledVertexLeftStretchAnglesAndMagnitudes", &Parametrizer::upsampledVertexLeftStretchAnglesAndMagnitudes, py::arg("nsubdiv"), py::arg("agreementThreshold") = M_PI / 8)
        .def("upsampledUV",                                   &Parametrizer::upsampledUV,                                   py::arg("nsubdiv"))
        .def("numFlips", &Parametrizer::numFlips)
        ;

    ////////////////////////////////////////////////////////////////////////////////
    // Local-global solver
    ////////////////////////////////////////////////////////////////////////////////
    py::class_<LocalGlobalParametrizer, Parametrizer>(m, "LocalGlobalParametrizer")
        .def(py::init<const std::shared_ptr<Mesh> &, const UVMap &>())
        .def_property("alphaMin", [](const LocalGlobalParametrizer &lg) { return lg.alphaMin(); },
                                  [](      LocalGlobalParametrizer &lg, Real val) { return lg.setAlphaMin(val); })
        .def_property("alphaMax", [](const LocalGlobalParametrizer &lg) { return lg.alphaMax(); },
                                  [](      LocalGlobalParametrizer &lg, Real val) { return lg.setAlphaMax(val); })
        .def("energy", &LocalGlobalParametrizer::energy)
        .def("getAlphas", &LocalGlobalParametrizer::getAlphas)
        .def("runIteration", &LocalGlobalParametrizer::runIteration)
        ;

    ////////////////////////////////////////////////////////////////////////////////
    // Local-global generic solver
    ////////////////////////////////////////////////////////////////////////////////
    py::class_<LocalGlobalGenericParametrizer, Parametrizer>(m, "LocalGlobalGenericParametrizer")
        .def(py::init<const std::shared_ptr<Mesh> &, const UVMap &>())
        .def_property("alphaMin", [](const LocalGlobalGenericParametrizer &lg) { return lg.alphaMin(); },
                                  [](      LocalGlobalGenericParametrizer &lg, Real val) { return lg.setAlphaMin(val); })
        .def_property("alphaMax", [](const LocalGlobalGenericParametrizer &lg) { return lg.alphaMax(); },
                                  [](      LocalGlobalGenericParametrizer &lg, Real val) { return lg.setAlphaMax(val); })
        .def_property("betaMin", [](const LocalGlobalGenericParametrizer &lg) { return lg.betaMin(); },
                                  [](      LocalGlobalGenericParametrizer &lg, Real val) { return lg.setBetaMin(val); })
        .def_property("betaMax", [](const LocalGlobalGenericParametrizer &lg) { return lg.betaMax(); },
                                  [](      LocalGlobalGenericParametrizer &lg, Real val) { return lg.setBetaMax(val); })
        .def("energy", &LocalGlobalGenericParametrizer::energy)
        .def("getAlphas", &LocalGlobalGenericParametrizer::getAlphas)
        .def("getBetas", &LocalGlobalGenericParametrizer::getBetas)
        .def("runIteration", &LocalGlobalGenericParametrizer::runIteration)
        .def("getLines", &LocalGlobalGenericParametrizer::getLines)
        .def("setLines", &LocalGlobalGenericParametrizer::setLines, py::arg("lines"))
        .def("projectPointInLines", &LocalGlobalGenericParametrizer::projectPointInLines, py::arg("pt"))
        .def("isPointInLines", &LocalGlobalGenericParametrizer::isPointInLines, py::arg("pt"), py::arg("eps") = 0.0)
        ;


    ////////////////////////////////////////////////////////////////////////////////
    // Laplacian regularization stencil
    ////////////////////////////////////////////////////////////////////////////////
    using DLS = DualLaplacianStencil<Mesh>;
    py::class_<DLS> pyDualLapStencil(m, "DualLaplacianStencil");

    py::enum_<DLS::Type>(pyDualLapStencil, "Type")
        .value("DualGraph",   DLS::Type::DualGraph)
        .value("DualMeshIDT", DLS::Type::DualMeshIDT)
        ;

    pyDualLapStencil
        .def(py::init<const Mesh &>(), py::arg("mesh"))
        .def_readwrite("type", &DLS::type)
        .def_readwrite("useUniformGraphWeights", &DLS::useUniformGraphWeights)
        .def("visit", [](const DLS &stencil, const size_t i, const std::function<void(size_t, size_t, Real)> &visitor) { stencil.visit(i, visitor); }, py::arg("i"), py::arg("visitor"))
        // .def("visit", &DLS::visit, py::arg("i"), py::arg("visitor"))
        ;

    ////////////////////////////////////////////////////////////////////////////////
    // Regularized parametrization energy
    ////////////////////////////////////////////////////////////////////////////////
    py::class_<RegularizedParametrizer, Parametrizer> pyRegParam(m, "RegularizedParametrizer");

    py::enum_<RegularizedParametrizer::EnergyType>(pyRegParam, "EnergyType")
        .value("Full"               , RegularizedParametrizer::EnergyType::Full)
        .value("Fitting"            , RegularizedParametrizer::EnergyType::Fitting)
        .value("PhiRegularization"  , RegularizedParametrizer::EnergyType::PhiRegularization)
        .value("AlphaRegularization", RegularizedParametrizer::EnergyType::AlphaRegularization)
        ;

    pyRegParam
        .def(py::init<LocalGlobalParametrizer &>())

        .def(   "uvOffset", &RegularizedParametrizer::   uvOffset)
        .def(    "uOffset", &RegularizedParametrizer::    uOffset)
        .def(    "vOffset", &RegularizedParametrizer::    vOffset)
        .def(  "phiOffset", &RegularizedParametrizer::  phiOffset)
        .def(  "psiOffset", &RegularizedParametrizer::  psiOffset)
        .def("alphaOffset", &RegularizedParametrizer::alphaOffset)

        .def_property("variableAlpha", [](const RegularizedParametrizer &rp)           { return rp.   variableAlpha(); },
                                       [](      RegularizedParametrizer &rp, bool val) { return rp.setVariableAlpha(val); })
        .def_property("alphaMin",      [](const RegularizedParametrizer &rp)           { return rp.   alphaMin(); },
                                       [](      RegularizedParametrizer &rp, Real val) { return rp.setAlphaMin(val); })
        .def_property("alphaMax",      [](const RegularizedParametrizer &rp)           { return rp.   alphaMax(); },
                                       [](      RegularizedParametrizer &rp, Real val) { return rp.setAlphaMax(val); })

        .def("numVars", &RegularizedParametrizer::numVars)
        .def("getVars", &RegularizedParametrizer::getVars)
        .def("setVars", &RegularizedParametrizer::setVars)

        .def("getAlphas", &RegularizedParametrizer::getAlphas)

        .def("energy", &RegularizedParametrizer::energy)

        .def("gradient", &RegularizedParametrizer::gradient, py::arg("energyType") = RegularizedParametrizer::EnergyType::Full)

        .def("hessian", [](const RegularizedParametrizer &rparam, RegularizedParametrizer::EnergyType et) { return rparam.hessian(et); }, py::arg("energyType") = RegularizedParametrizer::EnergyType::Full)
        .def("hessianSparsityPattern", &RegularizedParametrizer::hessianSparsityPattern, py::arg("val"))

        .def_property("alphaRegW", [](const RegularizedParametrizer &rp)           { return rp.   alphaRegW(); },
                                   [](      RegularizedParametrizer &rp, Real val) { return rp.setAlphaRegW(val); })
        .def_property("alphaRegP", [](const RegularizedParametrizer &rp)           { return rp.   alphaRegP(); },
                                   [](      RegularizedParametrizer &rp, Real val) { return rp.setAlphaRegP(val); })
        .def_property(  "phiRegW", [](const RegularizedParametrizer &rp)           { return rp.     phiRegW(); },
                                   [](      RegularizedParametrizer &rp, Real val) { return rp.  setPhiRegW(val); })
        .def_property(  "phiRegP", [](const RegularizedParametrizer &rp)           { return rp.     phiRegP(); },
                                   [](      RegularizedParametrizer &rp, Real val) { return rp.  setPhiRegP(val); })
        ;

    ////////////////////////////////////////////////////////////////////////////////
    // Regularized generic parametrization energy - both stretch factors as variables
    ////////////////////////////////////////////////////////////////////////////////
    py::class_<RegularizedGenericParametrizer, Parametrizer> pyRegGenParam(m, "RegularizedGenericParametrizer");

    py::enum_<RegularizedGenericParametrizer::EnergyType>(pyRegGenParam, "EnergyType")
        .value("Full"                    , RegularizedGenericParametrizer::EnergyType::Full)
        .value("Fitting"                 , RegularizedGenericParametrizer::EnergyType::Fitting)
        .value("PhiRegularization"       , RegularizedGenericParametrizer::EnergyType::PhiRegularization)
        .value("StretchRegularization"   , RegularizedGenericParametrizer::EnergyType::StretchRegularization)
        .value("DiffRegularization"      , RegularizedGenericParametrizer::EnergyType::DiffRegularization)
        .value("BendingRegularization"   , RegularizedGenericParametrizer::EnergyType::BendingRegularization)
        ;

    pyRegGenParam
        .def(py::init<LocalGlobalGenericParametrizer &>(), py::arg("lgparam"))
        .def(   "uvOffset", &RegularizedGenericParametrizer::   uvOffset)
        .def(    "uOffset", &RegularizedGenericParametrizer::    uOffset)
        .def(    "vOffset", &RegularizedGenericParametrizer::    vOffset)
        .def(  "phiOffset", &RegularizedGenericParametrizer::  phiOffset)
        .def(  "psiOffset", &RegularizedGenericParametrizer::  psiOffset)
        .def("alphaOffset", &RegularizedGenericParametrizer::alphaOffset)
        .def("betaOffset", &RegularizedGenericParametrizer::betaOffset)
        .def("stretchOffset", &RegularizedGenericParametrizer::stretchOffset)

        .def("jacobian_M", &RegularizedGenericParametrizer::jacobian_M, py::arg("tri_idx"))

        .def_property("variableStretch", [](const RegularizedGenericParametrizer &rgp)           { return rgp.   variableStretch(); },
                                         [](      RegularizedGenericParametrizer &rgp, bool val) { return rgp.setVariableStretch(val); })
        .def_property("useBarrier",      [](const RegularizedGenericParametrizer &rgp)           { return rgp.   useBarrier(); },
                                         [](      RegularizedGenericParametrizer &rgp, bool val) { return rgp.setUseBarrier(val); })
        .def_property("alphaMin",        [](const RegularizedGenericParametrizer &rgp)           { return rgp.   alphaMin(); },
                                         [](      RegularizedGenericParametrizer &rgp, Real val) { return rgp.setAlphaMin(val); })
        .def_property("alphaMax",        [](const RegularizedGenericParametrizer &rgp)           { return rgp.   alphaMax(); },
                                         [](      RegularizedGenericParametrizer &rgp, Real val) { return rgp.setAlphaMax(val); })
        .def_property("betaMin",         [](const RegularizedGenericParametrizer &rgp)           { return rgp.   betaMin(); },
                                         [](      RegularizedGenericParametrizer &rgp, Real val) { return rgp.setBetaMin(val); })
        .def_property("betaMax",         [](const RegularizedGenericParametrizer &rgp)           { return rgp.   betaMax(); },
                                         [](      RegularizedGenericParametrizer &rgp, Real val) { return rgp.setBetaMax(val); })


        .def("numVars", &RegularizedGenericParametrizer::numVars)
        .def("getVars", &RegularizedGenericParametrizer::getVars)
        .def("setVars", &RegularizedGenericParametrizer::setVars)

        .def("setLines", &RegularizedGenericParametrizer::setLines, py::arg("lines"))

        .def("getAlphas", &RegularizedGenericParametrizer::getAlphas)
        .def("getBetas", &RegularizedGenericParametrizer::getBetas)
        .def("getBarriers", &RegularizedGenericParametrizer::getBarriers)

        .def("kappa", &RegularizedGenericParametrizer::kappa, py::arg("tri"), py::arg("i"))
        .def("kappaAngle", &RegularizedGenericParametrizer::kappaAngle, py::arg("tri"), py::arg("i"))
        .def("getKappa", &RegularizedGenericParametrizer::getKappa, py::arg("i"))
        .def("getKappaAngle", &RegularizedGenericParametrizer::getKappaAngle)
        .def("curvature3d", &RegularizedGenericParametrizer::curvature3d, py::arg("tri"), py::arg("i"))

        // .def("energy", &RegularizedGenericParametrizer::energy)
        .def("energy", py::overload_cast<RegularizedGenericParametrizer::EnergyType>(&RegularizedGenericParametrizer::energy, py::const_), py::arg("energyType") = RegularizedGenericParametrizer::EnergyType::Full)

        .def("gradient", &RegularizedGenericParametrizer::gradient, py::arg("energyType") = RegularizedGenericParametrizer::EnergyType::Full)

        .def("hessian", [](const RegularizedGenericParametrizer &rparam, RegularizedGenericParametrizer::EnergyType et) { return rparam.hessian(et); }, py::arg("energyType") = RegularizedGenericParametrizer::EnergyType::Full)
        .def("hessianSparsityPattern", &RegularizedGenericParametrizer::hessianSparsityPattern, py::arg("val"))

        .def_property("stretchRegW", [](const RegularizedGenericParametrizer &rgp)           { return rgp.   stretchRegW(); },
                                     [](      RegularizedGenericParametrizer &rgp, Real val) { return rgp.setStretchRegW(val); })
        .def_property("stretchRegP", [](const RegularizedGenericParametrizer &rgp)           { return rgp.   stretchRegP(); },
                                     [](      RegularizedGenericParametrizer &rgp, Real val) { return rgp.setStretchRegP(val); })
        .def_property(  "phiRegW",   [](const RegularizedGenericParametrizer &rgp)           { return rgp.     phiRegW(); },
                                     [](      RegularizedGenericParametrizer &rgp, Real val) { return rgp.  setPhiRegW(val); })
        .def_property(  "phiRegP",   [](const RegularizedGenericParametrizer &rgp)           { return rgp.     phiRegP(); },
                                     [](      RegularizedGenericParametrizer &rgp, Real val) { return rgp.  setPhiRegP(val); })
        .def_property("diffRegW",    [](const RegularizedGenericParametrizer &rgp)           { return rgp.   diffRegW(); },
                                     [](      RegularizedGenericParametrizer &rgp, Real val) { return rgp.setDiffRegW(val); })
        .def_property("barrierA",    [](const RegularizedGenericParametrizer &rgp)           { return rgp.   barrierA(); },
                                     [](      RegularizedGenericParametrizer &rgp, Real val) { return rgp.setBarrierA(val); })
        .def_property("barrierB",    [](const RegularizedGenericParametrizer &rgp)           { return rgp.   barrierB(); },
                                     [](      RegularizedGenericParametrizer &rgp, Real val) { return rgp.setBarrierB(val); })                                   
        .def_readonly("dualLaplacianStencil",         &RegularizedGenericParametrizer::dualLaplacianStencil)                                   
        ;


    ///////////////////////////////////////////////////////////////////////////////////
    // Regularized parametrization energy - using pattern parameters as variables
    ///////////////////////////////////////////////////////////////////////////////////
    using RPP = RegularizedPatternParametrizer;

    py::class_<RPP, RegularizedGenericParametrizer> pyRegPatParam(m, "RegularizedPatternParametrizer");


    py::enum_<RPP::PatternEnergyType>(pyRegPatParam, "PatternEnergyType")
        .value("Full",                     RPP::PatternEnergyType::Full)
        .value("RGP",                      RPP::PatternEnergyType::RGP)
        .value("Bending",                  RPP::PatternEnergyType::Bending)
        .value("PatternRegularization", RPP::PatternEnergyType::PatternRegularization)
        .value("PatternBoundConstraint", RPP::PatternEnergyType::PatternBoundConstraint)
        .value("DEBUG_Fitting", RPP::PatternEnergyType::DEBUG_Fitting)
        .value("DEBUG_PhiRegularization", RPP::PatternEnergyType::DEBUG_PhiRegularization)
        .value("DEBUG_StretchRegularization", RPP::PatternEnergyType::DEBUG_StretchRegularization)
        ;

    pyRegPatParam
        .def(py::init<LocalGlobalGenericParametrizer &, RPP::MatInfoFunction, Eigen::VectorXd, size_t>(), py::arg("lgparam"), py::arg("newMatInfo"), py::arg("default_pattern_params"), py::arg("input_num_pattern_vars"))

        .def("numVars", &RPP::numVars)
        .def("getVars", &RPP::getVars)
        .def("setVars", &RPP::setVars)
        .def("getPatternParams", [](const RPP &rparam) { return rparam.getPatternParams(); })
        .def("patternOffset", &RPP::patternOffset)

        .def_property("patternRegW", [](const RPP &rparam)           { return rparam.   patternRegW(); },
                                     [](      RPP &rparam, Real val) { return rparam.setPatternRegW(val); })
        .def_property("patternRegP", [](const RPP &rparam)           { return rparam.   patternRegP(); },
                                     [](      RPP &rparam, Real val) { return rparam.setPatternRegP(val); })
        .def_property("patternParamBounds", [](const RPP &rparam)           { return rparam.   patternParamBounds(); },
                                             [](      RPP &rparam, Eigen::MatrixXd val) { return rparam.setPatternParamBounds(val); })
        .def_property("patternParamNormalizationFactors", [](const RPP &rparam)           { return rparam.   patternParamNormalizationFactors(); },
                                                          [](      RPP &rparam, Eigen::VectorXd val) { return rparam.setPatternParamNormalizationFactors(val); })
        .def_property("bendRegW", [](const RPP &rparam)           { return rparam.   bendRegW(); },
                                  [](      RPP &rparam, Real val) { return rparam.setBendRegW(val); })
        .def_property("phiRegW",  [](const RPP &rparam)           { return rparam.    phiRegW(); },
                                  [](      RPP &rparam, Real val) { return rparam. setPhiRegW(val); })
        
        .def("energy", py::overload_cast<RPP::PatternEnergyType>(&RPP::energy, py::const_), py::arg("energyType") = RPP::PatternEnergyType::Full)

        .def("gradient", &RPP::gradient, py::arg("energyType") = RPP::PatternEnergyType::Full)

        .def("hessian", [](const RPP &rparam, RPP::PatternEnergyType et) { return rparam.hessian(et); }, py::arg("energyType") = RPP::PatternEnergyType::Full)
        .def("hessianSparsityPattern", &RPP::hessianSparsityPattern, py::arg("val"))    
        .def("perVertexPatternParams", &RPP::perVertexPatternParams)
        .def("upsampledVertexLeftStretchAnglesAndPatternParameters", &RPP::upsampledVertexLeftStretchAnglesAndPatternParameters, py::arg("nsubdiv"), py::arg("agreementThreshold") = M_PI / 8)
        .def("getMatInfoArgs", &RPP::getMatInfoArgs)
        .def("get_stretch_angle_offset_from_pattern_params", &RPP::get_stretch_angle_offset_from_pattern_params, py::arg("query_pattern_params"))
        .def_property("useBarrier", [](const RPP &rparam)           { return rparam.getConstraintBarrier(); },
                                    [](      RPP &rparam, bool val) { return rparam.setConstraintBarrier(val); })
        ;


    ////////////////////////////////////////////////////////////////////////////////
    // Regularized parametrization energy, SVD version
    ////////////////////////////////////////////////////////////////////////////////
    py::class_<RegularizedParametrizerSVD, Parametrizer> pyRegParamSVD(m, "RegularizedParametrizerSVD");

    py::enum_<RegularizedParametrizerSVD::EnergyType>(pyRegParamSVD, "EnergyType")
        .value("Full",                  RegularizedParametrizerSVD::EnergyType::Full)
        .value("Fitting",               RegularizedParametrizerSVD::EnergyType::Fitting)
        .value("PhiRegularization",     RegularizedParametrizerSVD::EnergyType::PhiRegularization)
        .value("AlphaRegularization",   RegularizedParametrizerSVD::EnergyType::AlphaRegularization)
        .value("BendingRegularization", RegularizedParametrizerSVD::EnergyType::BendingRegularization)
        ;

    pyRegParamSVD
        .def(py::init<const std::shared_ptr<Mesh> &, const UVMap &, Real, Real, bool>(), py::arg("mesh"), py::arg("uv"), py::arg("alphaMin") = 1.0, py::arg("alphaMax") = 1.0, py::arg("transformForRigidMotionConstraint") = true)
        .def(py::init<LocalGlobalParametrizer &>(), py::arg("lgparam"))

        .def(   "uvOffset", &RegularizedParametrizerSVD::uvOffset)
        .def(    "uOffset", &RegularizedParametrizerSVD:: uOffset)
        .def(    "vOffset", &RegularizedParametrizerSVD:: vOffset)

        .def_property_readonly("rigidMotionPinVars", [&](const RegularizedParametrizerSVD &rsvd) { return rsvd.rigidMotionPinVars(); })

        .def_property("alphaMin",      [](const RegularizedParametrizerSVD &rsvd)           { return rsvd.   alphaMin(); },
                                       [](      RegularizedParametrizerSVD &rsvd, Real val) { return rsvd.setAlphaMin(val); })
        .def_property("alphaMax",      [](const RegularizedParametrizerSVD &rsvd)           { return rsvd.   alphaMax(); },
                                       [](      RegularizedParametrizerSVD &rsvd, Real val) { return rsvd.setAlphaMax(val); })

        .def("numVars", &RegularizedParametrizerSVD::numVars)
        .def("getVars", &RegularizedParametrizerSVD::getVars)
        .def("setVars", &RegularizedParametrizerSVD::setVars)

        .def("getAlphas",            &RegularizedParametrizerSVD::getAlphas)
        .def("getMinSingularValues", &RegularizedParametrizerSVD::getMinSingularValues)

        .def("tubeDirections", &RegularizedParametrizerSVD::tubeDirections)

        .def("curvature3d", &RegularizedParametrizerSVD::curvature3d, py::arg("tri"), py::arg("i"))

        .def("energy", py::overload_cast<RegularizedParametrizerSVD::EnergyType>(&RegularizedParametrizerSVD::energy, py::const_), py::arg("energyType") = RegularizedParametrizerSVD::EnergyType::Full)

        .def("gradient", &RegularizedParametrizerSVD::gradient, py::arg("energyType") = RegularizedParametrizerSVD::EnergyType::Full)

        .def("hessian", [](const RegularizedParametrizerSVD &rparam, RegularizedParametrizerSVD::EnergyType et, bool projectionMask) { return rparam.hessian(et, projectionMask); }, py::arg("energyType") = RegularizedParametrizerSVD::EnergyType::Full, py::arg("projectionMask") = false)
        .def("hessianSparsityPattern", &RegularizedParametrizerSVD::hessianSparsityPattern, py::arg("val"))

        .def_property("alphaRegW", [](const RegularizedParametrizerSVD &rsvd)           { return rsvd.   alphaRegW(); },
                                   [](      RegularizedParametrizerSVD &rsvd, Real val) { return rsvd.setAlphaRegW(val); })
        .def_property("alphaRegP", [](const RegularizedParametrizerSVD &rsvd)           { return rsvd.   alphaRegP(); },
                                   [](      RegularizedParametrizerSVD &rsvd, Real val) { return rsvd.setAlphaRegP(val); })
        .def_property(  "phiRegW", [](const RegularizedParametrizerSVD &rsvd)           { return rsvd.     phiRegW(); },
                                   [](      RegularizedParametrizerSVD &rsvd, Real val) { return rsvd.  setPhiRegW(val); })
        .def_property(  "phiRegP", [](const RegularizedParametrizerSVD &rsvd)           { return rsvd.     phiRegP(); },
                                   [](      RegularizedParametrizerSVD &rsvd, Real val) { return rsvd.  setPhiRegP(val); })
        .def_property( "bendRegW", [](const RegularizedParametrizerSVD &rsvd)           { return rsvd.    bendRegW(); },
                                   [](      RegularizedParametrizerSVD &rsvd, Real val) { return rsvd. setBendRegW(val); })
        .def_property( "stretchDeviationP", [](const RegularizedParametrizerSVD &rsvd)           { return rsvd.    stretchDeviationP(); },
                                            [](      RegularizedParametrizerSVD &rsvd, Real val) { return rsvd. setStretchDeviationP(val); })
        .def_readonly("dualLaplacianStencil",         &RegularizedParametrizerSVD::dualLaplacianStencil)
        .def_readwrite("scaleInvariantFittingEnergy", &RegularizedParametrizerSVD::scaleInvariantFittingEnergy)

        .def("setAlphas", &RegularizedParametrizerSVD::setAlphas, py::arg("newAlphas")) // for debugging/analysis only!
        ;

    ////////////////////////////////////////////////////////////////////////////////
    // Regularized generic parametrization energy, SVD version
    ////////////////////////////////////////////////////////////////////////////////
    py::class_<RegularizedGenericParametrizerSVD, Parametrizer> pyRegGenParamSVD(m, "RegularizedGenericParametrizerSVD");

    py::enum_<RegularizedGenericParametrizerSVD::EnergyType>(pyRegGenParamSVD, "EnergyType")
        .value("Full",                  RegularizedGenericParametrizerSVD::EnergyType::Full)
        .value("Fitting",               RegularizedGenericParametrizerSVD::EnergyType::Fitting)
        .value("PhiRegularization",     RegularizedGenericParametrizerSVD::EnergyType::PhiRegularization)
        .value("AlphaRegularization",   RegularizedGenericParametrizerSVD::EnergyType::AlphaRegularization)
        .value("BetaRegularization",   RegularizedGenericParametrizerSVD::EnergyType::BetaRegularization)
        .value("BendingRegularization", RegularizedGenericParametrizerSVD::EnergyType::BendingRegularization)
        ;

    pyRegGenParamSVD
        .def(py::init<const std::shared_ptr<Mesh> &, const UVMap &, Real, Real, Real, Real, bool>(), py::arg("mesh"), py::arg("uv"), py::arg("alphaMin") = 1.0, py::arg("alphaMax") = 1.0, py::arg("betaMin") = 1.0, py::arg("betaMax") = 1.0, py::arg("transformForRigidMotionConstraint") = true)
        .def(py::init<LocalGlobalGenericParametrizer &>(), py::arg("lgparam"))

        .def(   "uvOffset", &RegularizedGenericParametrizerSVD::uvOffset)
        .def(    "uOffset", &RegularizedGenericParametrizerSVD:: uOffset)
        .def(    "vOffset", &RegularizedGenericParametrizerSVD:: vOffset)

        .def_property_readonly("rigidMotionPinVars", [&](const RegularizedGenericParametrizerSVD &rsvd) { return rsvd.rigidMotionPinVars(); })

        .def_property("alphaMin",      [](const RegularizedGenericParametrizerSVD &rsvd)           { return rsvd.   alphaMin(); },
                                       [](      RegularizedGenericParametrizerSVD &rsvd, Real val) { return rsvd.setAlphaMin(val); })
        .def_property("alphaMax",      [](const RegularizedGenericParametrizerSVD &rsvd)           { return rsvd.   alphaMax(); },
                                       [](      RegularizedGenericParametrizerSVD &rsvd, Real val) { return rsvd.setAlphaMax(val); })
        .def_property("betaMin",      [](const RegularizedGenericParametrizerSVD &rsvd)           { return rsvd.   betaMin(); },
                                       [](      RegularizedGenericParametrizerSVD &rsvd, Real val) { return rsvd.setBetaMin(val); })
        .def_property("betaMax",      [](const RegularizedGenericParametrizerSVD &rsvd)           { return rsvd.   betaMax(); },
                                       [](      RegularizedGenericParametrizerSVD &rsvd, Real val) { return rsvd.setBetaMax(val); })

        .def("numVars", &RegularizedGenericParametrizerSVD::numVars)
        .def("getVars", &RegularizedGenericParametrizerSVD::getVars)
        .def("setVars", &RegularizedGenericParametrizerSVD::setVars)

        .def("getAlphas",            &RegularizedGenericParametrizerSVD::getAlphas)
        .def("getBetas",            &RegularizedGenericParametrizerSVD::getBetas)
        .def("getMinSingularValues", &RegularizedGenericParametrizerSVD::getMinSingularValues)

        .def("tubeDirections", &RegularizedGenericParametrizerSVD::tubeDirections)

        .def("curvature3d", &RegularizedGenericParametrizerSVD::curvature3d, py::arg("tri"), py::arg("i"))

        .def("energy", py::overload_cast<RegularizedGenericParametrizerSVD::EnergyType>(&RegularizedGenericParametrizerSVD::energy, py::const_), py::arg("energyType") = RegularizedGenericParametrizerSVD::EnergyType::Full)

        .def("gradient", &RegularizedGenericParametrizerSVD::gradient, py::arg("energyType") = RegularizedGenericParametrizerSVD::EnergyType::Full)

        .def("hessian", [](const RegularizedGenericParametrizerSVD &rparam, RegularizedGenericParametrizerSVD::EnergyType et, bool projectionMask) { return rparam.hessian(et, projectionMask); }, py::arg("energyType") = RegularizedGenericParametrizerSVD::EnergyType::Full, py::arg("projectionMask") = false)
        .def("hessianSparsityPattern", &RegularizedGenericParametrizerSVD::hessianSparsityPattern, py::arg("val"))

        .def_property("stretchRegW", [](const RegularizedGenericParametrizerSVD &rsvd)           { return rsvd.   stretchRegW(); },
                                   [](      RegularizedGenericParametrizerSVD &rsvd, Real val) { return rsvd.setStretchRegW(val); })
        .def_property("stretchRegP", [](const RegularizedGenericParametrizerSVD &rsvd)           { return rsvd.   stretchRegP(); },
                                   [](      RegularizedGenericParametrizerSVD &rsvd, Real val) { return rsvd.setStretchRegP(val); })
        .def_property(  "phiRegW", [](const RegularizedGenericParametrizerSVD &rsvd)           { return rsvd.     phiRegW(); },
                                   [](      RegularizedGenericParametrizerSVD &rsvd, Real val) { return rsvd.  setPhiRegW(val); })
        .def_property(  "phiRegP", [](const RegularizedGenericParametrizerSVD &rsvd)           { return rsvd.     phiRegP(); },
                                   [](      RegularizedGenericParametrizerSVD &rsvd, Real val) { return rsvd.  setPhiRegP(val); })
        .def_property( "bendRegW", [](const RegularizedGenericParametrizerSVD &rsvd)           { return rsvd.    bendRegW(); },
                                   [](      RegularizedGenericParametrizerSVD &rsvd, Real val) { return rsvd. setBendRegW(val); })
        .def_property( "stretchDeviationP", [](const RegularizedGenericParametrizerSVD &rsvd)           { return rsvd.    stretchDeviationP(); },
                                            [](      RegularizedGenericParametrizerSVD &rsvd, Real val) { return rsvd. setStretchDeviationP(val); })
        .def_readonly("dualLaplacianStencil",         &RegularizedGenericParametrizerSVD::dualLaplacianStencil)
        .def_readwrite("scaleInvariantFittingEnergy", &RegularizedGenericParametrizerSVD::scaleInvariantFittingEnergy)

        .def("setAlphas", &RegularizedGenericParametrizerSVD::setAlphas, py::arg("newAlphas")) // for debugging/analysis only!
        .def("setBetas", &RegularizedGenericParametrizerSVD::setBetas, py::arg("newBetas")) // for debugging/analysis only!
        ;

    ////////////////////////////////////////////////////////////////////////////////
    // Enable output redirection from Python side
    ////////////////////////////////////////////////////////////////////////////////
    py::add_ostream_redirect(m, "ostream_redirect");
}
