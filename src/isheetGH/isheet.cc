#include <MeshFEM/GlobalBenchmark.hh>
#include <MeshFEM/MeshIO.hh>
#include <MeshFEM/Future.hh>
#include <MeshFEM/newton_optimizer/newton_optimizer.hh>

#include "InflatableSheet.hh"
#include "inflation_newton.hh"
#include "InflatablePeriodicUnit.hh"
#include "InflatableMidSurfacePeriodicUnit.hh"
#include "periodic_stiffness_analysis.hh"
#include "TargetAttractedInflation.hh"
#include "MultilayerInflatable.hh"

extern "C"
{
#include "isheet.h"
}

namespace InflatableSheetGH {
    
    std::shared_ptr<InflatableSheet::Mesh> buildMesh(int numVertices, int numTrias, double *inCoords, int *inTrias){
        // Vertices
        std::vector<MeshIO::IOVertex> vertices;
        for (int i = 0; i < numVertices; i++)
        {
            vertices.emplace_back(inCoords[3 * i], inCoords[3 * i + 1], inCoords[3 * i + 2]);
        }

        // Triangles
        std::vector<MeshIO::IOElement> triangles;
        for (int i = 0; i < numTrias; i++)
        {
            triangles.emplace_back(inTrias[3 * i], inTrias[3 * i + 1], inTrias[3 * i + 2]);
        }

        return std::make_shared<InflatableSheet::Mesh>(triangles, vertices);
    }

    void getVisualizationMesh(InflatableSheet::Mesh mesh, double **outCoords, int **outElements, size_t *numCoords, size_t *numElements){
        const size_t nv = mesh.numVertices();
        const size_t nt = mesh.numElements();

        // Parsing Vertices
        *numCoords = nv * 3;
        auto sizeCoords = (*numCoords) * sizeof(double);
        *outCoords = static_cast<double *>(malloc(sizeCoords));
        std::vector<double> coords;
        for (const auto vi : mesh.vertices()){
            auto node = vi.node()->p.cast<Real>();
            coords.push_back(node[0]);
            coords.push_back(node[1]);
            coords.push_back(node[2]);
        }
        std::memcpy(*outCoords, coords.data(), sizeCoords);

        // Parsing Triangles
        *numElements = nt * 3;
        auto sizeElements = (*numElements) * sizeof(int);
        *outElements = static_cast<int *>(malloc(sizeElements));
        std::vector<int> trias;
        for (const auto tri : mesh.elements()) {
            trias.push_back(tri.vertex(0).index());
            trias.push_back(tri.vertex(1).index());
            trias.push_back(tri.vertex(2).index());
        }
        std::memcpy(*outElements, trias.data(), sizeElements);  
    }

    void getConvergenceReport(ConvergenceReport report, int numIterations, double **outReport){
        std::vector<double> flatReport;
        flatReport.push_back(report.success);
        flatReport.push_back(report.backtracking_failure);
        flatReport.insert(flatReport.end(), report.energy.begin(), report.energy.end());
        flatReport.insert(flatReport.end(), report.freeGradientNorm.begin(), report.freeGradientNorm.end());
        flatReport.insert(flatReport.end(), report.stepLength.begin(), report.stepLength.end());
        flatReport.insert(flatReport.end(), report.indefinite.begin(), report.indefinite.end());

        auto sizeReport = (numIterations * 4 + 2) * sizeof(double);
        *outReport = static_cast<double *>(malloc(sizeReport));
        std::memcpy(*outReport, flatReport.data(), sizeReport);
    }

    ///////////////////////////////////////////////////////////////////////////
    // Multi layer inflatable Sheet
    ///////////////////////////////////////////////////////////////////////////
    ISHEET_API MultilayerInflatable *isheet_multilayerInflatableSheet_build(int numVertices, int numTrias, Real *inCoords, int *inTrias, int numSheets, Real *inPressures, int *inReducedVarIdxForVertexOnSheet, const char **errorMessage){    
        try{
            const auto m = buildMesh(numVertices, numTrias, inCoords, inTrias);

            std::vector<Real> pressures(numSheets-1);
            Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> reducedVarIdxForVertexOnSheet(numVertices, numSheets);
            for (int i = 0; i < numSheets; i++){
                if(i<numSheets-1) pressures[i] = inPressures[i];
                for (int j = 0; j < numVertices; j++){
                    reducedVarIdxForVertexOnSheet(j, i) = inReducedVarIdxForVertexOnSheet[i * numVertices + j];
                }
            }

            return new MultilayerInflatable(m, numSheets,  pressures, reducedVarIdxForVertexOnSheet);
        }
        catch (const std::runtime_error &error)
        {
            *errorMessage = error.what();
            return nullptr;
        }
        catch (const std::out_of_range &error)
        {
            *errorMessage = error.what();
            return nullptr;
        }
        catch (...)
        {
            *errorMessage = "Unknown Error from the Unmanaged Code.";
            return nullptr;
        }
    }

    // Set methods
    ISHEET_API void isheet_multilayerInflatableSheet_setPressure(MultilayerInflatable *sheet, double *inPressures, size_t numPressures){
        std::vector<Real> pressures(numPressures);
        for (size_t i = 0; i < numPressures; i++) pressures[i] = inPressures[i];
        sheet->setPressure(pressures);
    }

    ISHEET_API void isheet_multilayerInflatableSheet_setUseTensionFieldEnergy(MultilayerInflatable *sheet, int useTensionFieldEnergy){
        sheet->setUseTensionFieldEnergy(useTensionFieldEnergy);
    }

    ISHEET_API void isheet_multilayerInflatableSheet_setUseHessianProjectedEnergy(MultilayerInflatable *sheet, int useHessianProjectedEnergy){
        sheet->setUseHessianProjectedEnergy(useHessianProjectedEnergy);
    }

    ISHEET_API void isheet_multilayerInflatableSheet_setVars(MultilayerInflatable *sheet, double *inVars, size_t numVars){
        Eigen::VectorXd vars = Eigen::Map<Eigen::VectorXd>(inVars, numVars, 1);
        sheet->setVars(vars);
    }
    
    ISHEET_API void isheet_multilayerInflatableSheet_setGravity(MultilayerInflatable *sheet, double* inVector){
        const InflatableSheet::V3d gravity = InflatableSheet::V3d(inVector[0], inVector[1], inVector[2]);
        sheet->setGravity(gravity);
    }

    ISHEET_API void isheet_multilayerInflatableSheet_setMassDensity(MultilayerInflatable *sheet, double rho){
        sheet->setRho(rho);
    }

    ISHEET_API void isheet_multilayerInflatableSheet_setThickness(MultilayerInflatable *sheet, double thickness){
        sheet->setThickness(thickness);
    }

    ISHEET_API void isheet_multilayerInflatableSheet_setYoungModulus(MultilayerInflatable *sheet, double youngModulus){
        sheet->setYoungModulus(youngModulus);
    }

    ISHEET_API void isheet_multilayerInflatableSheet_setReferenceVolume(MultilayerInflatable *sheet, double *inVolumes, size_t numVolumes)
    {
        std::vector<Real> volumes(numVolumes);
        for (size_t i = 0; i < numVolumes; i++) volumes[i] = inVolumes[i];
        sheet->setReferenceVolume(volumes);
    }

    // Get methods
    ISHEET_API void isheet_multilayerInflatableSheet_getMeshVisualization(MultilayerInflatable *sheet, double **outCoords, int **outElements, size_t *numCoords, size_t *numElements)
    {
        const auto m = *sheet->visualizationMesh();

        const size_t nv = m.numVertices();
        const size_t nt = m.numElements();

        //////////////////////////

        *numCoords = nv * 3;
        auto sizeCoords = (*numCoords) * sizeof(double);
        *outCoords = static_cast<double *>(malloc(sizeCoords));
        std::vector<double> coords;
        for (const auto vi : m.vertices()){
            auto node = vi.node()->p.cast<Real>();
            coords.push_back(node[0]);
            coords.push_back(node[1]);
            coords.push_back(node[2]);
        }
        std::memcpy(*outCoords, coords.data(), sizeCoords);

        *numElements = nt * 3;
        auto sizeElements = (*numElements) * sizeof(int);
        *outElements = static_cast<int *>(malloc(sizeElements));
        std::vector<int> trias;
        for (const auto tri : m.elements()) {
            trias.push_back(tri.vertex(0).index());
            trias.push_back(tri.vertex(1).index());
            trias.push_back(tri.vertex(2).index());
        }
        std::memcpy(*outElements, trias.data(), sizeElements);
    }

    ISHEET_API double isheet_multilayerInflatableSheet_getMassDensity(MultilayerInflatable *sheet){
        return sheet->getRho();
    }

    ISHEET_API void isheet_multilayerInflatableSheet_getVolumes(MultilayerInflatable *sheet,  double **outVolumes, size_t *numVolumes)
    {
        const auto volumes = sheet->volume();
        *numVolumes = volumes.size();
        auto sizeVolumes = (*numVolumes) * sizeof(double);
        *outVolumes = static_cast<double *>(malloc(sizeVolumes));
        std::memcpy(*outVolumes, volumes.data(), sizeVolumes);
    }

    ISHEET_API double isheet_multilayerInflatableSheet_getYoungModulus(MultilayerInflatable *sheet){
        return sheet->getYoungModulus();
    }

    ISHEET_API double isheet_multilayerInflatableSheet_getThickness(MultilayerInflatable *sheet){
        return sheet->getThickness();
    }

    ISHEET_API void isheet_multilayerInflatableSheet_getCenterFixedVars(MultilayerInflatable *sheet, int **outCenterFixedVars, size_t *numCenterFixedVars)
    {
        const auto fixedVxIdx = sheet->center_non_fused_vx_idx();
        const size_t startVarIdx = fixedVxIdx * 3;
        const size_t endVarIdx = 3 + fixedVxIdx * 3;
        *numCenterFixedVars = endVarIdx - startVarIdx;

        auto sizeCenterFixedVars = (*numCenterFixedVars) * sizeof(int);
        *outCenterFixedVars = static_cast<int *>(malloc(sizeCenterFixedVars));
        std::vector<int> centerFixedVars;
        for (size_t i = startVarIdx; i < endVarIdx; i++) centerFixedVars.push_back(i);
        std::memcpy(*outCenterFixedVars, centerFixedVars.data(), sizeCenterFixedVars);
    }

    ISHEET_API void isheet_multilayerInflatableSheet_getVertexVars(MultilayerInflatable *sheet, int sheetIndex, int vertexIndex, int **outVars, size_t *numVars){
        const auto fixedVxIdx = sheet->varIdx(sheetIndex, vertexIndex,0)/3;

        const size_t startVarIdx = fixedVxIdx * 3;
        const size_t endVarIdx = 3 + fixedVxIdx * 3;
        *numVars = endVarIdx - startVarIdx;

        auto sizeFixedVars = (*numVars) * sizeof(int);
        *outVars = static_cast<int *>(malloc(sizeFixedVars));
        std::vector<int> fixedVars;
        for (size_t i = startVarIdx; i < endVarIdx; i++) fixedVars.push_back(i);
        std::memcpy(*outVars, fixedVars.data(), sizeFixedVars);
    }

    ISHEET_API int isheet_multilayerInflatableSheet_getCenterNonFusedVertexIdx(MultilayerInflatable *sheet)
    {
        return sheet->center_non_fused_vx_idx();
    }

    ISHEET_API int isheet_multilayerInflatableSheet_getNumVars(MultilayerInflatable *sheet)
    {
        return sheet->numVars();
    }

    ISHEET_API void isheet_multilayerInflatableSheet_getVars(MultilayerInflatable *sheet, double **outVars, size_t *numVars)
    {
        const auto vars = sheet->getVars();
        *numVars = vars.size();
        auto sizeVars = (*numVars) * sizeof(double);
        *outVars = static_cast<double *>(malloc(sizeVars));
        std::memcpy(*outVars, vars.data(), sizeVars);
    }

    ISHEET_API void isheet_multilayerInflatableSheet_getStrains(MultilayerInflatable *sheet, double **outStrains, size_t *numStrains)
    {
        const auto ted = sheet->triEnergyDensities();
        *numStrains = ted.size() * 2;
        auto sizeStrains = (*numStrains) * sizeof(double);
        *outStrains = static_cast<double *>(malloc(sizeStrains));
        std::vector<double> strains;
        for (const auto data : ted){
            const auto vec = data.eigSensitivities().Lambda();
            strains.push_back(sqrt(vec(0))-1);
            strains.push_back(sqrt(vec(1))-1);
        }
        std::memcpy(*outStrains, strains.data(), sizeStrains);
    }

    ISHEET_API double isheet_multilayerInflatableSheet_getEnergy(MultilayerInflatable *sheet){
        return sheet->energy();
    }

    ISHEET_API double isheet_multilayerInflatableSheet_getGradientNorm(MultilayerInflatable *sheet){
        return sheet->gradient().norm();
    }

    ISHEET_API void isheet_multilayerInflatableSheet_getGradient(MultilayerInflatable *sheet, double **outVars, size_t *numVars){
        const auto grad = sheet->gradient();
        *numVars = grad.size();
        auto sizeGrad = (*numVars) * sizeof(double);
        *outVars = static_cast<double *>(malloc(sizeGrad));
        std::memcpy(*outVars, grad.data(), sizeGrad);
    }

    ISHEET_API void isheet_multilayerInflatableSheet_getGravity(MultilayerInflatable *sheet, double **outVars, size_t *numVars){
        const auto gravity = sheet->getGravity();
        *numVars = gravity.size();
        auto sizeGravity = (*numVars) * sizeof(double);
        *outVars = static_cast<double *>(malloc(sizeGravity));
        std::memcpy(*outVars, gravity.data(), sizeGravity);
    }

    ISHEET_API void isheet_multilayerInflatableSheet_getReferenceVolume(MultilayerInflatable *sheet, double **outVolumes, size_t *numVolumes){
        const auto volumes = sheet->referenceVolume();
        *numVolumes = volumes.size();
        auto sizeVolumes = (*numVolumes) * sizeof(double);
        *outVolumes = static_cast<double *>(malloc(sizeVolumes));
        std::memcpy(*outVolumes, volumes.data(), sizeVolumes);
    }

    ISHEET_API void isheet_multilayerInflatableSheet_getPressure(MultilayerInflatable *sheet, double **outPressures, size_t *numPressures){
        const auto pressures = sheet->getPressure();
        *numPressures = pressures.size();
        auto sizePressures = (*numPressures) * sizeof(double);
        *outPressures = static_cast<double *>(malloc(sizePressures));
        std::memcpy(*outPressures, pressures.data(), sizePressures);
    }
    // Solver
    ISHEET_API int isheet_multilayerInflatableSheet_newtonSolver(MultilayerInflatable *sheet, int numSupports, int *supports, int numIterations, double gradTol, int writeReport, double **outReport, const char **errorMessage)
    {
        try
        {
            std::vector<size_t> fixedVars;
            for (int i = 0; i < numSupports; i++) fixedVars.push_back(supports[i]);

            // Set real options: [0] numIterations, [1] gradTol, [2] hessianShift
            NewtonOptimizerOptions options;
            options.niter = numIterations;
            options.gradTol = gradTol;

            const auto report = inflation_newton<MultilayerInflatable>(*sheet, fixedVars, options);

            if (writeReport) getConvergenceReport(report, numIterations, outReport);

            *errorMessage = "";

            return report.success;
        }
        catch (const std::runtime_error &error)
        {
            *errorMessage = error.what();
            return -1;
        }
        catch (const std::out_of_range &error)
        {
            *errorMessage = error.what();
            return -1;
        }
        catch (...)
        {
            *errorMessage = "Unknown Error from the Unmanaged Code.";
            return -1;
        }
    }
    
    ///////////////////////////////////////////////////////////////////////////
    // Inflatable Sheet
    ///////////////////////////////////////////////////////////////////////////
    ISHEET_API InflatableSheet *isheet_inflatableSheet_build(int numVertices, int numTrias, Real *inCoords, int *inTrias, int *inFusedVertices)
    {
        auto m = buildMesh(numVertices, numTrias, inCoords, inTrias);

        std::vector<bool> fusedVertices(numVertices);
        for (int i = 0; i < numVertices; i++) fusedVertices[i] = inFusedVertices[i];

        return new InflatableSheet(m, fusedVertices);
    }

    // Set methods
    ISHEET_API void isheet_inflatableSheet_setPressure(InflatableSheet *sheet, double pressure){
        sheet->setPressure(pressure);
    }

    ISHEET_API void isheet_inflatableSheet_setUseTensionFieldEnergy(InflatableSheet *sheet, int useTensionFieldEnergy){
        sheet->setUseTensionFieldEnergy(useTensionFieldEnergy);
    }

    ISHEET_API void isheet_inflatableSheet_setUseHessianProjectedEnergy(InflatableSheet *sheet, int useHessianProjectedEnergy){
        sheet->setUseHessianProjectedEnergy(useHessianProjectedEnergy);
    }

    ISHEET_API void isheet_inflatableSheet_disableFusedRegionTensionFieldTheory(InflatableSheet *sheet, int disableFusedRegionTensionFieldTheory){
        sheet->disableFusedRegionTensionFieldTheory(disableFusedRegionTensionFieldTheory);
    }

    ISHEET_API void isheet_inflatableSheet_setVars(InflatableSheet *sheet, double *inVars, size_t numVars){
        Eigen::VectorXd vars = Eigen::Map<Eigen::VectorXd>(inVars, numVars, 1);
        sheet->setVars(vars);
    }
    
    ISHEET_API void isheet_inflatableSheet_setGravity(InflatableSheet *sheet, double* inVector){
        const InflatableSheet::V3d gravity = InflatableSheet::V3d(inVector[0], inVector[1], inVector[2]);
        sheet->setGravity(gravity);
    }

    ISHEET_API void isheet_inflatableSheet_setMassDensity(InflatableSheet *sheet, double rho){
        sheet->setRho(rho);
    }

    ISHEET_API void isheet_inflatableSheet_setThickness(InflatableSheet *sheet, double thickness){
        sheet->setThickness(thickness);
    }

    ISHEET_API void isheet_inflatableSheet_setYoungModulus(InflatableSheet *sheet, double youngModulus){
        sheet->setYoungModulus(youngModulus);
    }

    ISHEET_API void isheet_inflatableSheet_setReferenceVolume(InflatableSheet *sheet, double volume){
        sheet->setReferenceVolume(volume);
    }

    // Get methods
    ISHEET_API double isheet_inflatableSheet_getMassDensity(InflatableSheet *sheet){
        return sheet->getRho();
    }

    ISHEET_API void isheet_inflatableSheet_getMeshVisualization(InflatableSheet *sheet, double **outCoords, int **outElements, size_t *numCoords, size_t *numElements)
    {
        const auto m = *sheet->visualizationMesh();

        const size_t nv = m.numVertices();
        const size_t nt = m.numElements();

        //////////////////////////

        *numCoords = nv * 3;
        auto sizeCoords = (*numCoords) * sizeof(double);
        *outCoords = static_cast<double *>(malloc(sizeCoords));
        std::vector<double> coords;
        for (const auto vi : m.vertices()){
            auto node = vi.node()->p.cast<Real>();
            coords.push_back(node[0]);
            coords.push_back(node[1]);
            coords.push_back(node[2]);
        }
        std::memcpy(*outCoords, coords.data(), sizeCoords);

        *numElements = nt * 3;
        auto sizeElements = (*numElements) * sizeof(int);
        *outElements = static_cast<int *>(malloc(sizeElements));
        std::vector<int> trias;
        for (const auto tri : m.elements()) {
            trias.push_back(tri.vertex(0).index());
            trias.push_back(tri.vertex(1).index());
            trias.push_back(tri.vertex(2).index());
        }
        std::memcpy(*outElements, trias.data(), sizeElements);
    }

    ISHEET_API double isheet_inflatableSheet_getVolume(InflatableSheet *sheet){
        return sheet->volume();
    }

    ISHEET_API double isheet_inflatableSheet_getYoungModulus(InflatableSheet *sheet){
        return sheet->getYoungModulus();
    }

    ISHEET_API double isheet_inflatableSheet_getThickness(InflatableSheet *sheet){
        return sheet->getThickness();
    }

    ISHEET_API void isheet_inflatableSheet_getCenterFixedVars(InflatableSheet *sheet, int **outCenterFixedVars, size_t *numCenterFixedVars)
    {
        const auto fixedVxIdx = sheet->center_non_fused_vx_idx();
        const size_t startVarIdx = fixedVxIdx * 3;
        const size_t endVarIdx = 3 + fixedVxIdx * 3;
        *numCenterFixedVars = endVarIdx - startVarIdx;

        auto sizeCenterFixedVars = (*numCenterFixedVars) * sizeof(int);
        *outCenterFixedVars = static_cast<int *>(malloc(sizeCenterFixedVars));
        std::vector<int> centerFixedVars;
        for (size_t i = startVarIdx; i < endVarIdx; i++) centerFixedVars.push_back(i);
        std::memcpy(*outCenterFixedVars, centerFixedVars.data(), sizeCenterFixedVars);
    }

    ISHEET_API void isheet_inflatableSheet_getVertexVars(InflatableSheet *sheet, int sheetIndex, int vertexIndex, int **outVars, size_t *numVars){
        const auto fixedVxIdx = sheet->varIdx(sheetIndex, vertexIndex,0)/3;

        const size_t startVarIdx = fixedVxIdx * 3;
        const size_t endVarIdx = 3 + fixedVxIdx * 3;
        *numVars = endVarIdx - startVarIdx;

        auto sizeFixedVars = (*numVars) * sizeof(int);
        *outVars = static_cast<int *>(malloc(sizeFixedVars));
        std::vector<int> fixedVars;
        for (size_t i = startVarIdx; i < endVarIdx; i++) fixedVars.push_back(i);
        std::memcpy(*outVars, fixedVars.data(), sizeFixedVars);
    }

    ISHEET_API int isheet_inflatableSheet_getCenterNonFusedVertexIdx(InflatableSheet *sheet)
    {
        return sheet->center_non_fused_vx_idx();
    }

    ISHEET_API int isheet_inflatableSheet_getNumVars(InflatableSheet *sheet)
    {
        return sheet->numVars();
    }

    ISHEET_API void isheet_inflatableSheet_getVars(InflatableSheet *sheet, double **outVars, size_t *numVars)
    {
        const auto vars = sheet->getVars();
        *numVars = vars.size();
        auto sizeVars = (*numVars) * sizeof(double);
        *outVars = static_cast<double *>(malloc(sizeVars));
        std::memcpy(*outVars, vars.data(), sizeVars);
    }

    ISHEET_API void isheet_inflatableSheet_getStrains(InflatableSheet *sheet, double **outStrains, size_t *numStrains)
    {
        const auto ted = sheet->triEnergyDensities();
        *numStrains = ted.size() * 2;
        auto sizeStrains = (*numStrains) * sizeof(double);
        *outStrains = static_cast<double *>(malloc(sizeStrains));
        std::vector<double> strains;
        for (const auto data : ted){
            const auto vec = data.eigSensitivities().Lambda();
            strains.push_back(sqrt(vec(0))-1);
            strains.push_back(sqrt(vec(1))-1);
        }
        std::memcpy(*outStrains, strains.data(), sizeStrains);
    }

    ISHEET_API double isheet_inflatableSheet_getEnergy(InflatableSheet *sheet){
        return sheet->energy();
    }

    ISHEET_API double isheet_inflatableSheet_getGradientNorm(InflatableSheet *sheet){
        return sheet->gradient().norm();
    }

    ISHEET_API void isheet_inflatableSheet_getGradient(InflatableSheet *sheet, double **outVars, size_t *numVars){
        const auto grad = sheet->gradient();
        *numVars = grad.size();
        auto sizeGrad = (*numVars) * sizeof(double);
        *outVars = static_cast<double *>(malloc(sizeGrad));
        std::memcpy(*outVars, grad.data(), sizeGrad);
    }

    ISHEET_API void isheet_inflatableSheet_getGravity(InflatableSheet *sheet, double **outVars, size_t *numVars){
        const auto gravity = sheet->getGravity();
        *numVars = gravity.size();
        auto sizeGravity = (*numVars) * sizeof(double);
        *outVars = static_cast<double *>(malloc(sizeGravity));
        std::memcpy(*outVars, gravity.data(), sizeGravity);
    }

    ISHEET_API double isheet_inflatableSheet_getReferenceVolume(InflatableSheet *sheet){
        return sheet->referenceVolume();
    }
    
    ISHEET_API double isheet_inflatableSheet_getPressure(InflatableSheet *sheet){
        return sheet->getPressure();
    }

    // Solver
    ISHEET_API int isheet_inflatableSheet_newtonSolver(InflatableSheet *sheet, int numSupports, int *supports, int numIterations, double gradTol, int writeReport, double **outReport, const char **errorMessage)
    {
        try
        {
            std::vector<size_t> fixedVars;
            for (int i = 0; i < numSupports; i++) fixedVars.push_back(supports[i]);

            // Set real options: [0] numIterations, [1] gradTol, [2] hessianShift
            NewtonOptimizerOptions options;
            options.niter = numIterations;
            options.gradTol = gradTol;

            const auto report = inflation_newton<InflatableSheet>(*sheet, fixedVars, options);

            if (writeReport) getConvergenceReport(report, numIterations, outReport);

            *errorMessage = "";

            return report.success;
        }
        catch (const std::runtime_error &error)
        {
            *errorMessage = error.what();
            return -1;
        }
        catch (const std::out_of_range &error)
        {
            *errorMessage = error.what();
            return -1;
        }
        catch (...)
        {
            *errorMessage = "Unknown Error from the Unmanaged Code.";
            return -1;
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // Periodic Unit
    ///////////////////////////////////////////////////////////////////////////
    ISHEET_API InflatablePeriodicUnit *isheet_periodicUnit_build(int numVertices, int numTrias, double *inCoords, int *inTrias, int *inFusedVertices, double epsilon){
        auto m = buildMesh(numVertices, numTrias, inCoords, inTrias);

        std::vector<bool> fusedVertices(numVertices);
        for (int i = 0; i < numVertices; i++) fusedVertices[i] = inFusedVertices[i];

        return new InflatablePeriodicUnit(m, fusedVertices, epsilon);
    }

    // Void methods
    ISHEET_API void isheet_periodicUnit_reparametrizeVerticalOffset(InflatablePeriodicUnit *ipu)
    {
        ipu->reparametrize_vertical_offset();
    }   

    // Set methods
    ISHEET_API void isheet_periodicUnit_setPressure(InflatablePeriodicUnit *ipu, double pressure){
        ipu->sheet.setPressure(pressure);
    }

    ISHEET_API void isheet_periodicUnit_setUseTensionFieldEnergy(InflatablePeriodicUnit *ipu, int useTensionFieldEnergy){
        ipu->sheet.setUseTensionFieldEnergy(useTensionFieldEnergy);
    }

    ISHEET_API void isheet_periodicUnit_setUseHessianProjectedEnergy(InflatablePeriodicUnit *ipu, int useHessianProjectedEnergy){
        ipu->sheet.setUseHessianProjectedEnergy(useHessianProjectedEnergy);
    }

    ISHEET_API void isheet_periodicUnit_disableFusedRegionTensionFieldTheory(InflatablePeriodicUnit *ipu, int disableFusedRegionTensionFieldTheory){
        ipu->sheet.disableFusedRegionTensionFieldTheory(disableFusedRegionTensionFieldTheory);
    }

    ISHEET_API void isheet_periodicUnit_setVars(InflatablePeriodicUnit *ipu, double *inVars, size_t numVars){
        Eigen::VectorXd vars = Eigen::Map<Eigen::VectorXd>(inVars, numVars, 1);
        ipu->sheet.setVars(vars);
    }
    
    ISHEET_API void isheet_periodicUnit_setMassDensity(InflatablePeriodicUnit *ipu, double rho){
        ipu->sheet.setRho(rho);
    }

    ISHEET_API void isheet_periodicUnit_setGravity(InflatablePeriodicUnit *ipu, double* inVector){
        const InflatableSheet::V3d gravity = InflatableSheet::V3d(inVector[0], inVector[1], inVector[2]);
        ipu->sheet.setGravity(gravity);
    }

    ISHEET_API void isheet_periodicUnit_setThickness(InflatablePeriodicUnit *ipu, double thickness){
        ipu->sheet.setThickness(thickness);
    }

    ISHEET_API void isheet_periodicUnit_setYoungModulus(InflatablePeriodicUnit *ipu, double youngModulus){
        ipu->sheet.setYoungModulus(youngModulus);
    }

    ISHEET_API void isheet_periodicUnit_setReferenceVolume(InflatablePeriodicUnit *ipu, double volume){
        ipu->sheet.setReferenceVolume(volume);
    }

    ISHEET_API double isheet_periodicUnit_getPressure(InflatablePeriodicUnit *ipu){
        return ipu->sheet.getPressure();
    }

    // Get methods
    ISHEET_API double isheet_periodicUnit_getMassDensity(InflatablePeriodicUnit *ipu){
        return ipu->sheet.getRho();
    }

    ISHEET_API void isheet_periodicUnit_getVertexVars(InflatablePeriodicUnit *ipu,  int sheetIndex, int vertexIndex, int **outVars, size_t *numVars){
        const auto fixedVxIdx = ipu->get_IPU_vidx_for_inflatable_vidx(ipu->sheet.varIdx(sheetIndex, vertexIndex,0)/3);

        const size_t startVarIdx = 3 + fixedVxIdx * 3;
        const size_t endVarIdx = 6 + fixedVxIdx * 3;
        *numVars = endVarIdx - startVarIdx;

        auto sizeFixedVars = (*numVars) * sizeof(int);
        *outVars = static_cast<int *>(malloc(sizeFixedVars));
        std::vector<int> fixedVars;
        for (size_t i = startVarIdx; i < endVarIdx; i++) fixedVars.push_back(i);
        std::memcpy(*outVars, fixedVars.data(), sizeFixedVars);
    }

    ISHEET_API void isheet_periodicUnit_getMeshVisualization(InflatablePeriodicUnit *ipu, double **outCoords, int **outElements, size_t *numCoords, size_t *numElements){
        const auto m = *ipu->visualizationMesh();
        getVisualizationMesh(m, outCoords, outElements, numCoords, numElements);
    }

    ISHEET_API void isheet_periodicUnit_getCenterFixedVars(InflatablePeriodicUnit *ipu, int **outCenterFixedVars, size_t *numCenterFixedVars)
    {
        const auto fixedVxIdx = ipu->get_IPU_vidx_for_inflatable_vidx(ipu->sheet.center_non_fused_vx_idx());
        const size_t startVarIdx = 3 + fixedVxIdx * 3;
        const size_t endVarIdx = 3 + 3 + fixedVxIdx * 3;
        *numCenterFixedVars = endVarIdx - startVarIdx;

        auto sizeCenterFixedVars = (*numCenterFixedVars) * sizeof(int);
        *outCenterFixedVars = static_cast<int *>(malloc(sizeCenterFixedVars));
        std::vector<int> centerFixedVars;
        for (size_t i = startVarIdx; i < endVarIdx; i++) centerFixedVars.push_back(i);
        std::memcpy(*outCenterFixedVars, centerFixedVars.data(), sizeCenterFixedVars);
    }

    ISHEET_API int isheet_periodicUnit_getCenterNonFusedVertexIdx(InflatablePeriodicUnit *ipu)
    {
        return ipu->get_IPU_vidx_for_inflatable_vidx(ipu->sheet.center_non_fused_vx_idx());
    }

    ISHEET_API int isheet_periodicUnit_getNumVars(InflatablePeriodicUnit *ipu)
    {
        return ipu->numVars();
    }

    ISHEET_API void isheet_periodicUnit_getVars(InflatablePeriodicUnit *ipu, double **outVars, size_t *numVars)
    {
        const auto vars = ipu->getVars();
        *numVars = vars.size();
        auto sizeVars = (*numVars) * sizeof(double);
        *outVars = static_cast<double *>(malloc(sizeVars));
        std::memcpy(*outVars, vars.data(), sizeVars);
    }

    ISHEET_API void isheet_periodicUnit_getStrains(InflatablePeriodicUnit *ipu, double **outStrains, size_t *numStrains)
    {
        const auto ted = ipu->sheet.triEnergyDensities();
        *numStrains = ted.size() * 2;
        auto sizeStrains = (*numStrains) * sizeof(double);
        *outStrains = static_cast<double *>(malloc(sizeStrains));
        std::vector<double> strains;
        for (const auto data : ted){
            const auto vec = data.eigSensitivities().Lambda();
            strains.push_back(sqrt(vec(0))-1);
            strains.push_back(sqrt(vec(1))-1);
        }
        std::memcpy(*outStrains, strains.data(), sizeStrains);
    }
    
    ISHEET_API double isheet_periodicUnit_getEnergy(InflatablePeriodicUnit *ipu){
        return ipu->energy();
    }

    ISHEET_API void isheet_periodicUnit_getBendingStiffnessFixedVars(InflatablePeriodicUnit *ipu, int **outFixedVars, size_t *numFixedVars){
        const auto fixedVars = ipu->getBendingStiffnessFixedVars();
        *numFixedVars = fixedVars.size();
        auto sizeFixedVars = (*numFixedVars) * sizeof(int);
        *outFixedVars = static_cast<int *>(malloc(sizeFixedVars));
        std::memcpy(*outFixedVars, fixedVars.data(), sizeFixedVars);
    }

    ISHEET_API void isheet_periodicUnit_getStretchingStiffnessFixedVars(InflatablePeriodicUnit *ipu, int **outFixedVars, size_t *numFixedVars){
        const auto fixedVars = ipu->getStretchingStiffnessFixedVars();
        *numFixedVars = fixedVars.size();
        auto sizeFixedVars = (*numFixedVars) * sizeof(int);
        *outFixedVars = static_cast<int *>(malloc(sizeFixedVars));
        std::memcpy(*outFixedVars, fixedVars.data(), sizeFixedVars);
    }

    ISHEET_API double isheet_periodicUnit_getGradientNorm(InflatablePeriodicUnit *ipu){
        return ipu->gradient().norm();
    }

    ISHEET_API void isheet_periodicUnit_getGradient(InflatablePeriodicUnit *ipu, double **outVars, size_t *numVars){
        const auto grad = ipu->gradient();
        *numVars = grad.size();
        auto sizeGrad = (*numVars) * sizeof(double);
        *outVars = static_cast<double *>(malloc(sizeGrad));
        std::memcpy(*outVars, grad.data(), sizeGrad);
    }

    ISHEET_API void isheet_periodicUnit_getGravity(InflatablePeriodicUnit *ipu, double **outVars, size_t *numVars){
        const auto gravity = ipu->sheet.getGravity();
        *numVars = gravity.size();
        auto sizeGravity = (*numVars) * sizeof(double);
        *outVars = static_cast<double *>(malloc(sizeGravity));
        std::memcpy(*outVars, gravity.data(), sizeGravity);
    }

    ISHEET_API double isheet_periodicUnit_getReferenceVolume(InflatablePeriodicUnit *ipu){
        return ipu->sheet.referenceVolume();
    }

    ISHEET_API double isheet_periodicUnit_getYoungModulus(InflatablePeriodicUnit *ipu){
            return ipu->sheet.getYoungModulus();
    }

    ISHEET_API double isheet_periodicUnit_getThickness(InflatablePeriodicUnit *ipu){
        return ipu->sheet.getThickness();
    }

    ISHEET_API double isheet_periodicUnit_getVolume(InflatablePeriodicUnit *ipu){
            return ipu->sheet.volume();
    }
    // Solver
    ISHEET_API int isheet_periodicUnit_newtonSolver(InflatablePeriodicUnit *ipu, int numSupports, int *supports, int numIterations, double gradTol, double hessianShift, int writeReport, double **outReport, const char **errorMessage)
    {
        try
        {
            std::vector<size_t> fixedVars;
            for (int i = 0; i < numSupports; i++) fixedVars.push_back(supports[i]);

            NewtonOptimizerOptions options;
            options.niter = numIterations;
            options.gradTol = gradTol;

            const auto report = inflation_newton<InflatablePeriodicUnit>(*ipu, fixedVars, options, nullptr, hessianShift);

            if (writeReport) getConvergenceReport(report, numIterations, outReport);

            *errorMessage = "";

            return report.success;
        }
        catch (const std::runtime_error &error)
        {
            *errorMessage = error.what();
            return -1;
        }
        catch (const std::out_of_range &error)
        {
            *errorMessage = error.what();
            return -1;
        }
        catch (...)
        {
            *errorMessage = "Unknown Error from the Unmanaged Code.";
            return -1;
        }
    }

    ISHEET_API int isheet_periodicUnit_newtonStepSolver(InflatablePeriodicUnit *ipu, double pressure, int numSupports, int *supports, int numIterations, double gradTol, double hessianShift, int writeReport, double **outReport, const char **errorMessage)
    {
        try
        {
            ipu->setPressure(pressure);

            std::vector<size_t> fixedVars;
            for (int i = 0; i < numSupports; i++) fixedVars.push_back(supports[i]);

            NewtonOptimizerOptions options;
            options.niter = numIterations;
            options.gradTol = gradTol;

            const auto report = inflation_newton<InflatablePeriodicUnit>(*ipu, fixedVars, options, nullptr, hessianShift);

            if (writeReport) getConvergenceReport(report, numIterations, outReport);

            *errorMessage = "";

            return report.success;
        }
        catch (const std::runtime_error &error)
        {
            *errorMessage = error.what();
            return -1;
        }
        catch (const std::out_of_range &error)
        {
            *errorMessage = error.what();
            return -1;
        }
        catch (...)
        {
            *errorMessage = "Unknown Error from the Unmanaged Code.";
            return -1;
        }
    }

    ISHEET_API int isheet_periodicUnit_computeBendingStiffness(InflatablePeriodicUnit *ipu, int numSupports, int *supports, int numAlphas, double* alphas, int numIterations, double gradTol, double hessianShift, int useBases, int* outNumBendingStiffness, double** outBendingStiffness, int* outNumStiffnessCoefficient, double** outStiffnessCoefficient, const char **errorMessage){
        try
        {
            std::vector<size_t> fixedVars;
            for (int i = 0; i < numSupports; i++) fixedVars.push_back(supports[i]);
            
            Eigen::VectorXd alphaVars(numAlphas);
            for (int i = 0; i < numAlphas; i++) alphaVars(i) = alphas[i];

            NewtonOptimizerOptions options;
            options.niter = numIterations;
            options.gradTol = gradTol;

            const auto optimizer = get_inflation_optimizer<InflatablePeriodicUnit>(*ipu, fixedVars, options, nullptr, hessianShift);
            const auto stiffnessPair = getBendingStiffnessUsingBases(*ipu, alphaVars, *optimizer, hessianShift, fixedVars);

            *outNumBendingStiffness = stiffnessPair.first.size();
            auto sizeStiffness = (*outNumBendingStiffness) * sizeof(double);
            *outBendingStiffness = static_cast<double *>(malloc(sizeStiffness));
            std::memcpy(*outBendingStiffness, stiffnessPair.first.data(), sizeStiffness);

            if(useBases){
                *outNumStiffnessCoefficient = stiffnessPair.second.size();
                auto sizeStiffnessCoef = (*outNumStiffnessCoefficient) * sizeof(double);
                *outStiffnessCoefficient = static_cast<double *>(malloc(sizeStiffnessCoef));
                std::memcpy(*outStiffnessCoefficient, stiffnessPair.second.data(), sizeStiffnessCoef);
            }

            *errorMessage = "";
            return 1;
        }
        catch (const std::runtime_error &error)
        {
            *errorMessage = error.what();
            return -1;
        }
        catch (const std::out_of_range &error)
        {
            *errorMessage = error.what();
            return -1;
        }
        catch (...)
        {
            *errorMessage = "Unknown Error from the Unmanaged Code.";
            return -1;
        }
    }

    ISHEET_API int isheet_periodicUnit_computeStretchingStiffness(InflatablePeriodicUnit *ipu, int numSupports, int *supports, int numAlphas, double* alphas, int numIterations, double gradTol, double hessianShift, int* outNumStretchingStiffness, double** outStretchingStiffness, const char **errorMessage)
    {
        try
        {
            std::vector<size_t> fixedVars;
            for (int i = 0; i < numSupports; i++) fixedVars.push_back(supports[i]);
            
            Eigen::VectorXd alphaVars(numAlphas);
            for (int i = 0; i < numAlphas; i++) alphaVars(i) = alphas[i];

            NewtonOptimizerOptions options;
            options.niter = numIterations;
            options.gradTol = gradTol;

            const auto optimizer = get_inflation_optimizer<InflatablePeriodicUnit>(*ipu, fixedVars, options, nullptr, hessianShift);
            const auto stiffnessPair = getStretchingStiffness(*ipu, alphaVars, *optimizer, hessianShift, fixedVars);

            *outNumStretchingStiffness = stiffnessPair.size();
            auto sizeStiffness = (*outNumStretchingStiffness) * sizeof(double);
            *outStretchingStiffness = static_cast<double *>(malloc(sizeStiffness));
            std::memcpy(*outStretchingStiffness, stiffnessPair.data(), sizeStiffness);

            *errorMessage = "";
            return 1;
        }
        catch (const std::runtime_error &error)
        {
            *errorMessage = error.what();
            return -1;
        }
        catch (const std::out_of_range &error)
        {
            *errorMessage = error.what();
            return -1;
        }
        catch (...)
        {
            *errorMessage = "Unknown Error from the Unmanaged Code.";
            return -1;
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // Target attracted inflation
    ///////////////////////////////////////////////////////////////////////////

    ISHEET_API TargetAttractedInflation *isheet_targetAttractedInflation_build(int numVertices, int numTrias, Real *inCoords, int *inTrias, int *inFusedVertices, int numTargetVertices, int numTargetTrias, double *inTargetCoords, int *inTargetTrias, const char **errorMessage)
    {
        try
        {
            // Target surface
            std::vector<MeshIO::IOVertex> vertices;
            for (int i = 0; i < numTargetVertices; i++) vertices.emplace_back(inTargetCoords[3 * i], inTargetCoords[3 * i + 1], inTargetCoords[3 * i + 2]);
            // Triangles
            std::vector<MeshIO::IOElement> triangles;
            for (int i = 0; i < numTargetTrias; i++) triangles.emplace_back(inTargetTrias[3 * i], inTargetTrias[3 * i + 1], inTargetTrias[3 * i + 2]);
            
            const TargetAttractedInflation::Mesh targetSrf(triangles, vertices);

            // Inflatable sheet
            auto m = buildMesh(numVertices, numTrias, inCoords, inTrias);

            std::vector<bool> fusedVertices(numVertices);
            for (int i = 0; i < numVertices; i++) fusedVertices[i] = inFusedVertices[i];

            auto sheet = std::make_shared<InflatableSheet>(m, fusedVertices);

            return new TargetAttractedInflation(sheet, targetSrf);
        }
        catch (const std::runtime_error &error)
        {
            *errorMessage = error.what();
            return nullptr;
        }
        catch (const std::out_of_range &error)
        {
            *errorMessage = error.what();
            return nullptr;
        }
        catch (...)
        {
            *errorMessage = "Unknown Error from the Unmanaged Code.";
            return nullptr;
        }
    }

    ISHEET_API TargetAttractedInflation *isheet_targetAttractedInflation_build_new(InflatableSheet *sheet, int numTargetVertices, int numTargetTrias, double *inTargetCoords, int *inTargetTrias, const char **errorMessage)
    {
        try
        {
            // Target surface
            std::vector<MeshIO::IOVertex> vertices;
            for (int i = 0; i < numTargetVertices; i++) vertices.emplace_back(inTargetCoords[3 * i], inTargetCoords[3 * i + 1], inTargetCoords[3 * i + 2]);
            // Triangles
            std::vector<MeshIO::IOElement> triangles;
            for (int i = 0; i < numTargetTrias; i++) triangles.emplace_back(inTargetTrias[3 * i], inTargetTrias[3 * i + 1], inTargetTrias[3 * i + 2]);
            
            const TargetAttractedInflation::Mesh targetSrf(triangles, vertices);

            return new TargetAttractedInflation(std::shared_ptr<InflatableSheet>(sheet), targetSrf);
        }
        catch (const std::runtime_error &error)
        {
            *errorMessage = error.what();
            return nullptr;
        }
        catch (const std::out_of_range &error)
        {
            *errorMessage = error.what();
            return nullptr;
        }
        catch (...)
        {
            *errorMessage = "Unknown Error from the Unmanaged Code.";
            return nullptr;
        }
    }

    // Get methods
    ISHEET_API void isheet_targetAttractedInflation_getMeshVisualization(TargetAttractedInflation *attractedSheet, double **outCoords, int **outElements, size_t *numCoords, size_t *numElements){
        const auto m = *attractedSheet->sheet().visualizationMesh();
        getVisualizationMesh(m, outCoords, outElements, numCoords, numElements);
    }
 
    ISHEET_API double isheet_targetAttractedInflation_getMassDensity(TargetAttractedInflation *attractedSheet){
        return attractedSheet->sheet().getRho();
    }

    ISHEET_API double isheet_targetAttractedInflation_getVolume(TargetAttractedInflation *attractedSheet){
        return attractedSheet->sheet().volume();
    }

    ISHEET_API double isheet_targetAttractedInflation_getYoungModulus(TargetAttractedInflation *attractedSheet){
        return attractedSheet->sheet().getYoungModulus();
    }

    ISHEET_API double isheet_targetAttractedInflation_getThickness(TargetAttractedInflation *attractedSheet){
        return attractedSheet->sheet().getThickness();
    }

    ISHEET_API void isheet_targetAttractedInflation_getCenterFixedVars(TargetAttractedInflation *attractedSheet, int **outCenterFixedVars, size_t *numCenterFixedVars)
    {
        const auto fixedVxIdx = attractedSheet->sheet().center_non_fused_vx_idx();
        const size_t startVarIdx = fixedVxIdx * 3;
        const size_t endVarIdx = 3 + fixedVxIdx * 3;
        *numCenterFixedVars = endVarIdx - startVarIdx;

        auto sizeCenterFixedVars = (*numCenterFixedVars) * sizeof(int);
        *outCenterFixedVars = static_cast<int *>(malloc(sizeCenterFixedVars));
        std::vector<int> centerFixedVars;
        for (size_t i = startVarIdx; i < endVarIdx; i++) centerFixedVars.push_back(i);
        std::memcpy(*outCenterFixedVars, centerFixedVars.data(), sizeCenterFixedVars);
    }

    ISHEET_API void isheet_targetAttractedInflation_getVertexVars(TargetAttractedInflation *attractedSheet, int sheetIndex, int vertexIndex, int **outVars, size_t *numVars){
        const auto fixedVxIdx = attractedSheet->sheet().varIdx(sheetIndex, vertexIndex,0)/3;

        const size_t startVarIdx = fixedVxIdx * 3;
        const size_t endVarIdx = 3 + fixedVxIdx * 3;
        *numVars = endVarIdx - startVarIdx;

        auto sizeFixedVars = (*numVars) * sizeof(int);
        *outVars = static_cast<int *>(malloc(sizeFixedVars));
        std::vector<int> fixedVars;
        for (size_t i = startVarIdx; i < endVarIdx; i++) fixedVars.push_back(i);
        std::memcpy(*outVars, fixedVars.data(), sizeFixedVars);
    }

    ISHEET_API int isheet_targetAttractedInflation_getCenterNonFusedVertexIdx(TargetAttractedInflation *attractedSheet)
    {
        return attractedSheet->sheet().center_non_fused_vx_idx();
    }

    ISHEET_API int isheet_targetAttractedInflation_getNumVars(TargetAttractedInflation *attractedSheet)
    {
        return attractedSheet->sheet().numVars();
    }

    ISHEET_API void isheet_targetAttractedInflation_getVars(TargetAttractedInflation *attractedSheet, double **outVars, size_t *numVars)
    {
        const auto vars = attractedSheet->sheet().getVars();
        *numVars = vars.size();
        auto sizeVars = (*numVars) * sizeof(double);
        *outVars = static_cast<double *>(malloc(sizeVars));
        std::memcpy(*outVars, vars.data(), sizeVars);
    }

    ISHEET_API void isheet_targetAttractedInflation_getStrains(TargetAttractedInflation *attractedSheet, double **outStrains, size_t *numStrains)
    {
        const auto ted = attractedSheet->sheet().triEnergyDensities();
        *numStrains = ted.size() * 2;
        auto sizeStrains = (*numStrains) * sizeof(double);
        *outStrains = static_cast<double *>(malloc(sizeStrains));
        std::vector<double> strains;
        for (const auto data : ted){
            const auto vec = data.eigSensitivities().Lambda();
            strains.push_back(sqrt(vec(0))-1);
            strains.push_back(sqrt(vec(1))-1);
        }
        std::memcpy(*outStrains, strains.data(), sizeStrains);
    }

    ISHEET_API double isheet_targetAttractedInflation_getEnergy(TargetAttractedInflation *attractedSheet){
        return attractedSheet->sheet().energy();
    }

    ISHEET_API double isheet_targetAttractedInflation_getGradientNorm(TargetAttractedInflation *attractedSheet){
        return attractedSheet->sheet().gradient().norm();
    }

    ISHEET_API void isheet_targetAttractedInflation_getGradient(TargetAttractedInflation *attractedSheet, double **outVars, size_t *numVars){
        const auto grad = attractedSheet->sheet().gradient();
        *numVars = grad.size();
        auto sizeGrad = (*numVars) * sizeof(double);
        *outVars = static_cast<double *>(malloc(sizeGrad));
        std::memcpy(*outVars, grad.data(), sizeGrad);
    }

    ISHEET_API void isheet_targetAttractedInflation_getGravity(TargetAttractedInflation *attractedSheet, double **outVars, size_t *numVars){
        const auto gravity = attractedSheet->sheet().getGravity();
        *numVars = gravity.size();
        auto sizeGravity = (*numVars) * sizeof(double);
        *outVars = static_cast<double *>(malloc(sizeGravity));
        std::memcpy(*outVars, gravity.data(), sizeGravity);
    }

    ISHEET_API double isheet_targetAttractedInflation_getReferenceVolume(TargetAttractedInflation *attractedSheet){
        return attractedSheet->sheet().referenceVolume();
    }

    ISHEET_API double isheet_targetAttractedInflation_getPressure(TargetAttractedInflation *attractedSheet){
        return attractedSheet->sheet().getPressure();
    }

    // Set methods
    ISHEET_API void isheet_targetAttractedInflation_setPressure(TargetAttractedInflation *attractedSheet, double pressure){
        attractedSheet->sheet().setPressure(pressure);
    }

    ISHEET_API void isheet_targetAttractedInflation_setUseTensionFieldEnergy(TargetAttractedInflation *attractedSheet, int useTensionFieldEnergy){
        attractedSheet->sheet().setUseTensionFieldEnergy(useTensionFieldEnergy);
    }

    ISHEET_API void isheet_targetAttractedInflation_setUseHessianProjectedEnergy(TargetAttractedInflation *attractedSheet, int useHessianProjectedEnergy){
        attractedSheet->sheet().setUseHessianProjectedEnergy(useHessianProjectedEnergy);
    }

    ISHEET_API void isheet_targetAttractedInflation_disableFusedRegionTensionFieldTheory(TargetAttractedInflation *attractedSheet, int disableFusedRegionTensionFieldTheory){
        attractedSheet->sheet().disableFusedRegionTensionFieldTheory(disableFusedRegionTensionFieldTheory);
    }

    ISHEET_API void isheet_targetAttractedInflation_setVars(TargetAttractedInflation *attractedSheet, double *inVars, size_t numVars){
        Eigen::VectorXd vars = Eigen::Map<Eigen::VectorXd>(inVars, numVars, 1);
        attractedSheet->sheet().setVars(vars);
    }
    
    ISHEET_API void isheet_targetAttractedInflation_setMassDensity(TargetAttractedInflation *attractedSheet, double rho){
        attractedSheet->sheet().setRho(rho);
    }

    ISHEET_API void isheet_targetAttractedInflation_setGravity(TargetAttractedInflation *attractedSheet, double* inVector){
        const InflatableSheet::V3d gravity = InflatableSheet::V3d(inVector[0], inVector[1], inVector[2]);
        attractedSheet->sheet().setGravity(gravity);
    }

    ISHEET_API void isheet_targetAttractedInflation_setThickness(TargetAttractedInflation *attractedSheet, double thickness){
        attractedSheet->sheet().setThickness(thickness);
    }

    ISHEET_API void isheet_targetAttractedInflation_setYoungModulus(TargetAttractedInflation *attractedSheet, double youngModulus){
        attractedSheet->sheet().setYoungModulus(youngModulus);
    }

    ISHEET_API void isheet_targetAttractedInflation_setReferenceVolume(TargetAttractedInflation *attractedSheet, double volume){
        attractedSheet->sheet().setReferenceVolume(volume);
    }

    // Solver
    ISHEET_API int isheet_targetAttractedInflation_newtonSolver(TargetAttractedInflation *sheet, int numSupports, int *supports, int numIterations, double gradTol, int writeReport, double **outReport, const char **errorMessage)
    {
        try
        {
            std::vector<size_t> fixedVars;
            for (int i = 0; i < numSupports; i++) fixedVars.push_back(supports[i]);

            // Set real options: [0] numIterations, [1] gradTol, [2] hessianShift
            NewtonOptimizerOptions options;
            options.niter = numIterations;
            options.gradTol = gradTol;

            const auto report = inflation_newton<TargetAttractedInflation>(*sheet, fixedVars, options);

            if (writeReport) getConvergenceReport(report, numIterations, outReport);

            *errorMessage = "";

            return report.success;
        }
        catch (const std::runtime_error &error)
        {
            *errorMessage = error.what();
            return -1;
        }
        catch (const std::out_of_range &error)
        {
            *errorMessage = error.what();
            return -1;
        }
        catch (...)
        {
            *errorMessage = "Unknown Error from the Unmanaged Code.";
            return -1;
        }
    }
}
