#include <MeshFEM/GlobalBenchmark.hh>
#include <MeshFEM/MeshIO.hh>
#include <MeshFEM/Future.hh>
#include <MeshFEM/newton_optimizer/newton_optimizer.hh>

#include "InflatableSheet.hh"
#include "inflation_newton.hh"
#include "InflatablePeriodicUnit.hh"
#include "TargetAttractedInflation.hh"
#include "MultilayerInflatable.hh"

#ifndef ISHEET_H
#define ISHEET_H

#define ISHEET_API_VERSION "isheet"

#if defined(ISHEET_DLL)
#if defined(ISHEET_DLL_EXPORT)
#define ISHEET_API __declspec(dllexport)
#else
#define ISHEET_API __declspec(dllimport)
#endif
#else
#define ISHEET_API __attribute__((visibility("default")))
#endif

namespace InflatableSheetGH {
    ///////////////////////////////////////////////////////////////////////////
    // Multi layer inflatable Sheet
    ///////////////////////////////////////////////////////////////////////////
    ISHEET_API MultilayerInflatable *isheet_multilayerInflatableSheet_build(int numVertices, int numTrias, Real *inCoords, int *inTrias, int numSheets, Real *inPressures, int *inReducedVarIdxForVertexOnSheet, const char **errorMessage);

    // Solver
    ISHEET_API int isheet_multilayerInflatableSheet_newtonSolver(MultilayerInflatable *sheet, int numSupports, int *supports, int numIterations, double gradTol, int writeReport, double **outReport, const char **errorMessage);

    // Set methods
    ISHEET_API void isheet_multilayerInflatableSheet_setPressure(MultilayerInflatable *sheet, double *inPressures, size_t numPressures);

    ISHEET_API void isheet_multilayerInflatableSheet_setUseTensionFieldEnergy(MultilayerInflatable *sheet, int useTensionFieldEnergy);

    ISHEET_API void isheet_multilayerInflatableSheet_setUseHessianProjectedEnergy(MultilayerInflatable *sheet, int useHessianProjectedEnergy);

    ISHEET_API void isheet_multilayerInflatableSheet_setVars(MultilayerInflatable *sheet, double *inVars, size_t numVars);
    
    ISHEET_API void isheet_multilayerInflatableSheet_setGravity(MultilayerInflatable *sheet, double* inVector);

    ISHEET_API void isheet_multilayerInflatableSheet_setMassDensity(MultilayerInflatable *sheet, double rho);

    ISHEET_API void isheet_multilayerInflatableSheet_setThickness(MultilayerInflatable *sheet, double thickness);

    ISHEET_API void isheet_multilayerInflatableSheet_setYoungModulus(MultilayerInflatable *sheet, double youngModulus);

    ISHEET_API void isheet_multilayerInflatableSheet_setReferenceVolume(MultilayerInflatable *sheet, double *inVolumes, size_t numVolumes);

    // Get methods
    ISHEET_API void isheet_multilayerInflatableSheet_getMeshVisualization(MultilayerInflatable *sheet, double **outCoords, int **outElements, size_t *numCoords, size_t *numElements);

    ISHEET_API double isheet_multilayerInflatableSheet_getMassDensity(MultilayerInflatable *sheet);
    
    ISHEET_API void isheet_multilayerInflatableSheet_getVolumes(MultilayerInflatable *sheet,  double **outVolumes, size_t *numVolumes);

    ISHEET_API double isheet_multilayerInflatableSheet_getYoungModulus(MultilayerInflatable *sheet);

    ISHEET_API double isheet_multilayerInflatableSheet_getThickness(MultilayerInflatable *sheet);

    ISHEET_API void isheet_multilayerInflatableSheet_getCenterFixedVars(MultilayerInflatable *sheet, int **outCenterFixedVars, size_t *numCenterFixedVars);

    ISHEET_API int isheet_multilayerInflatableSheet_getCenterNonFusedVertexIdx(MultilayerInflatable *sheet);

    ISHEET_API int isheet_multilayerInflatableSheet_getNumVars(MultilayerInflatable *sheet);

    ISHEET_API void isheet_multilayerInflatableSheet_getVars(MultilayerInflatable *sheet, double **outVars, size_t *numVars);

    ISHEET_API void isheet_multilayerInflatableSheet_getStrains(MultilayerInflatable *sheet, double **outStrains, size_t *numStrains);

    ISHEET_API double isheet_multilayerInflatableSheet_getEnergy(MultilayerInflatable *sheet);

    ISHEET_API double isheet_multilayerInflatableSheet_getGradientNorm(MultilayerInflatable *sheet);

    ISHEET_API void isheet_multilayerInflatableSheet_getGradient(MultilayerInflatable *sheet, double **outVars, size_t *numVars);

    ISHEET_API void isheet_multilayerInflatableSheet_getVertexVars(MultilayerInflatable *sheet, int sheetIndex, int vertexIndex, int **outVars, size_t *numVars);

    ISHEET_API void isheet_multilayerInflatableSheet_getGravity(MultilayerInflatable *sheet, double **outVars, size_t *numVars);

    ISHEET_API void isheet_multilayerInflatableSheet_getReferenceVolume(MultilayerInflatable *sheet, double **outVolumes, size_t *numVolumes);

    ISHEET_API void isheet_multilayerInflatableSheet_getPressure(MultilayerInflatable *sheet, double **outPressures, size_t *numPressures);

    ///////////////////////////////////////////////////////////////////////////
    // Inflatable Sheet
    ///////////////////////////////////////////////////////////////////////////
    ISHEET_API InflatableSheet *isheet_inflatableSheet_build(int numVertices, int numElements, Real *inCoords, int *inElements, int *inFusedVertices);

    // Set methods
    ISHEET_API void isheet_inflatableSheet_setPressure(InflatableSheet *sheet, double pressure);

    ISHEET_API void isheet_inflatableSheet_setUseTensionFieldEnergy(InflatableSheet *sheet, int useTensionFieldEnergy);

    ISHEET_API void isheet_inflatableSheet_setUseHessianProjectedEnergy(InflatableSheet *sheet, int useHessianProjectedEnergy);

    ISHEET_API void isheet_inflatableSheet_disableFusedRegionTensionFieldTheory(InflatableSheet *sheet, int disableFusedRegionTensionFieldTheory);

    ISHEET_API void isheet_inflatableSheet_setVars(InflatableSheet *sheet, double *inVars, size_t numVars);

    ISHEET_API void isheet_inflatableSheet_setGravity(InflatableSheet *sheet, double* inVector);

    ISHEET_API void isheet_inflatableSheet_setMassDensity(InflatableSheet *sheet, double rho);

    ISHEET_API void isheet_inflatableSheet_setThickness(InflatableSheet *sheet, double thickness);

    ISHEET_API void isheet_inflatableSheet_setYoungModulus(InflatableSheet *sheet, double youngModulus);

    ISHEET_API void isheet_inflatableSheet_setReferenceVolume(InflatableSheet *sheet, double volume);

    // Get methods
    ISHEET_API double isheet_inflatableSheet_getMassDensity(InflatableSheet *sheet);

    ISHEET_API void isheet_inflatableSheet_getMeshVisualization(InflatableSheet *sheet, double **outCoords, int **outElements, size_t *numCoords, size_t *numElements);
    
    ISHEET_API double isheet_inflatableSheet_getVolume(InflatableSheet *sheet);

    ISHEET_API double isheet_inflatableSheet_getYoungModulus(InflatableSheet *sheet);

    ISHEET_API double isheet_inflatableSheet_getThickness(InflatableSheet *sheet);

    ISHEET_API void isheet_inflatableSheet_getCenterFixedVars(InflatableSheet *sheet, int **outCenterFixedVars, size_t *numCenterFixedVars);

    ISHEET_API int isheet_inflatableSheet_getCenterNonFusedVertexIdx(InflatableSheet *sheet);

    ISHEET_API int isheet_inflatableSheet_getNumVars(InflatableSheet *sheet);

    ISHEET_API void isheet_inflatableSheet_getVars(InflatableSheet *sheet, double **outVars, size_t *numVars);

    ISHEET_API void isheet_inflatableSheet_getStrains(InflatableSheet *sheet, double **outStrains, size_t *numStrains);

    ISHEET_API double isheet_inflatableSheet_getEnergy(InflatableSheet *sheet);

    ISHEET_API double isheet_inflatableSheet_getGradientNorm(InflatableSheet *sheet);

    ISHEET_API void isheet_inflatableSheet_getGradient(InflatableSheet *sheet, double **outVars, size_t *numVars);

    ISHEET_API void isheet_inflatableSheet_getVertexVars(InflatableSheet *sheet, int sheetIndex, int vertexIndex, int **outVars, size_t *numVars);

    ISHEET_API void isheet_inflatableSheet_getGravity(InflatableSheet *sheet, double **outVars, size_t *numVars);

    ISHEET_API double isheet_inflatableSheet_getReferenceVolume(InflatableSheet *sheet);

    ISHEET_API double isheet_inflatableSheet_getPressure(InflatableSheet *sheet);

    // Solver
    ISHEET_API int isheet_inflatableSheet_newtonSolver(InflatableSheet *sheet, int numSupports, int *supports, int numIterations, double gradTol, int writeReport, double **outReport, const char **errorMessage);

    ///////////////////////////////////////////////////////////////////////////
    // Periodic Unit
    ///////////////////////////////////////////////////////////////////////////
    ISHEET_API InflatablePeriodicUnit *isheet_periodicUnit_build(int numVertices, int numElements, double *inCoords, int *inElements, int *inFusedVertices, double epsilon);

    // Void methods
    ISHEET_API void isheet_periodicUnit_reparametrizeVerticalOffset(InflatablePeriodicUnit *ipu);

    // Set methods
    ISHEET_API void isheet_periodicUnit_setPressure(InflatablePeriodicUnit *ipu, double pressure);

    ISHEET_API void isheet_periodicUnit_setUseTensionFieldEnergy(InflatablePeriodicUnit *ipu, int useTensionFieldEnergy);

    ISHEET_API void isheet_periodicUnit_setUseHessianProjectedEnergy(InflatablePeriodicUnit *ipu, int useHessianProjectedEnergy);

    ISHEET_API void isheet_periodicUnit_disableFusedRegionTensionFieldTheory(InflatablePeriodicUnit *ipu, int disableFusedRegionTensionFieldTheory);

    ISHEET_API void isheet_periodicUnit_setVars(InflatablePeriodicUnit *ipu, double *inVars, size_t numVars);

    ISHEET_API void isheet_periodicUnit_setGravity(InflatablePeriodicUnit *ipu, double* inVector);

    ISHEET_API void isheet_periodicUnit_setMassDensity(InflatablePeriodicUnit *ipu, double rho);

    ISHEET_API void isheet_periodicUnit_setThickness(InflatablePeriodicUnit *ipu, double thickness);

    ISHEET_API void isheet_periodicUnit_setYoungModulus(InflatablePeriodicUnit *ipu, double youngModulus);

    ISHEET_API void isheet_periodicUnit_setReferenceVolume(InflatablePeriodicUnit *ipu, double volume);

    // Get methods
    ISHEET_API double isheet_periodicUnit_getMassDensity(InflatablePeriodicUnit *ipu);

    ISHEET_API void isheet_periodicUnit_getVertexVars(InflatablePeriodicUnit *ipu, int sheetIndex, int vertexIndex, int **outVars, size_t *numVars);

    ISHEET_API void isheet_periodicUnit_getMeshVisualization(InflatablePeriodicUnit *ipu, double **outCoords, int **outElements, size_t *numCoords, size_t *numElements);

    ISHEET_API void isheet_periodicUnit_getCenterFixedVars(InflatablePeriodicUnit *ipu, int **outCenterFixedVars, size_t *numCenterFixedVars);

    ISHEET_API int isheet_periodicUnit_getCenterNonFusedVertexIdx(InflatablePeriodicUnit *ipu);

    ISHEET_API int isheet_periodicUnit_getNumVars(InflatablePeriodicUnit *ipu);

    ISHEET_API void isheet_periodicUnit_getVars(InflatablePeriodicUnit *ipu, double **outVars, size_t *numVars);

    ISHEET_API void isheet_periodicUnit_getStrains(InflatablePeriodicUnit *ipu, double **outStrains, size_t *numStrains);

    ISHEET_API double isheet_periodicUnit_getEnergy(InflatablePeriodicUnit *ipu);

    ISHEET_API void isheet_periodicUnit_getBendingStiffnessFixedVars(InflatablePeriodicUnit *ipu, int **outFixedVars, size_t *numFixedVars);

    ISHEET_API void isheet_periodicUnit_getStretchingStiffnessFixedVars(InflatablePeriodicUnit *ipu, int **outFixedVars, size_t *numFixedVars);

    ISHEET_API double isheet_periodicUnit_getGradientNorm(InflatablePeriodicUnit *ipu);

    ISHEET_API void isheet_periodicUnit_getGradient(InflatablePeriodicUnit *ipu, double **outVars, size_t *numVars);

    ISHEET_API void isheet_periodicUnit_getGravity(InflatablePeriodicUnit *ipu, double **outVars, size_t *numVars);

    ISHEET_API double isheet_periodicUnit_getReferenceVolume(InflatablePeriodicUnit *ipu);

    ISHEET_API double isheet_periodicUnit_getYoungModulus(InflatablePeriodicUnit *ipu);

    ISHEET_API double isheet_periodicUnit_getThickness(InflatablePeriodicUnit *ipu);

    ISHEET_API double isheet_periodicUnit_getVolume(InflatablePeriodicUnit *ipu);

    ISHEET_API double isheet_periodicUnit_getPressure(InflatablePeriodicUnit *ipu);

    // Solver
    ISHEET_API int isheet_periodicUnit_newtonSolver(InflatablePeriodicUnit *ipu, int numSupports, int *supports, int numIterations, double gradTol, double hessianShift, int writeReport, double **outReport, const char **errorMessage);

    ISHEET_API int isheet_periodicUnit_newtonStepSolver(InflatablePeriodicUnit *ipu, double pressure, int numSupports, int *supports, int numIterations, double gradTol, double hessianShift, int writeReport, double **outReport, const char **errorMessage);

    ISHEET_API int isheet_periodicUnit_computeBendingStiffness(InflatablePeriodicUnit *ipu, int numSupports, int *supports, int numAlphas, double* alphas, int numIterations, double gradTol, double hessianShift, int useBases, int* outNumBendingStiffness, double** outBendingStiffness, int* outNumStiffnessCoefficient, double** outStiffnessCoefficient, const char **errorMessage);

    ISHEET_API int isheet_periodicUnit_computeStretchingStiffness(InflatablePeriodicUnit *ipu, int numSupports, int *supports, int numAlphas, double* alphas, int numIterations, double gradTol, double hessianShift, int* outNumStretchingStiffness, double** outStretchingStiffness, const char **errorMessage);

    ///////////////////////////////////////////////////////////////////////////
    // Target attracted inflation
    ///////////////////////////////////////////////////////////////////////////

    ISHEET_API TargetAttractedInflation *isheet_targetAttractedInflation_build(int numVertices, int numTrias, Real *inCoords, int *inTrias, int *inFusedVertices, int numTargetVertices, int numTargetTrias, double *inTargetCoords, int *inTargetTrias, const char **errorMessage);

    ISHEET_API TargetAttractedInflation *isheet_targetAttractedInflation_build_new(InflatableSheet *sheet, int numTargetVertices, int numTargetTrias, double *inTargetCoords, int *inTargetTrias, const char **errorMessage);

    // Get methods
    ISHEET_API void isheet_targetAttractedInflation_getMeshVisualization(TargetAttractedInflation *attractedSheet, double **outCoords, int **outElements, size_t *numCoords, size_t *numElements);

    ISHEET_API double isheet_targetAttractedInflation_getMassDensity(TargetAttractedInflation *attractedSheet);

    ISHEET_API double isheet_targetAttractedInflation_getVolume(TargetAttractedInflation *attractedSheet);

    ISHEET_API double isheet_targetAttractedInflation_getYoungModulus(TargetAttractedInflation *attractedSheet);

    ISHEET_API double isheet_targetAttractedInflation_getThickness(TargetAttractedInflation *attractedSheet);

    ISHEET_API void isheet_targetAttractedInflation_getCenterFixedVars(TargetAttractedInflation *attractedSheet, int **outCenterFixedVars, size_t *numCenterFixedVars);

    ISHEET_API void isheet_targetAttractedInflation_getVertexVars(TargetAttractedInflation *attractedSheet, int sheetIndex, int vertexIndex, int **outVars, size_t *numVars);

    ISHEET_API int isheet_targetAttractedInflation_getCenterNonFusedVertexIdx(TargetAttractedInflation *attractedSheet);

    ISHEET_API int isheet_targetAttractedInflation_getNumVars(TargetAttractedInflation *attractedSheet);

    ISHEET_API void isheet_targetAttractedInflation_getVars(TargetAttractedInflation *attractedSheet, double **outVars, size_t *numVars);

    ISHEET_API void isheet_targetAttractedInflation_getStrains(TargetAttractedInflation *attractedSheet, double **outStrains, size_t *numStrains);

    ISHEET_API double isheet_targetAttractedInflation_getEnergy(TargetAttractedInflation *attractedSheet);

    ISHEET_API double isheet_targetAttractedInflation_getGradientNorm(TargetAttractedInflation *attractedSheet);

    ISHEET_API void isheet_targetAttractedInflation_getGradient(TargetAttractedInflation *attractedSheet, double **outVars, size_t *numVars);

    ISHEET_API void isheet_targetAttractedInflation_getGravity(TargetAttractedInflation *attractedSheet, double **outVars, size_t *numVars);

    ISHEET_API double isheet_targetAttractedInflation_getReferenceVolume(TargetAttractedInflation *attractedSheet);

    ISHEET_API double isheet_targetAttractedInflation_getPressure(TargetAttractedInflation *attractedSheet);

    // Set methods
    ISHEET_API void isheet_targetAttractedInflation_setPressure(TargetAttractedInflation *attractedSheet, double pressure);

    ISHEET_API void isheet_targetAttractedInflation_setUseTensionFieldEnergy(TargetAttractedInflation *attractedSheet, int useTensionFieldEnergy);

    ISHEET_API void isheet_targetAttractedInflation_setUseHessianProjectedEnergy(TargetAttractedInflation *attractedSheet, int useHessianProjectedEnergy);

    ISHEET_API void isheet_targetAttractedInflation_disableFusedRegionTensionFieldTheory(TargetAttractedInflation *attractedSheet, int disableFusedRegionTensionFieldTheory);

    ISHEET_API void isheet_targetAttractedInflation_setVars(TargetAttractedInflation *attractedSheet, double *inVars, size_t numVars);
    
    ISHEET_API void isheet_targetAttractedInflation_setMassDensity(TargetAttractedInflation *attractedSheet, double rho);

    ISHEET_API void isheet_targetAttractedInflation_setGravity(TargetAttractedInflation *attractedSheet, double* inVector);

    ISHEET_API void isheet_targetAttractedInflation_setThickness(TargetAttractedInflation *attractedSheet, double thickness);

    ISHEET_API void isheet_targetAttractedInflation_setYoungModulus(TargetAttractedInflation *attractedSheet, double youngModulus);

    ISHEET_API void isheet_targetAttractedInflation_setReferenceVolume(TargetAttractedInflation *attractedSheet, double volume);

    // Solver
    ISHEET_API int isheet_targetAttractedInflation_newtonSolver(TargetAttractedInflation *sheet, int numSupports, int *supports, int numIterations, double gradTol, int writeReport, double **outReport, const char **errorMessage);
}
#endif