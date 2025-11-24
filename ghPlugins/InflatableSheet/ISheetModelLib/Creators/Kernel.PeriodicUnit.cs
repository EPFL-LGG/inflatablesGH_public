using System;
using System.Runtime.InteropServices;
using System.Security;

namespace ISheetModelLib.Creators
{
    public static partial class Kernel
    {
        private const string periodicUnit = "isheet_periodicUnit_";

        public static class PeriodicUnit
        {
            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = periodicUnit + "build")]
            public static extern IntPtr Build(int numVertices, int numElements, [In] double[] inCoords, [In] int[] inElements, [In] int[] inFusedVertices, [In] double epsilon);

            #region VoidMethods
            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = "isheet_periodicUnit_reparametrizeVerticalOffset")]
            public static extern void ReparametrizeVerticalOffset(IntPtr ipu);
            #endregion

            #region GetMethods
            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = periodicUnit + "getMeshVisualization")]
            public static extern IntPtr GetMeshVisualization(IntPtr ipu, out IntPtr outCoords, out IntPtr outElements, out int numCoords, out int numElements);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = periodicUnit + "getCenterFixedVars")]
            public static extern void GetCenterFixedVars(IntPtr ipu, out IntPtr outCenterFixedVars, out int numCenterFixedVars);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = periodicUnit + "getCenterNonFusedVertexIdx")]
            public static extern int GetCenterNonFusedVertexIdx(IntPtr ipu);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = periodicUnit + "getNumVars")]
            public static extern int GetNumVars(IntPtr ipu);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = periodicUnit + "getVars")]
            public static extern int GetVars(IntPtr ipu, out IntPtr outVars, out int numVars);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = periodicUnit + "getStrains")]
            public static extern void GetStrains(IntPtr ipu, out IntPtr outStrains, out int numStrains);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = periodicUnit + "getEnergy")]
            public static extern double GetEnergy(IntPtr ipu);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = periodicUnit + "getBendingStiffnessFixedVars")]
            public static extern void GetBendingStiffnessFixedVars(IntPtr ipu, out IntPtr outFixedVars, out int numFixedVars);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = periodicUnit + "getStretchingStiffnessFixedVars")]
            public static extern void GetStretchingStiffnessFixedVars(IntPtr ipu, out IntPtr outFixedVars, out int numFixedVars);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = periodicUnit + "getGradientNorm")]
            public static extern double GetGradientNorm(IntPtr ipu);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = periodicUnit + "getGradient")]
            public static extern void GetGradient(IntPtr ipu, out IntPtr outFixedVars, out int numFixedVars);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = periodicUnit + "getVertexVars")]
            public static extern void GetVertexVars(IntPtr ipu, [In] int sheetIdx, [In] int vertexIdx, out IntPtr outCenterFixedVars, out int numCenterFixedVars);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = periodicUnit + "getGravity")]
            public static extern void GetGravity(IntPtr ipu, out IntPtr outCoords, out int numCoords);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = periodicUnit + "getMassDensity")]
            public static extern double GetMaterialDensity(IntPtr ipu);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = periodicUnit + "getReferenceVolume")]
            public static extern double GetReferenceVolume(IntPtr ipu);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = periodicUnit + "getYoungModulus")]
            public static extern double GetYoungModulus(IntPtr sheet);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = periodicUnit + "getThickness")]
            public static extern double GetThickness(IntPtr sheet);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = periodicUnit + "getVolume")]
            public static extern double GetVolume(IntPtr sheet);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = periodicUnit + "getPressure")]
            public static extern double GetPressure(IntPtr sheet);
            #endregion

            #region SetMethods
            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = periodicUnit + "setPressure")]
            public static extern void SetPressure(IntPtr ipu, double pressure);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = periodicUnit + "setUseTensionFieldEnergy")]
            public static extern void SetUseTensionFieldEnergy(IntPtr ipu, int useTensionFieldEnergy);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = periodicUnit + "setUseHessianProjectedEnergy")]
            public static extern void SetUseHessianProjectedEnergy(IntPtr ipu, int useHessianProjectedEnergy);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = periodicUnit + "disableFusedRegionTensionFieldTheory")]
            public static extern void DisableFusedRegionTensionFieldTheory(IntPtr ipu, int disableFusedRegionTensionFieldTheory);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = periodicUnit + "setVars")]
            public static extern void SetVars(IntPtr ipu, [In] double[] inVars, int numVars);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = periodicUnit + "setGravity")]
            public static extern void SetGravity(IntPtr ipu, [In] double[] inCoords);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = periodicUnit + "setMassDensity")]
            public static extern void SetMaterialDensity(IntPtr ipu, double rho);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = periodicUnit + "setThickness")]
            public static extern void SetThickness(IntPtr ipu, double thickness);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = periodicUnit + "setYoungModulus")]
            public static extern void SetYoungModulus(IntPtr ipu, double youngModulus);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = periodicUnit + "setReferenceVolume")]
            public static extern void SetReferenceVolume(IntPtr ipu, double inVolume);

            #endregion


            #region Solver
            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = periodicUnit + "newtonSolver")]
            public static extern int NewtonSolver(IntPtr ipu, int numSupports, [In] int[] supports, int numIterations, double gradTol, double hessianShit, int writeReport, out IntPtr outReport, out IntPtr errorMessage);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = periodicUnit + "newtonStepSolver")]
            public static extern int NewtonStepSolver(IntPtr ipu, double pressure, int numSupports, [In] int[] supports, int numIterations, double gradTol, double hessianShit, int writeReport, out IntPtr outReport, out IntPtr errorMessage);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = periodicUnit + "computeBendingStiffness")]
            public static extern int ComputeBendingStiffness(IntPtr ipu, int numSupports, [In]int[] supports, int numAlphas, [In]double[] alphas, int numIterations, double gradTol, double hessianShift, int useBases, out int outNumBendingStiffness, out IntPtr outBendingStiffness, out int outNumStiffnessCoefficient, out IntPtr outStiffnessCoefficient, out IntPtr errorMessage);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = periodicUnit + "computeStretchingStiffness")]
            public static extern int ComputeStretchingStiffness(IntPtr ipu, int numSupports, [In] int[] supports, int numAlphas, [In] double[] alphas, int numIterations, double gradTol, double hessianShift, out int outNumStretchingStiffness, out IntPtr outStretchingStiffness, out IntPtr errorMessage);

            #endregion
        }
    }
}
