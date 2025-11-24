using System;
using System.Runtime.InteropServices;
using System.Security;

namespace ISheetModelLib.Creators
{
    public static partial class Kernel
    {
        private const string isheet_dylib = "libisheet.dylib";
        private const string inflatableSheet = "isheet_inflatableSheet_";

        public static class InflatableSheet
        {
            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = inflatableSheet + "build")]
            public static extern IntPtr Build(int numVertices, int numElements, [In] double[] inCoords, [In] int[] inElements, [In] int[] inFusedVertices);

            #region SetMethods
            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = inflatableSheet + "setPressure")]
            public static extern void SetPressure(IntPtr sheet, double pressure);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = inflatableSheet + "setUseTensionFieldEnergy")]
            public static extern void SetUseTensionFieldEnergy(IntPtr sheet, int useTensionFieldEnergy);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = inflatableSheet + "setUseHessianProjectedEnergy")]
            public static extern void SetUseHessianProjectedEnergy(IntPtr sheet, int useHessianProjectedEnergy);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = inflatableSheet + "disableFusedRegionTensionFieldTheory")]
            public static extern void DisableFusedRegionTensionFieldTheory(IntPtr sheet, int disableFusedRegionTensionFieldTheory);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = inflatableSheet + "setVars")]
            public static extern void SetVars(IntPtr sheet, [In] double[] inVars, int numVars);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = inflatableSheet + "setGravity")]
            public static extern void SetGravity(IntPtr sheet, [In] double[] inCoords);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = inflatableSheet + "setMassDensity")]
            public static extern void SetMaterialDensity(IntPtr sheet, double rho);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = inflatableSheet + "setThickness")]
            public static extern void SetThickness(IntPtr sheet, double thickness);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = inflatableSheet + "setYoungModulus")]
            public static extern void SetYoungModulus(IntPtr sheet, double youngModulus);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = inflatableSheet + "setReferenceVolume")]
            public static extern void SetReferenceVolume(IntPtr sheet, double inVolume);

            #endregion

            #region GetMethods
            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = inflatableSheet + "getMeshVisualization")]
            public static extern IntPtr GetMeshVisualization(IntPtr sheet, out IntPtr outCoords, out IntPtr outElements, out int numCoords, out int numElements);


            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = inflatableSheet + "getVolume")]
            public static extern double GetVolume(IntPtr sheet);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = inflatableSheet + "getYoungModulus")]
            public static extern double GetYoungModulus(IntPtr sheet);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = inflatableSheet + "getThickness")]
            public static extern double GetThickness(IntPtr sheet);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = inflatableSheet + "getCenterFixedVars")]
            public static extern void GetCenterFixedVars(IntPtr sheet, out IntPtr outCenterFixedVars, out int numCenterFixedVars);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = inflatableSheet + "getCenterNonFusedVertexIdx")]
            public static extern int GetCenterNonFusedVertexIdx(IntPtr sheet);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = inflatableSheet + "getNumVars")]
            public static extern int GetNumVars(IntPtr sheet);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = inflatableSheet + "getVars")]
            public static extern int GetVars(IntPtr ipu, out IntPtr outVars, out int numVars);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = inflatableSheet + "getStrains")]
            public static extern void GetStrains(IntPtr ipu, out IntPtr outStrains, out int numStrains);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = inflatableSheet + "getEnergy")]
            public static extern double GetEnergy(IntPtr ipu);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = inflatableSheet + "getGradientNorm")]
            public static extern double GetGradientNorm(IntPtr ipu);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = inflatableSheet + "getGradient")]
            public static extern void GetGradient(IntPtr ipu, out IntPtr outFixedVars, out int numFixedVars);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = inflatableSheet + "getVertexVars")]
            public static extern void GetVertexVars(IntPtr sheet, [In] int sheetIdx, [In] int vertexIdx, out IntPtr outCenterFixedVars, out int numCenterFixedVars);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = inflatableSheet + "getGravity")]
            public static extern void GetGravity(IntPtr sheet, out IntPtr outCoords, out int numCoords);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = inflatableSheet + "getMassDensity")]
            public static extern double GetMaterialDensity(IntPtr sheet);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = inflatableSheet + "getReferenceVolume")]
            public static extern double GetReferenceVolume(IntPtr sheet);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = inflatableSheet + "getPressure")]
            public static extern double GetPressure(IntPtr sheet);
            #endregion

            #region Solver
            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = inflatableSheet + "newtonSolver")]
            public static extern int NewtonSolver(IntPtr sheet, int numSupports, [In] int[] supports, int numIterations, double gradTol, int writeReport, out IntPtr outReport, out IntPtr errorMessage);

            #endregion
        }
    }
}
