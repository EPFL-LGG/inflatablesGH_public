using System;
using System.Runtime.InteropServices;
using System.Security;

namespace ISheetModelLib.Creators
{
	public partial class Kernel
	{
        private const string multiSheet = "isheet_multilayerInflatableSheet_";

        public static class MultipleSheets
        {
            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = multiSheet + "build")]
            public static extern IntPtr Build(int numVertices, int numElements, [In] double[] inCoords, [In] int[] inElements, [In]int numSheets, [In] double[] inPressures, [In] int[] inFusedVertices, out IntPtr errorMessage);

            #region SetMethods
            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = multiSheet + "setPressure")]
            public static extern void SetPressure(IntPtr sheet, [In] double[] inPressures, int numPressures);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = multiSheet + "setUseTensionFieldEnergy")]
            public static extern void SetUseTensionFieldEnergy(IntPtr sheet, int useTensionFieldEnergy);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = multiSheet + "setUseHessianProjectedEnergy")]
            public static extern void SetUseHessianProjectedEnergy(IntPtr sheet, int useHessianProjectedEnergy);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = multiSheet + "setVars")]
            public static extern void SetVars(IntPtr sheet, [In] double[] inVars, int numVars);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = multiSheet + "setGravity")]
            public static extern void SetGravity(IntPtr sheet, [In] double[] inCoords);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = multiSheet + "setMassDensity")]
            public static extern void SetMaterialDensity(IntPtr sheet, double rho);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = multiSheet + "setThickness")]
            public static extern void SetThickness(IntPtr sheet, double thickness);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = multiSheet + "setYoungModulus")]
            public static extern void SetYoungModulus(IntPtr sheet, double youngModulus);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = multiSheet + "setReferenceVolume")]
            public static extern void SetReferenceVolumes(IntPtr sheet, [In] double[] inVolumes, int numVolumes);

            #endregion

            #region GetMethods
            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = multiSheet + "getMeshVisualization")]
            public static extern IntPtr GetMeshVisualization(IntPtr sheet, out IntPtr outCoords, out IntPtr outElements, out int numCoords, out int numElements);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = multiSheet + "getVolumes")]
            public static extern void GetVolumes(IntPtr sheet, out IntPtr outVolumes, out int numVolumes);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = multiSheet + "getYoungModulus")]
            public static extern double GetYoungModulus(IntPtr sheet);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = multiSheet + "getThickness")]
            public static extern double GetThickness(IntPtr sheet);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = multiSheet + "getCenterFixedVars")]
            public static extern void GetCenterFixedVars(IntPtr sheet, out IntPtr outCenterFixedVars, out int numCenterFixedVars);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = multiSheet + "getCenterNonFusedVertexIdx")]
            public static extern int GetCenterNonFusedVertexIdx(IntPtr sheet);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = multiSheet + "getNumVars")]
            public static extern int GetNumVars(IntPtr sheet);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = multiSheet + "getVars")]
            public static extern int GetVars(IntPtr ipu, out IntPtr outVars, out int numVars);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = multiSheet + "getStrains")]
            public static extern void GetStrains(IntPtr ipu, out IntPtr outStrains, out int numStrains);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = multiSheet + "getEnergy")]
            public static extern double GetEnergy(IntPtr ipu);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = multiSheet + "getGradientNorm")]
            public static extern double GetGradientNorm(IntPtr ipu);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = multiSheet + "getGradient")]
            public static extern void GetGradient(IntPtr ipu, out IntPtr outFixedVars, out int numFixedVars);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = multiSheet + "getVertexVars")]
            public static extern void GetVertexVars(IntPtr sheet, [In] int sheetIdx, [In] int vertexIdx, out IntPtr outCenterFixedVars, out int numCenterFixedVars);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = multiSheet + "getGravity")]
            public static extern void GetGravity(IntPtr sheet, out IntPtr outCoords, out int numCoords);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = multiSheet + "getMassDensity")]
            public static extern double GetMaterialDensity(IntPtr sheet);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = multiSheet + "getReferenceVolume")]
            public static extern void GetReferenceVolumes(IntPtr sheet, out IntPtr outVolumes, out int numVolumes);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = multiSheet + "getPressure")]
            public static extern void GetPressures(IntPtr sheet, out IntPtr outPressures, out int numPressures);

            #endregion

            #region Solver
            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = multiSheet + "newtonSolver")]
            public static extern int NewtonSolver(IntPtr sheet, int numSupports, [In] int[] supports, int numIterations, double gradTol, int writeReport, out IntPtr outReport, out IntPtr errorMessage);

            #endregion
        }
    }
}

