using System;
using System.Runtime.InteropServices;
using System.Security;

namespace ISheetModelLib.Creators
{
	public static partial class Kernel
	{
        private const string attractedSheet = "isheet_targetAttractedInflation_";

        public static class TargetAttractedInflation
		{
            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = attractedSheet + "build")]
            public static extern IntPtr Build(int numVertices, int numElements, [In] double[] inCoords, [In] int[] inElements, [In] int[] inFusedVertices, int numTargetVertices, int numTargetTrias, [In] double[] inTargetCoords, [In] int[] inTargetTrias, out IntPtr error);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = attractedSheet + "build_new")]
            public static extern IntPtr Build2(IntPtr sheet, int numTargetVertices, int numTargetTrias, [In] double[] inTargetCoords, [In] int[] inTargetTrias, out IntPtr error);

            #region GetMethods

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = attractedSheet + "getMeshVisualization")]
            public static extern void GetMeshVisualization(IntPtr attractedSheet, out IntPtr outCoords, out IntPtr outElements, out int numCoords, out int numElements);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = attractedSheet + "getMassDensity")]
            public static extern double GetMaterialDensity(IntPtr attractedSheet);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = attractedSheet + "getVolume")]
            public static extern double GetVolume(IntPtr attractedSheet);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = attractedSheet + "getYoungModulus")]
            public static extern double GetYoungModulus(IntPtr attractedSheet);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = attractedSheet + "getThickness")]
            public static extern double GetThickness(IntPtr attractedSheet);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = attractedSheet + "getCenterFixedVars")]
            public static extern void GetCenterFixedVars(IntPtr attractedSheet, out IntPtr outCenterFixedVars, out int numCenterFixedVars);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = attractedSheet + "getVertexVars")]
            public static extern void GetVertexVars(IntPtr attractedSheet, int sheetIndex, int vertexIndex, out IntPtr outVars, out int numVars);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = attractedSheet + "getCenterNonFusedVertexIdx")]
            public static extern int GetCenterNonFusedVertexIdx(IntPtr attractedSheet);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = attractedSheet + "getNumVars")]
            public static extern int GetNumVars(IntPtr attractedSheet);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = attractedSheet + "getVars")]
            public static extern void GetVars(IntPtr attractedSheet, out IntPtr outVars, out int numVars);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = attractedSheet + "getStrains")]
            public static extern void GetStrains(IntPtr attractedSheet, out IntPtr outStrains, out int numStrains);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = attractedSheet + "getEnergy")]
            public static extern double GetEnergy(IntPtr attractedSheet);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = attractedSheet + "getGradientNorm")]
            public static extern double GetGradientNorm(IntPtr attractedSheet);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = attractedSheet + "getGradient")]
            public static extern void GetGradient(IntPtr attractedSheet, out IntPtr outVars, out int numVars);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = attractedSheet + "getGravity")]
            public static extern void GetGravity(IntPtr attractedSheet, out IntPtr outVars, out int numVars);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = attractedSheet + "getReferenceVolume")]
            public static extern double GetReferenceVolume(IntPtr attractedSheet);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = attractedSheet + "getPressure")]
            public static extern double GetPressure(IntPtr attractedSheet);
            #endregion

            #region SetMethods

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = attractedSheet + "setPressure")]
            public static extern void SetPressure(IntPtr attractedSheet, double pressure);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = attractedSheet + "setUseTensionFieldEnergy")]
            public static extern void SetUseTensionFieldEnergy(IntPtr attractedSheet, int useTensionFieldEnergy);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = attractedSheet + "setUseHessianProjectedEnergy")]
            public static extern void SetUseHessianProjectedEnergy(IntPtr attractedSheet, int useHessianProjectedEnergy);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = attractedSheet + "disableFusedRegionTensionFieldTheory")]
            public static extern void DisableFusedRegionTensionFieldTheory(IntPtr attractedSheet, int disableFusedRegionTensionFieldTheory);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = attractedSheet + "setVars")]
            public static extern void SetVars(IntPtr attractedSheet, [In] double[] inVars, int numVars);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = attractedSheet + "setMassDensity")]
            public static extern void SetMaterialDensity(IntPtr attractedSheet, double rho);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = attractedSheet + "setGravity")]
            public static extern void SetGravity(IntPtr attractedSheet, [In] double[] inVector);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = attractedSheet + "setThickness")]
            public static extern void SetThickness(IntPtr ipu, double thickness);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = attractedSheet + "setYoungModulus")]
            public static extern void SetYoungModulus(IntPtr ipu, double youngModulus);

            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = attractedSheet + "setReferenceVolume")]
            public static extern void SetReferenceVolume(IntPtr ipu, double inVolume);
            #endregion

            #region Solver
            [SuppressUnmanagedCodeSecurity]
            [DllImport(isheet_dylib, CallingConvention = CallingConvention.StdCall, EntryPoint = attractedSheet + "newtonSolver")]
            public static extern int NewtonSolver(IntPtr attractedSheet, int numSupports, [In] int[] supports, int numIterations, double gradTol, int writeReport, out IntPtr outReport, out IntPtr errorMessage);

            #endregion
        }
    }
}

