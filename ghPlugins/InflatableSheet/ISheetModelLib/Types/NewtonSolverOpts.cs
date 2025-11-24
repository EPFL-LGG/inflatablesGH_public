using System;
using System.Collections.Generic;
using GH_IO.Serialization;
using Grasshopper.Kernel.Types;

namespace ISheetModelLib.Types
{
    public struct NewtonSolverOpts : IGH_Goo
    {
        public int NumIterationsWithoutTempSupports { get; set; }
        public int NumIterationsFull { get; set; }
        public int NumSubIterationsPerTimeStep { get; set; }
        public double GradTol { get; set; }
        public double HessianShiftForRigidMotion { get; set; }
        public double HessianShiftForAlphaInPlanar { get; set; }

        public bool UseTensionFieldEnergy { get; set; }
        public bool UseHessianProjectedEnergy { get; set; }
        public bool DisableFusedRegionTensionFieldTheory { get; set; }
        public bool UseBendingStiffnessWithBases { get; set; }

        public bool TwoStepsPeriodicUnitSolver { get; set; }

        // Report types: (0)"No Report", (1) "Last-Step Report", (3) "Step-By-Step Report"
        public int ConvergenceReportType { get; set; }

        public NewtonSolverOpts(int numIterations = 20, int numSubIterationsPerTimeStep = 5, int numIterationsWithoutTempSupports = 0)
        {
            GradTol = 1e-10;
            NumIterationsFull = numIterations;
            NumIterationsWithoutTempSupports = numIterationsWithoutTempSupports;
            NumSubIterationsPerTimeStep = numSubIterationsPerTimeStep;
            HessianShiftForRigidMotion = 1e-10;
            HessianShiftForAlphaInPlanar = 1e-12;

            UseTensionFieldEnergy = true;
            UseHessianProjectedEnergy = false;
            DisableFusedRegionTensionFieldTheory = false;
            ConvergenceReportType = 0;

            UseBendingStiffnessWithBases = true;
            TwoStepsPeriodicUnitSolver = false;
        }

        public NewtonSolverOpts(NewtonSolverOpts opts)
        {
            GradTol = opts.GradTol;
            NumIterationsFull = opts.NumIterationsFull;
            NumIterationsWithoutTempSupports = opts.NumIterationsWithoutTempSupports;
            NumSubIterationsPerTimeStep = opts.NumSubIterationsPerTimeStep;
            HessianShiftForRigidMotion = opts.HessianShiftForRigidMotion;
            HessianShiftForAlphaInPlanar = opts.HessianShiftForAlphaInPlanar;

            UseTensionFieldEnergy = opts.UseTensionFieldEnergy;
            UseHessianProjectedEnergy = opts.UseHessianProjectedEnergy;
            DisableFusedRegionTensionFieldTheory = opts.DisableFusedRegionTensionFieldTheory;
            ConvergenceReportType = opts.ConvergenceReportType;

            UseBendingStiffnessWithBases = opts.UseBendingStiffnessWithBases;
            TwoStepsPeriodicUnitSolver = opts.TwoStepsPeriodicUnitSolver;
        }

        public int GetCombinedNumberIterationsWithoutTemporarySupports()
        {
            return NumIterationsWithoutTempSupports * NumSubIterationsPerTimeStep;
        }

        public int GetCombinedNumberIterationsFull()
        {
            return NumIterationsFull * NumSubIterationsPerTimeStep;
        }

        public bool IsConvergenceReportActivated(bool isLastIteration = false)
        {
            return (ConvergenceReportType == 1 && isLastIteration) || ConvergenceReportType == 2 ? true : false;
        }

        public override string ToString()
        {
            return "SolverOptions";
        }

        #region GH methods
        public bool IsValid => true;

        public string IsValidWhyNot => "Not enough data has been provided";

        public string TypeName => ToString();

        public string TypeDescription => ToString();

        public IGH_Goo Duplicate()
        {
            return (IGH_Goo)this.MemberwiseClone();
        }

        public IGH_GooProxy EmitProxy()
        {
            return null;
        }

        public bool CastFrom(object source)
        {
            return false;
        }

        public bool CastTo<T>(out T target)
        {
            target = default(T);
            return false;
        }

        public object ScriptVariable()
        {
            return null;
        }

        public bool Write(GH_IWriter writer)
        {
            return false;
        }

        public bool Read(GH_IReader reader)
        {
            return false;
        }
        #endregion
    }
}

