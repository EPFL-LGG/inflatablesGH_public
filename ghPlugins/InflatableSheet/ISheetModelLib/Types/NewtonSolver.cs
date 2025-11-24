using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using ISheetModelLib.Creators;
using ISheetDataLib.Types;
using ISheetDataLib.Utils;
using System.Linq;

namespace ISheetModelLib.Types
{
    public static class NewtonSolver
    {
        public enum SolverSteps { One_Step=0, Two_Steps=1 }

        public static bool InflationSolver(ElasticModel isheet, NewtonSolverOpts options, out ConvergenceReport report, bool updateMesh = true, bool isLastIteration=false, bool includeTemporary=true, double[] pressure=default)
        {
            int writeReport = options.IsConvergenceReportActivated(isLastIteration) ? 1 : 0;
            IntPtr ptrReport;

            int numSubIter = options.NumSubIterationsPerTimeStep;
            if (isLastIteration) numSubIter = options.NumSubIterationsPerTimeStep * 2;

            int[] supports = includeTemporary ? isheet.ModelIO.Supports.GetSupportsDoFsIndices() : isheet.ModelIO.Supports.GetNonTemporarySupportsDoFsIndices();
            int errorCode;
            if(pressure!=default) isheet.SetPressure(pressure);

            switch (isheet.ModelIO.ModelType)
            {
                case ElasticBodyType.Single_layer_inflatable:
                    errorCode = Kernel.InflatableSheet.NewtonSolver(isheet.ModelPtr, supports.Length, supports, numSubIter, options.GradTol, writeReport, out ptrReport, out isheet.ModelError);
                    break;
                case ElasticBodyType.Multi_layer_inflatable:
                    errorCode = Kernel.MultipleSheets.NewtonSolver(isheet.ModelPtr, supports.Length, supports, numSubIter, options.GradTol, writeReport, out ptrReport, out isheet.ModelError);
                    break;
                case ElasticBodyType.AttractedTargetInflation:
                    errorCode = Kernel.TargetAttractedInflation.NewtonSolver(isheet.ModelPtr, supports.Length, supports, numSubIter, options.GradTol, writeReport, out ptrReport, out isheet.ModelError);
                    break;
                default:
                    throw new Exception("Incompatible elastic model provided for Newton solver");
            }
            //if (isheet.ModelIO.ModelType == ElasticBodyType.Single_layer_inflatable) errorCode = Kernel.InflatableSheet.NewtonSolver(isheet.ModelPtr, supports.Length, supports, numSubIter, options.GradTol, writeReport, out ptrReport, out isheet.ModelError);
            //else if (isheet.ModelIO.ModelType == ElasticBodyType.Multi_layer_inflatable) errorCode = Kernel.MultipleSheets.NewtonSolver(isheet.ModelPtr, supports.Length, supports, numSubIter, options.GradTol, writeReport, out ptrReport, out isheet.ModelError);
            //else throw new Exception("Incompatible elastic model provided for Newton solver");

            if (errorCode == -1)
            {
                string errorMsg = isheet.ModelError.ToString();
                throw new Exception(errorMsg);
            }
            else
            {
                if (updateMesh) isheet.UpdateMeshVisualization();

                report = new ConvergenceReport();
                if (options.IsConvergenceReportActivated(isLastIteration))
                {
                    int numIter = includeTemporary ? options.NumIterationsFull : options.NumIterationsWithoutTempSupports;
                    int size = numIter * 4 + 2;
                    double[] data = new double[size];
                    Marshal.Copy(ptrReport, data, 0, size);

                    report = new ConvergenceReport(data, numIter);
                    Marshal.FreeCoTaskMem(ptrReport);
                }

                // Return true if the model converged
                if (errorCode == 1) return true;
                else return false;
            }
        }

        public static bool PeriodicUnitSolver(ElasticModel ipu, double pressure, NewtonSolverOpts options, out ExperimentLog log)
        {
            if (ipu.ModelIO.ModelType != ElasticBodyType.Periodic_Unit) throw new Exception("Invalid input: The inflatable model is not a periodic unit.");

            int writeReport = options.IsConvergenceReportActivated(true) ? 1 : 0;
            IntPtr ptrReport;

            ipu.SetUseTensionFieldEnergy(options.UseTensionFieldEnergy);
            ipu.SetUseHessianProjectedEnergy(options.UseHessianProjectedEnergy);
            ipu.DisableFusedRegionTensionFieldTheory(options.DisableFusedRegionTensionFieldTheory);
            int numVars = ipu.GetNumVars();
            int kappaIdx = numVars - 2;

            // Two-stage equilibrium solver before computing stiffness (fixing kappa)
            // First stage
            List<int> fixedVars = new List<int>() { kappaIdx };
            int errorCode = Kernel.PeriodicUnit.NewtonStepSolver(ipu.ModelPtr, pressure, fixedVars.Count, fixedVars.ToArray(), options.NumSubIterationsPerTimeStep, options.GradTol, options.HessianShiftForRigidMotion, writeReport, out ptrReport, out ipu.ModelError);
            if (errorCode == -1) throw new Exception(ipu.ModelError.ToString());
            double energyFirstStage = ipu.GetEnergy();

            //// Second stage 
            fixedVars = new List<int>(ipu.ModelIO.Supports.GetSupportsDoFsIndices()) { kappaIdx };
            if (options.TwoStepsPeriodicUnitSolver)
            {
                errorCode = Kernel.PeriodicUnit.NewtonStepSolver(ipu.ModelPtr, pressure, fixedVars.Count, fixedVars.ToArray(), options.NumSubIterationsPerTimeStep, options.GradTol, options.HessianShiftForAlphaInPlanar, writeReport, out ptrReport, out ipu.ModelError);
                if (errorCode == -1) throw new Exception(ipu.ModelError.ToString());
                double energySecondStage = ipu.GetEnergy();


                //// Check gradientNorm
                double gradNorm = ipu.GetGradientNorm();
                double kappaValue = ipu.GetVars()[kappaIdx];

                //// Log 
                log = new ExperimentLog(errorCode, kappaValue);
                log.GradientNorm = gradNorm;
                log.EnergyFirstStage = energyFirstStage;
                log.EnergySecondStage = energySecondStage;
                ipu.HasPlanarEquilibrium = log.PlanarEquilibrium;

                if (gradNorm != 0)
                {
                    // Planar equilibrium is not possible. Bending and stretching stiffness cannot be computed.
                    // Run an additional solver stage and terminate.
                    fixedVars = new List<int>(ipu.ModelIO.Supports.GetSupportsDoFsIndices());
                    errorCode = Kernel.PeriodicUnit.NewtonSolver(ipu.ModelPtr, fixedVars.Count, fixedVars.ToArray(), options.NumSubIterationsPerTimeStep, options.GradTol, options.HessianShiftForAlphaInPlanar, writeReport, out ptrReport, out ipu.ModelError);
                    if (errorCode == -1) throw new Exception(ipu.ModelError.ToString());
                    double energyThirdStage = ipu.GetEnergy();
                    log.EnergyThirdStage = energyThirdStage;
                }
            }
            else
            {
                //// Check gradientNorm
                double gradNorm = ipu.GetGradientNorm();
                double kappaValue = ipu.GetVars()[kappaIdx];

                log = new ExperimentLog(errorCode, kappaValue);
                log.GradientNorm = gradNorm;
                log.EnergyFirstStage = energyFirstStage;
            }

            ipu.UpdateMeshVisualization();

            if (errorCode == 1) return true;
            else return false;
        }

        public static bool PeriodicUnitComputeStiffness(PeriodicUnit ipu, NewtonSolverOpts options, out StiffnessResults results)
        {
            if (!ipu.HasPlanarEquilibrium)
            {
                results = new StiffnessResults();
                return false;
            }

            int numIterations = options.NumIterationsFull * options.NumSubIterationsPerTimeStep;
            results = new StiffnessResults();

            /////////////////////////////////////////////////////////////////////////////////////
            // First stage: Compute bending stiffness
            ipu.ReparametrizeVerticalOffset();
            int[] bendingStiffnessFixedVars = ipu.GetBendingStiffnessFixedVars();
            double[] bendingThetaSamples = Helpers.GetLinearSpaceDistribution(0, Math.PI, 1000);

            int numBendingStiffness, numStiffnessCoefficient;
            IntPtr bendingStiffnessPtr, stiffnessCoefficientPtr;
            int errorCode = Kernel.PeriodicUnit.ComputeBendingStiffness(ipu.ModelPtr, bendingStiffnessFixedVars.Length, bendingStiffnessFixedVars, bendingThetaSamples.Length, bendingThetaSamples, numIterations, options.GradTol, 0.0, options.UseBendingStiffnessWithBases ? 1 : 0, out numBendingStiffness, out bendingStiffnessPtr, out numStiffnessCoefficient, out stiffnessCoefficientPtr, out ipu.ModelError);
            if (errorCode == -1) throw new Exception(ipu.ModelError.ToString());

            // Parse results
            double[] bendingStiffness = new double[numBendingStiffness];
            Marshal.Copy(bendingStiffnessPtr, bendingStiffness, 0, numBendingStiffness);
            Marshal.FreeCoTaskMem(bendingStiffnessPtr);
            results.SetBendingStiffness(bendingStiffness, bendingThetaSamples);

            if (options.UseBendingStiffnessWithBases)
            {
                double[] stiffnessCoefficient = new double[numStiffnessCoefficient];
                Marshal.Copy(stiffnessCoefficientPtr, stiffnessCoefficient, 0, numStiffnessCoefficient);
                Marshal.FreeCoTaskMem(stiffnessCoefficientPtr);
                results.SetBendingCoefficient(stiffnessCoefficient);
            }

            /////////////////////////////////////////////////////////////////////////////////////
            // Second stage: Compute stretching stiffness (Only if no negative bending stiffness)
            if (bendingStiffness.Min() >= 0)
            {
                ipu.ReparametrizeVerticalOffset();
                int[] stretchingStiffnessFixedVars = ipu.GetStretchingStiffnessFixedVars();
                double[] stretchingThetaSamples = Helpers.GetLinearSpaceDistribution(0, 2*Math.PI, 1000);

                int numStretchingStiffness;
                IntPtr stretchingStiffnessPtr;
                errorCode = Kernel.PeriodicUnit.ComputeStretchingStiffness(ipu.ModelPtr, stretchingStiffnessFixedVars.Length, stretchingStiffnessFixedVars, stretchingThetaSamples.Length, stretchingThetaSamples, numIterations, options.GradTol, 0.0, out numStretchingStiffness, out stretchingStiffnessPtr, out ipu.ModelError);
                if (errorCode == -1) throw new Exception(ipu.ModelError.ToString());

                // Parse results
                double[] stretchingStiffness = new double[numStretchingStiffness];
                Marshal.Copy(stretchingStiffnessPtr, stretchingStiffness, 0, numStretchingStiffness);
                Marshal.FreeCoTaskMem(stretchingStiffnessPtr);
                results.SetStretchingStiffness(stretchingStiffness, stretchingThetaSamples);
            }
            

            if (errorCode == 1) return true;
            else return false;
        }
    }
}
