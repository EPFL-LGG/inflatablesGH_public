using System;
using System.Collections.Generic;

using Grasshopper;
using Grasshopper.Kernel;
using ISheetDataLib.Types;
using ISheetModelLib.Types;
using Rhino.Geometry;
using Rhino.Render;
using Rhino.UI;

namespace ISheetModel.Model
{
    public class NewtonSolverGH : GH_Component
    {
        private bool run;
        private ElasticModel copy;
        private ConvergenceReport report;
        private NewtonSolverOpts opts;
        private ExperimentLog log;

        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public NewtonSolverGH()
          : base("Inflation", "Inflation",
                "Newton solver for inflatable models.",
                "iSheet", "Model")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Model", "Model", "Inflatable sheet model", GH_ParamAccess.item);
            pManager.AddGenericParameter("Opts", "Opts", "Newton solver options.", GH_ParamAccess.item);
            pManager.AddBooleanParameter("Run", "Run", "Compute equilibrium.", GH_ParamAccess.item);
            pManager.AddBooleanParameter("Reset", "Reset", "Restart computation.", GH_ParamAccess.item);
            pManager[1].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Model", "Model", "Inflatable Model.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Report", "Report", "Convergence report", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object can be used to retrieve data from input parameters and 
        /// to store data in output parameters.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            ElasticModel model = null;
            bool reset = false;
            run = false;
            if (!DA.GetData(0, ref model)) return;
            if (!DA.GetData(1, ref opts)) opts = new NewtonSolverOpts(1, 100);
            DA.GetData(2, ref run);
            DA.GetData(3, ref reset);

            if (model == null)
            {
                this.AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Invalid Inflatable Sheet Model.");
                return;
            }
            else if (model.ModelIO.ModelType == ElasticBodyType.Periodic_Unit)
            {
                if (reset || copy == null)
                {
                    copy = (ElasticModel) model.Clone();
                    log = new ExperimentLog();
                    opts.NumSubIterationsPerTimeStep *= opts.NumIterationsFull;
                }

                if (run)
                {
                    NewtonSolver.PeriodicUnitSolver(copy, copy.ModelIO.Pressures[0], opts, out log);

                    if (opts.TwoStepsPeriodicUnitSolver == false) this.AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "[WARNING] Can not compute stiffness due to solver configuration. Please enable two-steps periodic unit solver!]");
                    else if (opts.TwoStepsPeriodicUnitSolver == true && !log.PlanarEquilibrium) this.AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "[WARNING] Can not compute stiffness due to non-planar equilibrium!]");
                }

                DA.SetData(0, copy);
                DA.SetData(1, log);
            }
            else
            {

                if (reset || copy == null)
                {
                    this.Message = "Reset";
                    copy = (ElasticModel) model.Clone();
                    report = new ConvergenceReport();
                }

                if (run)
                {
                    double[] dofs = copy.GetVars();
                    // Update positions of supports
                    foreach (Support sp in copy.ModelIO.Supports)
                    {
                        for (int i=0; i<sp.IndicesDoFs.Length; i++)
                        {
                            int idx = sp.IndicesDoFs[i];
                            if (idx!=-1)  dofs[idx] = sp.TargetPosition[i];
                        }
                    }
                    copy.SetVars(dofs);

                    NewtonSolverOpts opts1 = new NewtonSolverOpts(opts);
                    opts1.NumSubIterationsPerTimeStep = opts.GetCombinedNumberIterationsFull();
                    bool equilibrium = NewtonSolver.InflationSolver(copy, opts1, out report, true, false);

                    // Release temporary supports
                    if (opts.NumIterationsWithoutTempSupports > 0)
                    {
                        dofs = copy.GetVars();
                        // Update positions of supports
                        foreach (Support sp in copy.ModelIO.Supports)
                        {
                            if (sp.IsTemporary) continue;
                            for (int i = 0; i < sp.IndicesDoFs.Length; i++)
                            {
                                int idx = sp.IndicesDoFs[i];
                                if (idx != -1) dofs[idx] = sp.TargetPosition[i];
                            }
                        }
                        copy.SetVars(dofs);

                        NewtonSolverOpts opts2 = new NewtonSolverOpts(opts);
                        opts2.NumSubIterationsPerTimeStep = opts.GetCombinedNumberIterationsWithoutTemporarySupports();
                        equilibrium = NewtonSolver.InflationSolver(copy, opts2, out report, true, false, false);
                    }

                    if (equilibrium) this.Message = "Equilibrium";
                    else this.Message = "Equilibrium not found";
                }
                else this.Message = "Disabled";

                DA.SetData(0, copy);
                DA.SetData(1, report);
            }
        }

        public override GH_Exposure Exposure
        {
            get { return GH_Exposure.secondary; }
        }

        /// <summary>
        /// Provides an Icon for every component that will be visible in the User Interface.
        /// Icons need to be 24x24 pixels.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                // You can add image files to your project resources and access them like this:
                //return Resources.IconForThisComponent;
                return null;
            }
        }

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("40277234-743a-49ef-9c64-5d0a69991c30"); }
        }
    }
}
