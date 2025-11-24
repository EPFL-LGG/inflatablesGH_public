using System;
using System.Collections.Generic;

using Grasshopper;
using Grasshopper.Kernel;
using ISheetDataLib.Types;
using ISheetModelLib.Types;
using Rhino.Geometry;
using Rhino.Render;

namespace ISheetModel.Model
{
    public class AttractedInflatableSheetSolverGH : GH_Component
    {
        private bool run;
        private TargetAttractedInflation copy;
        private ConvergenceReport report;
        private NewtonSolverOpts opts;

        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public AttractedInflatableSheetSolverGH()
          : base("AttractedInflatableSolver", "AttractedSolver",
                "Newton solver for inflatable sheets.",
                "Isheet", "Model")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Model", "Model", "Attracted Inflatable sheet model", GH_ParamAccess.item);
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
            TargetAttractedInflation model = null;
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
            else if (model.ModelIO.ModelType == ElasticBodyType.PeriodicUnit)
            {
                this.AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Invalid Inflatable Sheet Model. The input model is a Periodic Unit.");
                return;
            }
            else
            {

                if (reset || copy == null)
                {
                    this.Message = "Reset";
                    copy = (TargetAttractedInflation) model.Clone();
                    report = new ConvergenceReport();
                }

                if (run)
                {
                    double[] dofs = copy.GetVars();
                    // Update positions of supports
                    foreach (Support sp in copy.ModelIO.Supports)
                    {
                        for (int i = 0; i < sp.IndicesDoFs.Length; i++)
                        {
                            int idx = sp.IndicesDoFs[i];
                            if (idx != -1) dofs[idx] = sp.TargetPosition[i];
                        }
                    }
                    copy.SetVars(dofs);

                    NewtonSolverOpts opts1 = new NewtonSolverOpts(opts);
                    opts1.NumSubIterationsPerTimeStep = opts.GetCombinedNumberIterationsFull();
                    NewtonSolver.InflatableSheetOptimizeMultipleSteps(copy, opts1, out report);
                }

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
            get { return new Guid("be93c471-54b9-40fa-8f29-448569217901"); }
        }
    }
}
