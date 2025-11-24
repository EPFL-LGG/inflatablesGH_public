using System;
using System.Collections.Generic;
using System.Linq;
using Grasshopper;
using Grasshopper.Kernel;
using ISheetDataLib.Types;
using ISheetModelLib.Types;
using ISheetModelLib.Utils;

namespace ISheetModel.Analysis
{
    public class ComputeStiffnessGH : GH_Component
    {
        NewtonSolverOpts opts;
        StiffnessResults results;

        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public ComputeStiffnessGH()
          : base("Stiffness", "Stiffness",
            "Compute bending and stretching stiffness for periodic units with planar equilibrium states.",
            "iSheet", "Analysis")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Model", "Model", "Periodic unit model", GH_ParamAccess.item);
            pManager.AddGenericParameter("Opts", "Opts", "Newton solver options.", GH_ParamAccess.item);
            pManager.AddBooleanParameter("ShowPlots", "ShowPlots", "Generate graph plots", GH_ParamAccess.item, false);
            pManager[1].Optional = true;
            pManager[2].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddNumberParameter("Bending", "Bending", "Bending stiffness values", GH_ParamAccess.list);
            pManager.AddNumberParameter("BendingCoef", "BendingCoef", "Bending coefficients", GH_ParamAccess.list);
            pManager.AddNumberParameter("Stretching", "Stretching", "Stretching stiffness values (Only if no negative bending stiffness)", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object can be used to retrieve data from input parameters and 
        /// to store data in output parameters.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            PeriodicUnit ipu = null;
            bool show = false;
            if (!DA.GetData(0, ref ipu)) return;
            if (!DA.GetData(1, ref opts)) opts = new NewtonSolverOpts(100, 5);
            DA.GetData(2, ref show);

            if (ipu == null)
            {
                this.AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Invalid Periodic Unit Model.");
                return;
            }

            if (ipu.ModelIO.ModelType == ElasticBodyType.Single_layer_inflatable)
            {
                this.AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Invalid Periodic Unit Model. The input model is an single-layer inflatable.");
                return;
            }

            if (ipu.ModelIO.ModelType == ElasticBodyType.Multi_layer_inflatable)
            {
                this.AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Invalid Periodic Unit Model. The input model is an multi-layer inflatable.");
                return;
            }

            if (!ipu.HasPlanarEquilibrium)
            {
                this.AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "No planar equilibrium state found for periodic unit.");
                return;
            }

            NewtonSolver.PeriodicUnitComputeStiffness(ipu, opts, out results);
            bool negativeBendingStiffness = false;
            if (results.BendingStiffness.Min() < 0)
            {
                this.AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "Negative bending stiffness values were found. Stretching stiffness cannot be computed.");
                negativeBendingStiffness = true;
            }

            if (show)
            {
                List<double> bendingStiffnessResults = new List<double>(results.BendingStiffness);
                bendingStiffnessResults.AddRange(results.BendingStiffness);
                List<double> bendingSamplesResults = new List<double>(results.BendingThetaSamples);
                bendingSamplesResults.AddRange(ISheetDataLib.Utils.Helpers.GetLinearSpaceDistribution(Math.PI,2*Math.PI, 1000));

                GraphPlotter.PolarChart("Bending Stiffness", bendingStiffnessResults.ToArray(), bendingSamplesResults.ToArray());
                if(!negativeBendingStiffness) GraphPlotter.PolarChart("Stretching Stiffness", results.StretchingStiffness, results.StretchingThetaSamples);
            }

            DA.SetDataList(0, results.BendingStiffness);
            DA.SetDataList(1, results.BendingCoefficient);
            DA.SetDataList(2, results.StretchingStiffness);
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
            get { return new Guid("cb619b79-d4f4-4381-be15-02d026b90623"); }
        }
    }
}
