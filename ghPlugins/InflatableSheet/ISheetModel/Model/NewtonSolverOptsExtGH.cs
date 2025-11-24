using System;
using System.Collections.Generic;
using GH_IO.Serialization;
using Grasshopper.Kernel;
using ISheetModel.Helpers;
using ISheetModelLib.Types;

namespace ISheetModel.Model
{
    public class NewtonSolverOptsExtGH : GH_Component
    {
        int reportType;
        List<List<string>> reportAttributes;
        List<string> selection;
        bool buildAttributes = true;
        readonly List<string> categories = new List<string>(new string[] { "Convergence Report" });
        readonly List<string> reportContent = new List<string>(new string[]
        {
            "No Report",
            "Last-Step Report",
            "Step-By-Step Report"
        });

        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public NewtonSolverOptsExtGH()
      : base("SolverOptionsExt", "SolverOptionsExt",
            "Full solver options for inflatable sheets and periodic units.",
            "iSheet", "Model")
        {
        }

        public override void CreateAttributes()
        {
            if (buildAttributes)
            {
                FunctionToSetSelectedContent(0, 0);
                buildAttributes = false;
            }
            m_attributes = new DropDownAttributesGH(this, FunctionToSetSelectedContent, reportAttributes, selection, categories);
        }

        public void FunctionToSetSelectedContent(int dropdownListId, int selectedItemId)
        {
            if (reportAttributes == null)
            {
                reportAttributes = new List<List<string>>();
                selection = new List<string>();
                reportAttributes.Add(reportContent);
                selection.Add(reportContent[reportType]);
            }

            if (dropdownListId == 0)
            {
                reportType = selectedItemId;
                selection[0] = reportAttributes[0][selectedItemId];
            }

            Params.OnParametersChanged();
            ExpireSolution(true);
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddIntegerParameter("Iterations0", "Iter0", "Number of iterations with all supports.", GH_ParamAccess.item, 20);
            pManager.AddIntegerParameter("Iterations1", "Iter1", "Number of iterations without temporarys supports.", GH_ParamAccess.item, 0);
            pManager.AddNumberParameter("GradTol", "GradTol", "Set the gradient tolerance.", GH_ParamAccess.item, 1e-10);
            pManager.AddNumberParameter("HessianShiftForRigidMotion", "HessianShiftForRigidMotion", "Set Hessian shift for rigid motion.", GH_ParamAccess.item, 1e-10);
            pManager.AddNumberParameter("HessianShiftForAlphaInPlanar", "HessianShiftForAlphaInPlanar", "Set Hessian shift for alpha in planar.", GH_ParamAccess.item, 1e-12);
            pManager.AddBooleanParameter("UseTensionFieldEnergy", "UseTensionFieldEnergy", "Use Tension Field Energy.", GH_ParamAccess.item, true);
            pManager.AddBooleanParameter("UseHessianProjectedEnergy", "UseHessianProjectedEnergy", "Use Hessian Projected Energy.", GH_ParamAccess.item, false);
            pManager.AddBooleanParameter("DisableFusedRegionTensionFieldTheory", "DisableFusedRegionTensionFieldTheory", "Disable Fused Region Tension Field Theory", GH_ParamAccess.item, false);
            pManager.AddBooleanParameter("TwoStepsPeridicUnitSolve", "TwoStepsPeridicUnitSolve", "Enable two steps for the peridic unit solver", GH_ParamAccess.item, true);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("SolverOpts", "SolverOpts", "Newton solver options for inflatable sheets and periodic units.", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object can be used to retrieve data from input parameters and 
        /// to store data in output parameters.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            bool useTensionFieldEnergy = true, useHessianProjectedEnergy = false, disableFusedRegionTensionFieldTheory = false, twoSteps = true;
            double gradTol = 1e-8, hessianShiftForRigidMotion = 1e-10, hessianShiftForAlphaInPlanar = 1e-12;
            int subIter = 1, iter0 = 100, iter1=0;

            DA.GetData(0, ref iter0);
            DA.GetData(1, ref iter1);
            DA.GetData(2, ref gradTol);
            DA.GetData(3, ref hessianShiftForRigidMotion);
            DA.GetData(4, ref hessianShiftForAlphaInPlanar);
            DA.GetData(5, ref useTensionFieldEnergy);
            DA.GetData(6, ref useHessianProjectedEnergy);
            DA.GetData(7, ref disableFusedRegionTensionFieldTheory);
            DA.GetData(8, ref twoSteps);

            NewtonSolverOpts opts = new NewtonSolverOpts(iter0, subIter, iter1);
            opts.GradTol = gradTol;
            opts.HessianShiftForRigidMotion = hessianShiftForRigidMotion;
            opts.HessianShiftForAlphaInPlanar = hessianShiftForAlphaInPlanar;
            opts.UseTensionFieldEnergy = useTensionFieldEnergy;
            opts.UseHessianProjectedEnergy = useHessianProjectedEnergy;
            opts.DisableFusedRegionTensionFieldTheory = disableFusedRegionTensionFieldTheory;
            opts.ConvergenceReportType = reportType;
            opts.TwoStepsPeriodicUnitSolver = twoSteps;

            DA.SetData(0, opts);
        }

        public override GH_Exposure Exposure
        {
            get { return GH_Exposure.tertiary; }
        }

        public override bool Write(GH_IWriter writer)
        {
            writer.SetInt32("reportType", reportType);
            return base.Write(writer);
        }

        public override bool Read(GH_IReader reader)
        {
            if (reader.TryGetInt32("reportType", ref reportType))
            {
                FunctionToSetSelectedContent(0, reportType);
                m_attributes = new DropDownAttributesGH(this, FunctionToSetSelectedContent, reportAttributes, selection, categories);
            }
            return base.Read(reader);
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
            get { return new Guid("c567dc6f-f5e9-4909-a272-463b293ed04f"); }
        }
    }
}
