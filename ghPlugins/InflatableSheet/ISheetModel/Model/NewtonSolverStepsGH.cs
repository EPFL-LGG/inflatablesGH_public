using System;
using System.Collections.Generic;
using GH_IO.Serialization;
using Grasshopper;
using Grasshopper.Kernel;
using ISheetDataLib.Types;
using ISheetModelLib.Types;
using Rhino.Geometry;
using Rhino.UI;

namespace ISheetModel.Model
{
    public class InflatableSheetDynamicSolverGH : GH_Component
    {
        private bool run;
        private int steps = 1;
        private ElasticModel copy;
        private ConvergenceReport report;
        private NewtonSolverOpts opts;
        private ExperimentLog log;
        private int numIterations;
        double refStep = 0;
        bool includeTemporary;
        bool equilibrium;
        double[] refPressure, sumPressure;
        private List<Mesh> copies;

        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public InflatableSheetDynamicSolverGH()
          : base("Step Inflation", "Step Inflation",
            "Newton solver with multiple inflation steps.",
            "Isheet", "Model")
        {
            copies = new List<Mesh>();
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
            pManager.AddBooleanParameter("Save Steps", "Save Steps", "Save a copy of the deformed mesh at each time step.", GH_ParamAccess.item, false);
            pManager[1].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Model", "Model", "Inflatable Model.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Report", "Report", "Convergence report", GH_ParamAccess.item);
            pManager.AddMeshParameter("Meshes", "Meshes", "Copies of the inflatable at each time step.", GH_ParamAccess.list);
        }

        protected override void AfterSolveInstance()
        {
            if (run && (steps < numIterations))
            {
                GH_Document document = base.OnPingDocument();
                if (document != null)
                {
                    GH_Document.GH_ScheduleDelegate callback = new GH_Document.GH_ScheduleDelegate(this.ScheduleCallback);
                    document.ScheduleSolution(1, callback);
                }
            }
            else if (steps == numIterations)
            {
                if(equilibrium) this.Message = "Equilibrium Found";
                else this.Message = "Inflation Ended";
            }
            else this.Message = "Stop at step " + steps;
        }

        private void ScheduleCallback(GH_Document doc)
        {
            this.ExpireSolution(false);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object can be used to retrieve data from input parameters and 
        /// to store data in output parameters.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            ElasticModel model = null;
            bool reset = false, saveSteps = false;
            run = false;
            if(!DA.GetData(0, ref model)) return;
            if (!DA.GetData(1, ref opts)) opts = new NewtonSolverOpts(20,5);
            DA.GetData(2, ref run);
            DA.GetData(3, ref reset);
            DA.GetData(4, ref saveSteps);

            if (model == null)
            {
                this.AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Invalid Inflatable Sheet Model.");
                return;
            }
            else if(model.ModelIO.ModelType == ElasticBodyType.Periodic_Unit)
            {
                if (reset || copy == null)
                {
                    this.Message = "Reset";
                    copy = (ElasticModel)model.Clone();
                    copies = new List<Mesh>();
                    log = new ExperimentLog();
                    numIterations = opts.NumIterationsFull;

                    sumPressure = new double[] { 0 };
                    refPressure = new double[] { copy.ModelIO.Pressures[0] / (numIterations - 1) };
                    if(saveSteps) copies.Add(copy.VisualizationMesh.DuplicateMesh());
                    steps = 1;
                }

                if (run)
                {
                    if (steps <= numIterations)
                    {
                        this.Message = "Inflation Step " + steps;
                        sumPressure[0] += refPressure[0];
                        equilibrium = NewtonSolver.PeriodicUnitSolver(copy, sumPressure[0], opts, out log);
                        if (saveSteps) copies.Add(copy.VisualizationMesh.DuplicateMesh());
                        steps++;
                    }

                    if (opts.TwoStepsPeriodicUnitSolver == false) this.AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "[WARNING] Can not compute stiffness due to solver configuration. Please use two-steps solver!]");
                    else if (opts.TwoStepsPeriodicUnitSolver == true && !log.PlanarEquilibrium) this.AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "[WARNING] Can not compute stiffness due to non-planar equilibrium!]");
                }

                DA.SetData(0, copy);
                DA.SetData(1, log);
                DA.SetDataList(2, copies);
            }
            else
            {

                if (reset || copy == null)
                {
                    this.Message = "Reset";
                    copy = (ElasticModel) model.Clone();
                    report = new ConvergenceReport();
                    copies = new List<Mesh>();

                    numIterations = opts.NumIterationsFull + opts.NumIterationsWithoutTempSupports;
                    refStep = 1.0 / (opts.NumIterationsFull - 1);
                    steps = 1;
                    includeTemporary = true;
                    //equilibrium = false;

                    sumPressure = new double[copy.ModelIO.NumberSheets-1];
                    refPressure = new double[copy.ModelIO.NumberSheets-1];
                    for (int i=0; i< copy.ModelIO.NumberSheets-1; i++)
                    {
                        sumPressure[i] = 0.0;
                        refPressure[i] = copy.ModelIO.Pressures[i] / (numIterations - 1);
                    }
                    if (saveSteps) copies.Add(copy.VisualizationMesh.DuplicateMesh());

                }

                if (run)
                {
                
                    if (steps < numIterations)
                    {
                        if (steps < opts.NumIterationsFull) includeTemporary = true;
                        else includeTemporary = false;

                        this.Message = "Inflation Step " + steps;

                        double[] dofs = copy.GetVars();
                        // Update positions of supports
                        foreach(Support sp in copy.ModelIO.Supports)
                        {
                            if (sp.IsTemporary && !includeTemporary) continue;

                            // Compute linear interpolation between initial position and target position
                            Line ln = new Line(sp.InitialPosition, sp.TargetPosition);
                            double t = refStep * (steps - 1);
                            var tPos = ln.PointAt(t>1 ? 1 : t);
                            for (int i = 0; i < sp.IndicesDoFs.Length; i++)
                            {
                                int idx = sp.IndicesDoFs[i];
                                if (idx != -1) dofs[idx] = tPos[i];
                            }
                        }
                        copy.SetVars(dofs);

                        for (int i = 0; i < copy.ModelIO.NumberSheets-1; i++) sumPressure[i] += refPressure[i];

                        equilibrium = NewtonSolver.InflationSolver(copy, opts, out report, true, true, includeTemporary, sumPressure);
                        if (saveSteps) copies.Add(copy.VisualizationMesh.DuplicateMesh());

                        report.IterationStep = steps;
                        steps++;
                    }
                }
                else this.Message = "Disabled";

                DA.SetData(0, copy);
                DA.SetData(1, report);
                DA.SetDataList(2, copies);
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
            get { return new Guid("9f9cea46-4fb0-4390-b9cf-4508bb2448e0"); }
        }
    }
}
