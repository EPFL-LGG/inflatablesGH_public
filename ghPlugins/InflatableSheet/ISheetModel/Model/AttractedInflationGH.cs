using System;
using System.Collections.Generic;

using Grasshopper;
using Grasshopper.Kernel;
using ISheetDataLib.Types;
using ISheetModelLib.Types;
using Rhino.Geometry;

namespace ISheetModel.Model
{
    public class AttractedInflationGH : GH_Component
    {
        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public AttractedInflationGH()
          : base("Attracted Inflation", "Attracted Inflation",
            "Build an inflatable model using a target surface to guide the inflation.",
            "iSheet", "Model")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("M", "M", "Inflatable single layer model.", GH_ParamAccess.item);
            pManager.AddMeshParameter("Target Surface", "Target surface", "Target surface for attracted inflation.", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Model", "Model", "Inflatable model with target surface for attracted inflation.", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object can be used to retrieve data from input parameters and 
        /// to store data in output parameters.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            ElasticModel model = null;
            Mesh mesh = new Mesh();
            DA.GetData(0, ref model);
            DA.GetData(1, ref mesh);

            InflatableData io = (InflatableData) model.ModelIO.Clone();
            io.SetTargetSurface(mesh);

            TargetAttractedInflation aModel = new TargetAttractedInflation(io);

            DA.SetData(0, aModel);
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
            get { return new Guid("14090ff1-2ec8-4982-8a9a-09eca78fcb00"); }
        }
    }
}
