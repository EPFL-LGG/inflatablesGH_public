using System;
using System.Collections.Generic;
using Grasshopper;
using Grasshopper.Kernel;
using ISheetDataLib.Types;
using ISheetModelLib.Types;
using Rhino.Geometry;

namespace ISheetModel.Model
{
    public class EditInflatableModelGH : GH_Component
    {
        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public EditInflatableModelGH()
          : base("EditInflatable", "EditInflatable",
            "Modify an inflatable model.",
            "Isheet", "Model")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Linkage", "Linkage", "Linkage model to modify.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Supports", "Supports", "Set of support conditions.", GH_ParamAccess.list);
            pManager.AddBooleanParameter("CleanSupports", "CleanSupports", "Remove previous support conditions", GH_ParamAccess.item, true);
            pManager[1].Optional = true;
            pManager[2].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Linkage", "Linkage", "Elastic linkage model.", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object can be used to retrieve data from input parameters and 
        /// to store data in output parameters.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            ElasticModel model = null;
            List<Support> supports = new List<Support>();
            bool cleanSp = true;
            DA.GetData(0, ref model);
            DA.GetDataList(1, supports);
            DA.GetData(2, ref cleanSp);

            ElasticModel copy = (ElasticModel) model.Clone();

            if (cleanSp) copy.ModelIO.CleanSupports();
            copy.ModelIO.AddSupports(supports);
            copy.InitSupports();

            DA.SetData(0, copy);
        }

        public override GH_Exposure Exposure
        {
            get { return GH_Exposure.tertiary; }
        }

        /// <summary>
        /// Provides an Icon for every component that will be visible in the User Interface.
        /// Icons need to be 24x24 pixels.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
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
            get { return new Guid("a1a46587-02f2-4d90-ad1b-1cd6f24043f4"); }
        }
    }
}
