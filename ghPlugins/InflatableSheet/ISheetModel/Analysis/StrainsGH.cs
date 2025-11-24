using System;
using System.Collections.Generic;
using System.Linq;
using GH_IO.Serialization;
using Grasshopper;
using Grasshopper.Kernel;
using ISheetDataLib.Types;
using ISheetModel.Helpers;
using ISheetModelLib.Types;
using Rhino.Geometry;
using static ISheetModelLib.Utils.ColorMaps;

namespace ISheetModel.Analysis
{
    public class StrainsGH : GH_Component
    {
        int cmapIdx;
        List<List<string>> menuAttributes;
        List<string> selection;
        bool buildAttributes = true;

        #region dropdownmenu content
        readonly List<string> categories = new List<string>(new string[] { "ColorMaps" });
        readonly List<string> cmapTypes = ((ColorMapTypes[])Enum.GetValues(typeof(ColorMapTypes))).Select(t => t.ToString()).ToList();
        #endregion

        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public StrainsGH()
          : base("Strains", "Strains",
                "Compute strains on a inflatable model.",
                "iSheet", "Analysis")
        {
        }

        public override void CreateAttributes()
        {
            if (buildAttributes)
            {
                FunctionToSetSelectedContent(0, 0);
                buildAttributes = false;
            }
            m_attributes = new DropDownAttributesGH(this, FunctionToSetSelectedContent, menuAttributes, selection, categories);
        }

        public void FunctionToSetSelectedContent(int dropdownListId, int selectedItemId)
        {
            if (menuAttributes == null)
            {
                menuAttributes = new List<List<string>>();
                selection = new List<string>();
                menuAttributes.Add(cmapTypes);
                selection.Add(cmapTypes[cmapIdx]);
            }

            if (dropdownListId == 0)
            {
                cmapIdx = selectedItemId;
                selection[0] = menuAttributes[0][selectedItemId];
            }

            Params.OnParametersChanged();
            ExpireSolution(true);
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Model", "Model", "Inflatable model", GH_ParamAccess.item);
            pManager.AddNumberParameter("Alpha", "Alpha", "Set the alpha value (from 0.0 to 1.0) to control the transparency of the visualization", GH_ParamAccess.item, 0.3);
            pManager.AddNumberParameter("LowerBound", "LowerBound", "Lower bound of the data set. If no value is explicitly provided, the minimum value of the data set is assumed.", GH_ParamAccess.item);
            pManager.AddNumberParameter("UpperBound", "UpperBound", "Upper bound of the data set. If no value is explicitly provided, the maximum value of the data set is assumed.", GH_ParamAccess.item);
            pManager[1].Optional = true;
            pManager[2].Optional = true;
            pManager[3].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("Mesh", "M", "Inflatable mesh.", GH_ParamAccess.item);
            pManager.AddNumberParameter("PerFace", "PerFace", "Strains computed per face", GH_ParamAccess.list);
            pManager.AddNumberParameter("PerVertex", "PerVertex", "Strains per vertex computed as the weigthed average strains per face.", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object can be used to retrieve data from input parameters and 
        /// to store data in output parameters.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            ElasticModel model = null;
            double alpha = 0.3, lowerBound = default, upperBound = default;
            if (!DA.GetData(0, ref model)) return;
            DA.GetData(1, ref alpha);
            DA.GetData(2, ref lowerBound);
            DA.GetData(3, ref upperBound);

            ColorMapTypes cmapType = ((ColorMapTypes[])Enum.GetValues(typeof(ColorMapTypes)))[cmapIdx];

            if (model == null)
            {
                this.AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Invalid Model.");
                return;
            }

            InflatableFields inflatable = new InflatableFields(model, cmapType, (int)(alpha * 255), lowerBound, upperBound);

            DA.SetData(0, inflatable.Mesh);
            DA.SetDataList(1, model.GetStrainsPerMeshFace());
            DA.SetDataList(2, model.GetStrainsPerMeshVertex());
        }

        public override bool Write(GH_IWriter writer)
        {
            writer.SetInt32("cmapIdx", cmapIdx);
            return base.Write(writer);
        }

        public override bool Read(GH_IReader reader)
        {
            if (reader.TryGetInt32("cmapIdx", ref cmapIdx))
            {
                FunctionToSetSelectedContent(0, cmapIdx);
                m_attributes = new DropDownAttributesGH(this, FunctionToSetSelectedContent, menuAttributes, selection, categories);
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
            get { return new Guid("6276a650-d727-4084-8bc1-aff39a97ef84"); }
        }
    }
}
