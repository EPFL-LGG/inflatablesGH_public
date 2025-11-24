using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using Eto.Forms;
using GH_IO.Serialization;
using Grasshopper.Kernel;
using ISheetModel.Helpers;
using ISheetModelLib.Utils;

namespace ISheetModel.Analysis
{
    public class GraphPlotterGH : GH_Component
    {
        int graphIdx;
        List<List<string>> graphAttributes;
        List<string> selection;
        bool buildAttributes = true;
        readonly List<string> categories = new List<string>(new string[] { "Graph" });
        readonly List<string> graphType = new List<string>(new string[]
        {
            "Bending",
            "Stretching",
        });

        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public GraphPlotterGH()
          : base("GraphPlotter", "GraphPlotter",
            "Graph plotter using matplot",
            "Isheet", "Analysis")
        {
        }

        public override void CreateAttributes()
        {
            if (buildAttributes)
            {
                FunctionToSetSelectedContent(0, 0);
                buildAttributes = false;
            }
            m_attributes = new DropDownAttributesGH(this, FunctionToSetSelectedContent, graphAttributes, selection, categories);
        }

        public void FunctionToSetSelectedContent(int dropdownListId, int selectedItemId)
        {
            if (graphAttributes == null)
            {
                graphAttributes = new List<List<string>>();
                selection = new List<string>();
                graphAttributes.Add(graphType);
                selection.Add(graphType[graphIdx]);
            }

            if (dropdownListId == 0)
            {
                graphIdx = selectedItemId;
                selection[0] = graphAttributes[0][selectedItemId];
            }

            Params.OnParametersChanged();
            ExpireSolution(true);
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("Data", "Data", "Data to plot", GH_ParamAccess.list);
            pManager.AddTextParameter("PyExec","PyExec","Python exectuable",GH_ParamAccess.item);
            pManager.AddTextParameter("PyFile","PyFile","Python file to run.",GH_ParamAccess.item);
            pManager.AddTextParameter("Path", "Path", "Path where the graph is going to be saved", GH_ParamAccess.item);
            pManager.AddTextParameter("Filename", "Filename", "Filename of the graph.", GH_ParamAccess.item);
            pManager.AddBooleanParameter("ShowGraph", "ShowGraph", "Show graph", GH_ParamAccess.item, false);
            pManager.AddBooleanParameter("Run", "Run", "Run graph.", GH_ParamAccess.item);
            pManager[5].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddTextParameter("File","File","Saved file.",GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object can be used to retrieve data from input parameters and 
        /// to store data in output parameters.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            List<double> data = new List<double>();
            bool run = false, show = false;
            string pyExec = "", pyFile="", path = "", filename = "";
            if(!DA.GetDataList(0, data)) return;
            DA.GetData(1, ref pyExec);
            DA.GetData(2, ref pyFile);
            DA.GetData(3, ref path);
            DA.GetData(4, ref filename);
            DA.GetData(5, ref show);
            DA.GetData(6, ref run);

            string dataString = PyCallback.StringFromDoubleArray(data);

            string[] txt = new string[0];
            if (run)
            {
                //GraphPlotter.PolarChart("stretching", data.ToArray());
                //PyCallback.CallScript(pyExec, pyFile, dataString, selection[0], path + filename, show.ToString());
            }

            string result = "File saved\n" + filename;
            DA.SetData(0, result);
        }

        public override bool Write(GH_IWriter writer)
        {
            writer.SetInt32("graphType", graphIdx);
            return base.Write(writer);
        }

        public override bool Read(GH_IReader reader)
        {
            if (reader.TryGetInt32("graphType", ref graphIdx))
            {
                FunctionToSetSelectedContent(0, graphIdx);
                m_attributes = new DropDownAttributesGH(this, FunctionToSetSelectedContent, graphAttributes, selection, categories);
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
            get { return new Guid("d4a48934-816e-4c4b-b3b8-a69e9763ea8f"); }
        }
    }
}
