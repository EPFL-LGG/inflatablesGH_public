using System;
using System.Collections.Generic;

using Grasshopper;
using Grasshopper.Kernel;
using ISheetDataLib.Types;
using Rhino.Geometry;

namespace ISheetData.Data
{
  public class ISheetWriteDataGH : GH_Component
  {
        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public ISheetWriteDataGH()
          : base("Write Data", "Write Data",
            "Write a Json file.",
            "iSheet", "IO")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Data", "Data", "Data to write.", GH_ParamAccess.item);
            pManager.AddTextParameter("Path", "Path", "Directory path.", GH_ParamAccess.item);
            pManager.AddTextParameter("Filename", "Filename", "Name of the file.", GH_ParamAccess.item);
            pManager.AddBooleanParameter("Write", "Write", "Write file.", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddTextParameter("File", "File", "Json file.", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object can be used to retrieve data from input parameters and 
        /// to store data in output parameters.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            string path = "";
            string filename = "";
            bool write = false;
            InflatableData data = null;
            if(!DA.GetData(0, ref data)) return;
            DA.GetData(1, ref path);
            DA.GetData(2, ref filename);
            DA.GetData(3, ref write);

            if (data == null) this.AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Null data");

            if (write) data.WriteJsonFile(path, filename);

            path += filename;

            DA.SetData(0, path);
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
      get { return new Guid("41e0e342-4d4b-4dac-9302-6b560aef3629"); }
    }
  }
}
