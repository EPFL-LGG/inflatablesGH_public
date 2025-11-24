using System;
using System.Collections.Generic;

using Grasshopper;
using Grasshopper.Kernel;
using ISheetDataLib.Types;
using Rhino.Geometry;

namespace ISheetData.Data
{
    public class SheetMaterialGH : GH_Component
    {
        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public SheetMaterialGH()
          : base("Material", "Mat",
            "Set material properties for the inflatable model.",
            "iSheet", "IO")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("E", "E", "Set the Young Modulus.", GH_ParamAccess.item);
            pManager.AddNumberParameter("Thickness", "Thickness", "Set the material thickness.", GH_ParamAccess.item);
            pManager.AddNumberParameter("Rho", "Rho", "Set the material density.", GH_ParamAccess.item, 0.00);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Material", "Material", "Material properties of the inflatable model.", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object can be used to retrieve data from input parameters and 
        /// to store data in output parameters.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            double E=0, thickness=0, rho=0;
            DA.GetData(0, ref E);
            DA.GetData(1, ref thickness);
            DA.GetData(2, ref rho);

            MaterialData mat = new MaterialData(E, thickness, rho);

            DA.SetData(0, mat);
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
            get { return new Guid("908c9dd2-5e8a-44df-8389-89274a93e747"); }
        }
    }
}
