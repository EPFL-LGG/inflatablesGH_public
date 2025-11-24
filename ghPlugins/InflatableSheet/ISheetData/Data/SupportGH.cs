using System;
using System.Collections.Generic;

using Grasshopper;
using Grasshopper.Kernel;
using ISheetDataLib.Types;
using Rhino.Geometry;

namespace ISheetData.Data
{
    public class SupportGH : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the SupportGH class.
        /// </summary>
        public SupportGH()
          : base("Support", "Support",
              "Set the support conditions for the inflatable model.",
              "iSheet", "IO")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("Point", "Pt", "Initial position of the support.", GH_ParamAccess.item);
            pManager.AddBooleanParameter("X", "X", "Fix translation along the X-axis.", GH_ParamAccess.item, true);
            pManager.AddBooleanParameter("Y", "Y", "Fix translation along the Y-axis.", GH_ParamAccess.item, true);
            pManager.AddBooleanParameter("Z", "Z", "Fix translation along the Z-axis.", GH_ParamAccess.item, true);
            pManager.AddBooleanParameter("IsTopSheet", "IsTopSheet", "Set the support on the top or bottom sheet.", GH_ParamAccess.item, true);
            pManager.AddPointParameter("Target", "Target", "Target position of the support.", GH_ParamAccess.item);
            pManager.AddBooleanParameter("IsTemporary", "IsTemporary", "Set a temporary support.", GH_ParamAccess.item, false);
            pManager[5].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Support", "S", "Inflatable sheet support data.", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Point3d pos = new Point3d();
            Point3d target = new Point3d();
            bool[] lockedDoFs = new bool[3];
            bool isTopSheet = true, isTemporary=false;
            DA.GetData(0, ref pos);
            DA.GetData(1, ref lockedDoFs[0]);
            DA.GetData(2, ref lockedDoFs[1]);
            DA.GetData(3, ref lockedDoFs[2]);
            DA.GetData(4, ref isTopSheet);
            if(!DA.GetData(5, ref target)) target = pos;
            DA.GetData(6, ref isTemporary);

            Support support = new Support(pos, lockedDoFs, isTopSheet, -1, isTemporary);
            support.TargetPosition = target;

            DA.SetData(0, support);
        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return null;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("f9760659-b060-4cb0-88d0-0aec561a1c79"); }
        }
    }
}
