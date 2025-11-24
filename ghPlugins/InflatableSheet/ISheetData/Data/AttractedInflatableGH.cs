using System;
using System.Collections.Generic;

using Grasshopper;
using Grasshopper.Kernel;
using ISheetDataLib.Types;
using Rhino.FileIO;
using Rhino.Geometry;

namespace ISheetData.Data
{
    public class AttractedInflatableDataGH : GH_Component
    {
        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public AttractedInflatableDataGH()
          : base("AttractedInflatableData", "AttractedInflatableData",
            "Data set for initializing an inflatable model.",
            "Isheet", "Data")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("Mesh", "M", "Input triangular mesh.", GH_ParamAccess.item);
            pManager.AddPointParameter("Pts", "Pts", "Points where inflation sheet is fused.", GH_ParamAccess.list);
            pManager.AddGenericParameter("Supports", "Supports", "Support conditions.", GH_ParamAccess.list);
            pManager.AddNumberParameter("Pressure", "Pressure", "Set pressure", GH_ParamAccess.item);
            pManager.AddMeshParameter("Target","Target","Target mesh", GH_ParamAccess.item);
            pManager[2].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Data", "Data", "Attracted inflatable Data.", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object can be used to retrieve data from input parameters and 
        /// to store data in output parameters.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Mesh m = null, target= null;
            double pressure = 0;
            List<Point3d> fusedVertices = new List<Point3d>();
            List<Support> supports = new List<Support>();

            if (!DA.GetData(0, ref m)) return;
            DA.GetDataList(1, fusedVertices);
            DA.GetDataList(2, supports);
            DA.GetData(3, ref pressure);
            DA.GetData(4, ref target);

            if (m == null) return;

            TextLog log = new TextLog();
            MeshCheckParameters parameters = new MeshCheckParameters();
            parameters.CheckForBadNormals = true;
            parameters.CheckForDegenerateFaces = true;
            parameters.CheckForDisjointMeshes = true;
            parameters.CheckForDuplicateFaces = true;
            parameters.CheckForExtremelyShortEdges = true;
            parameters.CheckForInvalidNgons = true;
            parameters.CheckForNonManifoldEdges = true;
            parameters.CheckForRandomFaceNormals = true;
            parameters.CheckForSelfIntersection = true;
            parameters.CheckForUnusedVertices = true;
            m.Check(log, ref parameters);

            if (parameters.RandomFaceNormalCount > 0)
            {
                this.AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Mesh has " + parameters.RandomFaceNormalCount + " random face normals.");
                return;
            }
            if (parameters.InvalidNgonCount > 0)
            {
                this.AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Mesh has " + parameters.InvalidNgonCount + " invalid ngons.");
                return;
            }
            if (parameters.NonManifoldEdgeCount > 0)
            {
                this.AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Mesh has " + parameters.NonManifoldEdgeCount + " non manifold edges.");
                return;
            }
            if (parameters.NonUnitVectorNormalCount > 0)
            {
                this.AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Mesh has " + parameters.NonUnitVectorNormalCount + " non unit vector normals.");
                return;
            }
            if (parameters.ZeroLengthNormalCount > 0)
            {
                this.AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Mesh has " + parameters.ZeroLengthNormalCount + " zero length normals.");
                return;
            }
            if (parameters.UnusedVertexCount > 0)
            {
                this.AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Mesh has " + parameters.UnusedVertexCount + " unused vertex count.");
                return;
            }
            if (parameters.VertexFaceNormalsDifferCount > 0)
            {
                this.AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Mesh has " + parameters.VertexFaceNormalsDifferCount + " vertex normals that differ greatly from face normals.");
                return;
            }
            if (parameters.ExtremelyShortEdgeCount > 0)
            {
                this.AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Mesh has " + parameters.ExtremelyShortEdgeCount + " extremely short edges.");
                return;
            }
            if (parameters.DegenerateFaceCount > 0)
            {
                this.AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Mesh has " + parameters.DegenerateFaceCount + " degenerate faces.");
                return;
            }
            if (parameters.SelfIntersectingPairsCount > 0)
            {
                this.AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Mesh has " + parameters.SelfIntersectingPairsCount + " self-intersections.");
                return;
            }

            InflatableData sheet = new InflatableData(m, fusedVertices, false, target);
            sheet.AddSupports(supports);
            sheet.Pressure = pressure;


            DA.SetData(0, sheet);
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
            get { return new Guid("2836cc5a-aae4-4de7-a685-4b0c70904824"); }
        }
    }
}
