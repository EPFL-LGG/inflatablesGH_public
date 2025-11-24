using System;
using System.Collections.Generic;

using Grasshopper;
using Grasshopper.Kernel;
using ISheetModelLib.Types;
using ISheetModelLib.Utils;
using Rhino.Geometry;

namespace ISheetModel.Analysis
{
    public class SplitInflatableGH : GH_Component
    {
        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public SplitInflatableGH()
          : base("SplitLayers", "SplitLayers",
            "Separates the inflatable model into meshes corresponding to each individual sheet layer.",
            "iSheet", "Analysis")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Model", "Model", "Inflated model to be split into individual sheet-layer meshes.", GH_ParamAccess.item);

        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("Sheets", "Sheets", "List of meshes defining the geometry of the individual sheets in the inflatable model.", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object can be used to retrieve data from input parameters and 
        /// to store data in output parameters.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            ElasticModel model = null;
            if (!DA.GetData(0, ref model)) return;

            int numSheets = model.ModelIO.NumberSheets;
            Mesh[] sheetMeshes = new Mesh[numSheets];
            Mesh mesh = model.VisualizationMesh;
            mesh.UnifyNormals();
            mesh.Normals.ComputeNormals();

            int numVertices = mesh.Vertices.Count;
            int numRefVertices = numVertices / numSheets;
            int numRefFaces = mesh.Faces.Count / numSheets;
            int idxV = 0;
            for (int i = 0; i < numSheets; i++)
            {
                Mesh m = new Mesh();
                // Vertices
                for (int j = 0; j < numRefVertices; j++) m.Vertices.Add(mesh.Vertices[idxV + j]);
                idxV += numRefVertices;

                // Faces
                for (int j = 0; j < numRefFaces; j++)
                {
                    MeshFace f1 = mesh.Faces[j];
                    m.Faces.AddFace(new MeshFace(f1.A, f1.B, f1.C, f1.D));
                }

                m.UnifyNormals();
                m.Normals.ComputeNormals();
                sheetMeshes[i] = m;
            }

            DA.SetDataList(0, sheetMeshes);
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
            get { return new Guid("99766385-95ea-4b24-8118-b90e6f2cdf26"); }
        }
    }
}
