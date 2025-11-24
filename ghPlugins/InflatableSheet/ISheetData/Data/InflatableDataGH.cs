using System;
using System.Collections.Generic;
using System.Linq;
using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using ISheetDataLib.Types;
using Rhino.FileIO;
using Rhino.Geometry;
using ISheetDataLib.Utils;
using System.Security.Policy;
using ISheetData.Helpers;
using GH_IO.Serialization;

namespace ISheetData.Data
{
    public class InflatableDataGH : GH_Component
    {
        int dataType;
        List<List<string>> dataAttributes;
        List<string> selection;
        bool buildAttributes = true;
        readonly List<string> categories = new List<string>(new string[] { "Inflatable Type" });
        readonly List<string> dataContent = new List<string>(new string[]
        {
            "Single-Layer",
            "Multi-layer",
            "Periodic-Unit",
        });

        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public InflatableDataGH()
          : base("Inflatable Data", "Inflatable Data",
            "Data for initializing the inflatable model.",
            "iSheet", "IO")
        {
        }

        public override void CreateAttributes()
        {
            if (buildAttributes)
            {
                FunctionToSetSelectedContent(0, 0);
                buildAttributes = false;
            }
            m_attributes = new DropDownAttributesGH(this, FunctionToSetSelectedContent, dataAttributes, selection, categories);
        }

        public void FunctionToSetSelectedContent(int dropdownListId, int selectedItemId)
        {
            if (dataAttributes == null)
            {
                dataAttributes = new List<List<string>>();
                selection = new List<string>();
                dataAttributes.Add(dataContent);
                selection.Add(dataContent[dataType]);
            }

            if (dropdownListId == 0)
            {
                dataType = selectedItemId;
                selection[0] = dataAttributes[0][selectedItemId];
            }

            Params.OnParametersChanged();
            ExpireSolution(true);
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("Mesh", "M", "Input triangular mesh.", GH_ParamAccess.item);
            pManager.AddPointParameter("Pts", "Pts", "Defines the points where the sheets are fused to create air channels. Multilayer inflatable structures are initialized using a tree structure, where each branch represents an interface between two consecutive sheets.The total number of sheets is equal to the number of interfaces plus one.", GH_ParamAccess.tree);
            pManager.AddGenericParameter("Supports", "Supports", "Support conditions.", GH_ParamAccess.list);
            pManager.AddNumberParameter("Pressure", "Pressure", "Inflation pressure applied to air channels. Multilayer sheets require a list of pressures, one per interface. If only one value is provided, it will be used for all channels.", GH_ParamAccess.list);
            pManager.AddGenericParameter("Material", "Material", "Material properties of the inflatable model.", GH_ParamAccess.item);
            pManager[2].Optional = true;
            pManager[4].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Data", "Data", "Multilayer Inflatable Data.", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object can be used to retrieve data from input parameters and 
        /// to store data in output parameters.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Mesh m = null;
            List<double> pressures = new List<double>();
            GH_Structure<GH_Point> fusingInterfaces = new GH_Structure<GH_Point>();
            List<Support> supports = new List<Support>();
            MaterialData mat = new MaterialData();

            if (!DA.GetData(0, ref m)) return;
            DA.GetDataTree(1, out fusingInterfaces);
            DA.GetDataList(2, supports);
            DA.GetDataList(3, pressures);
            DA.GetData(4, ref mat);

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

            if (supports.Count > 0)
            {
                bool tempSuppOnly = true;
                foreach (var sp in supports)
                {
                    if (!sp.IsTemporary) tempSuppOnly = false;
                }


                if (tempSuppOnly) this.AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Only temporary supports have been found. Please define at least one permanent.");
            }

            InflatableData sheet;

            // Check if it's single layered inflatble or a multi-layer inflatable
            int numSheets = fusingInterfaces.Branches.Count + 1;
            bool isMultilayer = numSheets > 2 ? true : false;



            if (dataType == 1)
            {
                if (!isMultilayer)
                {
                    this.AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Invalid input: A single branch of fused vertices detected for multiple-layers of air channels.");
                    return;
                }
                else
                {
                    // Fused vertices by interface (number of sheets is the number of interfaces + 1)
                    PointCloud cloud = new PointCloud(m.Vertices.ToPoint3dArray());
                    int numVertices = m.Vertices.Count;
                    int[][] fusedVerticesIndices = new int[numSheets][];
                    fusedVerticesIndices[0] = Enumerable.Range(0, numVertices).ToArray();

                    int nextIdx = numVertices;
                    for (int i = 1; i < numSheets; i++)
                    {
                        fusedVerticesIndices[i] = new int[numVertices];
                        PointCloud layerVertices = new PointCloud(fusingInterfaces.Branches[i - 1].Select(pt => pt.Value));

                        for (int j = 0; j < numVertices; j++)
                        {
                            var vtx = cloud[j].Location;
                            int fusedIdx = layerVertices.ClosestPoint(vtx);

                            // Check if the the current mesh vertex is a fused vertex.
                            // If it is, save the vertex index of the previous sheet layer. This means that the two sheets are fused.
                            // Otherwise, create a new vertex index. 
                            if (layerVertices[fusedIdx].Location.DistanceTo(vtx) < 1e-6) fusedVerticesIndices[i][j] = fusedVerticesIndices[i - 1][j];
                            else
                            {
                                fusedVerticesIndices[i][j] = nextIdx;
                                nextIdx++;
                            }
                        }
                    }

                    sheet = new InflatableData(m, numSheets, fusedVerticesIndices);
                    sheet.AddSupports(supports);
                    sheet.SetPressures(pressures);
                }
            }
            else
            {
                if (isMultilayer)
                {
                    this.AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Invalid input: Multiple branches of fused vertices detected for a single-layer air channel.");
                    return;
                }
                else
                {
                    sheet = new InflatableData(m, fusingInterfaces[0].Select(pt => pt.Value), dataType == 2 ? true : false);
                    sheet.AddSupports(supports);
                    sheet.SetPressure(pressures[0]);
                }
            }
            sheet.Material = mat;

            DA.SetData(0, sheet);
        }

        public override bool Write(GH_IWriter writer)
        {
            writer.SetInt32("dataType", dataType);
            return base.Write(writer);
        }

        public override bool Read(GH_IReader reader)
        {
            if (reader.TryGetInt32("dataType", ref dataType))
            {
                FunctionToSetSelectedContent(0, dataType);
                m_attributes = new DropDownAttributesGH(this, FunctionToSetSelectedContent, dataAttributes, selection, categories);
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
            get { return new Guid("9575b671-2212-497f-94cf-4e98bee211f7"); }
        }
    }
}