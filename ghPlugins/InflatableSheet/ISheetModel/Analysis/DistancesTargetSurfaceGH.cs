using System;
using System.Collections.Generic;

using Grasshopper;
using Grasshopper.Kernel;
using Rhino.Geometry;
using System.Linq;
using System.Windows.Forms;
using System.Drawing;
using ISheetModelLib.Utils;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using ISheetModelLib.Types;
using static ISheetModelLib.Types.NewtonSolver;
using ISheetModel.Helpers;
using GH_IO.Serialization;
using static ISheetModelLib.Utils.ColorMaps;

namespace ISheetModel.Analysis
{
    public class DistancesTargetSurfaceGH : GH_Component
    {

        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public DistancesTargetSurfaceGH()
          : base("SurfaceDeviation", "SurfaceDeviation",
            "Calculate the distances to the target surface using the diagonal length of the bounding box.",
            "iSheet", "Analysis")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("InflatedMesh", "InflatedMesh", "Mesh of the inflated model. If multiple meshes are used then the target deviation will be computed using all the meshes.", GH_ParamAccess.list);
            pManager.AddMeshParameter("TargetMesh", "TargetMesh", "Target surface represented as a mesh.If the number of inflated meshes doesn't match the number of target meshes then the first target mesh will be used for running all calculations.", GH_ParamAccess.list);
            pManager.AddBooleanParameter("UseSqrDist", "UseSqrDist", "If true, it utilizes squared distances for fast computation; otherwise, it calculates exact distances.", GH_ParamAccess.item, true);
            pManager.AddBooleanParameter("UseMid", "UseMid", "If true, it utilizes mid surface to calculate distances.", GH_ParamAccess.item, true);
            pManager[2].Optional = true;
            pManager[3].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("Distance Mesh", "DistMesh", "Color-coded mesh with target surface deviation.", GH_ParamAccess.list);
            pManager.AddNumberParameter("Normalized Distances", "NormDistances", "Normalized distances to the target surface.", GH_ParamAccess.tree);
            pManager.AddNumberParameter("Distances", "Distances", "Distances to the target surface", GH_ParamAccess.tree);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object can be used to retrieve data from input parameters and 
        /// to store data in output parameters.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            List<Mesh> meshes = new List<Mesh>(), targets = new List<Mesh>();
            bool useSqr = true, useMid = true;
            if (!DA.GetDataList(0, meshes)) return;
            if (!DA.GetDataList(1, targets)) return;
            DA.GetData(2, ref useSqr);
            DA.GetData(3, ref useMid);

            Mesh mesh, target;
            int numMeshes = meshes.Count;
            int numTargets = targets.Count;

            List<Mesh> coloredMeshes = new List<Mesh>();
            GH_Structure<GH_Number> distancesTree = new GH_Structure<GH_Number>();
            GH_Structure<GH_Number> normDistancesTree = new GH_Structure<GH_Number>();
            List<double> allDistances = new List<double>();

            // Compute distances
            for (int i = 0; i < numMeshes; i++)
            {
                mesh = meshes[i];
                target = numMeshes == numTargets ? targets[i] : targets[0];
                mesh.UnifyNormals();
                mesh.Normals.ComputeNormals();

                int numVertices = mesh.Vertices.Count;
                int numHalfVertices = numVertices / 2;
                BoundingBox bb = mesh.GetBoundingBox(true);
                double diagLength = bb.Diagonal.Length;

                double[] distances = new double[numVertices];

                for (int j = 0; j < numHalfVertices; j++)
                {
                    Point3d p1 = mesh.Vertices[j];
                    Point3d p2 = mesh.Vertices[numHalfVertices + j];

                    if (useMid)
                    {
                        Point3d midP = ((p1 + p2) * 0.5);
                        distances[j] = useSqr ? target.ClosestPoint(midP).DistanceToSquared(midP) / diagLength : target.ClosestPoint(midP).DistanceTo(midP) / diagLength;
                        distances[numHalfVertices + j] = distances[j];
                    }
                    else
                    {
                        distances[j] = useSqr ? target.ClosestPoint(p1).DistanceToSquared(p1) / diagLength : target.ClosestPoint(p1).DistanceTo(p1) / diagLength;
                        distances[numHalfVertices + j] = useSqr ? target.ClosestPoint(p2).DistanceToSquared(p2) / diagLength : target.ClosestPoint(p2).DistanceTo(p2) / diagLength;
                    }
                }

                allDistances.AddRange(distances);
                distancesTree.AppendRange(distances.Select(d => new GH_Number(d)), new GH_Path(i));
            }


            // Compute normalized values and colormaps
            for (int i = 0; i < numMeshes; i++)
            {
                mesh = meshes[i];
                int numVertices = mesh.Vertices.Count;
                int numHalfVertices = numVertices / 2;

                // Normalized values
                double min = allDistances.Min();
                double max = allDistances.Max();
                double range = max - min;

                Color[] colorMap = ColorMaps.GetColorMap(ColorMapTypes.Viridis_Sequential, 255);
                double[] normalizedDistances = new double[numVertices];

                for (int j = 0; j < numVertices; j++)
                {
                    GH_Number dist = (GH_Number)distancesTree.get_Branch(i)[j];
                    double normalizedDist = (dist.Value - min) / range;
                    normalizedDistances[j] = normalizedDist;

                    int colorIndex = (int)(normalizedDist * (colorMap.Length - 1));
                    Color color = colorMap[colorIndex];

                    mesh.VertexColors.Add(color);
                }

                normDistancesTree.AppendRange(normalizedDistances.Select(d => new GH_Number(d)), new GH_Path(i));
                coloredMeshes.Add(mesh);
            }

            DA.SetDataList(0, coloredMeshes);
            DA.SetDataTree(1, normDistancesTree);
            DA.SetDataTree(2, distancesTree);
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
            get { return new Guid("c455c714-df26-4753-9aed-3c3e103871d7"); }
        }
    }
}
