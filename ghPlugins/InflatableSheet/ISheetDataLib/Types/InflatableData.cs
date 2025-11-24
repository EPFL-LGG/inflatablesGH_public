using System;
using System.Collections.Generic;
using System.Linq;
using Rhino.Geometry;
using Newtonsoft.Json;
using System.IO;
using Grasshopper.Kernel.Types;
using GH_IO.Serialization;
using ISheetDataLib.Utils;
using static System.Windows.Forms.VisualStyles.VisualStyleElement.TreeView;

namespace ISheetDataLib.Types
{
    public enum ElasticBodyType { Single_layer_inflatable, Periodic_Unit, AttractedTargetInflation, Multi_layer_inflatable }

    public partial class InflatableData : ICloneable
    {
        public double[] Pressures { get; private set; }
        public int[] FusedVertices { get; set; }
        public double[][] Vertices { get; set; }
        public int[][] Faces { get; set; }
        public TargetSurfaceData TargetSurface { get; private set; }
        public ElasticBodyType ModelType { get; protected set; }
        public double Epsilon = 1e-9;
        public SupportCollection Supports { get; private set; }
        public int NumberSheets { get; private set; }
        public MaterialData Material { get; set; }

        private PointCloud _cloud { get; set;  }

        public InflatableData(InflatableData data)
        {
            Material = data.Material;
            Pressures = data.Pressures;
            FusedVertices = data.FusedVertices.ToArray();
            Vertices = data.Vertices.ToArray();
            Faces = data.Faces.ToArray();
            ModelType = data.ModelType;
            Supports = (SupportCollection) data.Supports.Clone();
            _cloud = (PointCloud) data._cloud.Duplicate();
            TargetSurface = data.TargetSurface==null ? null : (TargetSurfaceData) data.TargetSurface.Clone();
            NumberSheets = data.NumberSheets;
        }

        public InflatableData(Mesh mesh, IEnumerable<Point3d> fusedVertices = default, bool periodic=false, Mesh targetSurface=default)
		{
            // Default values
            Material = new MaterialData(300, 0.075, 0.0);
            Pressures = new double[] { 0.0 };
            NumberSheets = 2;

            // Store mesh values
            Vertices = new double[mesh.Vertices.Count][];
            for (int i=0; i<mesh.Vertices.Count; i++)
            {
                var p = mesh.Vertices[i];
                Vertices[i] = new double[] { p.X, p.Y, p.Z };
            }
            Faces = new int[mesh.Faces.Count][];
            for (int i = 0; i < mesh.Faces.Count; i++)
            {
                var f = mesh.Faces[i];
                Faces[i] = new int[] { f.A, f.B, f.C };
            }

            // Generate mesh data
            int numVertices, numTrias;
            double[] coords;
            int[] trias;
            Helpers.GetFlattenMeshData(mesh, out numVertices, out numTrias, out coords, out trias);

            FusedVertices = Enumerable.Repeat(0, numVertices).ToArray();
            _cloud = new PointCloud(mesh.Vertices.ToPoint3dArray());

            // Fused vertices  data
            // Step 1: naked boundaries
            ModelType = ElasticBodyType.Periodic_Unit;
            if (!periodic)
            {
                ModelType = ElasticBodyType.Single_layer_inflatable;

                var nakedBoundaries = mesh.GetNakedEdges();
                for (int i = 0; i < nakedBoundaries.Length; i++)
                {
                    var pts = nakedBoundaries[i];
                    foreach (var p in pts)
                    {
                        int idx = _cloud.ClosestPoint(p);
                        FusedVertices[idx] = 1;
                    }
                }
            }

            // Step 2: input vertices
            if (fusedVertices != default)
            {
                for (int i = 0; i < fusedVertices.Count(); i++)
                {
                    var p = fusedVertices.ElementAt(i);
                    int idx = _cloud.ClosestPoint(p);
                    if (_cloud[idx].Location.DistanceTo(p) < 1e-3)
                    {
                        FusedVertices[idx] = 1;
                    }
                }
            }

            IsValid = true;

            Supports = new SupportCollection();

            if (targetSurface == default) TargetSurface = null;
            else
            {
                ModelType = ElasticBodyType.AttractedTargetInflation;
                TargetSurface = new TargetSurfaceData(targetSurface);
            }
        }

        // Constructor for multilayer inflatbles
        public InflatableData(Mesh mesh, int numSheets, int[][] fusedVerticesIndices)
        {
            if (numSheets < 3) throw new Exception("Invalid number of sheets: a multilayer inflatable must have at least three sheets.");
            // Default values
            Material = new MaterialData(300, 0.075, 0.0);
            NumberSheets = numSheets;
            Pressures = new double[NumberSheets-1];

            // Store mesh values
            Vertices = new double[mesh.Vertices.Count][];
            for (int i = 0; i < mesh.Vertices.Count; i++)
            {
                var p = mesh.Vertices[i];
                Vertices[i] = new double[] { p.X, p.Y, p.Z };
            }
            Faces = new int[mesh.Faces.Count][];
            for (int i = 0; i < mesh.Faces.Count; i++)
            {
                var f = mesh.Faces[i];
                Faces[i] = new int[] { f.A, f.B, f.C };
            }

            // Generate mesh data
            int numVertices, numTrias;
            double[] coords;
            int[] trias;
            Helpers.GetFlattenMeshData(mesh, out numVertices, out numTrias, out coords, out trias);

            _cloud = new PointCloud(mesh.Vertices.ToPoint3dArray());

            // Fused vertices  data
            ModelType = ElasticBodyType.Multi_layer_inflatable;
            FusedVertices = Helpers.FlattenIntArray(fusedVerticesIndices);

            IsValid = true;

            Supports = new SupportCollection();
        }

        public object Clone()
        {
            return new InflatableData(this);
        }

        public Point3d GetVertex(int index)
        {
            return _cloud[index].Location;
        }

        public void SetTargetSurface(Mesh targetMesh, double weight = 0.0001)
        {
            if (ModelType != ElasticBodyType.Single_layer_inflatable) throw new Exception("Target surface assignment supported only for single-layer inflatables.");
            TargetSurface = new TargetSurfaceData(targetMesh, weight);
            ModelType = ElasticBodyType.AttractedTargetInflation;
        }

        public int GetCenterVertexIndex()
        {
            BoundingBox bb = new BoundingBox(_cloud.GetPoints());
            return _cloud.ClosestPoint(bb.Center);
        }

        public void AddSupport(Support support)
        {
            Point3d p = support.InitialPosition;

            int idx = _cloud.ClosestPoint(p);
            if (idx == -1) throw new Exception("Invalid inflatable sheet. Found 0 vertices.");
            else{
                support.Index = idx;
                Supports.Add(support);
            }
        }

        public void AddSupports(IEnumerable<Support> supports)
        {
            foreach (var sp in supports) AddSupport(sp);
        }

        public void CleanSupports()
        {
            Supports.Clear();
        }

        public void AddCentralSupport()
        {
            var p = _cloud.GetBoundingBox(false).Center;
            int idx = _cloud.ClosestPoint(p);
            Supports.Add(new Support(_cloud[idx].Location, new bool[] { true, true, true }, true, idx));
        }

        public void SetPressure(double pressure, int sheetLayer=0)
        {
            Pressures[sheetLayer] = pressure;
        }

        public void SetPressures(IEnumerable<double> pressures)
        {
            for (int i=0; i<NumberSheets-1; i++) Pressures[i] = i < pressures.Count() ? pressures.ElementAt(i) : pressures.ElementAt(0);
        }
    }
}

