using System;
using System.Collections.Generic;
using Grasshopper.Kernel.Data;
using Rhino.Geometry;

namespace ISheetDataLib.Utils
{
    public static class Helpers
    {
        public static double[] GetLinearSpaceDistribution(double min, double max, int resolution)
        {
            Interval domain = new Interval(min, max);
            double step = 1.0 / (resolution - 1);
            double[] result = new double[resolution];
            for (int i=0; i<resolution; i++) result[i] = domain.ParameterAt(step * i);

            return result;
        }

        public static Mesh GetQuadMesh(double[] coords, int[] quads)
        {
            Mesh m = new Mesh();
            int vCount = coords.Length / 3;
            for (int i = 0; i < vCount; i++)
            {
                m.Vertices.Add(new Point3d(coords[i * 3], coords[i * 3 + 1], coords[i * 3 + 2]));
            }

            int eCount = quads.Length / 4;
            for (int i = 0; i < eCount; i++)
            {
                m.Faces.AddFace(new MeshFace(quads[i * 4], quads[i * 4 + 1], quads[i * 4 + 2], quads[i * 4 + 3]));
            }
            m.UnifyNormals();

            return m;
        }

        public static Mesh GetTriasMesh(double[] coords, int[] trias)
        {
            Mesh m = new Mesh();
            int vCount = coords.Length / 3;
            for (int i = 0; i < vCount; i++)
            {
                m.Vertices.Add(new Point3d(coords[i * 3], coords[i * 3 + 1], coords[i * 3 + 2]));
            }

            int eCount = trias.Length / 3;
            for (int i = 0; i < eCount; i++)
            {
                m.Faces.AddFace(new MeshFace(trias[i * 3], trias[i * 3 + 1], trias[i * 3 + 2]));
            }
            m.Normals.ComputeNormals();

            return m;
        }

        public static void GetFlattenMeshData(Mesh m, out int numVertices, out int numTrias, out double[] outCoords, out int[] outTrias)
        {
            m.Faces.ConvertQuadsToTriangles();

            numVertices = m.Vertices.Count;
            outCoords = new double[numVertices * 3];
            for (int i = 0; i < numVertices; i++)
            {
                Point3d p = m.Vertices[i];
                outCoords[i * 3] = p.X;
                outCoords[i * 3 + 1] = p.Y;
                outCoords[i * 3 + 2] = p.Z;
            }

            numTrias = m.Faces.Count;
            outTrias = new int[numTrias * 3];
            for (int i = 0; i < numTrias; i++)
            {
                MeshFace f = m.Faces[i];
                outTrias[i * 3] = f.A;
                outTrias[i * 3 + 1] = f.B;
                outTrias[i * 3 + 2] = f.C;
            }
        }

        public static double[] FlattenDoubleArray(double[][] data)
        {
            List<double> list = new List<double>();
            for (int i = 0; i < data.Length; i++)
            {
                for (int j = 0; j < data[i].Length; j++)
                {
                    list.Add(data[i][j]);
                }
            }
            return list.ToArray();
        }

        public static int[] FlattenIntArray(int[][] data)
        {
            List<int> list = new List<int>();
            for (int i = 0; i < data.Length; i++)
            {
                for (int j = 0; j < data[i].Length; j++)
                {
                    list.Add(data[i][j]);
                }
            }
            return list.ToArray();
        }

        public static double[] FlattenPoint3dArray(Point3d[] data)
        {
            double[] list = new double[data.Length*3];
            for (int i = 0; i < data.Length; i++)
            {
                Point3d p = data[i];
                list[i * 3] = p.X;
                list[i * 3 + 1] = p.Y;
                list[i * 3 + 2] = p.Z;
            }
            return list;
        }

        public static int[] FlattenTriangularFacesArray(MeshFace[] trias)
        {
            int[] list = new int[trias.Length * 3];
            for (int i = 0; i < trias.Length; i++)
            {
                MeshFace f = trias[i];
                list[i * 3] = f.A;
                list[i * 3 + 1] = f.B;
                list[i * 3 + 2] = f.C;
            }
            return list;
        }

        public static T[][] Transpose<T>(T[][] input)
        {
            if (input == null || input.Length == 0) return Array.Empty<T[]>();

            int rowCount = input.Length;
            int colCount = input[0].Length;

            // Optional: Validate all rows have the same length
            for (int i = 1; i < rowCount; i++)
            {
                if (input[i].Length != colCount) throw new ArgumentException("All rows must have the same length.");
            }

            T[][] result = new T[colCount][];
            for (int i = 0; i < colCount; i++)
            {
                result[i] = new T[rowCount];
                for (int j = 0; j < rowCount; j++) result[i][j] = input[j][i];
            }

            return result;
        }

    }
}
