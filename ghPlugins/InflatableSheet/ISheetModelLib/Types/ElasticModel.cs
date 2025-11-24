using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using Eto.Drawing;
using GH_IO.Serialization;
using Grasshopper.Kernel.Types;
using ISheetDataLib.Types;
using ISheetDataLib.Utils;
using ISheetModelLib.Creators;
using Plotly.NET.TraceObjects;
using Rhino.Geometry;

namespace ISheetModelLib.Types
{
    public abstract partial class ElasticModel : IGH_Goo, ICloneable
    {
        public bool HasPlanarEquilibrium { get; set; } // Flag for periodic units
        public Mesh VisualizationMesh { get; protected set; }
        public InflatableData ModelIO { get; protected set; }
        public IntPtr ModelPtr, ModelError;

        private string MSG_ERROR_MODEL = "Unknown inflatable model";

        public abstract object Clone();

        protected void Init(InflatableData modelData)
        {
            ModelIO = new InflatableData(modelData);
            VisualizationMesh = new Mesh();

            // Parse Data
            int numVertices = ModelIO.Vertices.Length, numTrias = ModelIO.Faces.Length;
            double[] coords = Helpers.FlattenDoubleArray(ModelIO.Vertices);
            int[] trias = Helpers.FlattenIntArray(ModelIO.Faces);

            // Build Inflatable Model
            if (ModelIO.ModelType == ElasticBodyType.Single_layer_inflatable)
            {
                ModelPtr = Kernel.InflatableSheet.Build(numVertices, numTrias, coords, trias, ModelIO.FusedVertices);
            }
            else if (ModelIO.ModelType == ElasticBodyType.Periodic_Unit)
            {
                ModelPtr = Kernel.PeriodicUnit.Build(numVertices, numTrias, coords, trias, ModelIO.FusedVertices, ModelIO.Epsilon);
            }
            else if (ModelIO.ModelType == ElasticBodyType.Multi_layer_inflatable)
            {
                ModelPtr = Kernel.MultipleSheets.Build(numVertices, numTrias, coords, trias, ModelIO.NumberSheets, ModelIO.Pressures, ModelIO.FusedVertices, out ModelError);
            }
            else if (ModelIO.ModelType == ElasticBodyType.AttractedTargetInflation)
            {
                int numTargetVertices = ModelIO.TargetSurface.Vertices.Length, numTargetTrias = ModelIO.TargetSurface.Trias.Length;
                double[] targetCoords = Helpers.FlattenDoubleArray(ModelIO.TargetSurface.Vertices);
                int[] targetTrias = Helpers.FlattenIntArray(ModelIO.TargetSurface.Trias);
                ModelPtr = Kernel.TargetAttractedInflation.Build(numVertices, numTrias, coords, trias, ModelIO.FusedVertices, numTargetVertices, numTargetTrias, targetCoords, targetTrias, out ModelError);
                if (ModelPtr == IntPtr.Zero) throw new Exception(Marshal.PtrToStringAuto(ModelError));
            }
            else throw new Exception(MSG_ERROR_MODEL);

            if (ModelPtr == IntPtr.Zero) throw new Exception(Marshal.PtrToStringAuto(ModelError));

            SetPressure(ModelIO.Pressures);
            SetYoungModulus(ModelIO.Material.E);
            SetThickness(ModelIO.Material.Thickness);
            SetMaterialDensity(ModelIO.Material.MaterialDensity);

            // Define support conditions
            InitSupports();
            // Initialize meshes
            InitMeshVisualization();
        }

        public void InitSupports()
        {
            if (ModelIO.Supports.Count == 0)
            {
                var index = ModelIO.GetCenterVertexIndex();
                var pos = ModelIO.GetVertex(index);
                ModelIO.AddSupport(new Support(pos, new bool[] { true, true, true }, true, index));
            }
            foreach (var sp in ModelIO.Supports) sp.SetDoFsIndices(GetVertexVars(sp.IsTopSheet, sp.Index));
        }

        public int GetNumVars()
        {
            if (ModelIO.ModelType == ElasticBodyType.Single_layer_inflatable) return Kernel.InflatableSheet.GetNumVars(ModelPtr);
            else if (ModelIO.ModelType == ElasticBodyType.Periodic_Unit) return Kernel.PeriodicUnit.GetNumVars(ModelPtr);
            else if (ModelIO.ModelType == ElasticBodyType.Multi_layer_inflatable) return Kernel.MultipleSheets.GetNumVars(ModelPtr);
            else if (ModelIO.ModelType == ElasticBodyType.AttractedTargetInflation) return Kernel.TargetAttractedInflation.GetNumVars(ModelPtr);
            else throw new Exception(MSG_ERROR_MODEL);
        }

        public double[] GetVars()
        {
            int numVars;
            IntPtr varsPtr;

            if (ModelIO.ModelType == ElasticBodyType.Single_layer_inflatable) Kernel.InflatableSheet.GetVars(ModelPtr, out varsPtr, out numVars);
            else if (ModelIO.ModelType == ElasticBodyType.Periodic_Unit) Kernel.PeriodicUnit.GetVars(ModelPtr, out varsPtr, out numVars);
            else if (ModelIO.ModelType == ElasticBodyType.Multi_layer_inflatable) Kernel.MultipleSheets.GetVars(ModelPtr, out varsPtr, out numVars);
            else if (ModelIO.ModelType == ElasticBodyType.AttractedTargetInflation) Kernel.TargetAttractedInflation.GetVars(ModelPtr, out varsPtr, out numVars);
            else throw new Exception(MSG_ERROR_MODEL);

            double[] vars = new double[numVars];
            Marshal.Copy(varsPtr, vars, 0, numVars);
            Marshal.FreeCoTaskMem(varsPtr);

            return vars;
        }

        public void SetVars(double[] vars)
        {
            if (ModelIO.ModelType == ElasticBodyType.Single_layer_inflatable) Kernel.InflatableSheet.SetVars(ModelPtr, vars, vars.Length);
            else if (ModelIO.ModelType == ElasticBodyType.Periodic_Unit) Kernel.PeriodicUnit.SetVars(ModelPtr, vars, vars.Length);
            else if (ModelIO.ModelType == ElasticBodyType.Multi_layer_inflatable) Kernel.MultipleSheets.SetVars(ModelPtr, vars, vars.Length);
            else if (ModelIO.ModelType == ElasticBodyType.AttractedTargetInflation) Kernel.TargetAttractedInflation.SetVars(ModelPtr, vars, vars.Length);
            else throw new Exception(MSG_ERROR_MODEL);
        }

        public double GetEnergy()
        {
            if (ModelIO.ModelType == ElasticBodyType.Single_layer_inflatable) return Kernel.InflatableSheet.GetEnergy(ModelPtr);
            else if (ModelIO.ModelType == ElasticBodyType.Periodic_Unit) return Kernel.PeriodicUnit.GetEnergy(ModelPtr);
            else if (ModelIO.ModelType == ElasticBodyType.Multi_layer_inflatable) return Kernel.MultipleSheets.GetEnergy(ModelPtr);
            else if (ModelIO.ModelType == ElasticBodyType.AttractedTargetInflation) return Kernel.TargetAttractedInflation.GetEnergy(ModelPtr);
            else throw new Exception(MSG_ERROR_MODEL);
        }

        public double GetGradientNorm()
        {
            if (ModelIO.ModelType == ElasticBodyType.Single_layer_inflatable) return Kernel.InflatableSheet.GetGradientNorm(ModelPtr);
            else if (ModelIO.ModelType == ElasticBodyType.Periodic_Unit) return Kernel.PeriodicUnit.GetGradientNorm(ModelPtr);
            else if (ModelIO.ModelType == ElasticBodyType.Multi_layer_inflatable) return Kernel.MultipleSheets.GetGradientNorm(ModelPtr);
            else if (ModelIO.ModelType == ElasticBodyType.AttractedTargetInflation) return Kernel.TargetAttractedInflation.GetGradientNorm(ModelPtr);
            else throw new Exception(MSG_ERROR_MODEL);
        }

        public double[] GetGradient()
        {
            int numVars;
            IntPtr varsPtr;

            if (ModelIO.ModelType == ElasticBodyType.Single_layer_inflatable) Kernel.InflatableSheet.GetGradient(ModelPtr, out varsPtr, out numVars);
            else if (ModelIO.ModelType == ElasticBodyType.Periodic_Unit) Kernel.PeriodicUnit.GetGradient(ModelPtr, out varsPtr, out numVars);
            else if (ModelIO.ModelType == ElasticBodyType.Multi_layer_inflatable) Kernel.MultipleSheets.GetGradient(ModelPtr, out varsPtr, out numVars);
            else if (ModelIO.ModelType == ElasticBodyType.AttractedTargetInflation) Kernel.TargetAttractedInflation.GetGradient(ModelPtr, out varsPtr, out numVars);
            else throw new Exception(MSG_ERROR_MODEL);

            double[] gradient = new double[numVars];
            Marshal.Copy(varsPtr, gradient, 0, numVars);
            Marshal.FreeCoTaskMem(varsPtr);
            return gradient;
        }

        public double[] GetStrainsPerMeshFace()
        {
            int numStrains;
            IntPtr strainsPtr;

            if (ModelIO.ModelType == ElasticBodyType.Single_layer_inflatable) Kernel.InflatableSheet.GetStrains(ModelPtr, out strainsPtr, out numStrains);
            else if (ModelIO.ModelType == ElasticBodyType.Periodic_Unit) Kernel.PeriodicUnit.GetStrains(ModelPtr, out strainsPtr, out numStrains);
            else if (ModelIO.ModelType == ElasticBodyType.Multi_layer_inflatable) Kernel.MultipleSheets.GetStrains(ModelPtr, out strainsPtr, out numStrains);
            else if (ModelIO.ModelType == ElasticBodyType.AttractedTargetInflation) Kernel.TargetAttractedInflation.GetStrains(ModelPtr, out strainsPtr, out numStrains);
            else throw new Exception(MSG_ERROR_MODEL);

            double[] ted = new double[numStrains];
            Marshal.Copy(strainsPtr, ted, 0, numStrains);
            Marshal.FreeCoTaskMem(strainsPtr);

            int count = numStrains / 2;
            double[] strains = new double[count];
            for (int i=0; i<count; i++)
            {
                strains[i] = ted[i * 2];
            }
            return strains;
        }

        public double[] GetStrainsPerMeshVertex()
        {
            double[] strainsPerFace = GetStrainsPerMeshFace();
            int numVertices = VisualizationMesh.Vertices.Count;
            Point3d[] pts = VisualizationMesh.Vertices.ToPoint3dArray();
            double[] strainsPerVertex = new double[numVertices];

            for(int i=0; i<numVertices; i++)
            {
                int[] faceIdx = VisualizationMesh.Vertices.GetVertexFaces(i);

                double weightedStrain=0, totalFaceAreas = 0;
                for (int j=0; j<faceIdx.Length; j++)
                {
                    MeshFace face = VisualizationMesh.Faces.GetFace(faceIdx[j]);
                    double faceStrain = strainsPerFace[faceIdx[j]];
                    double faceArea = AreaMassProperties.Compute(new PolylineCurve(new Point3d[] { pts[face.A], pts[face.B], pts[face.C], pts[face.A] })).Area;
                    weightedStrain += faceStrain * faceArea;
                    totalFaceAreas += faceArea;
                }

                strainsPerVertex[i] = weightedStrain / totalFaceAreas;
            }

            return strainsPerVertex;
        }

        public void SetUseTensionFieldEnergy(bool useTensionFieldEnergy)
        {
            if (ModelIO.ModelType == ElasticBodyType.Single_layer_inflatable) Kernel.InflatableSheet.SetUseTensionFieldEnergy(ModelPtr, useTensionFieldEnergy ? 1 : 0);
            else if (ModelIO.ModelType == ElasticBodyType.Periodic_Unit) Kernel.PeriodicUnit.SetUseTensionFieldEnergy(ModelPtr, useTensionFieldEnergy ? 1 : 0);
            else if (ModelIO.ModelType == ElasticBodyType.Multi_layer_inflatable) Kernel.MultipleSheets.SetUseTensionFieldEnergy(ModelPtr, useTensionFieldEnergy ? 1 : 0);
            else if (ModelIO.ModelType == ElasticBodyType.AttractedTargetInflation) Kernel.TargetAttractedInflation.SetUseTensionFieldEnergy(ModelPtr, useTensionFieldEnergy ? 1 : 0);
            else throw new Exception(MSG_ERROR_MODEL);
        }

        public void SetUseHessianProjectedEnergy(bool useHessianProjectedEnergy)
        {
            if (ModelIO.ModelType == ElasticBodyType.Single_layer_inflatable) Kernel.InflatableSheet.SetUseHessianProjectedEnergy(ModelPtr, useHessianProjectedEnergy ? 1 : 0);
            else if (ModelIO.ModelType == ElasticBodyType.Periodic_Unit) Kernel.PeriodicUnit.SetUseHessianProjectedEnergy(ModelPtr, useHessianProjectedEnergy ? 1 : 0);
            else if (ModelIO.ModelType == ElasticBodyType.Multi_layer_inflatable) Kernel.MultipleSheets.SetUseHessianProjectedEnergy(ModelPtr, useHessianProjectedEnergy ? 1 : 0);
            else if (ModelIO.ModelType == ElasticBodyType.AttractedTargetInflation) Kernel.TargetAttractedInflation.SetUseHessianProjectedEnergy(ModelPtr, useHessianProjectedEnergy ? 1 : 0);
            else throw new Exception(MSG_ERROR_MODEL);
        }

        public void DisableFusedRegionTensionFieldTheory(bool disableFusedRegionTensionFieldTheory)
        {
            if (ModelIO.ModelType == ElasticBodyType.Single_layer_inflatable) Kernel.InflatableSheet.DisableFusedRegionTensionFieldTheory(ModelPtr, disableFusedRegionTensionFieldTheory ? 1 : 0);
            else if (ModelIO.ModelType == ElasticBodyType.Periodic_Unit) Kernel.PeriodicUnit.DisableFusedRegionTensionFieldTheory(ModelPtr, disableFusedRegionTensionFieldTheory ? 1 : 0);
            else if (ModelIO.ModelType == ElasticBodyType.AttractedTargetInflation) Kernel.TargetAttractedInflation.DisableFusedRegionTensionFieldTheory(ModelPtr, disableFusedRegionTensionFieldTheory ? 1 : 0);
            else throw new Exception(MSG_ERROR_MODEL);
        }

        public int[] GetCenterFixedVars()
        {
            int numVars;
            IntPtr varsPtr;
            if (ModelIO.ModelType == ElasticBodyType.Single_layer_inflatable) Kernel.InflatableSheet.GetCenterFixedVars(ModelPtr, out varsPtr, out numVars);
            else if (ModelIO.ModelType == ElasticBodyType.Periodic_Unit) Kernel.PeriodicUnit.GetCenterFixedVars(ModelPtr, out varsPtr, out numVars);
            else if (ModelIO.ModelType == ElasticBodyType.Multi_layer_inflatable) Kernel.MultipleSheets.GetCenterFixedVars(ModelPtr, out varsPtr, out numVars);
            else if (ModelIO.ModelType == ElasticBodyType.AttractedTargetInflation) Kernel.TargetAttractedInflation.GetCenterFixedVars(ModelPtr, out varsPtr, out numVars);
            else throw new Exception(MSG_ERROR_MODEL);

            int[] fixedVars = new int[numVars];
            Marshal.Copy(varsPtr, fixedVars, 0, numVars);
            Marshal.FreeCoTaskMem(varsPtr);

            return fixedVars;
        }

        public int[] GetVertexVars(bool isTopSheet, int vertexIdx)
        {
            int sheetIdx = isTopSheet ? 0 : 1;
            int numVars;
            IntPtr varsPtr;
            if (ModelIO.ModelType == ElasticBodyType.Single_layer_inflatable) Kernel.InflatableSheet.GetVertexVars(ModelPtr, sheetIdx, vertexIdx, out varsPtr, out numVars);
            else if (ModelIO.ModelType == ElasticBodyType.Periodic_Unit) Kernel.PeriodicUnit.GetVertexVars(ModelPtr, sheetIdx, vertexIdx, out varsPtr, out numVars);
            else if (ModelIO.ModelType == ElasticBodyType.Multi_layer_inflatable) Kernel.MultipleSheets.GetVertexVars(ModelPtr, sheetIdx, vertexIdx, out varsPtr, out numVars);
            else if (ModelIO.ModelType == ElasticBodyType.AttractedTargetInflation) Kernel.TargetAttractedInflation.GetVertexVars(ModelPtr, sheetIdx, vertexIdx, out varsPtr, out numVars);
            else throw new Exception(MSG_ERROR_MODEL);

            int[] fixedVars = new int[numVars];
            Marshal.Copy(varsPtr, fixedVars, 0, numVars);
            Marshal.FreeCoTaskMem(varsPtr);

            return fixedVars;
        }

        protected void GetMeshVisualization(out double[] outCoords, out int[] outElements)
        {
            int numCoords, numElements;
            IntPtr cPtr, qPtr;

            if (ModelIO.ModelType == ElasticBodyType.Single_layer_inflatable) Kernel.InflatableSheet.GetMeshVisualization(ModelPtr, out cPtr, out qPtr, out numCoords, out numElements);
            else if (ModelIO.ModelType == ElasticBodyType.Periodic_Unit) Kernel.PeriodicUnit.GetMeshVisualization(ModelPtr, out cPtr, out qPtr, out numCoords, out numElements);
            else if (ModelIO.ModelType == ElasticBodyType.Multi_layer_inflatable) Kernel.MultipleSheets.GetMeshVisualization(ModelPtr, out cPtr, out qPtr, out numCoords, out numElements);
            else if (ModelIO.ModelType == ElasticBodyType.AttractedTargetInflation) Kernel.TargetAttractedInflation.GetMeshVisualization(ModelPtr, out cPtr, out qPtr, out numCoords, out numElements);
            else throw new Exception(MSG_ERROR_MODEL);

            outCoords = new double[numCoords];
            outElements = new int[numElements];
            Marshal.Copy(cPtr, outCoords, 0, numCoords);
            Marshal.Copy(qPtr, outElements, 0, numElements);
            Marshal.FreeCoTaskMem(cPtr);
            Marshal.FreeCoTaskMem(qPtr);
        }

        private void InitMeshVisualization()
        {
            double[] outCoords;
            int[] outElements;
            GetMeshVisualization(out outCoords, out outElements);

            VisualizationMesh = new Mesh();
            VisualizationMesh = Helpers.GetTriasMesh(outCoords, outElements);
        }

        public void UpdateMeshVisualization()
        {
            double[] outCoords;
            int[] outElements;
            GetMeshVisualization(out outCoords, out outElements);

            int vCount = outCoords.Length / 3;
            for (int i = 0; i < vCount; i++)
            {
                VisualizationMesh.Vertices.SetVertex(i, outCoords[i * 3], outCoords[i * 3 + 1], outCoords[i * 3 + 2]);
            }
        }

        public void SetGravityVector(Vector3d vector)
        {
            double[] vec = new double[] { vector.X, vector.Y, vector.Z};
            if (ModelIO.ModelType == ElasticBodyType.Single_layer_inflatable) Kernel.InflatableSheet.SetGravity(ModelPtr, vec);
            else if (ModelIO.ModelType == ElasticBodyType.Periodic_Unit) Kernel.PeriodicUnit.SetGravity(ModelPtr, vec);
            else if (ModelIO.ModelType == ElasticBodyType.Multi_layer_inflatable) Kernel.MultipleSheets.SetGravity(ModelPtr, vec);
            else if (ModelIO.ModelType == ElasticBodyType.AttractedTargetInflation) Kernel.TargetAttractedInflation.SetGravity(ModelPtr, vec);
            else throw new Exception(MSG_ERROR_MODEL);
        }

        public Vector3d GetGravityVector()
        {
            int numVars;
            IntPtr varsPtr;

            if (ModelIO.ModelType == ElasticBodyType.Single_layer_inflatable) Kernel.InflatableSheet.GetGravity(ModelPtr, out varsPtr, out numVars);
            else if (ModelIO.ModelType == ElasticBodyType.Periodic_Unit) Kernel.PeriodicUnit.GetGravity(ModelPtr, out varsPtr, out numVars);
            else if (ModelIO.ModelType == ElasticBodyType.Multi_layer_inflatable) Kernel.MultipleSheets.GetGravity(ModelPtr, out varsPtr, out numVars);
            else if (ModelIO.ModelType == ElasticBodyType.AttractedTargetInflation) Kernel.TargetAttractedInflation.GetGravity(ModelPtr, out varsPtr, out numVars);
            else throw new Exception(MSG_ERROR_MODEL);

            double[] xyz = new double[numVars];
            Marshal.Copy(varsPtr, xyz, 0, numVars);
            Marshal.FreeCoTaskMem(varsPtr);

            return new Vector3d(xyz[0], xyz[1], xyz[2]);
        }

        public void SetMaterialDensity(double rho)
        {
            if (ModelIO.ModelType == ElasticBodyType.Single_layer_inflatable) Kernel.InflatableSheet.SetMaterialDensity(ModelPtr, rho);
            else if (ModelIO.ModelType == ElasticBodyType.Periodic_Unit) Kernel.PeriodicUnit.SetMaterialDensity(ModelPtr, rho);
            else if (ModelIO.ModelType == ElasticBodyType.Multi_layer_inflatable) Kernel.MultipleSheets.SetMaterialDensity(ModelPtr, rho);
            else if(ModelIO.ModelType == ElasticBodyType.AttractedTargetInflation) Kernel.TargetAttractedInflation.SetMaterialDensity(ModelPtr, rho);
            else throw new Exception(MSG_ERROR_MODEL);
        }

        public double GetMaterialDensity()
        {
            if (ModelIO.ModelType == ElasticBodyType.Single_layer_inflatable) return Kernel.InflatableSheet.GetMaterialDensity(ModelPtr);
            else if (ModelIO.ModelType == ElasticBodyType.Periodic_Unit) return Kernel.PeriodicUnit.GetMaterialDensity(ModelPtr);
            else if (ModelIO.ModelType == ElasticBodyType.Multi_layer_inflatable) return Kernel.MultipleSheets.GetMaterialDensity(ModelPtr);
            else if (ModelIO.ModelType == ElasticBodyType.AttractedTargetInflation) return Kernel.TargetAttractedInflation.GetMaterialDensity(ModelPtr);
            else throw new Exception(MSG_ERROR_MODEL);
        }

        public void SetYoungModulus(double youngModulus)
        {
            if (ModelIO.ModelType == ElasticBodyType.Single_layer_inflatable) Kernel.InflatableSheet.SetYoungModulus(ModelPtr, youngModulus);
            else if (ModelIO.ModelType == ElasticBodyType.Periodic_Unit) Kernel.PeriodicUnit.SetYoungModulus(ModelPtr, youngModulus);
            else if (ModelIO.ModelType == ElasticBodyType.Multi_layer_inflatable) Kernel.MultipleSheets.SetYoungModulus(ModelPtr, youngModulus);
            else if (ModelIO.ModelType == ElasticBodyType.AttractedTargetInflation) Kernel.TargetAttractedInflation.SetYoungModulus(ModelPtr, youngModulus);
            else throw new Exception(MSG_ERROR_MODEL);
        }

        public double GetYoungModulus()
        {
            if (ModelIO.ModelType == ElasticBodyType.Single_layer_inflatable) return Kernel.InflatableSheet.GetYoungModulus(ModelPtr);
            else if (ModelIO.ModelType == ElasticBodyType.Periodic_Unit) return Kernel.PeriodicUnit.GetYoungModulus(ModelPtr);
            else if (ModelIO.ModelType == ElasticBodyType.Multi_layer_inflatable) return Kernel.MultipleSheets.GetYoungModulus(ModelPtr);
            else if (ModelIO.ModelType == ElasticBodyType.AttractedTargetInflation) return Kernel.TargetAttractedInflation.GetYoungModulus(ModelPtr);
            else throw new Exception(MSG_ERROR_MODEL);
        }

        public void SetThickness(double thickness)
        {
            if (ModelIO.ModelType == ElasticBodyType.Single_layer_inflatable) Kernel.InflatableSheet.SetThickness(ModelPtr, thickness);
            else if (ModelIO.ModelType == ElasticBodyType.Periodic_Unit) Kernel.PeriodicUnit.SetThickness(ModelPtr, thickness);
            else if (ModelIO.ModelType == ElasticBodyType.Multi_layer_inflatable) Kernel.MultipleSheets.SetThickness(ModelPtr, thickness);
            else if (ModelIO.ModelType == ElasticBodyType.AttractedTargetInflation) Kernel.TargetAttractedInflation.SetThickness(ModelPtr, thickness);
            else throw new Exception(MSG_ERROR_MODEL);
        }

        public double GetThickness()
        {
            if (ModelIO.ModelType == ElasticBodyType.Single_layer_inflatable) return Kernel.InflatableSheet.GetThickness(ModelPtr);
            else if (ModelIO.ModelType == ElasticBodyType.Periodic_Unit) return Kernel.PeriodicUnit.GetThickness(ModelPtr);
            else if (ModelIO.ModelType == ElasticBodyType.Multi_layer_inflatable) return Kernel.MultipleSheets.GetThickness(ModelPtr);
            else if (ModelIO.ModelType == ElasticBodyType.AttractedTargetInflation) return Kernel.TargetAttractedInflation.GetThickness(ModelPtr);
            else throw new Exception(MSG_ERROR_MODEL);
        }

        public double[] GetVolume()
        {
            if (ModelIO.ModelType == ElasticBodyType.Single_layer_inflatable) return new double[] { Kernel.InflatableSheet.GetVolume(ModelPtr) };
            else if (ModelIO.ModelType == ElasticBodyType.Periodic_Unit) return new double[] { Kernel.PeriodicUnit.GetVolume(ModelPtr) };
            else if (ModelIO.ModelType == ElasticBodyType.Multi_layer_inflatable)
            {
                int numVolumes;
                IntPtr volumesPtr;
                Kernel.MultipleSheets.GetVolumes(ModelPtr, out volumesPtr, out numVolumes);
                double[] volumes = new double[numVolumes];
                Marshal.Copy(volumesPtr, volumes, 0, numVolumes);
                Marshal.FreeCoTaskMem(volumesPtr);
                return volumes;
            }
            else if (ModelIO.ModelType == ElasticBodyType.AttractedTargetInflation) return new double[] { Kernel.TargetAttractedInflation.GetVolume(ModelPtr) };
            else throw new Exception(MSG_ERROR_MODEL);
        }

        public void SetReferenceVolume(double[] volumes)
        {
            if (ModelIO.ModelType == ElasticBodyType.Single_layer_inflatable) Kernel.InflatableSheet.SetReferenceVolume(ModelPtr, volumes[0]);
            else if (ModelIO.ModelType == ElasticBodyType.Periodic_Unit) Kernel.PeriodicUnit.SetReferenceVolume(ModelPtr, volumes[0]);
            else if (ModelIO.ModelType == ElasticBodyType.Multi_layer_inflatable) Kernel.MultipleSheets.SetReferenceVolumes(ModelPtr, volumes, volumes.Length);
            else if (ModelIO.ModelType == ElasticBodyType.AttractedTargetInflation) Kernel.TargetAttractedInflation.SetReferenceVolume(ModelPtr, volumes[0]);
            else throw new Exception(MSG_ERROR_MODEL);
        }

        public double[] GetReferenceVolume()
        {
            if (ModelIO.ModelType == ElasticBodyType.Single_layer_inflatable) return new double[] { Kernel.InflatableSheet.GetReferenceVolume(ModelPtr) };
            else if (ModelIO.ModelType == ElasticBodyType.Periodic_Unit) return new double[] { Kernel.PeriodicUnit.GetReferenceVolume(ModelPtr) };
            else if (ModelIO.ModelType == ElasticBodyType.Multi_layer_inflatable)
            {
                int numVolumes;
                IntPtr volumesPtr;
                Kernel.MultipleSheets.GetReferenceVolumes(ModelPtr, out volumesPtr, out numVolumes);
                double[] volumes = new double[numVolumes];
                Marshal.Copy(volumesPtr, volumes, 0, numVolumes);
                Marshal.FreeCoTaskMem(volumesPtr);
                return volumes;
            }
            else if (ModelIO.ModelType == ElasticBodyType.AttractedTargetInflation) return new double[] { Kernel.TargetAttractedInflation.GetReferenceVolume(ModelPtr) };
            else throw new Exception(MSG_ERROR_MODEL);
        }

        public void SetPressure(double[] pressure)
        {
            if (ModelIO.ModelType == ElasticBodyType.Single_layer_inflatable) Kernel.InflatableSheet.SetPressure(ModelPtr, pressure[0]);
            else if (ModelIO.ModelType == ElasticBodyType.Periodic_Unit) Kernel.PeriodicUnit.SetPressure(ModelPtr, pressure[0]);
            else if (ModelIO.ModelType == ElasticBodyType.Multi_layer_inflatable) Kernel.MultipleSheets.SetPressure(ModelPtr, pressure, pressure.Length);
            else if (ModelIO.ModelType == ElasticBodyType.AttractedTargetInflation) Kernel.TargetAttractedInflation.SetPressure(ModelPtr, pressure[0]);
            else throw new Exception(MSG_ERROR_MODEL);
        }

        public double[] GetPressure()
        {
            if (ModelIO.ModelType == ElasticBodyType.Single_layer_inflatable) return new double[] { Kernel.InflatableSheet.GetPressure(ModelPtr) };
            else if (ModelIO.ModelType == ElasticBodyType.Periodic_Unit) return new double[] { Kernel.PeriodicUnit.GetPressure(ModelPtr) };
            else if (ModelIO.ModelType == ElasticBodyType.Multi_layer_inflatable)
            {
                int numPressures;
                IntPtr pressuresPtr;
                Kernel.MultipleSheets.GetPressures(ModelPtr, out pressuresPtr, out numPressures);
                double[] pressures = new double[numPressures];
                Marshal.Copy(pressuresPtr, pressures, 0, numPressures);
                Marshal.FreeCoTaskMem(pressuresPtr);
                return pressures;
            }
            else if (ModelIO.ModelType == ElasticBodyType.AttractedTargetInflation) return new double[] { Kernel.TargetAttractedInflation.GetPressure(ModelPtr) };
            else throw new Exception(MSG_ERROR_MODEL);
        }

        #region GH_Methods
        public bool IsValid
        {
            get
            {
                if (ModelPtr != null || ModelPtr != IntPtr.Zero) return true;
                else return false;
            }
        }

        public string IsValidWhyNot => "";

        public string TypeName => ToString();

        public string TypeDescription => "Elastic Body";

        public IGH_Goo Duplicate()
        {
            return (IGH_Goo)this.MemberwiseClone();
        }

        public IGH_GooProxy EmitProxy()
        {
            return null;
        }

        public bool CastFrom(object source)
        {
            return false;
        }

        public bool CastTo<T>(out T target)
        {
            if (typeof(T).Equals(typeof(GH_Mesh)))
            {
                UpdateMeshVisualization();
                target = (T)(object)new GH_Mesh(VisualizationMesh);
                return true;
            }

            target = default(T);
            return false;
        }

        public object ScriptVariable()
        {
            return null;
        }

        public bool Write(GH_IWriter writer)
        {
            return false;
        }

        public bool Read(GH_IReader reader)
        {
            return false;
        }
        #endregion
    }
}

