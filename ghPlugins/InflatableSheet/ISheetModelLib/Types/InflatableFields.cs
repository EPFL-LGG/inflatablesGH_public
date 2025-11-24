using System;
using GH_IO.Serialization;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Types;
using Rhino;
using Rhino.DocObjects;
using Rhino.Geometry;
using System.Drawing.Imaging;
using System.Linq;
using System.Drawing;
using static ISheetModelLib.Utils.ColorMaps;
using ISheetModelLib.Utils;

namespace ISheetModelLib.Types
{
	public class InflatableFields : IGH_PreviewData, IGH_Goo, IGH_BakeAwareData
    {
        //public double[] Data { get; private set; }
        public double[] DataPerVertex { get; private set; }
        public double[] NormalizedData { get; private set; }

        public Mesh Mesh { get; private set; }
        private Color _wireColor;

        public InflatableFields(ElasticModel model, ColorMapTypes cMap, int alpha, double lowerBound = default, double upperBound = default)
        {
            Mesh = model.VisualizationMesh.DuplicateMesh();

            //Data = rod.GetScalarFieldSqrtBendingEnergies();
            DataPerVertex = model.GetStrainsPerMeshVertex();

            double tol = +1.0e-8;
            double min = lowerBound == default ? DataPerVertex.Min() : lowerBound == 0 ? tol : lowerBound;
            double max = upperBound == default ? DataPerVertex.Max() : upperBound;
            double range = max - min + tol;

            Color[] colormap = ColorMaps.GetColorMap(cMap, alpha == 0 ? 1 : alpha);
            int lastColor = colormap.Length - 1;
            _wireColor = Color.FromArgb(alpha, 250, 235, 215); // AntiqueWhite
            NormalizedData = new double[DataPerVertex.Length];
            for (int i = 0; i < DataPerVertex.Length; i++)
            {
                NormalizedData[i] = (DataPerVertex[i] - min) / range;

                int colorIndex = (int)(NormalizedData[i] * lastColor);
                if (colorIndex < 0) colorIndex = 0;
                else if (colorIndex > lastColor) colorIndex = lastColor;
                Mesh.VertexColors.Add(colormap[colorIndex]);
            }
        }

        public override string ToString()
        {
            return "Inflatable Strains";
        }

        public bool BakeGeometry(RhinoDoc doc, ObjectAttributes att, out Guid obj_guid)
        {
            obj_guid = Guid.Empty;

            if (att == null) att = doc.CreateDefaultAttributes();

            string id = Guid.NewGuid().ToString();
            int idxGr = doc.Groups.Add(ToString() + id);

            ObjectAttributes att1 = att.Duplicate();
            att1.AddToGroup(idxGr);

            doc.Objects.AddMesh(Mesh);

            return true;
        }

        public BoundingBox ClippingBox
        {
            get
            {
                return new BoundingBox(Mesh.Vertices.ToPoint3dArray());
            }
        }

        public void DrawViewportMeshes(GH_PreviewMeshArgs args)
        {
            args.Pipeline.DrawMeshFalseColors(Mesh);
        }

        public void DrawViewportWires(GH_PreviewWireArgs args)
        {
            args.Pipeline.DrawMeshWires(Mesh, _wireColor);
        }

        #region GH_Methods
        public bool IsValid
        {
            get
            {
                return true;
            }
        }

        public string IsValidWhyNot => "";

        public string TypeName => "Inflatable Strains";

        public string TypeDescription => "";

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

