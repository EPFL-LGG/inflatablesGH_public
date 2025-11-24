using System;
using Grasshopper.Kernel;
using ISheetDataLib.Types;
using Rhino;
using Rhino.DocObjects;
using Rhino.Geometry;

namespace ISheetModelLib.Types
{
    public abstract partial class ElasticModel : IGH_PreviewData, IGH_BakeAwareData
    {
        public BoundingBox ClippingBox
        {
            get
            {
                return new BoundingBox(VisualizationMesh.Vertices.ToPoint3dArray());
            }
        }

        public bool BakeGeometry(RhinoDoc doc, ObjectAttributes att, out Guid obj_guid)
        {
            obj_guid = Guid.Empty;

            if (att == null) att = doc.CreateDefaultAttributes();

            string id = Guid.NewGuid().ToString();
            int idxGr = doc.Groups.Add(ToString() + id);

            ObjectAttributes att1 = att.Duplicate();
            att1.AddToGroup(idxGr);

            doc.Objects.AddMesh(VisualizationMesh);

            return true;
        }

        public void DrawViewportMeshes(GH_PreviewMeshArgs args)
        {
            if (ModelIO.ModelType == ElasticBodyType.Single_layer_inflatable) args.Pipeline.DrawPoints(ModelIO.Supports.GetSupportsAsPoint3dArray(), Rhino.Display.PointStyle.RoundControlPoint, (float)7, System.Drawing.Color.Blue);
            args.Pipeline.DrawMeshShaded(VisualizationMesh, args.Material);
        }

        public void DrawViewportWires(GH_PreviewWireArgs args)
        {
            if (ModelIO.ModelType == ElasticBodyType.Single_layer_inflatable) args.Pipeline.DrawPoints(ModelIO.Supports.GetSupportsAsPoint3dArray(), Rhino.Display.PointStyle.RoundControlPoint, (float)7, System.Drawing.Color.Blue);
            args.Pipeline.DrawMeshWires(VisualizationMesh, args.Color);
        }
    }
}
