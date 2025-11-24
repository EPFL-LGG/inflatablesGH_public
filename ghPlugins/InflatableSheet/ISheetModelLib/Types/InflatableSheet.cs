using System;
using GH_IO.Serialization;
using Grasshopper.Kernel.Types;
using System.Runtime.InteropServices;
using Rhino.Geometry;
using ISheetDataLib.Types;
using ISheetDataLib.Utils;
using ISheetModelLib.Creators;

namespace ISheetModelLib.Types
{
    public partial class InflatableSheet : ElasticModel
    {
        public InflatableSheet(InflatableData data)
        {
            Init(data);
        }

        public InflatableSheet(InflatableSheet model)
        {
            Init(model.ModelIO);
        }

        public override object Clone()
        {
            return new InflatableSheet(this);
        }
    }
}

