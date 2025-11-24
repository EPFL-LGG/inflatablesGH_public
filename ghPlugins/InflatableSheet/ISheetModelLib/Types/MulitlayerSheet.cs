using System;
using ISheetDataLib.Types;

namespace ISheetModelLib.Types
{
	public class MultipleSheets : ElasticModel
    {
        public MultipleSheets(InflatableData data)
        {
            Init(data);
        }

        public MultipleSheets(MultipleSheets model)
        {
            Init(model.ModelIO);
        }

        public override object Clone()
        {
            return new MultipleSheets(this);
        }
    }
}

