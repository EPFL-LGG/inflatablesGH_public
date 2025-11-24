using System;
using ISheetDataLib.Types;

namespace ISheetModelLib.Types
{
	public class TargetAttractedInflation : ElasticModel
    {
        public TargetAttractedInflation(InflatableData data)
        {
            Init(data);
        }

        public TargetAttractedInflation(TargetAttractedInflation model)
        {
            Init(model.ModelIO);
        }

        public override object Clone()
        {
            return new TargetAttractedInflation(this);
        }
    }
}

