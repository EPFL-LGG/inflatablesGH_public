using System;
using Newtonsoft.Json.Linq;

namespace ISheetDataLib.Types
{
	public class MaterialData
	{
		public double E;
		public double Thickness;
		public double MaterialDensity;

		public MaterialData(JToken data)
		{
            E = (double)data["E"];
            Thickness = (double)data["Thickness"];
            MaterialDensity = (double)data["MaterialDensity"];
        }

        public MaterialData(double e=300, double thickness=0.075, double rho=0.0)
        {
            E = e;
            Thickness = thickness;
            MaterialDensity = rho;
        }

        public MaterialData(MaterialData data)
        {
            E = data.E;
            Thickness = data.Thickness;
            MaterialDensity = data.MaterialDensity;
        }

        public override string ToString()
        {
            return "InflatableMaterial";
        }
    }
}

