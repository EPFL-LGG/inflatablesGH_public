using System;
namespace ISheetModelLib.Types
{
	public struct StiffnessResults
	{
		public double[] BendingStiffness { get; private set; }
        public double[] BendingThetaSamples { get; private set; }
        public double[] BendingCoefficient { get; private set; }
        public double[] StretchingStiffness { get; private set; }
        public double[] StretchingThetaSamples { get; private set; }


        public void SetBendingStiffness(double[] bendingStiffness, double[] bendingThetaSamples)
        {
            BendingStiffness = bendingStiffness;
            BendingThetaSamples = bendingThetaSamples;
        }

        public void SetBendingCoefficient(double[] bendingCoefficient)
        {
            BendingCoefficient = bendingCoefficient;
        }

        public void SetStretchingStiffness(double[] stretchingStiffness, double[] stretchingThetaSamples)
        {
            StretchingStiffness = stretchingStiffness;
            StretchingThetaSamples = stretchingThetaSamples;
        }
    }
}

