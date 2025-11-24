using System;
using Plotly.NET;

namespace ISheetModelLib.Utils
{
	public static class GraphPlotter
	{
        public static void PolarChart(string title, double[] r, double[] theta, bool toDegrees = true)
        {
            double[] thetaDeg = new double[theta.Length] ;
            if (toDegrees)
            {
                for (int i = 0; i < theta.Length; i++)
                {
                    thetaDeg[i] = theta[i] * 180 / Math.PI;
                }
            }
            else thetaDeg = theta;

            GenericChart.GenericChart chart = ChartPolar.Chart.SplinePolar<double, double, double>(r, thetaDeg);
            chart.WithTitle(title);
            chart.Show();
        }
	}
}

