using System;
namespace ISheetModelLib.Types
{
	public struct ExperimentLog
	{
        public bool Success { get; private set; }
        public double KappaValue { get; private set; }
        public bool PlanarEquilibrium { get; private set; }
        public double EnergyFirstStage { get; set; }
        public double EnergySecondStage { get; set; }
        public double EnergyThirdStage { get; set; }
        public double GradientNorm { get; set; }
        public bool ContainsInformation { get; set; }

        public ExperimentLog(int success, double kappaValue)
        {
            ContainsInformation = true;
            Success = Convert.ToBoolean(success);
            KappaValue = kappaValue;
            if (Math.Abs(KappaValue) > 1e-10) PlanarEquilibrium = false;
            else PlanarEquilibrium = true;

            EnergyFirstStage = -1;
            EnergySecondStage = -1;
            EnergyThirdStage = -1;
            GradientNorm = -1;
        }

        public override string ToString()
        {
            if (!ContainsInformation) return "--- Empty Experiment Log ---\"";

            string txt = "--- Experiment Log ---";
            txt += "\nInflatable periodic unit simulation succeed: " + Success;
            txt += "\nEnergies per solver stage: [1] " + EnergyFirstStage + " [2] " + EnergySecondStage;
            txt += "\nSimulation Kappa value: " + KappaValue;
            txt += "\nGradientNorm: " + GradientNorm;

            if(!PlanarEquilibrium) txt += "\nPlanar equilibrium: False - [WARNING - Can not compute stiffness due to non-planar equilibrium!]";
            else txt += "\nPlanar equilibrium: True";

            txt += "\n\n------------------------------\n\n";
            txt += "--- Additional solver stage [kappa not equal zero] ---";
            txt += "\nEnergy third solver stage: " + EnergyThirdStage;
            return txt;
        }
    }
}

