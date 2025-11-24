using System;
using System.Runtime.InteropServices;
using ISheetDataLib.Types;
using ISheetModelLib.Creators;

namespace ISheetModelLib.Types
{
	public class PeriodicUnit : ElasticModel
    {
        public PeriodicUnit(InflatableData data)
        {
            Init(data);
        }

        public PeriodicUnit(PeriodicUnit model)
        {
            Init(model.ModelIO);
        }

        public void ReparametrizeVerticalOffset()
        {
            Kernel.PeriodicUnit.ReparametrizeVerticalOffset(ModelPtr);
        }

        public int[] GetBendingStiffnessFixedVars()
        {
            int numFixedVars;
            IntPtr fixedVarsPtr;
            Kernel.PeriodicUnit.GetBendingStiffnessFixedVars(ModelPtr, out fixedVarsPtr, out numFixedVars);

            int[] fixedVars = new int[numFixedVars];
            Marshal.Copy(fixedVarsPtr, fixedVars, 0, numFixedVars);
            Marshal.FreeCoTaskMem(fixedVarsPtr);
            return fixedVars;
        }

        public int[] GetStretchingStiffnessFixedVars()
        {
            int numFixedVars;
            IntPtr fixedVarsPtr;
            Kernel.PeriodicUnit.GetStretchingStiffnessFixedVars(ModelPtr, out fixedVarsPtr, out numFixedVars);

            int[] fixedVars = new int[numFixedVars];
            Marshal.Copy(fixedVarsPtr, fixedVars, 0, numFixedVars);
            Marshal.FreeCoTaskMem(fixedVarsPtr);
            return fixedVars;
        }

        public override object Clone()
        {
            return new PeriodicUnit(this);
        }
    }
}

