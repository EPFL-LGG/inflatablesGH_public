using System;
using System.Collections.Generic;
using System.Linq;
using Newtonsoft.Json.Linq;
using Rhino.Geometry;

namespace ISheetDataLib.Types
{
	public class Support
	{
        public Point3d InitialPosition { get; set; }
        public Point3d TargetPosition { get; set; }
        public int Index { get; set; }
        public bool[] LockedDoFs { get; private set; }
        public bool IsTopSheet { get; private set; }
        public int[] IndicesDoFs { get; private set; }
        public bool IsTemporary { get; set; }

        public Support(JToken data)
        {
            // Indexes
            var token = data["Position"];
            InitialPosition = new Point3d((double)token[0], (double)token[1], (double)token[2]);

            token = data["TargetPosition"];
            TargetPosition = new Point3d((double)token[0], (double)token[1], (double)token[2]);

            // Indexes
            Index = (int) data["Index"];

            // Locked DOF
            token = data["LockedDOF"];
            LockedDoFs = new bool[3];
            for (int i = 0; i < 3; i++) LockedDoFs[i] = (bool)token[i];

            // Indices DOF
            token = data["IndicesDoFs"];
            IndicesDoFs = new int[data.Count()];
            for (int i = 0; i < data.Count(); i++) IndicesDoFs[i] = (int)token[i];

            // Sheet
            IsTopSheet = (bool) data["IsTopSheet"];

            IsTemporary = (bool)data["IsTemporary"];
        }

        public Support(Point3d p, bool[] DOF, bool isTopSheet, int index=-1, bool isTemporary=false)
        {
            InitialPosition = p;
            TargetPosition = p;
            Index = index;
            if (DOF.Length != 3) throw new Exception("Invalid input for support conditions");
            else LockedDoFs = DOF;
            IsTopSheet = isTopSheet;
            IndicesDoFs = new int[] { -1,-1,-1 };
            IsTemporary = isTemporary;
        }

        public Support(Support support)
        {
            InitialPosition = support.InitialPosition;
            TargetPosition = support.TargetPosition;
            Index = support.Index;
            LockedDoFs = support.LockedDoFs.ToArray();
            IndicesDoFs = support.IndicesDoFs.ToArray();
            IsTopSheet = support.IsTopSheet;
            IsTemporary = support.IsTemporary;
        }

        public void SetDoFsIndices(int[] indices)
        {
            for (int i = 0; i < LockedDoFs.Count(); i++)
            {
                if (LockedDoFs[i]) IndicesDoFs[i] = indices[i];
            }
        }

        public override string ToString()
        {
            return IsTemporary ?  "TempSupportData" : "SupportData";
        }
    }
}

