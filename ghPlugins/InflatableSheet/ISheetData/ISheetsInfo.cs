using System;
using System.Drawing;
using Grasshopper;
using Grasshopper.Kernel;

namespace ISheets
{
    public class ISheetsInfo : GH_AssemblyInfo
    {
        public override string Name => "ISheets Info";

        //Return a 24x24 pixel bitmap to represent this GHA library.
        public override Bitmap Icon => null;

        //Return a short string describing the purpose of this GHA library.
        public override string Description => "";

        public override Guid Id => new Guid("5A2B9DCB-FB61-4A3F-BF4F-FF4053C3896E");

        //Return a string identifying you or your company.
        public override string AuthorName => "";

        //Return a string representing your preferred contact details.
        public override string AuthorContact => "";
    }
}
