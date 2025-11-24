using System;
using System.Drawing;
using Grasshopper;
using Grasshopper.Kernel;

namespace ISheetModel
{
  public class ISheetModelInfo : GH_AssemblyInfo
  {
    public override string Name => "ISheetModel";

    //Return a 24x24 pixel bitmap to represent this GHA library.
    public override Bitmap Icon => null;

    //Return a short string describing the purpose of this GHA library.
    public override string Description => "";

    public override Guid Id => new Guid("b9c103d3-3c18-4a64-96be-d16c43539811");

    //Return a string identifying you or your company.
    public override string AuthorName => "";

    //Return a string representing your preferred contact details.
    public override string AuthorContact => "";
  }
}
