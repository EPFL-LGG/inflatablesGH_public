using System;
using System.IO;
using GH_IO.Serialization;
using Grasshopper.Kernel.Types;
using Newtonsoft.Json;

namespace ISheetDataLib.Types
{
    public partial class InflatableData : IGH_Goo
    {
        public override string ToString()
        {
            return ModelType.ToString();
        }

        public void WriteJsonFile(string path, string filename)
        {
            // Serialize JSON directly to a file
            using (StreamWriter file = File.CreateText(@path + filename + ".json"))
            {
                JsonSerializer serializer = new JsonSerializer();
                serializer.Serialize(file, this);
            }
        }

        #region GH_Methods
        public bool IsValid {get;set;}

        public string IsValidWhyNot => "";

        public string TypeName => ToString();

        public string TypeDescription => "Data set for initializing an inflatable model.";

        public bool CastFrom(object source)
        {
            throw new NotImplementedException();
        }

        public bool CastTo<T>(out T target)
        {
            throw new NotImplementedException();
        }

        public IGH_Goo Duplicate()
        {
            throw new NotImplementedException();
        }

        public IGH_GooProxy EmitProxy()
        {
            throw new NotImplementedException();
        }

        public bool Read(GH_IReader reader)
        {
            throw new NotImplementedException();
        }

        public object ScriptVariable()
        {
            throw new NotImplementedException();
        }

        public bool Write(GH_IWriter writer)
        {
            throw new NotImplementedException();
        }
        #endregion
    }
}

