using System;
using System.Collections;
using System.Collections.Generic;
using Rhino.Geometry;
using System.Linq;

namespace ISheetDataLib.Types
{
	public class SupportCollection : IList<Support>, ICloneable
    {
        private List<Support> _supports;

        public SupportCollection()
        {
            _supports = new List<Support>();
        }

        public SupportCollection(SupportCollection supports)
        {
            _supports = new List<Support>(supports._supports);
        }

        public Support this[int index] { get => _supports[index]; set => _supports[index] = value; }

        public int Count => _supports.Count;

        public bool IsReadOnly => false;

        public void Add(Support support)
        {
            _supports.Add(support);
        }

        public void Clear()
        {
            _supports.Clear();
        }

        public bool Contains(Support support)
        {
            return _supports.Contains(support);
        }

        public void CopyTo(Support[] array, int arrayIndex)
        {
            _supports.ToList().CopyTo(array, arrayIndex);
        }

        public IEnumerator<Support> GetEnumerator()
        {
            return _supports.GetEnumerator();
        }

        public int IndexOf(Support support)
        {
            return _supports.IndexOf(support);
        }

        public void Insert(int index, Support support)
        {
            _supports.Insert(index, support);
        }

        public bool Remove(Support support)
        {
            return _supports.Remove(support);
        }

        public void RemoveAt(int index)
        {
            _supports.RemoveAt(index);
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return _supports.GetEnumerator();
        }

        public Point3d[] GetSupportsAsPoint3dArray(bool targetPositions=true)
        {
            return _supports.Select(sp => targetPositions ? sp.TargetPosition : sp.InitialPosition).ToArray();
        }

        public int[] GetSupportsDoFsIndices()
        {
            return _supports.Select(sp => sp.IndicesDoFs).SelectMany(dof => dof).Where(idx => idx>-1).ToHashSet().ToArray();
        }

        public int[] GetNonTemporarySupportsDoFsIndices()
        {
            List<int> indices = new List<int>();
            foreach (Support sp in _supports)
            {
                if (sp.IsTemporary) continue;
                indices.AddRange(sp.IndicesDoFs.Select(dof => dof).Where(idx => idx > -1).ToHashSet().ToArray());
            }
            return indices.ToArray();
        }

        public object Clone()
        {
            return new SupportCollection(this);
        }
    }
}

