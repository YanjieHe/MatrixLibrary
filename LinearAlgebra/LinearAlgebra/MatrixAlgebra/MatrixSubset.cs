using System;
using System.Collections.Generic;
namespace LinearAlgebra.MatrixAlgebra
{
    internal sealed class MatrixSubset
    {
        private double[,] mat;
        internal MatrixSubset(double[,] Mat)
        {
            mat = Mat;
        }
        internal IEnumerable<double> GetRow(int index)
        {
            int nCols = mat.ColumnCount();
            for (int j = 0; j < nCols; j++)
                yield return mat[index, j];
        }
        internal IEnumerable<double> GetColumn(int index)
        {
            int nRows = mat.RowCount();
            for (int i = 0; i < nRows; i++)
                yield return mat[i, index];
        }
        internal IEnumerable<double> GetDiagonal(bool mainDiagonal = true)
        {
            int n = Math.Min(mat.RowCount(), mat.ColumnCount());
            if (mainDiagonal)
                for (int i = 0; i < n; i++)
                    yield return mat[i, i];
            else
            {
                int nCols = mat.ColumnCount();
                for (int i = 0; i < n; i++)
                    yield return mat[i, nCols - i - 1];
            }
        }
        internal void SetRow(int index, IEnumerable<double> data)
        {
            int nCols = mat.ColumnCount();
            using (var x = data.GetEnumerator())
            {
                for (int j = 0; j < nCols; j++)
                {
                    if (!x.MoveNext()) throw new ArgumentOutOfRangeException();
                    mat[index, j] = x.Current;
                }
                if (x.MoveNext()) throw new ArgumentOutOfRangeException();
            }
        }
        internal void SetColumn(int index, IEnumerable<double> data)
        {
            int nRows = mat.RowCount();
            using (var x = data.GetEnumerator())
            {
                for (int i = 0; i < nRows; i++)
                {
                    if (!x.MoveNext()) throw new ArgumentOutOfRangeException();
                    mat[i, index] = x.Current;
                }
                if (x.MoveNext()) throw new ArgumentOutOfRangeException();
            }
        }
        internal void SetDiagonal(IEnumerable<double> data, bool mainDiagonal = true)
        {
            int n = Math.Min(mat.RowCount(), mat.ColumnCount());
            int nCols = mat.ColumnCount();
            using (var x = data.GetEnumerator())
            {
                for (int i = 0; i < n; i++)
                {
                    if (!x.MoveNext()) throw new ArgumentOutOfRangeException();
                    if (mainDiagonal) mat[i, i] = x.Current;
                    else mat[i, nCols - i - 1] = x.Current;
                }
                if (x.MoveNext()) throw new ArgumentOutOfRangeException();
            }
        }
    }
}
