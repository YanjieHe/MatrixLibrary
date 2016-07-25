using System;
using System.Collections.Generic;
using System.Text;

namespace LinearAlgebra.MatrixAlgebra
{
    internal static class Utility
    {
        const int DOUBLE_SIZE = sizeof(double);
        /// <summary>
        /// 交换矩阵的两个指定位置的元素
        /// </summary>
        /// <param name="x">指向矩阵的指针</param>
        /// <param name="u"></param>
        /// <param name="v"></param>
        internal static unsafe void Swap(double* x, int u, int v)
        {
            double z = x[u];
            x[u] = x[v];
            x[v] = z;
        }

        /// <summary>
        /// 交换数组的两个指定位置的元素
        /// </summary>
        /// <param name="x"></param>
        /// <param name="i"></param>
        /// <param name="j"></param>
        internal static void Swap<T>(T[] x, int i, int j)
        {
            T z = x[i];
            x[i] = x[j];
            x[j] = z;
        }
        internal static unsafe void SwapRow(double* x, int i, int j, int nCols)
        {
            int u = i * nCols, v = j * nCols;
            for (int m = 0; m < nCols; m++)//交换行
            {
                Swap(x, u, v);
                u++;
                v++;
            }
        }
        internal static unsafe void SwapColumn(double* x, int i, int j, int nRows, int nCols)
        {
            int u = i, v = j;
            for (int m = 0; m < nRows; m++)//交换列
            {
                Swap(x, u, v);
                u += nCols;
                v += nCols;
            }
        }
        internal static int RowCount(this double[,] mat)
        {
            return mat.GetLength(0);
        }
        internal static int ColumnCount(this double[,] mat)
        {
            return mat.GetLength(1);
        }
        internal static double[,] ArrayToMatrix(int nRows, int nCols, double[] Elements)
        {
            int count = nRows * nCols;
            if (count != Elements.Length)
                throw new Exception("输入的元素与设置的矩阵大小不匹配");
            var mat = new Double[nRows, nCols];
            Buffer.BlockCopy(Elements, 0, mat, 0, DOUBLE_SIZE * count);
            return mat;
        }
        internal static double[,] IEnumerableToMatrix(int nRows, int nCols, IEnumerable<double> Elements)
        {
            int count = nRows * nCols;
            var M = new Double[nRows, nCols];
            unsafe
            {
                fixed (double* mat = M)
                using (var x = Elements.GetEnumerator())
                {
                    for (int i = 0; i < count; i++)
                    {
                        if (!x.MoveNext())
                            throw new Exception("输入的元素与设置的矩阵大小不匹配");
                        mat[i] = x.Current;
                    }
                    if (x.MoveNext())
                        throw new Exception("输入的元素与设置的矩阵大小不匹配");
                }
            }
            return M;
        }
        internal static string MatrixToString(double[,] mat)
        {
            int m = mat.RowCount();
            int n = mat.ColumnCount() - 1;
            var s = new StringBuilder();
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    s.Append(mat[i, j].ToString());
                    s.Append("\t");
                }
                s.AppendLine(mat[i, n].ToString());
            }
            return s.ToString();
        }
        internal static string MatrixToString(double[,] mat, string format, IFormatProvider formatProvider)
        {
            int m = mat.RowCount();
            int n = mat.ColumnCount() - 1;
            var s = new StringBuilder();
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    s.Append(mat[i, j].ToString(format, formatProvider));
                    s.Append("\t");
                }
                s.AppendLine(mat[i, n].ToString(format, formatProvider));
            }
            return s.ToString();
        }
        internal static double[,] CopyMatrix(this double[,] mat)
        {
            var NewMatrix = new Double[mat.RowCount(), mat.ColumnCount()];
            Array.Copy(mat, NewMatrix, mat.Length);
            return NewMatrix;
        }
        internal static double[] CopyVector(this double[] vector)
        {
            var NewVector = new Double[vector.Length];
            Array.Copy(vector, NewVector, vector.Length);
            return NewVector;
        }
        internal static Matrix rBind(Matrix[] Matrices)
        {
            int rowsSum = 0;
            int nCols = Matrices[0].ColumnCount;
            for (int i = 0; i < Matrices.Length; i++)
            {
                if (nCols != Matrices[i].ColumnCount)
                    throw new Exception("进行纵向拼接的矩阵的列数必须相同");
                rowsSum += Matrices[i].RowCount;
            }
            var result = new double[rowsSum, nCols];
            int position = 0;
            int currentCount = 0;
            for (int i = 0; i < Matrices.Length; i++)
            {
                currentCount = Matrices[i].Count;
                Array.Copy(Matrices[i].Elements, 0, result, position, currentCount);
                position += currentCount;
            }
            return result;
        }
        internal static Matrix cBind(Matrix[] Matrices)
        {
            int colsSum = 0;
            int nRows = Matrices[0].RowCount;
            for (int i = 0; i < Matrices.Length; i++)
            {
                if (nRows != Matrices[i].RowCount)
                    throw new Exception("进行横向拼接的矩阵的行数必须相同");
                colsSum += Matrices[i].ColumnCount;
            }
            var result = new double[nRows, colsSum];
            int currentCount = 0;
            for (int i = 0; i < Matrices.Length; i++)
            {
                int mat_cCount = Matrices[i].ColumnCount;
                for (int j = 0; j < nRows; j++)
                {
                    Array.Copy(Matrices[i].Elements, j * mat_cCount,
                        result, currentCount + j * colsSum, mat_cCount);
                }
                currentCount += mat_cCount;
            }
            return result;
        }
    }
}
