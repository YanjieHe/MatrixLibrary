using System;
using System.Threading.Tasks;

namespace LinearAlgebra.MatrixAlgebra
{
    internal static class MatrixComputation
    {
        const Double eps = 1e-9;
        internal static double[,] UnaryMinus(double[,] mat)
        {
            var ans = new double[mat.RowCount(), mat.ColumnCount()];
            unsafe
            {
                fixed (double* a = ans)
                fixed (double* m = mat)
                for (int i = mat.Length - 1; i >= 0; i--)
                    a[i] = -m[i];
            }
            return ans;
        }
        internal static double[,] Add(double[,] A, double[,] B)
        {
            int mA = A.RowCount(), mB = B.RowCount(),
                nA = A.ColumnCount(), nB = B.ColumnCount();
            if (mA != mB || nA != nB)
                throw new ArgumentException("进行矩阵加法运算的两个矩阵的形状必须相同");
            var ans = new double[mA, nA];
            unsafe
            {
                fixed (double* answer = ans)
                fixed (double* a = A)
                fixed (double* b = B)
                    for (int i = mA * nA - 1; i >= 0; i--)
                    answer[i] = a[i] + b[i];
            }
            return ans;
        }
        internal static double[,] Subtract(double[,] A, double[,] B)
        {
            int mA = A.RowCount(), mB = B.RowCount(),
                nA = A.ColumnCount(), nB = B.ColumnCount();
            if (mA != mB || nA != nB)
                throw new ArgumentException("进行矩阵减法运算的两个矩阵的形状必须相同");
            var ans = new double[mA, nA];
            unsafe
            {
                fixed (double* answer = ans)
                fixed (double* a = A)
                fixed (double* b = B)
                   for (int i = mA * nA - 1; i >= 0; i--)
                    answer[i] = a[i] - b[i];
            }
            return ans;
        }
        internal static double[,] Multiply(double[,] A, double[,] B)
        {
            int mA = A.RowCount(), mB = B.RowCount(),
                nA = A.ColumnCount(), nB = B.ColumnCount();
            if (nA != mB)
                throw new ArgumentException(
                    "进行矩阵乘法运算时，第一个矩阵的列数和第二个矩阵的行数必须相同");
            /*
            [A][B][C]   *   [G][H]      [A * G + B * I + C * K][A * H + B * J + C * L]
            [D][E][F]       [I][J]  =   [D * G + E * I + F * K][D * H + E * J + F * L]
                            [K][L]
            */
            var ans = new double[mA, nB];
            unsafe
            {
                Parallel.For(0, mA, delegate (int i)/*并行计算矩阵乘法*/
                {
                    fixed (double* mat = ans)
                    fixed (double* a = A)
                    fixed (double* b = B)
                    CacheAwareMultiply(mat, a, b, i, nA, nB);
                });
            }
            return ans;
        }
        private static unsafe void CacheAwareMultiply
         (double* mat, double* a, double* b, int i, int nA, int nB)
        {
            int inA = i * nA;
            int inB = i * nB;
            int knB = 0;//int knB = k * nB;
            for (int k = 0; k < nA; k++, knB += nB)
            {
                double inAk = a[inA + k];
                for (int j = 0; j < nB; j++)
                    mat[inB + j] += inAk * b[knB + j];
            }
        }
        internal static double[,] Divide(double[,] A, double[,] B)
        {
            return Multiply(A, Inverse(B));
        }
        internal static double[,] Inverse(double[,] Mat)
        {
            int nRows = Mat.RowCount();
            int nCols = Mat.ColumnCount();
            if (nRows != nCols) throw new ArgumentException("只有方阵才可以求逆");
            double[,] M = Mat.CopyMatrix();
            var pnRow = new int[nCols];
            var pnCol = new int[nCols];
            double d = 0.0, p = 0.0;
            int k, u, v;
            unsafe
            {
                fixed (double* mat = M)
                {
                    //消元
                    for (k = 0; k < nCols; k++)
                    {
                        d = 0.0;
                        for (int i = k; i < nCols; i++)
                            for (int j = k; j < nCols; j++)
                            {
                                u = i * nCols + j;
                                p = Math.Abs(mat[u]);
                                if (p > d)//选主元
                                {
                                    d = p;
                                    pnRow[k] = i;
                                    pnCol[k] = j;
                                }
                            }

                        if (d == 0.0)
                            throw new Exception("解不出啊！");

                        if (pnRow[k] != k)
                            Utility.SwapRow(mat, k, pnRow[k], nCols);
                        if (pnCol[k] != k)
                            Utility.SwapColumn(mat, k, pnCol[k], nRows, nCols);

                        v = k * nCols + k;
                        mat[v] = 1.0 / mat[v];
                        for (int j = 0; j < nCols; j++)
                            if (j != k)
                            {
                                u = k * nCols + j;
                                mat[u] *= mat[v];
                            }

                        for (int i = 0; i < nCols; i++)
                            if (i != k)
                                for (int j = 0; j < nCols; j++)
                                    if (j != k)
                                    {
                                        u = i * nCols + j;
                                        mat[u] -= mat[i * nCols + k] * mat[k * nCols + j];
                                    }

                        for (int i = 0; i < nCols; i++)
                            if (i != k)
                            {
                                u = i * nCols + k;
                                mat[u] *= -mat[v];
                            }
                    }

                    //恢复行列次序
                    for (k = nCols - 1; k >= 0; k--)
                    {
                        if (pnCol[k] != k)
                            Utility.SwapRow(mat, k, pnCol[k], nCols);
                        if (pnRow[k] != k)
                            Utility.SwapColumn(mat, k, pnRow[k], nRows, nCols);
                    }
                }
            }
            return M;
        }
        
        internal static double[,] Transpose(double[,] A)
        {
            int nRows = A.RowCount();
            int nCols = A.ColumnCount();
            var result = new double[nCols, nRows];
            unsafe
            {
                fixed (double* a = A)
                fixed (double* mat = result)
                for (int i = 0; i < nRows; i++)
                    for (int j = 0; j < nCols; j++)
                        mat[j * nRows + i] = a[i * nCols + j];
            }
            return result;
        }
        internal static int Rank(double[,] Mat)
        {
            int nRows = Mat.RowCount();
            int nCols = Mat.ColumnCount();
            double[,] M = Mat.CopyMatrix();
            int k = Math.Min(nRows, nCols);
            double d, p;
            int u, v;
            int nis = 0, js = 0;
            int rank = 0;
            unsafe
            {
                fixed (double* mat = M)
                for (int w = 0; w < k; w++)//消元
                {
                    d = 0.0;
                    for (int i = w; i < nRows; i++)
                        for (int j = w; j < nCols; j++)
                        {
                            u = i * nCols + j;
                            p = Math.Abs(mat[u]);
                            if (p > d)//选主元
                            {
                                d = p;
                                nis = i;
                                js = j;
                            }
                        }

                    if (Math.Abs(d) < eps)//由于浮点数的误差，这里需要与精度进行比较
                        return rank;

                    rank++;
                    if (nis != w)
                        for (int j = w; j < nCols; j++)
                        {
                            u = w * nCols + j;
                            v = nis * nCols + j;
                            Utility.Swap(mat, u, v);
                        }

                    if (js != w)
                        for (int i = w; i < nRows; i++)
                        {
                            u = i * nCols + js;
                            v = i * nCols + w;
                            Utility.Swap(mat, u, v);
                        }

                    v = w * nCols + w;
                    for (int i = w + 1; i < nRows; i++)//书上这里有bug，这里应该是行数
                    {
                        p = mat[i * nCols + w] / mat[v];
                        for (int j = w + 1; j < nCols; j++)
                        {
                            u = i * nCols + j;
                            mat[u] -= p * mat[w * nCols + j];
                        }
                    }
                }
            }
            return rank;
        }
    }
}
