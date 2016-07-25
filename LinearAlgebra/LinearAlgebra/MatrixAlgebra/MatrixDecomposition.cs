using System;
using System.Linq;

namespace LinearAlgebra.MatrixAlgebra
{
    internal sealed class MatrixDecomposition
    {
        const Double eps = 1e-9;
        double[,] mat;
        public MatrixDecomposition(double[,] Mat)
        {
            mat = Mat;
        }
        internal MatrixLU LU()
        {
            /*
            A = LU
            A = [1          ]   [*   *   ... *]
                |*  1       |   |    *   ... *|
                |...        |   |        * ...|   
                [*  * ...  1]   [            *]
            */
            int nRows = mat.RowCount();
            int nCols = mat.ColumnCount();
            if (nRows != nCols)
                throw new Exception("进行LU分解的矩阵需要是方阵");
            double[,] L = new double[nCols, nCols];
            double[,] U = mat;
            var pnRows = Enumerable.Range(0, nCols).ToArray();
            int w = 0, v;
            int Parity = 1;//用于记录奇偶性
            unsafe
            {
                fixed (double* l = L)
                fixed (double* u = U)
                for (int k = 0; k < nCols - 1; k++)
                {
                    double p = 0.0;
                    for (int i = k; i < nRows; i++)  //选主元
                    {
                        double d = Math.Abs(u[i * nCols + k]);
                        if (d > p)
                        {
                            p = d;
                            w = i;
                        }
                    }

                    if (p == 0)
                        throw new Exception("奇异矩阵解不出啊！");

                    Utility.Swap(pnRows, k, w);

                    for (int i = 0; i < k; i++)
                        Utility.Swap(l, k * nCols + i, w * nCols + i);

                    if (k != w) Parity *= -1;//行变换会改变奇偶性
                    Utility.SwapRow(u, k, w, nCols);
                    for (int i = k + 1; i < nRows; i++)
                    {
                        v = i * nCols + k;
                        l[v] = u[v] / u[k * nCols + k];
                        for (int j = k; j < nCols; j++)
                        {
                            v = i * nCols + j;
                            u[v] -= l[i * nCols + k] * u[k * nCols + j];
                        }
                    }
                }
            }
            return new MatrixLU
            {
                L = L,
                U = U,
                Parity = Parity
            };
        }
        internal MatrixQR QR()
        {
            int nRows = mat.RowCount();
            int nCols = mat.ColumnCount();
            if (nRows < nCols)
                throw new Exception("进行QR变换的矩阵的行数需要大于等于列数");
            double[,] Q = SpecialMatrices.Eye(nRows, nRows);
            double[,] R = mat;
            int n = nCols;
            if (nRows == nCols) n--;
            double d, w;
            int p, u, v;
            double alpha, t;
            unsafe
            {
                fixed (double* q = Q)
                fixed (double* r = R)
                {
                    for (int k = 0; k < n; k++)
                    {
                        d = 0.0;
                        u = k * nCols + k;
                        for (int i = k; i < nRows; i++)
                        {
                            v = i * nCols + k;
                            w = Math.Abs(r[v]);
                            if (w > d)
                                d = w;
                        }

                        alpha = 0.0;
                        for (int i = k; i < nRows; i++)
                        {
                            v = i * nCols + k;
                            t = r[v] / d;
                            alpha += t * t;
                        }

                        if (r[u] > 0.0) d = -d;

                        alpha = d * Math.Sqrt(alpha);
                        if (Math.Abs(alpha) < eps)
                            throw new Exception("解不出……");

                        d = Math.Sqrt(2.0 * alpha * (alpha - r[u]));
                        if (d > eps)
                        {
                            r[u] = (r[u] - alpha) / d;
                            for (int i = k + 1; i < nRows; i++)
                            {
                                p = i * nCols + k;
                                r[p] /= d;
                            }

                            for (int j = 0; j < nRows; j++)
                            {
                                t = 0.0;
                                for (int m = k; m < nRows; m++)
                                    t += r[m * nCols + k] * q[m * nRows + j];

                                for (int i = k; i < nRows; i++)
                                {
                                    p = i * nRows + j;
                                    q[p] -= 2.0 * t * r[i * nCols + k];
                                }
                            }

                            for (int j = k + 1; j < nCols; j++)
                            {
                                t = 0.0;
                                for (int m = k; m < nRows; m++)
                                    t += r[m * nCols + k] * r[m * nCols + j];

                                for (int i = k; i < nRows; i++)
                                {
                                    p = i * nCols + j;
                                    r[p] -= 2.0 * t * r[i * nCols + k];
                                }
                            }

                            r[u] = alpha;
                            for (int i = k + 1; i < nRows; i++)
                                r[i * nCols + k] = 0.0;
                        }

                    }

                    for (int i = 0; i < nRows - 1; i++)
                        for (int j = i + 1; j < nRows; j++)
                        {
                            p = i * nRows + j;
                            u = j * nRows + i;
                            Utility.Swap(q, u, p);
                        }
                }
            }
            return new MatrixQR
            {
                Q = Q,
                R = R
            };
        }
    }
}
