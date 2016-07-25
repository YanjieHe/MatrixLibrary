using System;
using System.Text;
using LinearAlgebra.MatrixAlgebra;
namespace LinearAlgebra.LinearEquations
{
    public class GaussElimination
    {
        const double eps = 1e-9;
        public double[] X { get; private set; }
        public GaussElimination(Matrix A, double[] B, bool CreateNewInstance = true)
        {
            if (CreateNewInstance)
                X = Solve(A.Elements.CopyMatrix(), B.CopyVector());
            else
                X = Solve(A.Elements, B);
        }
        private static double[] Solve(double[,] A, double[] B)
        {
            int nRows = A.RowCount();
            int nCols = A.ColumnCount();
            if (nRows != nCols) throw new ArgumentException("系数矩阵的行数与列数必须相等");
            int n = nCols;
            int[] pnCols = new int[nCols];
            int pnRow = 0;
            unsafe
            {
                fixed (double* a = A)
                fixed (double* b = B)
                {
                    for (int k = 0; k < n - 1; k++)
                    {
                        double d = 0.0;
                        for (int i = k; i < n; i++)
                        {
                            for (int j = k; j < n; j++)
                            {
                                double t = Math.Abs(a[i * n + j]);
                                if (t > d)
                                {
                                    d = t;
                                    pnCols[k] = j;
                                    pnRow = i;
                                }
                            }
                        }
                        if (Math.Abs(d) < eps) throw new Exception("求解失败……");

                        if (pnCols[k] != k)
                            Utility.SwapColumn(a, k, pnCols[k], nRows, nCols);

                        if (pnRow != k)
                        {
                            for (int j = k; j < n; j++)
                                Utility.Swap(a, k * n + j, pnRow * n + j);

                            Utility.Swap(b, k, pnRow);
                        }

                        d = a[k * n + k];
                        for (int j = k + 1; j < n; j++)
                        {
                            int u = k * n + j;
                            a[u] /= d;
                        }

                        b[k] /= d;
                        for (int i = k + 1; i < n; i++)
                        {
                            for (int j = k + 1; j < n; j++)
                            {
                                int u = i * n + j;
                                a[u] -= a[i * n + k] * a[k * n + j];
                            }
                            b[i] -= a[i * n + k] * b[k];
                        }
                    }

                    double q = a[(n - 1) * n + (n - 1)];
                    if (Math.Abs(q) < eps) throw new Exception("求解失败……");

                    b[n - 1] /= q;
                    for (int i = n - 2; i >= 0; i--)//回代
                    {
                        double sum = 0.0;
                        for (int j = i + 1; j < n; j++)
                        {
                            sum += a[i * n + j] * b[j];
                        }
                        b[i] -= sum;
                    }

                    pnCols[n - 1] = n - 1;
                    for (int k = n - 1; k >= 0; k--)
                        if (pnCols[k] != k)
                            Utility.Swap(b, k, pnCols[k]);
                }
            }
            return B;
        }
        public override string ToString()
        {
            var s = new StringBuilder();
            for (int i = 0; i < X.Length; i++)
                s.Append("X").Append(i).Append(" = ")
                    .AppendFormat("{0:F6}", X[i]).AppendLine();
            return s.ToString();
        }
    }
}
