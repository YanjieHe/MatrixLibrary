using System;
using LinearAlgebra.MatrixAlgebra;
namespace LinearAlgebra.LinearEquations
{
    public class ConjugateGradient
    {
        public double[] X { get; private set; }
        /// <summary>
        /// 求解对称正定方程组的共轭梯度法
        /// A * X = B
        /// </summary>
        /// <param name="A">系数矩阵</param>
        /// <param name="B">常数向量</param>
        public ConjugateGradient(Matrix A, double[] B)
        {
            X = Solve(A.Elements, B);
        }
        private static double[] Solve(double[,] A, double[] B)
        {
            const Double eps = 1e-9;
            int nRows = A.RowCount();
            int nCols = A.ColumnCount();
            if (nRows != nCols) throw new ArgumentException();
            int n = nCols;//未知数个数
            double[] x = new double[n];
            double[] p = new double[n];
            double[] r = new double[n];
            for (int i = 0; i < n; i++)
            {
                double t = B[i];
                p[i] = t;
                r[i] = t;
            }
            for (int i = 0; i < n; i++)
            {
                double[] s = MatrixMultiplyVector(A, p);
                double d = 0.0;
                double e = 0.0;
                for (int k = 0; k < n; k++)
                {
                    d += p[k] * B[k];
                    e += p[k] * s[k];
                }

                double alpha = d / e;
                for (int k = 0; k < n; k++)
                    x[k] += alpha * p[k];

                double[] q = MatrixMultiplyVector(A, x);
                d = 0.0;
                for (int k = 0; k < n; k++)
                {
                    r[k] = B[k] - q[k];
                    d += r[k] * s[k];
                }

                double beta = d / e;
                d = 0.0;
                for (int k = 0; k < n; k++)
                    d += r[k] * r[k];

                d = Math.Sqrt(d);
                if (d < eps) return x;
                for (int k = 0; k < n; k++)
                    p[k] = r[k] - beta * p[k];
            }
            return x;
        }
        private static double[] MatrixMultiplyVector(double[,] A, double[] B)
        {
            if (A.RowCount() != A.ColumnCount()) throw new ArgumentException();
            if (A.RowCount() != B.Length) throw new ArgumentException();
            int n = B.Length;
            var ans = new double[n];
            unsafe
            {
                fixed (double* a = A)
                fixed (double* b = B)
                for (int i = 0; i < n; i++)
                {
                    double sum = 0.0;
                    int i_n = i * n;
                    for (int j = 0; j < n; j++)
                        sum += a[i_n + j] * b[j];
                    ans[i] = sum;
                }
            }
            return ans;
        }
    }
}
