using System;
namespace LinearAlgebra.MatrixAlgebra
{
    internal sealed class SpecialMatrices
    {
        internal static double[,] Ones(int nRows, int nCols)
        {
            var mat = new double[nRows, nCols];
            unsafe
            {
                fixed (double* m = mat)
                for (int i = nRows * nCols - 1; i >= 0; i--)
                    m[i] = 1.0;
            }
            return mat;
        }
        internal static double[,] Eye(int nRows, int nCols)
        {
            int n = Math.Min(nRows, nCols);
            var mat = new double[nRows, nCols];
            for (int i = 0; i < n; i++)
                mat[i, i] = 1.0;
            return mat;
        }
        internal static double[,] RandomMatrix(int nRows, int nCols, int seed)
        {
            return RandomMatrix(nRows, nCols, new Random(seed));
        }
        internal static double[,] RandomMatrix(int nRows, int nCols)
        {
            return RandomMatrix(nRows, nCols, new Random());
        }
        private static double[,] RandomMatrix(int nRows, int nCols, Random rnd)
        {
            double[,] mat = new double[nRows, nCols];
            int n = nRows * nCols;
            unsafe
            {
                fixed (double* m = mat)
                    for (int i = 0; i < n; i++)
                    m[i] = rnd.NextDouble();
            }
            return mat;
        }
    }
}
