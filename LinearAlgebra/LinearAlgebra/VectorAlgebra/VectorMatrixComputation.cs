using System;
namespace LinearAlgebra.VectorAlgebra
{
    internal static class VectorMatrixComputation
    {
        internal static Vector RowVectorMultiplyMatrix(Vector X, Matrix Y)
        {
            if (X.vectorType == VectorType.Row)
            {
                if (X.Count != Y.RowCount) throw new ArgumentOutOfRangeException();
                int n = X.Count;
                int nCols = Y.ColumnCount;
                var Z = new double[nCols];
                unsafe
                {
                    fixed (double* x = X.Elements)
                    fixed (double* y = Y.Elements)
                    fixed (double* z = Z)
                    for (int j = 0; j < nCols; j++)
                    {
                        double sum = 0.0;
                        for (int i = 0; i < n; i++)
                            sum += x[i] * y[i * nCols + j];
                        z[j] = sum;
                    }
                }
                return new Vector(Z, VectorType.Row);
            }
            else throw new ArgumentOutOfRangeException();
        }
        internal static Matrix ColumnVectorMultiplyMatrix(Vector X, Matrix Y)
        {
            if (X.vectorType == VectorType.Column)
            {
                if (Y.RowCount != 1) throw new ArgumentOutOfRangeException();
                int nRows = X.Count;
                int nCols = Y.ColumnCount;
                double[,] Z = new double[nRows, nCols];
                unsafe
                {
                    fixed (double* x = X.Elements)
                    fixed (double* y = Y.Elements)
                    fixed (double* z = Z)
                    for (int j = 0; j < nCols; j++)
                        for (int i = 0; i < nRows; i++)
                            z[i * nCols + j] = x[i] * y[j];
                }
                return Z;
            }
            else throw new ArgumentOutOfRangeException();
        }
        internal static Matrix MatrixMultiplyRowVector(Matrix X, Vector Y)
        {
            if (Y.vectorType == VectorType.Row)
            {
                if (X.ColumnCount != 1) throw new ArgumentOutOfRangeException();
                int nRows = X.RowCount;
                int nCols = Y.Count;
                double[,] Z = new double[nRows, nCols];
                unsafe
                {
                    fixed (double* x = X.Elements)
                    fixed (double* y = Y.Elements)
                    fixed (double* z = Z)
                    for (int j = 0; j < nCols; j++)
                        for (int i = 0; i < nRows; i++)
                            z[i * nCols + j] = x[i] * y[j];
                }
                return Z;
            }
            else throw new ArgumentOutOfRangeException();
        }
        internal static Vector MatrixMultiplyColumnVector(Matrix X, Vector Y)
        {
            if (Y.vectorType == VectorType.Column)
            {
                if (X.ColumnCount != Y.Count) throw new ArgumentOutOfRangeException();
                int nRows = X.RowCount;
                int nCols = X.ColumnCount;
                double[] Z = new double[nRows];
                unsafe
                {
                    fixed (double* x = X.Elements)
                    fixed (double* y = Y.Elements)
                    fixed (double* z = Z)
                    for (int i = 0; i < nRows; i++)
                    {
                        double sum = 0.0;
                        for (int j = 0; j < nCols; j++)
                            sum += x[i * nCols + j] * y[j];
                        Z[i] = sum;
                    }
                }
                return new Vector(Z, VectorType.Column);
            }
            else throw new ArgumentOutOfRangeException();
        }
    }
}
