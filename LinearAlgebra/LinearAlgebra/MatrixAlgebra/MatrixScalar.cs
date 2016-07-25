namespace LinearAlgebra.MatrixAlgebra
{
    internal static class MatrixScalar
    {
        internal static double[,] AddScalar(double[,] Mat, double Scalar)
        {
            int m = Mat.RowCount(), n = Mat.ColumnCount();
            var result = new double[m, n];
            unsafe
            {
                fixed (double* mat = result)
                fixed (double* a = Mat)
                    for (int i = m * n - 1; i >= 0; i--)
                    mat[i] = a[i] + Scalar;
            }
            return result;
        }
        internal static double[,] AddScalar(double Scalar, double[,] Mat)
        {
            return AddScalar(Mat, Scalar);
        }
        internal static double[,] SubtractScalar(double[,] Mat, double Scalar)
        {
            int m = Mat.RowCount(), n = Mat.ColumnCount();
            var result = new double[m, n];
            unsafe
            {
                fixed (double* mat = result)
                fixed (double* a = Mat)
                    for (int i = m * n - 1; i >= 0; i--)
                    mat[i] = a[i] - Scalar;
            }
            return result;
        }
        internal static double[,] SubtractScalar(double Scalar, double[,] Mat)
        {
            int m = Mat.RowCount(), n = Mat.ColumnCount();
            var result = new double[m, n];
            unsafe
            {
                fixed (double* mat = result)
                fixed (double* b = Mat)
                    for (int i = m * n - 1; i >= 0; i--)
                    mat[i] = Scalar - b[i];
            }
            return result;
        }
        internal static double[,] MultiplyScalar(double[,] Mat, double Scalar)
        {
            int m = Mat.RowCount(), n = Mat.ColumnCount();

            var result = new double[m, n];
            unsafe
            {
                fixed (double* mat = result)
                fixed (double* a = Mat)
                    for (int i = m * n - 1; i >= 0; i--)
                    mat[i] = a[i] * Scalar;
            }
            return result;
        }
        internal static double[,] MultiplyScalar(double Scalar, double[,] Mat)
        {
            return MultiplyScalar(Mat, Scalar);
        }
        internal static double[,] DivideScalar(double[,] Mat, double Scalar)
        {
            return MultiplyScalar(Mat, 1.0 / Scalar);
        }
        internal static double[,] DivideScalar(double Scalar, double[,] InverseMat)
        {
            //注意这里的第二个参数需要是原矩阵的逆矩阵
            return MultiplyScalar(Scalar, InverseMat);
        }
    }
}
