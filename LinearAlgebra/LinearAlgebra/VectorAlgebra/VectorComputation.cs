using System;
namespace LinearAlgebra.VectorAlgebra
{
    internal static class VectorComputation
    {
        internal static Vector UnaryMinus(Vector X)
        {
            int n = X.Count;
            var Z = new double[n];
            unsafe
            {
                fixed (double* x = X.Elements)
                fixed (double* z = Z)
                for (int i = 0; i < n; i++)
                    z[i] = -x[i];
            }
            return new Vector(Z, X.vectorType);
        }
        internal static Vector Add(Vector X, Vector Y)
        {
            VectorCheck(X, Y);
            int n = X.Count;
            var Z = new double[n];
            unsafe
            {
                fixed (double* x = X.Elements)
                fixed (double* y = Y.Elements)
                fixed (double* z = Z)
                for (int i = 0; i < n; i++)
                    z[i] = x[i] + y[i];
            }
            return new Vector(Z, X.vectorType);
        }
        internal static Vector Subtract(Vector X, Vector Y)
        {
            VectorCheck(X, Y);
            int n = X.Count;
            var Z = new double[n];
            unsafe
            {
                fixed (double* x = X.Elements)
                fixed (double* y = Y.Elements)
                fixed (double* z = Z)
                for (int i = 0; i < n; i++)
                    z[i] = x[i] - y[i];
            }
            return new Vector(Z, X.vectorType);
        }
        internal static Vector Multiply(Vector X, Vector Y)
        {
            VectorCheck(X, Y);
            int n = X.Count;
            var Z = new double[n];
            unsafe
            {
                fixed (double* x = X.Elements)
                fixed (double* y = Y.Elements)
                fixed (double* z = Z)
                for (int i = 0; i < n; i++)
                    z[i] = x[i] * y[i];
            }
            return new Vector(Z, X.vectorType);
        }
        internal static Vector Divide(Vector X, Vector Y)
        {
            VectorCheck(X, Y);
            int n = X.Count;
            var Z = new double[n];
            unsafe
            {
                fixed (double* x = X.Elements)
                fixed (double* y = Y.Elements)
                fixed (double* z = Z)
                for (int i = 0; i < n; i++)
                    z[i] = x[i] / y[i];
            }
            return new Vector(Z, X.vectorType);
        }
        internal static Vector Map(Vector X, Func<double, double> func)
        {
            int n = X.Count;
            var Z = new double[n];
            unsafe
            {
                fixed (double* x = X.Elements)
                fixed (double* z = Z)
                for (int i = 0; i < n; i++)
                    z[i] = func(x[i]);
            }
            return new Vector(Z, X.vectorType);
        }
        private static void VectorCheck(Vector X, Vector Y)
        {
            if (X.Count != Y.Count)
                throw new ArgumentException("进行运算的两个向量的长度必须一致");
            if (X.vectorType != Y.vectorType)
                throw new ArgumentException("进行运算的两个向量的类型必须一致");
        }
        internal static double DotProduct(Vector X, Vector Y)
        {
            VectorCheck(X, Y);
            double sum = 0.0;
            unsafe
            {
                fixed (double* x = X.Elements)
                fixed (double* y = Y.Elements)
                for (int i = X.Count - 1; i >= 0; i--)
                    sum += x[i] * y[i];
            }
            return sum;
        }
    }
}
