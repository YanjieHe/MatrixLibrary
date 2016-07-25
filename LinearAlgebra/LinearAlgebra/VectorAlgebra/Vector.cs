using System;
using System.Collections;
using System.Collections.Generic;
using LinearAlgebra.MatrixAlgebra;
namespace LinearAlgebra.VectorAlgebra
{
    public enum VectorType
    {
        Row, Column
    };
    public sealed class Vector : IEquatable<Vector>, IEnumerable<double>
    {
        double[] elements;
        public double[] Elements { get { return elements; } set { elements = value; } }
        VectorType vType;
        public VectorType vectorType { get { return vType; } set { vType = value; } }
        public int Count { get { return elements.Length; } }
        public double this[int index] { get { return elements[index]; } set { elements[index] = value; } }
        public Vector(double[] _elements, VectorType _vectorType)
        {
            elements = _elements;
            vType = _vectorType;
        }
        public static Vector operator +(Vector X)
        {
            return new Vector(X.elements.CopyVector(), X.vType);
        }
        public static Vector operator -(Vector X)
        {
            return VectorComputation.UnaryMinus(X);
        }
        public static Vector operator +(Vector X, Vector Y)
        {
            return VectorComputation.Add(X, Y);
        }
        public static Vector operator -(Vector X, Vector Y)
        {
            return VectorComputation.Subtract(X, Y);
        }
        public static Vector operator *(Vector X, Vector Y)
        {
            return VectorComputation.Multiply(X, Y);
        }
        public static Vector operator /(Vector X, Vector Y)
        {
            return VectorComputation.Divide(X, Y);
        }
        public Vector Map(Func<double, double> func)
        {
            return VectorComputation.Map(this, func);
        }
        public static double DotProduct(Vector X, Vector Y)
        {
            return VectorComputation.DotProduct(X, Y);
        }
        public static Vector RowVecMulMat(Vector X, Matrix Y)
        {
            return VectorMatrixComputation.RowVectorMultiplyMatrix(X, Y);
        }
        public static Matrix ColVecMulMat(Vector X, Matrix Y)
        {
            return VectorMatrixComputation.ColumnVectorMultiplyMatrix(X, Y);
        }
        public static Matrix MatMulRowVec(Matrix X, Vector Y)
        {
            return VectorMatrixComputation.MatrixMultiplyRowVector(X, Y);
        }
        public static Vector MatMulColVec(Matrix X, Vector Y)
        {
            return VectorMatrixComputation.MatrixMultiplyColumnVector(X, Y);
        }
        public override string ToString()
        {
            if (vType == VectorType.Row) return string.Join("  ", elements);
            else return string.Join("\r\n", elements);
        }
        public bool Equals(Vector other)
        {
            const double eps = 1e-09;
            if (vType != other.vType) return false;
            int n = Count;
            if (n != other.Count) return false;
            unsafe
            {
                fixed (double* x = elements)
                fixed (double* y = other.elements)
                for (int i = 0; i < n; i++)
                    if (Math.Abs(x[i] - y[i]) > eps)
                        return false;
            }
            return true;
        }
        public static Vector Range(VectorType _vectorType, double start, int count, double step = 1.0)
        {
            var v = new double[count];
            for (int i = 0; i < count; i++, start += step)
                v[i] = start;
            return new Vector(v, _vectorType);
        }
        public IEnumerator<double> GetEnumerator()
        {
            foreach (var item in elements) yield return item;
        }
        IEnumerator IEnumerable.GetEnumerator()
        {
            foreach (var item in elements) yield return item;
        }
    }
}
