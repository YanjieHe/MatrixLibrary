using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Collections;
using System.Globalization;
using LinearAlgebra.MatrixAlgebra;
namespace LinearAlgebra
{
    public sealed class Matrix : IEnumerable<Double>, IEquatable<Matrix>, IFormattable
    {
        Double[,] elements;
        const Double eps = 1e-9;
        const int DOUBLE_SIZE = sizeof(double);
        private MatrixSubset SubMat;
        /// <summary>
        /// 利用二维Double数组初始化矩阵
        /// </summary>
        /// <param name="Elements"></param>
        public Matrix(Double[,] Elements)
        {
            Initialize(Elements);
        }
        public Matrix(Matrix M)
        {
            Initialize(Utility.CopyMatrix(M.elements));
        }
        private void Initialize(double[,] Elements)
        {
            elements = Elements;
            SubMat = new MatrixSubset(elements);
        }
        /// <summary>
        /// 获取矩阵的行数
        /// </summary>
        public int RowCount
        {
            get { return elements.GetLength(0); }
        }
        /// <summary>
        /// 获取矩阵的列数
        /// </summary>
        public int ColumnCount
        {
            get { return elements.GetLength(1); }
        }
        /// <summary>
        /// 获取矩阵的元素数量
        /// </summary>
        public int Count
        {
            get { return elements.Length; }
        }
        public Double this[int row, int col]
        {
            get { return elements[row, col]; }
            set { elements[row, col] = value; }
        }
        /// <summary>
        /// 返回矩阵元素
        /// </summary>
        /// <returns></returns>
        public Double[,] Elements
        {
            get { return elements; }
            set { elements = value; }
        }
        /// <summary>
        /// 初始化矩阵
        /// </summary>
        /// <param name="nRows">行数</param>
        /// <param name="nCols">列数</param>
        /// <param name="Elements">矩阵元素</param>
        /// <returns></returns>
        public static Matrix Create(int nRows, int nCols, Double[] Elements)
        {
            return Utility.ArrayToMatrix(nRows, nCols, Elements);
        }
        /// <summary>
        /// 初始化矩阵
        /// </summary>
        /// <param name="nRows">行数</param>
        /// <param name="nCols">列数</param>
        /// <param name="Elements">矩阵元素</param>
        /// <returns></returns>
        public static Matrix Create(int nRows, int nCols, IEnumerable<Double> Elements)
        {
            return Utility.IEnumerableToMatrix(nRows, nCols, Elements);
        }
        public static implicit operator Matrix(double[,] Elements)
        {
            return new Matrix(Elements);
        }
        #region Overload operators 
        public static Matrix operator +(Matrix M)
        {
            return new Matrix(M);
        }
        public static Matrix operator -(Matrix M)
        {
            return MatrixComputation.UnaryMinus(M.elements);
        }
        public static Matrix operator +(Matrix A, Matrix B)
        {
            return MatrixComputation.Add(A.elements, B.elements);
        }
        public static Matrix operator -(Matrix A, Matrix B)
        {
            return MatrixComputation.Subtract(A.elements, B.elements);
        }
        public static Matrix operator *(Matrix A, Matrix B)
        {
            return MatrixComputation.Multiply(A.elements, B.elements);
        }
        public static Matrix operator /(Matrix A, Matrix B)
        {
            return MatrixComputation.Divide(A.elements, B.elements);
        }
        public static Matrix operator +(Matrix Mat, double Scalar)
        {
            return MatrixScalar.AddScalar(Mat.elements, Scalar);
        }
        public static Matrix operator +(double Scalar, Matrix Mat)
        {
            return MatrixScalar.AddScalar(Scalar, Mat.elements);
        }
        public static Matrix operator -(Matrix Mat, double Scalar)
        {
            return MatrixScalar.SubtractScalar(Mat.elements, Scalar);
        }
        public static Matrix operator -(double Scalar, Matrix Mat)
        {
            return MatrixScalar.SubtractScalar(Scalar, Mat.elements);
        }
        public static Matrix operator *(Matrix Mat, double Scalar)
        {
            return MatrixScalar.MultiplyScalar(Mat.elements, Scalar);
        }
        public static Matrix operator *(double Scalar, Matrix Mat)
        {
            return MatrixScalar.MultiplyScalar(Scalar, Mat.elements);
        }
        public static Matrix operator /(Matrix Mat, double Scalar)
        {
            return MatrixScalar.DivideScalar(Mat.elements, Scalar);
        }
        public static Matrix operator /(double Scalar, Matrix Mat)
        {
            return MatrixScalar.DivideScalar(Scalar, Mat.Inverse().elements);
        }
        #endregion
        /// <summary>
        /// 全选主元高斯-约当法求逆矩阵
        /// </summary>
        /// <returns>逆矩阵</returns>
        public Matrix Inverse()
        {
            return MatrixComputation.Inverse(elements);
        }
        /// <summary>
        /// 矩阵的转置
        /// </summary>
        /// <returns>返回转置后的矩阵</returns>
        public Matrix Transpose()
        {
            return MatrixComputation.Transpose(elements);
        }
        /// <summary>
        /// 矩阵的幂
        /// </summary>
        /// <param name="Exponent">幂次</param>
        /// <returns></returns>
        public Matrix Power(int Exponent)
        {
            switch (Exponent)
            {
                case -1: return Inverse();
                case 0: return Eye(RowCount, ColumnCount);
                case 1: return new Matrix(this);
                case 2: return this * this;
                case 3: return this * this * this;
                default:
                    if (Exponent < 0)
                        return Inverse().Power(-Exponent);
                    var M = Power(Exponent >> 1);
                    M *= M;
                    return (Exponent & 1) == 0 ? M : this * M;
            }
        }
        /// <summary>
        /// 全选主元高斯消去法求解矩阵的秩
        /// </summary>
        /// <returns>矩阵的秩</returns>
        public int Rank()
        {
            return MatrixComputation.Rank(elements);
        }
        /// <summary>
        /// 矩阵的LU分解
        /// </summary>
        /// <param name="CreateNewInstance">是否在运算中创建新实例</param>
        /// <returns></returns>
        public MatrixLU LU(bool CreateNewInstance = true)
        {
            return new MatrixDecomposition(
                CreateNewInstance ? elements.CopyMatrix() : elements).LU();
        }
        /// <summary>
        /// 矩阵的QR分解（Householder方法）
        /// </summary>
        /// <param name="CreateNewInstance">是否在运算中创建新实例</param>
        /// <returns></returns>
        public MatrixQR QR(bool CreateNewInstance = true)
        {
            return new MatrixDecomposition(
                CreateNewInstance ? elements.CopyMatrix() : elements).QR();
        }
        /// <summary>
        /// 求行列式的值
        /// </summary>
        /// <param name="CreateNewInstance">是否在运算中创建新实例</param>
        /// <returns></returns>
        public double Det(bool CreateNewInstance = true)
        {
            return Det(LU(CreateNewInstance));
        }
        /// <summary>
        /// 求行列式的值
        /// </summary>
        /// <param name="Mat_LU">矩阵的LU分解结果</param>
        /// <returns></returns>
        public static double Det(MatrixLU Mat_LU)
        {
            double det = Mat_LU.Parity;
            var U = Mat_LU.U;
            int n = U.ColumnCount;
            for (int i = 0; i < n; i++)
                det *= U[i, i];
            return det;
        }
        /// <summary>
        /// 利用初等相似变换将一般实矩阵约化为Hessen Berg矩阵
        /// （注：存在不同的变换方式，因此结果并不唯一）
        /// </summary>
        /// <param name="CreateNewInstance">是否在运算中创建新实例</param>
        /// <returns></returns>
        public Matrix HessenBerg(bool CreateNewInstance = true)
        {
            double max, t;
            int i = 0, u, v;
            int nRows = RowCount;
            int nCols = ColumnCount;
            if (nRows != nCols)
                throw new Exception("约化为HessenBerg矩阵的矩阵必须是方阵");
            Matrix H = CreateNewInstance ? new Matrix(this) : this;
            unsafe
            {
                fixed (double* h = H.elements)
                for (int k = 1; k < nCols - 1; k++)
                {
                    max = 0.0;
                    for (int j = k; j < nCols; j++)
                    {
                        u = j * nCols + k - 1;
                        t = h[u];
                        if (Math.Abs(t) > Math.Abs(max))
                        {
                            max = t;
                            i = j;
                        }
                    }

                    if (max != 0.0)
                    {
                        if (i != k)
                        {
                            for (int j = k - 1; j < nCols; j++)
                            {
                                u = i * nCols + j;
                                v = k * nCols + j;
                                Swap(h, u, v);
                            }
                            SwapColumn(h, i, k, nRows, nCols);
                        }

                        for (i = k + 1; i < nCols; i++)
                        {
                            u = i * nCols + k - 1;
                            t = h[u] / max;
                            h[u] = 0.0;
                            for (int j = k; j < nCols; j++)
                            {
                                v = i * nCols + j;
                                h[v] -= t * h[k * nCols + j];
                            }

                            for (int j = 0; j < nCols; j++)
                            {
                                v = j * nCols + k;
                                h[v] += t * h[j * nCols + i];
                            }
                        }
                    }
                }
            }
            return H;
        }
        /// <summary>
        /// 求实矩阵全部特征值
        /// （方法：先变换为Hessen Berg矩阵，再应用带原点位移的双重步QR方法进行迭代）
        /// </summary>
        /// <param name="MaxTimes">最大迭代次数</param>
        /// <param name="Precision">迭代精度</param>
        /// <returns></returns>
        public MatrixEigenValue Eigen(int MaxTimes = 100, double Precision = 0.0001)
        {
            return Eigen(HessenBerg(), MaxTimes, Precision);
        }
        /// <summary>
        /// 求Hessen Berg矩阵的全部特征根
        /// （方法：QR迭代）
        /// </summary>
        /// <param name="matHB">Hessen Berg矩阵</param>
        /// <param name="MaxTimes">最大迭代次数</param>
        /// <param name="Precision">迭代精度</param>
        /// <returns></returns>
        private static MatrixEigenValue Eigen(Matrix matHB, int MaxTimes, double Precision)
        {
            int nRows = matHB.RowCount;
            int nCols = matHB.ColumnCount;
            var Real = new double[nCols];
            var Imaginary = new double[nCols];
            var IterationTimes = 0;
            int m = nCols;
            double p = 0.0, q = 0.0, r = 0.0, xy = 0.0;
            unsafe
            {
                fixed (double* mat = matHB.elements)
                while (m != 0)
                {
                    int t = m - 1;
                    while ((t > 0) &&
                        (Math.Abs(mat[t * nCols + t - 1])
                        > Precision * (Math.Abs(mat[(t - 1) * nCols + t - 1]) + Math.Abs(mat[t * nCols + t]))))
                        t--;

                    int ii = (m - 1) * nCols + (m - 1);
                    int jj = (m - 1) * nCols + (m - 2);
                    int kk = (m - 2) * nCols + (m - 1);
                    int tt = (m - 2) * nCols + (m - 2);

                    if (t == m - 1)
                    {
                        Real[m - 1] = mat[(m - 1) * nCols + (m - 1)];
                        Imaginary[m - 1] = 0.0;
                        m--;
                        IterationTimes = 0;
                    }
                    else if (t == m - 2)
                    {
                        double b = -(mat[ii] + mat[tt]);
                        double c = mat[ii] * mat[tt] - mat[jj] * mat[kk];
                        double w = b * b - 4.0 * c;
                        double y = Math.Sqrt(Math.Abs(w));

                        if (w > 0.0)
                        {
                            xy = 1.0;
                            if (b < 0.0)
                                xy = -1.0;
                            Real[m - 1] = (-b - xy * y) / 2.0;
                            Real[m - 2] = c / Real[m - 1];
                            Imaginary[m - 1] = 0.0;
                            Imaginary[m - 2] = 0.0;
                        }
                        else
                        {
                            Real[m - 1] = -b / 2.0;
                            Real[m - 2] = Real[m - 1];
                            Imaginary[m - 1] = y / 2.0;
                            Imaginary[m - 2] = -Imaginary[m - 1];
                        }

                        m -= 2;
                        IterationTimes = 0;
                    }
                    else
                    {
                        if (IterationTimes > MaxTimes)
                            throw new Exception(string.Format("迭代了{0}次，可还是没解出来……", IterationTimes));

                        IterationTimes++;

                        for (int j = t + 2; j < m; j++)
                            mat[j * nCols + j - 2] = 0.0;
                        for (int j = t + 3; j < m; j++)
                            mat[j * nCols + j - 3] = 0.0;

                        for (int k = t; k < m - 1; k++)
                        {
                            if (k != t)
                            {
                                p = mat[k * nCols + k - 1];
                                q = mat[(k + 1) * nCols + k - 1];
                                r = 0.0;
                                if (k != m - 2)
                                    r = mat[(k + 2) * nCols + k - 1];
                            }
                            else
                            {
                                double x = mat[ii] + mat[tt];
                                double y = mat[tt] * mat[ii] - mat[kk] * mat[jj];

                                ii = t * nCols + t;
                                jj = t * nCols + t + 1;
                                kk = (t + 1) * nCols + t;
                                tt = (t + 1) * nCols + t + 1;

                                p = mat[ii] * (mat[ii] - x)
                                    + mat[jj] * mat[kk] + y;
                                q = mat[kk] * (mat[ii] + mat[tt] - x);
                                r = mat[kk] * mat[(t + 2) * nCols + t + 1];
                            }

                            if ((Math.Abs(p) + Math.Abs(q) + Math.Abs(r)) != 0.0)
                            {
                                xy = 1.0;
                                if (p < 0.0)
                                    xy = -1.0;
                                double s = xy * Math.Sqrt(p * p + q * q + r * r);
                                if (k != t)
                                    mat[k * nCols + k - 1] = -s;

                                double e = -q / s;
                                double f = -r / s;
                                double x = -p / s;
                                double y = -x - f * r / (p + s);
                                double g = e * r / (p + s);
                                double z = -x - e * q / (p + s);

                                for (int j = k; j < m; j++)
                                {
                                    ii = k * nCols + j;
                                    jj = (k + 1) * nCols + j;
                                    p = x * mat[ii] + e * mat[jj];
                                    q = e * mat[ii] + y * mat[jj];
                                    r = f * mat[ii] + g * mat[jj];

                                    if (k != m - 2)
                                    {
                                        kk = (k + 2) * nCols + j;
                                        p += f * mat[kk];
                                        q += g * mat[kk];
                                        r += z * mat[kk];
                                        mat[kk] = r;
                                    }

                                    mat[jj] = q;
                                    mat[ii] = p;
                                }

                                int u = k + 3;
                                if (u >= m - 1)
                                    u = m - 1;

                                for (int i = t; i <= u; i++)
                                {
                                    ii = i * nCols + k;
                                    jj = i * nCols + k + 1;
                                    p = x * mat[ii] + e * mat[jj];
                                    q = e * mat[ii] + y * mat[jj];
                                    r = f * mat[ii] + g * mat[jj];

                                    if (k != m - 2)
                                    {
                                        kk = i * nCols + k + 2;
                                        p += f * mat[kk];
                                        q += g * mat[kk];
                                        r += z * mat[kk];
                                        mat[kk] = r;
                                    }

                                    mat[jj] = q;
                                    mat[ii] = p;
                                }
                            }
                        }
                    }
                }
            }
            return new MatrixEigenValue
            {
                Real = Real,
                Imaginary = Imaginary
            };
        }
        /// <summary>
        /// 零矩阵
        /// </summary>
        /// <param name="nRows">行数</param>
        /// <param name="nCols">列数</param>
        /// <returns></returns>
        public static Matrix Zeros(int nRows, int nCols)
        {
            return new double[nRows, nCols];
        }
        /// <summary>
        /// 零矩阵（方阵）
        /// </summary>
        /// <param name="n">行数与列数</param>
        /// <returns></returns>
        public static Matrix Zeros(int n)
        {
            return Zeros(n, n);
        }
        /// <summary>
        /// 元素全为一的矩阵
        /// </summary>
        /// <param name="nRows">行数</param>
        /// <param name="nCols">列数</param>
        /// <returns></returns>
        public static Matrix Ones(int nRows, int nCols)
        {
            return SpecialMatrices.Ones(nRows, nCols);
        }
        /// <summary>
        /// 元素全为一的矩阵（方阵）
        /// </summary>
        /// <param name="n">行数与列数</param>
        /// <returns></returns>
        public static Matrix Ones(int n)
        {
            return Ones(n, n);
        }
        /// <summary>
        /// 单位矩阵
        /// </summary>
        /// <param name="nRows">行数</param>
        /// <param name="nCols">列数</param>
        /// <returns></returns>
        public static Matrix Eye(int nRows, int nCols)
        {
            return SpecialMatrices.Eye(nRows, nCols);
        }
        /// <summary>
        /// 单位矩阵（方阵）
        /// </summary>
        /// <param name="n">行数与列数</param>
        /// <returns></returns>
        public static Matrix Eye(int n)
        {
            return Eye(n, n);
        }
        public static Matrix Random(int nRows, int nCols)
        {
            return SpecialMatrices.RandomMatrix(nRows, nCols);
        }
        public static Matrix Random(int nRows, int nCols, int seed)
        {
            return SpecialMatrices.RandomMatrix(nRows, nCols, seed);
        }
        /// <summary>
        /// 交换矩阵的两个指定位置的元素
        /// </summary>
        /// <param name="x">指向矩阵的指针</param>
        /// <param name="u"></param>
        /// <param name="v"></param>
        private static unsafe void Swap(double* x, int u, int v)
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
        private static void Swap<T>(T[] x, int i, int j)
        {
            T z = x[i];
            x[i] = x[j];
            x[j] = z;
        }
        private static unsafe void SwapRow(double* x, int i, int j, int nCols)
        {
            int u = i * nCols, v = j * nCols;
            for (int m = 0; m < nCols; m++)//交换行
            {
                Swap(x, u, v);
                u++;
                v++;
            }
        }
        private static unsafe void SwapColumn(double* x, int i, int j, int nRows, int nCols)
        {
            int u = i, v = j;
            for (int m = 0; m < nRows; m++)//交换列
            {
                Swap(x, u, v);
                u += nCols;
                v += nCols;
            }
        }
        /// <summary>
        /// 纵向拼接相同列数的矩阵
        /// </summary>
        /// <param name="Matrices"></param>
        /// <returns></returns>
        public static Matrix rBind(params Matrix[] Matrices)
        {
            return Utility.rBind(Matrices);
        }
        /// <summary>
        /// 横向拼接相同行数的矩阵
        /// </summary>
        /// <param name="Matrices"></param>
        /// <returns></returns>
        public static Matrix cBind(params Matrix[] Matrices)
        {
            return Utility.cBind(Matrices);
        }
        /// <summary>
        /// 根据函数，将原矩阵的值映射到新矩阵上
        /// </summary>
        /// <param name="f">函数</param>
        /// <returns></returns>
        public Matrix Map(Func<double, double> f)
        {
            int m = RowCount;
            int n = ColumnCount;
            int count = Count;
            var res = new double[m, n];
            unsafe
            {
                fixed (double* mat = elements)
                fixed (double* result = res)
                    for (int i = 0; i < count; i++)
                    result[i] = f(mat[i]);
            }
            return res;
        }
        IEnumerator<Double> IEnumerable<Double>.GetEnumerator()
        {
            foreach (var item in elements)
                yield return item;
        }
        IEnumerator IEnumerable.GetEnumerator()
        {
            yield return this.AsEnumerable();
        }
        public IEnumerable<double> GetRow(int index)
        {
            return SubMat.GetRow(index);
        }
        public IEnumerable<double> GetColumn(int index)
        {
            return SubMat.GetColumn(index);
        }
        public IEnumerable<double> GetDiagonal(bool mainDiagonal = true)
        {
            return SubMat.GetDiagonal(mainDiagonal);
        }
        public void SetRow(int index, IEnumerable<double> data)
        {
            SubMat.SetRow(index, data);
        }
        public void SetColumn(int index, IEnumerable<double> data)
        {
            SubMat.SetColumn(index, data);
        }
        public void SetDiagonal(IEnumerable<double> data, bool mainDiagonal = true)
        {
            SubMat.SetDiagonal(data, mainDiagonal);
        }
        /// <summary>
        /// 矩阵判等
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public bool Equals(Matrix other)
        {
            int mA = RowCount, nA = ColumnCount;
            int mB = other.RowCount, nB = other.ColumnCount;
            if (mA == mB && nA == nB)
            {
                unsafe
                {
                    fixed (double* a = elements)
                    fixed (double* b = other.elements)
                        for (int i = Count - 1; i >= 0; i--)
                        if (Math.Abs(a[i] - b[i]) > eps)//考虑浮点数的误差
                            return false;
                }
                return true;
            }
            else return false;
        }
        /// <summary>
        /// 将矩阵以字符串的形式输出
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            return Utility.MatrixToString(elements);
        }
        public string ToString(string format, IFormatProvider formatProvider)
        {
            return Utility.MatrixToString(elements, format, formatProvider);
        }
        public string ToString(string format)
        {
            return Utility.MatrixToString(elements, format, CultureInfo.CurrentCulture);
        }
        // override object.Equals
        public override bool Equals(object obj)
        {
            if (obj == null || GetType() != obj.GetType())
                return false;
            return Equals((Matrix)obj);
        }
        // override object.GetHashCode
        public override int GetHashCode()
        {
            double sum = 0.0;
            foreach (var item in elements)
                sum += Math.Abs(item);
            return (int)Math.Sqrt(sum);
        }
        /// <summary>
        /// 从文本文件中加载矩阵
        /// </summary>
        /// <param name="filePath">文件路径</param>
        /// <param name="encoding">编码</param>
        /// <returns></returns>
        public static Matrix Load(string filePath, Encoding encoding)
        {
            return MatrixIO.Load(filePath, encoding);
        }
        /// <summary>
        /// 从文本文件中加载矩阵
        /// </summary>
        /// <param name="filePath">文件路径</param>
        /// <returns></returns>
        public static Matrix Load(string filePath)
        {
            return Load(filePath, Encoding.Default);
        }
        /// <summary>
        /// 保存矩阵
        /// </summary>
        /// <param name="filePath">文件路径</param>
        public void Save(string filePath)
        {
            MatrixIO.Save(elements, filePath);
        }
        public double[][] ToJaggedArray()
        {
            int nRows = RowCount;
            int nCols = ColumnCount;
            double[][] arr = new double[nRows][];
            for (int i = 0; i < nRows; i++)
            {
                arr[i] = new double[nCols];
                Buffer.BlockCopy(elements, i * nCols * DOUBLE_SIZE, arr[i], 0, nCols * DOUBLE_SIZE);
            }
            return arr;
        }
    }
}
