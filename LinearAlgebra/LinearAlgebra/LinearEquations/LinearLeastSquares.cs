using System;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using LinearAlgebra.MatrixAlgebra;
namespace LinearAlgebra.LinearEquations
{
    public class LinearLeastSquares
    {
        public double[] Beta { get; private set; }
        public double TSS { get; private set; }
        public double ESS { get; private set; }
        public double RSS { get; private set; }
        public int VariableCount { get; private set; }
        public int Obs { get; private set; }
        public double F
        {
            get
            {
                return (RSS * (Obs - VariableCount - 1)) / (ESS * VariableCount);
            }
        }
        public double Rsq
        {
            get
            {
                return RSS / TSS;
            }
        }
        public double AdjRsq
        {
            get
            {
                int n = Obs;
                int k = VariableCount;
                return 1.0 - (1.0 - Rsq) * ((n - 1.0) / (n - k - 1.0));
            }
        }
        public LinearLeastSquares(Matrix X, double[] Y)
        {
            int xRows = X.RowCount;
            int xCols = X.ColumnCount;
            int yLength = Y.Length;
            if (xRows != Y.Length) throw new ArgumentException("X矩阵与Y矩阵的行数需要相同");
            VariableCount = xCols;
            Obs = xRows;
            TSS = SST(Y);
            var XY = BindXY(X, Y);
            var R = QR(XY);
            Beta = BackSolve(R, xCols + 1);//加上截距项
            ESS = SSE(R, xCols);
            RSS = TSS - ESS;
        }
        private static double[] BackSolve(double[,] R, int n)//回代
        {
            var beta = new double[n];
            beta[n - 1] = R[n - 1, n] / R[n - 1, n - 1];//Y的值在R矩阵的第n+1列上，下标为n
            for (int k = n - 2; k >= 0; k--)
            {
                double sum = 0.0;
                for (int j = k + 1; j < n; j++)
                    sum += R[k, j] * beta[j];
                beta[k] = (R[k, n] - sum) / R[k, k];
            }
            return beta;
        }
        private static double[,] BindXY(Matrix X, double[] Y)//拼接X矩阵与Y矩阵
        {
            int xRows = X.RowCount;
            int xCols = X.ColumnCount;
            double[,] XY = new double[xRows, xCols + 2];
            unsafe
            {
                Parallel.For(0, xRows, i =>
                {
                    fixed (double* x = X.Elements)
                    fixed (double* xy = XY)
                    fixed (double* y = Y)
                    {
                        int xyStart = i * (xCols + 2);
                        int xStart = i * xCols;
                        int xEnd = xStart + xCols;
                        xy[xyStart] = 1.0;/*添加截距项，X矩阵的第一列全为1*/
                        xyStart++;
                        for (int k = xStart; k < xEnd; k++, xyStart++)
                            xy[xyStart] = x[k];
                        xy[xyStart] = y[i];
                    }
                });
            }
            return XY;
        }
        public override string ToString()
        {
            var s = new StringBuilder(string.Concat("Y = ", Beta[0].ToString("F6")));
            for (int i = 1; i < Beta.Length; i++)//使输出结果更清晰
                s.Append(" + ").Append(Beta[i].ToString("F6")).Append(" X").Append(i);
            s
                .AppendLine()
                .AppendLine("--------------------------------------------")
                .AppendFormat("Obs = {0}\r\n", Obs)
                .AppendFormat("TSS = {0:F6}\t", TSS).AppendFormat("RSS = {0:F6}\t", RSS).AppendFormat("ESS = {0:F6}\r\n", ESS)
                .AppendFormat("F = {0:F6}\r\n", F)
                .AppendFormat("R^2 = {0:F6}\tAdjusted R^2 = {1:F6}", Rsq, AdjRsq);
            return s.ToString();
        }
        private static double[,] QR(double[,] mat)//QR分解
        {
            const Double eps = 1e-9;
            int nRows = mat.RowCount();
            int nCols = mat.ColumnCount();
            if (nRows < nCols)
                throw new Exception("进行QR变换的矩阵的行数需要大于等于列数");
            double[,] R = mat;
            int n = nCols;
            if (nRows == nCols) n--;
            unsafe
            {
                fixed (double* r = R)
                for (int k = 0; k < n; k++)
                {
                    double d = 0.0;
                    int u = k * nCols + k;
                    for (int i = k; i < nRows; i++)
                    {
                        int v = i * nCols + k;
                        double w = Math.Abs(r[v]);
                        if (w > d) d = w;
                    }

                    double alpha = 0.0;
                    for (int i = k; i < nRows; i++)
                    {
                        int v = i * nCols + k;
                        double t = r[v] / d;
                        alpha += t * t;
                    }

                    if (r[u] > 0.0) d = -d;

                    alpha = d * Math.Sqrt(alpha);
                    if (Math.Abs(alpha) < eps) throw new Exception("解不出……");

                    d = Math.Sqrt(2.0 * alpha * (alpha - r[u]));
                    if (d > eps)
                    {
                        r[u] = (r[u] - alpha) / d;
                        for (int i = k + 1; i < nRows; i++)
                            r[i * nCols + k] /= d;

                        for (int j = k + 1; j < nCols; j++)
                        {
                            double t = 0.0;
                            for (int m = k; m < nRows; m++)
                                t += r[m * nCols + k] * r[m * nCols + j];

                            for (int i = k; i < nRows; i++)
                                r[i * nCols + j] -= 2.0 * t * r[i * nCols + k];
                        }
                        r[u] = alpha;
                        for (int i = k + 1; i < nRows; i++)
                            r[i * nCols + k] = 0.0;
                    }
                }
            }
            return R;//返回R矩阵
        }
        private static double SST(double[] Y)
        {
            double yBar = Y.Average();
            double sum = 0.0;
            unsafe
            {
                fixed (double* y = &Y[0])
                for (int i = Y.Length - 1; i >= 0; i--)
                {
                    double t = (y[i] - yBar);
                    sum += t * t;
                }
            }
            return sum;
        }
        private static double SSE(double[,] R, int xCols)
        {
            double sum = 0.0;
            int nRows = R.RowCount();
            int nCols = R.ColumnCount();
            unsafe
            {
                fixed (double* r = R)
                for (int i = xCols + 1; i < nRows; i++)
                {
                    int u = i * nCols + (nCols - 1);
                    double t = r[u];
                    sum += t * t;
                }
            }
            return sum;
        }
    }
}
