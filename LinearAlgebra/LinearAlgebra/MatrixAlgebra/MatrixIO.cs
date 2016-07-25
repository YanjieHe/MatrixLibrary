using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
namespace LinearAlgebra.MatrixAlgebra
{
    internal static class MatrixIO
    {
        internal static Matrix Load(string filePath, Encoding encoding)
        {
            using (var sr = new StreamReader(filePath, encoding))
            {
                string[] arr = sr.ReadLine().Split('\t', ',');
                int cols = arr.Length;//以第一行的列数为准
                List<double> list = new List<double>(arr.Select(double.Parse));
                int rows = 1;
                while (sr.Peek() != -1)
                {
                    string line = sr.ReadLine();
                    if (string.IsNullOrWhiteSpace(line)) continue;
                    arr = line.Split('\t', ',');
                    if (arr.Length != cols) throw new Exception("读取的矩阵列数不一致");
                    list.AddRange(arr.Select(double.Parse));
                    rows++;
                }
                return Utility.IEnumerableToMatrix(rows, cols, list);
            }
        }
        internal static void Save(double[,] mat, string filePath)
        {
            using (var sw = new StreamWriter(filePath))
            {
                int m = mat.RowCount();
                int n = mat.ColumnCount() - 1;
                for (int i = 0; i < m; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        sw.Write(mat[i, j]);
                        sw.Write('\t');
                    }
                    sw.WriteLine(mat[i, n]);
                }
            }
        }
    }
}
