using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GeoclawNetCore._1D
{
    public static class JaggedArrayExtension
    {
        public static string WriteLine<T>(this T[][] arr)
        {
            int rowCount = arr.Length;
            int columnCount = arr[0].Length;
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < rowCount; i++)
            {
                for (int j = 0; j < columnCount; j++)
                {
                    sb.Append($" {arr[i][j]} ");
                }
                sb.AppendLine();
            }
            return sb.ToString();
        }

        public static string WriteLine<T>(this T[] arr)
        {
            int rowCount = arr.Length;
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < rowCount; i++)
            {
                sb.Append($" {arr[i]} ");
            }
            return sb.ToString();
        }

        public static T[][] TransposeRowsAndColumns<T>(this T[][] arr)
        {
            int rowCount = arr.Length;
            int columnCount = arr[0].Length;
            T[][] transposed = new T[columnCount][];
            if (rowCount == columnCount)
            {
                transposed = (T[][])arr.Clone();
                for (int i = 1; i < rowCount; i++)
                {
                    for (int j = 0; j < i; j++)
                    {
                        T temp = transposed[i][j];
                        transposed[i][j] = transposed[j][i];
                        transposed[j][i] = temp;
                    }
                }
            }
            else
            {
                for (int column = 0; column < columnCount; column++)
                {
                    transposed[column] = new T[rowCount];
                    for (int row = 0; row < rowCount; row++)
                    {
                        transposed[column][row] = arr[row][column];
                    }
                }
            }
            return transposed;
        }

        /// <summary>
        /// This transformation assumes that the array length
        /// are same for every sub elements
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="arr"></param>
        /// <returns></returns>
        public static T[][][] TransposeRowsAndColumns<T>(this T[][][] arr)
        {
            T[][][] newArray;

            int mx = arr.Length;
            int mv = arr[0].Length;
            int m = arr[0][0].Length;

            newArray = new T[m][][];
            for (int i = 0; i < m; i++)
            {
                newArray[i] = new T[mv][];
                for (int j = 0; j < mv; j++)
                {
                    newArray[i][j] = new T[mx];
                    for (int k = 0; k < mx; k++)
                    {
                        newArray[i][j][k] = arr[k][j][i];
                    }
                }
            }
            return newArray;
        }
    }
}
