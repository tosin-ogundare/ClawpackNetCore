namespace ClassicClawPackNetCore._1D
{
    /// <summary>
    ///     copy the contents of q1 into q2
    /// </summary>
    public class copyq1
    {

        public void Run(ref double[][] q1, ref double[][] q2)
        {
            for (int i = 0; i < q1.Length; i++) for (int j = 0; j < q1[i].Length; j++) q2[i][j] = q1[i][j];
        }
    }
}
