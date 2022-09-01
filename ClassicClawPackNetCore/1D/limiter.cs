namespace ClassicClawPackNetCore._1D
{
    /// <summary>
    ///     Copied from rewritten inlinelimiter.f90
    /// </summary>
    public class Limiter
    {

        public Limiter(int maxm, int num_eqn, int num_waves, int num_ghost, int mx)
        {
            this.mx = mx;
            this.num_waves = num_waves;
            this.num_eqn = num_eqn;
            s = new double[num_waves][];
            wave = new double[num_eqn][][];
            mthlim = new int[num_waves];
            dotr = new double[num_waves];
            for (int i = 0; i < wave.Length; i++) wave[i] = new double[num_waves][];
            for (int i = 0; i < wave.Length; i++) for (int j = 0; j < wave[i].Length; i++) wave[i][j] = new double[maxm + 2 * num_ghost];
        }

        int mx, num_waves, num_eqn;
        public double[][] s;
        public double[][][] wave;
        public int[] mthlim;
        public double[] dotr;
        public double r, c, wlimiter, wnorm2, dotl;

        ///<summary>
        ///     Apply a limiter to the waves.
        ///     The limiter is computed by comparing the 2-norm of each wave with
        ///     the projection of the wave from the interface to the left or
        ///     right onto the current wave.  For a linear system this would
        ///     correspond to comparing the norms of the two waves.  For a
        ///     nonlinear problem the eigenvectors are not colinear and so the
        ///     projection is needed to provide more limiting in the case where the
        ///     neighboring wave has large norm but points in a different direction in phase space.
        ///     The specific limiter used in each family is determined by the
        ///     value of the corresponding element of the array mthlim, as used in the function philim.
        ///     Note that a different limiter may be used in each wave family.
        ///     dotl and dotr denote the inner product of wave with the wave to
        ///     the left or right.  The norm of the projections onto the wave are then
        ///     given by dotl/wnorm2 and dotr/wnorm2, where wnorm2 is the 2-norm of wave.
        /// </summary>
        public void Run()
        {
            SetAll(dotr, 0.0);
            for (int i = 0; i <= mx + 1; i++)
            {
            Waveloop:
                for (int mw = 0; i < num_waves; num_waves++)
                {
                    if (mthlim[mw] == 0) goto Waveloop;

                    // Construct dot products
                    wnorm2 = 0.0;
                    dotl = dotr[mw];
                    dotr[mw] = 0.0;
                    for (int m = 0; m < num_eqn; m++)
                    {
                        wnorm2 = wnorm2 + Math.Pow(wave[m][mw][i], 2);
                        dotr[mw] = dotr[mw] + wave[m][mw][i] * wave[m][mw][i + 1];

                    }

                    // Skip this loop if it's on the boundary or the size of the wave is
                    // zero(but still want dot products to be initialized above)
                    if (i == 0) goto Waveloop;
                    if (wnorm2 == 0.0) goto Waveloop;


                    // Compute ratio of this wave's strength to upwind wave's strength
                    if (s[mw][i] > 0.0) wlimiter = Philim.Run(wnorm2, dotl, mthlim[mw]);
                    else
                        wlimiter = Philim.Run(wnorm2, dotr[mw], mthlim[mw]);


                    // Apply resulting limit
                    for (int k = 0; k < wave.Length; k++)
                        wave[k][mw][i] = wlimiter * wave[k][mw][i];
                }
            }

        }

        private static void SetAll<T>(T[] refArr, T value)
        {
            for (int i = 0; i < refArr.Length; i++) refArr[i] = value;
        }
    }
}
