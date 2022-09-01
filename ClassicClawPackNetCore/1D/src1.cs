using static GeoclawNetCore._1D.Setprob;

namespace ClassicClawPackNetCore._1D
{
    public class Src1
    {
        public Src1(int mx, int mbc, double dx, double dt, double[][] q)
        {
            this.mx = mx;
            this.mbc = mbc;
            this.dx = dx;
            this.q = q;
        }

        int mx, mbc;
        double dx, dt;
        double[][] q;

        ///<remarks>
        ///     source terms for radial symmetry
        ///     ndim should be set in setprob.f
        ///     ndim = 2  for cylindrical symmetry
        ///     ndim = 3  for spherical symmetry
        /// </remarks>
        public void Run()
        {
            for (int i = 0; i < mx + mbc; i++)
            {
                var xcell = xlower + (i - 0.50) * dx;
                q[0][i] = q[0][i] - dt * (ndim - 1) * bulk / xcell * q[1][i];
            }
        }
    }
}
