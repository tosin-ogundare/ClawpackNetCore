using static GeoclawNetCore._1D.Setprob;

namespace ClassicClawPackNetCore._1D
{
    public class Src1
    {
        public Src1(int mx, int mbc, double dx, double dt, double[][] q, double xlower_local)
        {
            this.mx = mx;
            this.mbc = mbc;
            this.dx = dx;
            this.dt = dt;
            this.q = q;
            this.xlower_local = xlower_local == 0 ? xlower: xlower_local;
        }

        int mx, mbc;
        double dx, dt;
        public double[][] q;
        double xlower_local;
        ///<remarks>
        ///     source terms for radial symmetry
        ///     ndim should be set in setprob.f
        ///     ndim = 2  for cylindrical symmetry
        ///     ndim = 3  for spherical symmetry
        /// </remarks>
        public void Run()
        {
            for (int i = 0; i < mx + 2 * mbc; i++)
            {
                var xcell = xlower_local + (i - 0.50) * dx;
                q[0][i] = q[0][i] - dt * (ndim - 1) * bulk / xcell * q[1][i];
            }
        }
    }
}
