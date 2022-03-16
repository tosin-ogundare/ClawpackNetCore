namespace ClassicClawPackNetCore._1D
{
    /// <summary>
    ///     Standard boundary condition choices for claw2
    /// </summary>
    public class bc1
    {
        public bc1(int maux, int[] mthbc, int mx = 4, int meqn = 4, int mbc = 2)
        {
            q = new double[meqn][];
            aux = new double[maux][];
            for (int i = 0; i < meqn; i++) q[i] = new double[mx + 2 * mbc];
            for (int i = 0; i < maux; i++) aux[i] = new double[mx + 2 * mbc];
            this.mthbc = mthbc;
            this.mbc = mbc;
            this.mx = mx;
            this.meqn = meqn;
        }
        int mbc, mx, meqn;
        public double[][] q;
        public double[][] aux;
        public int[] mthbc = new int[2];

        /// <remarks>
        /// At each boundary  k = 1 (left),  2 (right):
        ///   mthbc(k) =  0  for user-supplied BC's (must be inserted!)
        ///            =  1  for zero-order extrapolation
        ///            =  2  for periodic boundary coniditions
        ///            =  3  for solid walls, assuming this can be implemented
        ///                  by reflecting the data about the boundary and then
        ///                  negating the 2'nd component of q.
        /// Extend the data from the computational region
        ///      i = 1, 2, ..., mx2
        /// to the virtual cells outside the region, with
        ///      i = 1-ibc  and   i = mx+ibc   for ibc=1,...,mbc
        /// </remarks>
        public void Run()
        {
            int primaryTypeId = mthbc[0] + 1;
            int secondaryTypeId = mthbc[1] + 1;

            switch (primaryTypeId)
            {
                default: InvalidCase(); break;
                case 2: ZeroOrderExtrapolation(mbc, meqn, q); break;
                case 3: Periodic(mbc, meqn, mx, q); break;
                case 4: SolidWall(mbc, meqn, q); break;

            };

            switch (secondaryTypeId)
            {
                default: InvalidCase(); break;
                case 2: ZeroOrderExtrapolationSec(mbc, meqn, mx, q); break;
                case 3: PeriodicSec(mbc, meqn, mx, q); break;
                case 4: SolidWallSec(mbc, meqn, mx, q); break;

            };
        }

        private static readonly Action<int, int, double[][]> ZeroOrderExtrapolation = (mbc, meqn, q) => {
            for (int ibc = 0; ibc < mbc; ibc++) for (int m = 0; m < meqn; m++) q[m][1 - ibc] = q[m][2];
        };

        private static readonly Action<int, int, int, double[][]> ZeroOrderExtrapolationSec = (mbc, meqn, mx, q) => {
            for (int ibc = 0; ibc < mbc; ibc++) for (int m = 0; m < meqn; m++) q[m][mx + ibc] = q[m][mx];
        };

        private static Action InvalidCase = () => throw new InvalidOperationException("*** ERROR *** mthbc(1)=0 || 1 and no BCs specified in bc1 or bc2");

        private static readonly Action<int, int, int, double[][]> Periodic = (mbc, meqn, mx, q) => {
            for (int ibc = 0; ibc < mbc; ibc++) for (int m = 0; m < meqn; m++) q[m][1 - ibc] = q[m][mx+1-ibc];
        };

        private static readonly Action<int, int, int, double[][]> PeriodicSec = (mbc, meqn, mx, q) => {
            for (int ibc = 0; ibc < mbc; ibc++) for (int m = 0; m < meqn; m++) q[m][mx + ibc] = q[m][ibc];
        };

        private static readonly Action<int, int, double[][]> SolidWall = (mbc, meqn, q) => {
            for (int ibc = 0; ibc < mbc; ibc++)
            {
                for (int m = 0; m < meqn; m++) q[m][1 - ibc] = q[m][ibc];
                q[1][1 - ibc] = -q[2][ibc];
            }
        };

        private static readonly Action<int, int, int, double[][]> SolidWallSec = (mbc, meqn, mx,  q) => {
            for (int ibc = 0; ibc < mbc; ibc++)
            {
                for (int m = 0; m < meqn; m++) q[m][mx + ibc] = q[m][mx+1-ibc];
                q[1][mx + ibc] = -q[2][mx + 1 - ibc];
            }
        };


    }
}
