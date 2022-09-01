namespace ClassicClawPackNetCore._1D
{
    public class claw1
    {
        /// <summary>
        ///     claw1 class constructor.
        /// </summary>
        /// <param name="meqn">
        ///     meqn is the number of equations in the system of conservation laws.
        /// </param>
        /// <param name="mwaves"> 
        ///     mwaves is the number of waves that result from the
        ///     solution of each Riemann problem.Often mwaves = meqn but
        ///     for some problems these may be different.
        /// </param>
        /// <param name="maux"></param>
        /// <param name="mbc">
        /// mbc is the number of "ghost cells" that must be added on to each
        /// side of the domain to handle boundary conditions.The cells
        /// actually in the physical domain are labelled from 1 to mx in x.
        /// The arrays are dimensioned actually indexed from 1-mbc to mx+mbc.
        /// For the methods currently implemented, mbc = 2 should be used.
        /// If the user implements another method that has a larger stencil and
        ///       hence requires more ghost cells, a larger value of mbc could be used.
        ///       q is extended from the physical domain to the ghost cells by the
        ///       user-supplied routine bc1.
        /// </param>
        /// <param name="mx">
        /// mx is the number of grid cells in the x-direction, in the
        /// physical domain.In addition there are mbc grid cells
        /// along each edge of the grid that are used for boundary conditions
        /// </param>
        /// <param name="xlower"></param>
        /// <param name="dx">
        /// grid spacing in x.
        /// (for a computation in ax <= x <= bx, set dx = (bx - ax) / mx.)
        /// </param>
        /// <param name="tstart">initial time</param>
        /// <param name="tend">
        /// Desired final time (on input).
        /// If tend<tstart, then claw1 returns after a single successful
        ///                 time step has been taken (single-step mode).
        /// Otherwise, as many steps are taken as needed to reach tend,
        ///                 up to a maximum of nv(1).
        ///         = Actual time reached(on output).
        /// </param>
        /// <param name="mthbc"></param>
        /// <param name="work">
        /// work(mwork) = double precision work array of length at least mwork
        /// </param>
        /// <param name="mwork">
        /// mwork = length of work array.  Must be at least
        /// (mx + 2*mbc) * (2 + 4*meqn + mwaves + meqn* mwaves)
        /// If mwork is too small then the program returns with info = 4
        /// and prints the necessary value of mwork to unit 6.
        /// </param>
        /// <param name="use_fwave"></param>
        /// <param name="mthlim">
        /// mthlim(1:mwaves) = array of values specifying the flux limiter to be used
        ///                     in each wave family mw.Often the same value will be used
        ///                     for each value of mw, but in some cases it may be
        ///                     desirable to use different limiters.For example,
        ///                     for the Euler equations the superbee limiter might be
        /// used for the contact discontinuity (mw= 2) while another
        ///                     limiter is used for the nonlinear waves.Several limiters
        /// are built in and others can be added by modifying the
        /// subroutine philim.
        ///        mthlim(mw) = 0 for no limiter
        ///                   = 1 for minmod
        ///                   = 2 for superbee
        ///                   = 3 for van Leer
        ///                   = 4 for monotonized centered
        /// </param>
        /// <param name="dtv">
        /// array of values related to the time step:
        ///(Note: method(1)=1 indicates variable size time steps)
        /// dtv(1) = value of dt to be used in all steps if method(1) = 0
        ///                = value of dt to use in first step if method(1) = 1
        /// dtv(2) = unused if method(1) = 0.
        ///                = maximum dt allowed if method(1) = 1.
        /// dtv(3) = smallest dt used(on output)
        /// dtv(4) = largest dt used(on output)
        /// dtv(5) = dt used in last step(on output)
        /// </param>
        /// <param name="cflv">
        /// Array of values related to Courant number:
        /// cflv(1) = maximum Courant number to be allowed.With variable
        ///                   time steps the step is repeated if the Courant
        /// number is larger than this value.With fixed time
        ///                   steps the routine aborts.Usually cflv(1) = 1.0 should work.
        ///         cflv(2) = unused if method(1) = 0.
        /// = desired Courant number if method(1) = 1.
        ///                   Should be somewhat less than cflv(1), e.g. 0.9
        ///         cflv(3) = largest Courant number observed(on output).
        ///         cflv(4) = Courant number in last step(on output).
        /// </param>
        /// <param name="nv">
        /// nv(1:2) = array of values related to the number of time steps:
        /// nv(1) = unused if method(1) = 0
        ///       = maximum number of time steps allowed if method(1) = 1
        /// nv(2) = number of time steps taken(on output).
        /// </param>
        /// <param name="method">
        /// method(1:7) = array of values specifying the numerical method to use
        /// method(1) = 0 if fixed size time steps are to be taken.
        ///     In this case, dt = dtv(1) in all steps.
        /// = 1 if variable time steps are to be used.
        ///     In this case, dt = dtv(1) in the first step and
        ///                       thereafter the value cflv(2) is used to choose the
        ///                       next time step based on the maximum wave speed seen
        ///                       in the previous step.Note that since this value
        ///                       comes from the previous step, the Courant number will
        ///                       not in general be exactly equal to the desired value
        ///                       If the actual Courant number in the next step is
        /// greater than 1, then this step is redone with a
        ///                       smaller dt.
        ///
        ///         method(2) = 1 if Godunov's method is to be used, with no 2nd order
        ///                       corrections.
        /// = 2 if second order correction terms are to be added, with
        ///                       a flux limiter as specified by mthlim.
        ///        method(3) is not used in one - dimension.
        ///        method(4) = 0 to suppress printing
        /// = 1 to print dt and Courant number every time step
        ///        method(5) = 0 if there is no source term psi.In this case
        ///        the subroutine src1 is never called so a dummy
        ///                       parameter can be given.
        /// = 1 if there is a source term.In this case
        ///     the subroutine src1 must be provided and a
        ///                       fractional step method is used.
        /// In each time step the following sequence is followed:
        ///     call bc to extend data to ghost cells
        ///                            call step1 to advance hyperbolic eqn by dt
        ///                            call src1 to advance source terms by dt
        /// = 2 if there is a source term and Strang splitting is to
        ///                       be used instead of the Godunov splitting above.
        ///                       In each time step the following sequence is followed:
        ///     call bc to extend data to ghost cells
        ///                            call src1 to advance source terms by dt / 2
        ///                            call step1 to advance hyperbolic equation by dt
        ///                            call src1 to advance source terms by dt / 2
        ///                       For most problems 1 is recommended rather than 2
        ///                       since it is less expensive and works essentially as
        ///                       well on most problems.
        ///     method(6) = 0 if there is no capacity function capa.
        /// = mcapa > 0 if there is a capacity function.In this case
        ///     aux(i, mcapa) is the capacity of the i'th cell and you
        ///                       must also specify method(7).ge.mcapa and set aux.
        ///         method(7) = 0 if there is no aux array used.
        /// = maux > 0  if there are maux auxiliary variables.

        ///     The recommended choice of methods for most problems is
        /// method(1) = 1, method(2) = 2.
        /// </param>
        ///<remarks>
        /// Reference: https://github.com/clawpack/classic/blob/2d0f3d2449d8bd9ea8361031c9d25f44bd1d11f9/src/1d/claw1.f
        ///
        ///     q(meqn, 1-mbc:mx+mbc) 
        /// On input:  initial data at time tstart.
        /// On output: final solution at time tend.
        /// q(m, i) = value of mth component in the i'th cell.
        /// Values within the physical domain are in q(m, i)
        ///                for i = 1,2,...,mx
        ///        mbc extra cells on each end are needed for boundary conditions
        ///        as specified in the routine bc1.
        ///
        ///     aux(maux, 1-mbc:mx+mbc)
        /// Array of auxiliary variables that are used in specifying the problem.
        /// If method(7) = 0 then there are no auxiliary variables and aux
        /// can be a dummy variable.
        ///        If method(7) = maux > 0 then there are maux auxiliary variables
        /// and aux must be dimensioned as above.
        ///
        ///        Capacity functions are one particular form of auxiliary variable.
        /// These arise in some applications, e.g.variable coefficients in
        /// advection or acoustics problems.
        /// See Clawpack Note # 5 for examples.
        ///
        ///        If method(6) = 0 then there is no capacity function.
        /// If method(6) = mcapa > 0  then there is a capacity function and
        /// capa(i), the "capacity" of the i'th cell, is assumed to be
        /// stored in aux(mcapa, i).
        /// In this case we require method(7).ge.mcapa.
        /// </remarks>
        public claw1(int meqn, int mwaves, int maux, int mbc, int mx,
                    double xlower, double dx, double tstart, double tend,
                    double[] mthbc, double[] work, int mwork,int[] mthlim,
                    double[] dtv, double [] cflv, int[] nv, int[] method, bool use_fwave = true
                    )
        {
            this.mthlim = mthlim;
            this.mbc = mbc;
            this.mx = mx;
            this.meqn = meqn;
            q = new double[meqn][];
            aux = new double[maux][];
            this.dx = dx;
            this.tstart = tstart;
            this.tend = tend;
            this.mwork = mwork;
            q = new double[meqn][];
            aux = new double[maux][];
            f = new double[meqn][];
            s = new double[mwaves][];
            wave = new double[mwaves][][];
            qwork = new double[meqn][];
            amdq = new double[meqn][];
            apdq = new double[meqn][];
            dtdx = new double[maxmx + 2 * mbc];
            method = new int[7];
            mthlim = new int[mwaves];
            for (int i = 0; i < amdq.Length; i++) amdq[i] = new double[maxmx + 2 * mbc];
            for (int i = 0; i < apdq.Length; i++) apdq[i] = new double[maxmx + 2 * mbc];
            for (int i = 0; i < s.Length; i++) s[i] = new double[maxmx + 2 * mbc];
            for (int i = 0; i < f.Length; i++) f[i] = new double[maxmx + 2 * mbc];
            for (int i = 0; i < aux.Length; i++) aux[i] = new double[maxmx + 2 * mbc];
            for (int i = 0; i < q.Length; i++) q[i] = new double[maxmx + 2 * mbc];
            for (int i = 0; i < wave.Length; i++) wave[i] = new double[meqn][];
            for (int i = 0; i < wave.Length; i++) for (int j = 0; j < wave[i].Length; j++) wave[i][j] = new double[maxmx + 2 * mbc];
            for (int i = 0; i < qwork.Length; i++) qwork[i] = new double[maxmx + 2 * mbc];

            dtv = new double[5]; // dimension dtv(5)
            cflv = new double[4]; // dimension cflv(4)

            mthlim = new int[mwaves]; // dimension mthlim(mwaves)
            method = new int[7]; // dimension method(7)
            nv = new int[2]; // dimension nv(2)
        }

        int maxmx, meqn, mwaves, mbc, maux, mx, mwork;
        double dx;
        public double[][] q;
        public double[][] aux;
        public double[][] f;
        public double[][] s;
        public double[][][] wave;
        public double[][] qwork;
        public double[][] amdq;
        public double[][] apdq;
        public double[] dtdx;
        public int[] method;
        public int[] mthlim;
        public double tstart, tend;
        public double[] dtv;
        public int[] nv;
        public double[] cflv; // dimension cflv(4)
        // bcL: boundary cells at the valve location(left boundary)
        // bcR: boundary cells at the other end of pipe(right boundary)
        public double[][] bcL;
        public double[][] bcR;


        public int Run()
        {
            var runMainLoop = true;
            var errorInfo = 0;
            var t = tstart;
            var maxn = nv[0];
            var dt = dtv[0]; //# initial dt
            var cflmax = 0.0F;
            var dtmin = dt;
            var dtmax = dt;
            nv[1] = 0;

            //     # partition work array into pieces for passing into step1:
            var i0f = 0; //var i0f = 1; changed to 0 since C3 array index starts from 0
            var i0wave = i0f + (maxmx + 2 * mbc) * meqn;
            var i0s = i0wave + (maxmx + 2 * mbc) * meqn * mwaves;
            var i0dtdx = i0s + (maxmx + 2 * mbc) * mwaves;
            var i0qwork = i0dtdx + (maxmx + 2 * mbc);
            var i0amdq = i0qwork + (maxmx + 2 * mbc) * meqn;
            var i0apdq = i0amdq + (maxmx + 2 * mbc) * meqn;
            i0dtdx = i0apdq + (maxmx + 2 * mbc) * meqn;
            var i0end = i0dtdx + (maxmx + 2 * mbc) - 1;

            // check for errors in data:
            errorInfo = CalculateStepAndCheckforErrors(ref maxn, dt, tstart, tend);
            if (errorInfo != 0)
            {
                runMainLoop = false;
            }
            else
            {
                if (mwork < i0end)
                {
                    errorInfo = 4;
                    runMainLoop = false;
                }
            }

            if (runMainLoop && maxn != 0)
            {
                copyq1 cp1 = new copyq1();
                for (int n = 0; n < maxn; n++)
                {
                    //time at beginning of time step.
                    //adjust dt to hit tend exactly if we're near end of computation
                    //(unless tend < tstart, which is a flag to take only a single step)
                    var told = t; 
                    if (told + dt > tend && tstart < tend) dt = tend - told;

                    if (method[0] == 1)
                    {
                        //save old q in case we need to retake step with smaller dt:
                        cp1.Run(ref q, ref qwork);
                        //copyq1(maxmx, meqn, mbc, mx, q, work(i0qwork));
                    }

                    // midpoint in time for Strang splitting
                    var dt2 = dt / 2.0F;
                    var thalf = t + dt2;
                    t = told + dt; 

                    // extend data from grid to bordering boundary cells:               
                    var bc1 = new bc1(maux, new int[] { 1, 2 }, mx, meqn, mbc); // bc1(maxmx, meqn, mbc, maux, mx, x, q, aux, t, dx)
                    bc1.q = q;
                    bc1.aux = aux;
                    bc1.Run();

                    if (method[4] == 1)
                    {
                        //with source term:   use Strang splitting
                        Src1 src1 = new Src1(mx, mbc, dx, dt, q);
                        src1.Run();
                    }

                    Step1 step1 = new Step1(mwaves, mbc, meqn, maux, mx, (int)dt, (int)dx, true);

                    // passed following parameter
                    // maxmx, cfl, &rp1
                    //remaining paramter
                    step1.q = q;
                    step1.aux = aux;
                    step1.method = method;
                    step1.mthlim = mthlim;
                    step1.f = f;
                    step1.wave = wave;
                    step1.s = s;
                    step1.amdq = amdq;
                    step1.apdq = apdq;
                    step1.dtdx = dtdx;
                    step1.Run();

                    if (method[4] == 1)
                    {
                        Src1 src1 = new Src1(mx, mbc, dx, dt, q);
                        src1.Run();
                    }

                    for (int ieqn = 0; ieqn < meqn; ieqn++)
                    {
                        // first two and last two elements
                        // in fortran code it was updating -1,0 and mx + 1 and mx + 2
                        // this transform in C# 1, 0 and (mx + mbc) + 1 and (mx + mbc) + 2
                        for (int ibc = 1; ibc <= mbc; ibc++)
                        {
                            q[ieqn][2 - ibc] = bcL[ieqn][ibc - 1];
                            q[ieqn][mx + mbc + ibc] = bcR[ieqn][ibc - 1];
                        }
                    }
                }
            }

            if (method[0] == 1 && t < tend && nv[1] == maxn) errorInfo = 11; //too many timesteps
            if (method[0] == 0 && cflmax > cflv[0]) errorInfo = 12; //Courant number too large with fixed dt
            tend = t;
            cflv[2] = cflmax;
            cflv[3] = 0;//cfl;
            dtv[2] = dtmin;
            dtv[3] = dtmax;
            dtv[4] = dt;

            return errorInfo;
        }

        private int CalculateStepAndCheckforErrors(ref int maxn, double dt, double tstart, double tend)
        {
            if (method[0] == 0)
            {
                //fixed size time steps.  Compute the number of steps:
                if (tend < tstart)
                {
                    //single step mode
                    maxn = 1;
                }
                else
                {
                    maxn = (int)Math.Round((tend - tstart) / dt); // Round to nearest int
                    if ((maxn * dt - (tend - tstart)) > (1.0 * Math.Pow(10, -5)) * (tend - tstart))
                    {
                        var dabs = (maxn * dt - (tend - tstart)) * 1.0F;
                        if (dabs > ((1.0 * Math.Pow(10, -5)) * (tend - tstart)))
                        {
                            //dt doesn't divide time interval integer number of times                            
                            return 2;
                        }
                    }
                }
            }
            if (method[0] == 1 && cflv[1] > cflv[0])
            {
                return 3;
            }
            return 0;
        }
    }
}