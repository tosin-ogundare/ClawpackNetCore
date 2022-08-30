using DotNumerics.LinearAlgebra.CSLapack;
using static GeoclawNetCore._1D.Setprob;

namespace GeoclawNetCore._1D
{
    /// <summary>
    ///     Riemann solver for the two layer shallow water equations.
    /// </summary>
    public class rp1_inundation
    {
        public rp1_inundation(double[][] ql, double[][] qr, double[][] auxl, double[][] auxr, int mx, int maxmx, int meqn = 4, int mwaves = 4, int mbc = 2)
        {
            fwave = new double[maxmx + 2 * mbc][][];
            s = new double[maxmx + 2 * mbc][];
            eig_vec = new double[meqn][];
            amdq = new double[maxmx + 2 * mbc][];
            apdq = new double[maxmx + 2 * mbc][];
            trans_wave = new int[mx + 2 * mbc];
            wave_correction = new double[mx + 2 * mbc];
            speeds = new double[mwaves][];
            this.ql = ql;
            this.qr = qr;
            this.auxl = auxl;
            this.auxr = auxr;
            this.mbc = mbc;
            this.mx = mx;
            this.meqn = meqn;
            this.mwaves = mwaves;
            upperbound1 = mx + 2 * mbc;
            for (int i = 0; i < amdq.Length; i++) amdq[i] = new double[meqn];
            for (int i = 0; i < apdq.Length; i++) apdq[i] = new double[meqn];
            for (int i = 0; i < s.Length; i++) s[i] = new double[mwaves];
            for (int i = 0; i < eig_vec.Length; i++) eig_vec[i] = new double[mwaves];
            for (int i = 0; i < speeds.Length; i++) speeds[i] = new double[mbc];
            for (int i = 0; i < fwave.Length; i++) fwave[i] = new double[meqn][];
            for (int i = 0; i < fwave.Length; i++) for (int j = 0; j < fwave[i].Length; j++) fwave[i][j] = new double[mwaves];
        }

        public int[] trans_wave;
        public double[][][] fwave;
        public double[][] s;
        public double[] wave_correction;
        public double[][] amdq;
        public double[][] apdq;
        public double[][] ql;
        public double[][] qr;
        public double[][] auxl;
        public double[][] auxr;
        public double[][] eig_vec;
        public double[][] speeds;

        private int info, rare, layer_index, mbc, mx, meqn, mwaves, upperbound1;
        public int[] ipiv = new int[4];
        public double[] alpha = new double[4];
        public double[] beta = new double[4];
        public double[] delta = new double[4];
        public double[][] A = new double[4][];
        public double[] flux_l = new double[4];
        public double[] flux_r = new double[4];
        public double[] h_ave = new double[2];
        public double[] inundation_height = new double[2];
        public double[] h_l = new double[2];
        public double[] u_l = new double[2];
        public double[] hu_l = new double[2];
        public double[] h_r = new double[2];
        public double[] u_r = new double[2];
        public double[] hu_r = new double[2];
        public double[] h_hat_l = new double[2];
        public double[] h_hat_r = new double[2];
        public double[] momentum_transfer = new double[2];
        private bool[] dry_state_l = new bool[2];
        private bool[] dry_state_r = new bool[2];
        public bool inundation;
        public double b_l, tau, b_r, w_l, w_r, flux_transfer_r, flux_transfer_l, wind_speed;

        Func<double[], double> Product = (arr) => {
            double holder = 1.0;
            foreach (double item in arr) holder *= item;
            return holder;
        };

        public void Run()
        {
            for (int i = 1; i < mx + 2 * mbc; i++)
            {
                inundation = false;
                SetAll(dry_state_l, false);
                SetAll(dry_state_r, false);

                for (int j = 0; j < 2; j++)
                {

                    layer_index = 2 * j;
                    h_l[j] = qr[i - 1][layer_index] / rho[j];
                    h_r[j] = ql[i][layer_index] / rho[j];
                    hu_l[j] = qr[i - 1][layer_index + 1] / rho[j];
                    hu_r[j] = ql[i][layer_index + 1] / rho[j];

                    h_hat_l[j] = auxr[i - 1][j + 2];
                    h_hat_r[j] = auxl[i][j + 2];

                    // Check for dry states in this layer

                    if (h_l[j] < dry_tolerance)
                    {
                        dry_state_l[j] = true;
                        h_l[j] = 0.0;
                        u_l[j] = 0.0;
                    }
                    else
                        u_l[j] = qr[i - 1][layer_index + 2] / qr[i - 1][layer_index + 1];

                    if (h_r[j] < dry_tolerance)
                    {
                        dry_state_r[j] = true;
                        h_r[j] = 0.0;
                        u_r[j] = 0.0;
                    }
                    else
                        u_r[j] = ql[i][layer_index + 2] / ql[i][layer_index + 1]; ;
                }

                h_ave = h_l.Select((x, i) => 0.5 * (x + h_r[i])).ToArray();
                b_l = auxr[i - 1][0];
                b_r = auxl[i][0];
                w_l = auxr[i - 1][1];
                w_r = auxl[i][1];

                // Solve Single layer problem seperately

                if (dry_state_r[1] && dry_state_l[1])
                {
                    var ret = Rp1.SingleLayerEigen(h_l, h_r, u_l, u_r);
                    double[] lambda = ret.Item1;
                    eig_vec = ret.Item2;
                    trans_wave[i] = 0; wave_correction[i] = 0.0;
                    s[i] = lambda;

                    delta[1] = rho[1] * (hu_r[1] - hu_l[1]);
                    flux_r[2] = rho[1] * (Math.Pow(h_r[1] * u_r[1], 2) + 0.50 * Math.Pow(g * h_r[1], 2));
                    flux_l[2] = rho[1] * (Math.Pow(h_l[1] * u_l[1], 2) + 0.50 * Math.Pow(g * h_l[1], 2));
                    delta[2] = flux_r[2] - flux_l[2] + g * rho[1] * h_ave[1] * (b_r - b_l);

                    beta[1] = (delta[1] * s[i][4] - delta[2]) / (s[i][4] - s[i][1]);
                    beta[2] = 0.0;
                    beta[3] = 0.0;
                    beta[4] = (delta[2] - s[i][1] * beta[1]) / s[i][4];

                    // Calculate waves
                    for (int k = 0; k < fwave[i].Length; k++)
                        for (int mw = 0; mw < mwaves; mw++)
                            fwave[i][k][mw] = eig_vec[k][mw] * beta[mw];

                    break;
                }


                // Calculate eigen - space values
                //Inundation cases
                if (dry_state_r[2] && (!dry_state_l[2]) && (h_l[2] + b_l > b_r))
                {

                    Console.WriteLine("Right inundation problem");
                    inundation = true;

                    if (inundation_method == 0)
                        throw new Exception("Inundation not allowed.");

                    else if (inundation_method == 1)
                    {
                        // Linear eigensystem

                        inundation_height = new[] { h_r[0], 0.0 };

                        var ret = Linear_Eigen(h_l, inundation_height, u_l, u_r, ref trans_wave[i], ref wave_correction[i]);
                        double[] lambda = ret.Item1;
                        eig_vec = ret.Item2;
                        s[i] = lambda;

                        // Corrective wave

                        s[i][2] = u_l[1] + 2.0 * Math.Sqrt(g * (1.0 - r) * h_l[1]);
                        s[i][3] = u_r[0] + Math.Sqrt(g * h_r[0]);
                        alpha[2] = r * g * h_l[2] / (Math.Pow(s[i][2] - u_l[1], 2) - g * h_l[1]);
                        alpha[3] = 0.0;
                        eig_vec[0][2] = 1.0;
                        eig_vec[0][3] = 1.0;
                        eig_vec[1][2] = s[i][2];
                        eig_vec[1][3] = s[i][3];
                        eig_vec[2][2] = alpha[2];
                        eig_vec[2][3] = alpha[3];
                        eig_vec[3][2] = s[i][2] * alpha[2] + s[i][2] * alpha[3];
                        eig_vec[3][3] = s[i][3] * alpha[2] + s[i][3] * alpha[3];
                    }

                    else if (inundation_method == 2)
                    {
                        inundation_height = new[] { h_r[0], dry_tolerance };
                        h_r[1] = dry_tolerance;
                        var ret = Linear_Eigen(h_l, inundation_height, u_l, u_r, ref trans_wave[i], ref wave_correction[i]);
                        double[] lambda = ret.Item1;
                        eig_vec = ret.Item2;
                        s[i] = lambda;
                    }
                    else if (inundation_method == 3)
                    {
                        inundation_height = new[] { h_r[0], dry_tolerance };
                        var ret = Velocity_Eigen(h_l, inundation_height, u_l, u_r, ref trans_wave[i], ref wave_correction[i]);
                        double[] lambda = ret.Item1;
                        eig_vec = ret.Item2;
                        s[i] = lambda;

                        //Correction for fast wave

                        s[i][3] = u_r[1] + Math.Sqrt(g * h_r[0]);
                        eig_vec[0][3] = 1.0;
                        eig_vec[1][3] = s[i][3];
                        eig_vec[2][3] = 0.0;
                        eig_vec[3][3] = 0.0;
                    }
                    else if (inundation_method == 4)
                    {
                        inundation_height = new[] { h_r[0], dry_tolerance };

                        var ret = Lapack_Eigen(h_l, inundation_height, u_l, u_r, ref trans_wave[i], ref wave_correction[i]);
                        double[] lambda = ret.Item1;
                        eig_vec = ret.Item2;
                        s[i] = lambda;
                        // Correction for the fast waves

                        s[i][1] = u_r[0] + Math.Sqrt(g * h_r[0]);

                        eig_vec[0][1] = 1.0;
                        eig_vec[1][1] = s[i][1];
                        eig_vec[2][1] = 0.0;
                        eig_vec[3][1] = 0.0;
                    }

                    else if (inundation_method == 5)
                    {
                        inundation_height = new[] { h_r[0], 0.0 };
                        var ret = Lapack_Eigen(h_l, inundation_height, u_l, u_r, ref trans_wave[i], ref wave_correction[i]);
                        double[] lambda = ret.Item1;
                        eig_vec = ret.Item2;
                        s[i] = lambda;
                    }

                    else if (dry_state_l[1] && (!dry_state_r[2]) && (h_r[2] + b_r > b_l))
                    {

                        Console.WriteLine("Left inundation problem");
                        inundation = true;

                        // Inundation problem eigen
                        if (inundation_method == 0) throw new Exception("Inundation not allowed.");

                        else if (inundation_method == 1)
                        {
                            // Linear eigensystem

                            inundation_height = new[] { h_l[0], dry_tolerance };

                            var ret = Linear_Eigen(inundation_height, h_r, u_l, u_r, ref trans_wave[i], ref wave_correction[i]);
                            double[] lambda = ret.Item1;
                            eig_vec = ret.Item2;
                            s[i] = lambda;

                            // Corrections to internal wave
                            s[i][1] = u_r[1] - 2.0 * Math.Sqrt(g * (1.0 - r) * h_r[1]);
                            alpha[1] = r * g * h_r[2] / Math.Pow(s[i][1] - u_r[1], 2) - (g * h_r[1]);
                            eig_vec[0][1] = 1.0;
                            eig_vec[1][1] = s[i][1];
                            eig_vec[2][1] = alpha[1];
                            eig_vec[3][1] = alpha[1] * s[i][1];

                            // Correction for the fast waves
                            s[i][0] = u_l[0] - Math.Sqrt(g * h_l[1]);

                            eig_vec[0][0] = 1.0;
                            eig_vec[1][0] = s[i][0];
                            eig_vec[2][0] = 0.0;
                            eig_vec[3][0] = 0.0;
                        }
                        else if (inundation_method == 2)
                        {
                            // Use linearized eigensystem
                            inundation_height = new[] { h_l[0], dry_tolerance };

                            var ret = Linear_Eigen(inundation_height, h_r, u_l, u_r, ref trans_wave[i], ref wave_correction[i]);
                            double[] lambda = ret.Item1;
                            eig_vec = ret.Item2;
                            s[i] = lambda;

                            // Correction for the fast waves
                            s[i][0] = u_l[0] - Math.Sqrt(g * h_l[1]);
                            eig_vec[0][0] = 1.0;
                            eig_vec[1][0] = s[i][0];
                            eig_vec[2][0] = 0.0;
                            eig_vec[3][0] = 0.0;
                        }

                        else if (inundation_method == 3)
                        {
                            //Use velocity difference expansion eigensystems
                            inundation_height = new[] { h_l[0], dry_tolerance };

                            var ret = Velocity_Eigen(h_l, inundation_height, u_l, u_r, ref trans_wave[i], ref wave_correction[i]);
                            double[] lambda = ret.Item1;
                            eig_vec = ret.Item2;
                            s[i] = lambda;

                            // Correction for the fast waves
                            s[i][0] = u_r[0] - Math.Sqrt(g * h_r[0]);
                            eig_vec[0][0] = 1.0;
                            eig_vec[1][0] = s[i][0];
                            eig_vec[2][0] = 0.0;
                            eig_vec[3][0] = 0.0;
                        }

                        else if (inundation_method == 4)
                        {
                            // LAPACK solver with corrective wave and small wet layer
                            inundation_height = new[] { h_l[0], dry_tolerance };

                            var ret = Lapack_Eigen(h_l, inundation_height, u_l, u_r, ref trans_wave[i], ref wave_correction[i]);
                            double[] lambda = ret.Item1;
                            eig_vec = ret.Item2;
                            s[i] = lambda;

                            // Correction for the fast waves
                            s[i][0] = u_l[0] - Math.Sqrt(g * h_l[0]);
                            eig_vec[0][0] = 1.0;
                            eig_vec[1][0] = s[i][0];
                            eig_vec[2][0] = 0.0;
                            eig_vec[3][0] = 0.0;
                        }
                        else if (inundation_method == 5)
                        {
                            // Use the LAPACK solver with no correction
                            inundation_height = new[] { h_l[0], dry_tolerance };

                            var ret = Lapack_Eigen(h_l, inundation_height, u_l, u_r, ref trans_wave[i], ref wave_correction[i]);
                            double[] lambda = ret.Item1;
                            eig_vec = ret.Item2;
                            s[i] = lambda;
                        }

                    }
                    else
                    {
                        // Wall or wet case
                        // Wall dry state or completely wet case
                        if (eigen_method == 1)
                        {
                            var ret = Linear_Eigen(h_hat_l, h_hat_r, u_l, u_r, ref trans_wave[i], ref wave_correction[i]);
                            double[] lambda = ret.Item1;
                            eig_vec = ret.Item2;
                            s[i] = lambda;
                        }
                        else if (eigen_method == 2)
                        {
                            var ret = Linear_Eigen(h_l, h_r, u_l, u_r, ref trans_wave[i], ref wave_correction[i]);
                            double[] lambda = ret.Item1;
                            eig_vec = ret.Item2;
                            s[i] = lambda;
                        }
                        else if (eigen_method == 3)
                        {
                            var ret = Velocity_Eigen(h_l, h_r, u_l, u_r, ref trans_wave[i], ref wave_correction[i]);
                            double[] lambda = ret.Item1;
                            eig_vec = ret.Item2;
                            s[i] = lambda;
                        }
                        else if (eigen_method == 4)
                        {
                            if (dry_state_r[1] && (!dry_state_l[1]) || dry_state_l[1] && (!dry_state_r[1]))
                            {
                                var ret = Linear_Eigen(h_l, h_r, u_l, u_r, ref trans_wave[i], ref wave_correction[i]);
                                double[] lambda = ret.Item1;
                                s[i] = lambda;
                            }
                            else
                            {
                                var ret = Lapack_Eigen(h_l, h_r, u_l, u_r, ref trans_wave[i], ref wave_correction[i]);
                                double[] lambda = ret.Item1;
                                s[i] = lambda;
                            }
                        }
                        else
                            throw new Exception("Invalid eigensystem method requested, method = (1,4).");

                    }
                }

                //end of eigenspace calculation

                //Calculate flux vector to be projected onto e - space
                // Calculate jumps in fluxes
                if (dry_state_r[1] && (!dry_state_l[1]) && (!inundation))
                {
                    // Wall boundary conditions
                    h_r[1] = h_l[1];
                    hu_r[1] = -hu_l[1];
                    u_r[1] = -u_l[1];

                    // Top layer eta(2) = b_r - h_l(2) - b_l
                    momentum_transfer[0] = g * rho[0] * h_ave[0] * (b_r - h_l[1] - b_l);


                    momentum_transfer[1] = 0.0;
                    flux_transfer_r = 0;
                    flux_transfer_l = 0;
                    // Left state dry, right wet
                }
                else if (dry_state_l[1] && (!dry_state_r[1]) && (!inundation))
                {
                    // Wall boundary conditions
                    h_l[1] = h_r[1];
                    hu_l[1] = -hu_r[1];
                    u_l[1] = -u_r[1];

                    //Top layer eta(2) = h_r(2) + b_r - b_l
                    momentum_transfer[0] = g * rho[0] * h_ave[0] * (b_r + h_r[1] - b_l);
                    momentum_transfer[1] = 0.0;
                    flux_transfer_r = 0.0;
                    flux_transfer_l = 0.0;
                    // Fully wet bottom layer or inundation
                }
                else
                {
                    momentum_transfer[0] = g * rho[0] * h_ave[0] * (h_r[1] - h_l[1] + b_r - b_l);
                    momentum_transfer[1] = -g * rho[0] * h_ave[0] * (h_r[1] - h_l[1]) + g * rho[1] * h_ave[1] * (b_r - b_l);
                    // Bottom layer momentum transfer flux
                    flux_transfer_r = g * rho[0] * Product(h_r);
                    flux_transfer_l = g * rho[0] * Product(h_l);
                }

                // Flux jumps
                for (int j = 0; j < 1; j++) 
                {
                    layer_index = 2 * j;
                    flux_r[layer_index] = rho[j] * hu_r[j];
                    flux_r[layer_index + 1] = rho[j] * (h_r[j] * Math.Pow(u_r[j], 2) + 0.50 * g * Math.Pow(h_r[j], 2));

                    flux_l[layer_index] = rho[j] * hu_l[j];
                    flux_l[layer_index + 1] = rho[j] * (h_l[j] * Math.Pow(u_l[j],2 )+ 0.50 * g * Math.Pow(h_l[j], 2));
                }

                delta = flux_r.Select((x, i) => x - flux_l[i]).ToArray();

                // Bottom layer additional flux
                delta[3] += flux_transfer_r - flux_transfer_l;

                // Momentum source term from layers
                delta[1] = delta[1] + momentum_transfer[0];
                delta[3] = delta[1] + momentum_transfer[1];

                // Wind forcing
                wind_speed = 0.50 * (w_l + w_r);
                tau = Wind_Drag(wind_speed) * rho_air * wind_speed;
                delta[1] = delta[1] - tau * wind_speed;

                //Solve system, solution is stored in delta
                beta = Eval_Lapack_Solve(ref eig_vec, delta);

                // Calculate waves
                for (int k = 0; k < fwave[i].Length; k++)
                    for (int mw = 0; mw < mwaves; mw++)
                        fwave[i][k][mw] = eig_vec[k][mw] * beta[mw];
            }

            if (!entropy_fix)
            {
                for (int i = 1; i < mx + 2 * mbc; i++)
                {
                    for (int k = 0; k < meqn; k++)
                        for (int mw = 0; mw < mwaves; mw++)
                            if (s[i][mw] > 0.0) apdq[i][k] += fwave[i][k][mw];
                            else
                                amdq[i][k] += fwave[i][k][mw];
                }
            }
            else
            {
                for (int i = 1; i < mx + 2 * mbc; i++)
                {

                    if (trans_wave[i] != 0)
                    {
                        for (int k = 0; k < meqn; k++)
                            for (int mw = 0; mw < trans_wave[i] - 1; mw++)
                                amdq[i][k] += fwave[i][k][mw];
                        for (int k = 0; k < meqn; k++)
                        {
                            amdq[i][k] += wave_correction[i] * fwave[i][k][trans_wave[i]];
                            apdq[i][k] += (1.0 - wave_correction[i]) * fwave[i][k][trans_wave[i]];
                        }
                    }
                    else
                    {
                        for (int k = 0; k < meqn; k++)
                            for (int mw = 0; mw < mwaves; mw++)
                                if (s[i][mw] > 0.0) apdq[i][k] += fwave[i][k][mw];
                                else
                                    amdq[i][k] += fwave[i][k][mw];
                    }
                }
            }

        }
        private static void SetAll<T>(T[] refArr, T value)
        {
            for (int i = 0; i < refArr.Length; i++) refArr[i] = value;
        }

        internal (double[], double[][]) Lapack_Eigen(double[] h_l, double[] h_r, double[] u_l, double[] u_r, ref int transonic_wave, ref double wave_correction, int mwaves = 4)
        {

            double[] s = new double[mwaves];
            double[] s_l = new double[mwaves];
            double[] s_r = new double[mwaves];
            double[] h_ave = new double[2];
            double[] u_ave = new double[2];
            double[][] eig_vec = new double[mwaves][];

            // Solve eigenvalue problem
            h_ave = h_l.Select((x, i) => 0.50 * (x + h_r[i])).ToArray();
            u_ave = u_l.Select((x, i) => 0.50 * (x + u_r[i])).ToArray();

            Eval_Lapack_Eigen(h_ave, u_ave, ref s, ref eig_vec);

            // Determine wave speeds
            transonic_wave = 0;
            wave_correction = 0.0;

            if (entropy_fix)
            {
                // Check to see if we may be at a transonic rarefaction
                Eval_Lapack_Eigen(h_l, u_l, ref s_l, ref eig_vec);
                Eval_Lapack_Eigen(h_r, u_r, ref s_r, ref eig_vec);


                //Check each wave for a transonic problem
                for (int j = 0; j < 4; j++)
                {
                    if (s_l[j] < 0.0 && 0.0 < s_r[j])
                    {
                        Console.WriteLine($"Transonic wave detected in wave family {j}.");
                        transonic_wave = j;
                        wave_correction = (s_l[j] - s[j]) / (s_l[j] - s_r[j]);
                    }
                }
            }

            return (s, eig_vec);
        }

        internal (double[], double[][]) Linear_Eigen(double[] h_l, double[] h_r, double[] u_l, double[] u_r, ref int transonic_wave, ref double wave_correction)
        {
            double gamma_l, gamma_r;
            gamma_l = h_l[1] / h_l[0];
            gamma_r = h_r[1] / h_r[0];
            double[] s = new double[4];
            // Left state alphas
            alpha[0] = 0.50 * (gamma_l - 1.0 + Math.Sqrt(Math.Pow(gamma_l - 1.0, 2) + 4.0 * r * gamma_l));
            alpha[1] = 0.50 * (gamma_l - 1.0 - Math.Sqrt(Math.Pow(gamma_l - 1.0, 2) + 4.0 * r * gamma_l));
            // Right state alphas
            alpha[2] = 0.50 * (gamma_r - 1.0 - Math.Sqrt(Math.Pow(gamma_r - 1.0, 2) + 4.0 * r * gamma_r));
            alpha[3] = 0.50 * (gamma_r - 1.0 + Math.Sqrt(Math.Pow(gamma_r - 1.0, 2) + 4.0 * r * gamma_r));

            // Left state speeds
            speeds[0][0] = u_l[0] - Math.Sqrt(g * h_l[0] * (1 + alpha[0]));
            speeds[1][0] = u_l[1] - Math.Sqrt(g * h_l[0] * (1 + alpha[1]));
            speeds[2][0] = u_l[1] + Math.Sqrt(g * h_l[0] * (1 + alpha[1]));
            speeds[3][0] = u_l[0] + Math.Sqrt(g * h_l[0] * (1 + alpha[0]));

            // Right state speeds
            speeds[0][1] = u_l[0] - Math.Sqrt(g * h_r[0] * (1 + alpha[3]));
            speeds[1][1] = u_l[1] - Math.Sqrt(g * h_r[0] * (1 + alpha[2]));
            speeds[2][1] = u_l[1] + Math.Sqrt(g * h_r[0] * (1 + alpha[2]));
            speeds[3][1] = u_l[0] + Math.Sqrt(g * h_r[0] * (1 + alpha[3]));

            // Determine wave speeds
            transonic_wave = 0;
            wave_correction = 0.0;

            if (entropy_fix) 
            {
                transonic_wave = 0;
                wave_correction = 0.0;
                for (int mw = 0; mw < 4; mw++)
                {
                    // Both speeds right going
                    if (speeds[mw][0] > 0.0) s[mw] = speeds[mw][1];
                    // Transonic rarefaction
                    //s_l < 0.0 < s_r
                    else if (speeds[mw][1] > 0.0)
                    {
                        // Assign base speed for rarefaction, should approximate true speed

                        s[mw] = 0.50 * (speeds[mw][0] + speeds[mw][1]); ;
                        transonic_wave = mw;
                        wave_correction = Math.Abs(speeds[mw][0]) / (Math.Abs(speeds[mw][0]) + Math.Abs(speeds[mw][1]));
                    }
                    else
                        s[mw] = speeds[mw][0];
                }
            }
            else
                s = new[] { speeds[0][0], speeds[1][0], speeds[2][1], speeds[2][1] };

            SetAll(eig_vec[0], 1.0);
            eig_vec[1] = s;
            eig_vec[2] = alpha;
            eig_vec[3] = s.Select((x, i) => x * alpha[i]).ToArray();
            return (s, eig_vec);
        }

        internal (double[], double[][]) Velocity_Eigen(double[] h_l, double[] h_r, double[] u_l, double[] u_r, ref int transonic_wave, ref double wave_correction)
        {
            double[] s = new double[mwaves];
            double[][] eig_vec = new double[mwaves][];
            double r = Setprob.r;
            double g = Setprob.g;
            double one_minus_r = Setprob.one_minus_r;
            for (int i = 0; i < eig_vec.Length; i++) eig_vec[i] = new double[mwaves];

            double total_depth_l, total_depth_r, mult_depth_l = 1, mult_depth_r = 1;

            total_depth_l = h_l.Sum();
            total_depth_r = h_r.Sum();
            foreach (double item in h_l) mult_depth_l *= item;
            foreach (double item in h_r) mult_depth_r *= item;

            //Determine wave speeds
            transonic_wave = 0;
            wave_correction = 0.0;

            s[0] = (h_l[0] * u_l[0] + h_l[1] * u_l[1]) / total_depth_l - Math.Sqrt(g * total_depth_l);
            s[1] = (h_l[1] * u_l[0] + h_l[0] * u_l[1]) / total_depth_l - Math.Sqrt(g * one_minus_r * mult_depth_l / total_depth_l)
                    * (1 - Math.Pow(u_l[0] - u_l[1], 2) / (g * one_minus_r * total_depth_l));
            s[2] = (h_r[1] * u_r[0] + h_r[0] * u_r[1]) / total_depth_l - Math.Sqrt(g * one_minus_r * mult_depth_r / total_depth_r)
                    * (1 - Math.Pow(u_r[0] - u_r[1], 2) / (g * one_minus_r * total_depth_r));
            s[3] = (h_r[0] * u_r[0] + h_r[1] * u_r[1]) / total_depth_r - Math.Sqrt(g * total_depth_r);

            alpha[0] = (Math.Pow(s[0] - u_l[0], 2) - g * h_l[0]) / (g * h_l[0]);
            alpha[1] = (Math.Pow(s[1] - u_l[0], 2) - g * h_l[0]) / (g * h_l[0]);
            alpha[2] = (Math.Pow(s[2] - u_l[0], 2) - g * h_l[0]) / (g * h_l[0]);
            alpha[3] = (Math.Pow(s[3] - u_l[0], 2) - g * h_l[0]) / (g * h_l[0]);

            SetAll(eig_vec[0], 1.0);
            eig_vec[1] = s;
            eig_vec[2] = alpha;
            eig_vec[3] = s.Select((x, i) => x * alpha[i]).ToArray();

            return (s, eig_vec);
        }

        internal void Eval_Lapack_Eigen(double[] h, double[] u, ref double[] s, ref double[][] eig_vec)
        {
            double[] imag_s = new double[mwaves];
            double[] h_ave = new double[2];
            double[] u_ave = new double[2];
            double[][] A = new double[mwaves][];
            double r = Setprob.r;
            double g = Setprob.g;
            double[] empty = new double[mwaves * mwaves], work = new double[mwaves * mwaves];
            int info = default;
            int lwork = default;

            for (int i = 0; i < A.Length; i++) A[i] = new double[mwaves];

            A[0] = new double[] { 0.0, 1.0, 0.0, 0.0 };
            A[1] = new double[] { -Math.Pow(u_ave[0], 2) + g * h_ave[0], 2.0 * u_ave[0], g * h_ave[0], 0.0 };
            A[2] = new double[] { 0.0, 0.0, 0.0, 1.0 };
            A[3] = new double[] { g * r * h_ave[1], 0.0, -Math.Pow(u_ave[1], 2) + g * h_ave[1], 2.0 * u_ave[1] };
            double[] A_flattened = A.SelectMany(x => x).ToArray();
            double[] eig_vec_flattened = eig_vec.SelectMany(x => x).ToArray();

            var dgeev = new DGEEV();
            //Call LAPACK eigen solver
            dgeev.Run("N", "V", 4, ref A_flattened, 0, 4, ref s, 0, ref imag_s, 0, ref empty, 0, 1, ref eig_vec_flattened, 0, 4, ref work, 0, lwork, ref info);

            info = (info < 0) ? -info : info;
            for (int i = 0; i < mwaves; i++) eig_vec[i] = eig_vec_flattened.Skip(4 * i).Take(4).ToArray();
        }

        internal double[] Eval_Lapack_Solve(ref double[][] eig_vec, double[] delta)
        {
            double[] eig_vec_flattened = eig_vec.SelectMany(x => x).ToArray();
            double[] beta = delta;
            int[] ipiv = new int[mwaves];
            int info = 0;
            var dgesv = new DGESV();
            dgesv.Run(4, 1, ref eig_vec_flattened, 0, 4, ref ipiv, 0, ref beta, 0, delta.Length, ref info);

            for (int i = 0; i < mwaves; i++) eig_vec[i] = eig_vec_flattened.Skip(4 * i).Take(4).ToArray();
            return beta;
        }


    }
}
