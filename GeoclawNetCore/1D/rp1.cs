using DotNumerics.LinearAlgebra.CSLapack;
using static GeoclawNetCore._1D.Setprob;

namespace GeoclawNetCore._1D
{
    /// <summary>
    ///     Riemann solver for linearized multilayer shallow water equations
    /// </summary>
    public class Rp1
    {

        /// <summary>
        /// 
        /// </summary>
        /// <param name="maxmx"></param>
        /// <param name="meqn"></param>
        /// <param name="mwaves"></param>
        /// <param name="mbc"></param>
        /// <param name="mx"></param>
        /// <param name="ql"></param>
        /// <param name="qr"></param>
        /// <param name="auxl"></param>
        /// <param name="auxr"></param>
        public Rp1(double[][] ql, double[][] qr, double[][] auxl, double[][] auxr, int mx, int maxmx, int meqn = 4, int mwaves = 4, int mbc = 2)
        {
            fwave = new double[maxmx + 2 * mbc][][];
            s = new double[maxmx + 2 * mbc][];
            eig_vec = new double[meqn][];
            amdq = new double[maxmx + 2 * mbc][];
            apdq = new double[maxmx + 2 * mbc][];
            this.ql = ql;
            this.qr = qr;
            this.auxl = auxl;
            this.auxr = auxr;
            this.mbc = mbc;
            this.meqn = meqn;
            this.mwaves = mwaves;
            upperbound1 = mx + 2 * mbc;
            for (int i = 0; i < amdq.Length; i++) amdq[i] = new double[meqn];
            for (int i = 0; i < apdq.Length; i++) apdq[i] = new double[meqn];
            for (int i = 0; i < s.Length; i++) s[i] = new double[mwaves];
            for (int i = 0; i < eig_vec.Length; i++) eig_vec[i] = new double[mwaves];
            for (int i = 0; i < fwave.Length; i++) fwave[i] = new double[meqn][];
            for (int i = 0; i < fwave.Length; i++) for (int j = 0; j < fwave[i].Length; i++) fwave[i][j] = new double[mwaves];
        }

        public double[][][] fwave;
        public double[][] s;
        public double[][] amdq;
        public double[][] apdq;
        public double[][] ql;
        public double[][] qr;
        public double[][] auxl;
        public double[][] auxr;
        public double[][] eig_vec;

        private int info, rare, layer_index, mbc, meqn, mwaves, upperbound1;
        public int[] ipiv = new int[4];
        public double[] beta = new double[4];
        public double[] delta = new double[4];
        public double[][] A = new double[4][];
        public double[] flux_l = new double[4];
        public double[] flux_r = new double[4];
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

        public double b_l, tau, b_r, w_l, w_r, flux_transfer_r, flux_transfer_l, wind_speed;
        public double Run()
        {
            for (int i = 1; i < upperbound1; i++)
            {
                SetAll(dry_state_l, false);
                SetAll(dry_state_r, false);

                for (int j = 0; j < 2; j++)
                {
                    layer_index = 2 * (j - 1);
                    h_l[j] = qr[i - 1][layer_index + 1] / rho[j];
                    h_r[j] = ql[i][layer_index + 1] / rho[j];
                    hu_l[j] = qr[i - 1][layer_index + 2] / rho[j];
                    hu_r[j] = ql[i][layer_index + 2] / rho[j];


                    h_hat_l[j] = auxr[i - 1][j + 2];
                    h_hat_r[j] = auxl[i][j + 2];

                    // Check for dry states in this layer
                    if (h_l[j] < dry_tolerance)
                    {
                        dry_state_l[j] = true;
                        h_l[j] = 0;
                        u_l[j] = 0;
                    }
                    else u_l[j] = qr[i - 1][layer_index + 2] / qr[i - 1][layer_index + 1];


                    if (h_r[j] < dry_tolerance)
                    {
                        dry_state_r[j] = true;
                        h_r[j] = 0;
                        u_r[j] = 0;
                    }
                    else u_r[j] = ql[i][layer_index + 2] / ql[i][layer_index + 1];
                }

                double[] h_ave = h_l.Select((x, i) => (x + h_r[i]) * 0.5).ToArray();

                b_l = auxr[i - 1][i - 1 + mbc];
                b_r = auxl[i][i - 1 + mbc];
                w_l = auxr[i - 1][i + mbc];
                w_r = auxl[i][i + mbc];

                //Inundation test
                if (dry_state_r[2] && (!dry_state_l[2]) && h_l[2] + b_l > b_r)
                {
                    rare = 1;
                    throw new Exception("Right inundation problem not handled");
                }
                else if (dry_state_l[2] && (!dry_state_r[2]) && h_r[2] + b_r > b_l)
                {
                    rare = 2;
                    throw new Exception("Left inundation problem not handled");
                }
                else
                    rare = 0;

                //Solve Single layer problem seperately
                if (dry_state_r[2] && dry_state_l[2])
                {
                    var ret = SingleLayerEigen(h_l, h_r, u_l, u_r);
                    double[] lambda = ret.Item1;
                    eig_vec = ret.Item2;

                    s[i] = lambda;

                    delta[0] = rho[0] * (hu_r[0] - hu_l[0]);
                    flux_r[1] = rho[0] * Math.Pow(h_r[0] * u_r[0], 2.0) + 0.5 * g * Math.Pow(h_r[0], 2.0);
                    flux_l[1] = rho[0] * Math.Pow(h_r[0] * u_r[0], 2.0) + 0.5 * g * Math.Pow(h_l[1], 2.0);
                    delta[1] = flux_r[2] - flux_l[2] + g * rho[1] * h_ave[1] * (b_r - b_l);


                    beta[0] = (delta[0] * s[i][meqn - 1] - delta[1]) / (s[i][meqn - 1] - s[i][0]);
                    beta[1] = 0;
                    beta[2] = 0;
                    beta[3] = (delta[1] - s[i][0] * beta[0]) / s[i][meqn - 1];

                    // Calculate waves
                    for (int k = 0; k < meqn; k++)
                        for (int mw = 0; mw < mwaves; mw++)
                            fwave[i][k][mw] = eig_vec[k][mw] * beta[mw];

                }

                //Calculate eigen - space values
                int eigen_method = Setprob.eigen_method;

                if (eigen_method == 1)
                {
                    var ret = Linear_Eigen(h_hat_l, h_hat_r, u_l, u_r, rare);
                    double[] lambda = ret.Item1;
                    eig_vec = ret.Item2;
                    s[i] = lambda;
                }
                else if (eigen_method == 2)
                {
                    var ret = Linear_Eigen(h_l, h_r, u_l, u_r, rare);
                    double[] lambda = ret.Item1;
                    eig_vec = ret.Item2;
                    s[i] = lambda;
                }
                else if (eigen_method == 3)
                {
                    var ret = Velocity_Eigen(h_l, h_r, u_l, u_r, rare);
                    double[] lambda = ret.Item1;
                    eig_vec = ret.Item2;
                    s[i] = lambda;
                }
                else if (eigen_method == 4)
                {
                    if (dry_state_r[1] && (!dry_state_l[1]) || dry_state_l[1] && (!dry_state_r[1]))
                    {
                        var ret = Linear_Eigen(h_l, h_r, u_l, u_r, rare);
                        double[] lambda = ret.Item1;
                        eig_vec = ret.Item2;
                        s[i] = lambda;
                    }
                    else
                    {
                        var ret = Exact_Eigen(h_l, h_r, u_l, u_r, rare);
                        double[] lambda = ret.Item1;
                        eig_vec = ret.Item2;
                        s[i] = lambda;
                    }
                }
                else
                    throw new Exception("Invalid eigensystem method requested, method = (1,4).");

                //Calculate flux vector to be projected onto e - space
                // Right state dry, left wet
                if (dry_state_r[1] && (!dry_state_l[1]))
                {
                    // Inundation
                    if (rare == 1) throw new Exception("Inundation of right state not implemented.");
                    else
                    {
                        h_r[1] = h_l[1];
                        hu_r[1] = -hu_l[1];
                        u_r[1] = -u_l[1];
                        flux_transfer_r = 0.0;
                        flux_transfer_l = 0.0;
                        momentum_transfer[0] = g * rho[0] * h_ave[0] * (b_r - h_l[1] - b_l);
                        momentum_transfer[1] = 0.0;
                    }
                }
                // Left state dry, right wet
                else if (dry_state_l[1] && (!dry_state_r[1]))
                {
                    // Inundation
                    if (rare == 2) throw new Exception("Inundation of left state not implemented.");
                    else
                    {
                        h_l[1] = h_r[1];
                        hu_l[1] = -hu_r[1];
                        u_l[1] = -u_r[1];
                        flux_transfer_r = 0.0;
                        flux_transfer_l = 0.0;
                        momentum_transfer[0] = g * rho[0] * h_ave[0] * (b_r - h_r[1] - b_l);
                        momentum_transfer[1] = 0.0;
                    }
                }
                // Fully wet bottom layer
                else
                {
                    momentum_transfer[0] = g * rho[0] * h_ave[0] * (h_r[1] - h_l[1] + b_r - b_l);
                    momentum_transfer[1] = g * rho[0] * h_ave[0] * (h_r[1] - h_l[1]) + g * rho[1] * h_ave[1] * (b_r - b_l);
                    flux_transfer_r = g * rho[0] * h_r[0] * h_r[1];
                    flux_transfer_l = g * rho[0] * h_l[0] * h_l[1];
                }


                for (int j = 0; j < 2; j++) 
                {
                    layer_index = 2 * j;
                    flux_r[layer_index] = rho[j] * hu_r[j];
                    flux_r[layer_index + 1] = rho[j] * Math.Pow(Math.Pow(h_r[j] * u_r[j],2) + 0.50 * g * h_r[j] ,2);
                    flux_l[layer_index] = rho[j] * hu_l[j];
                    flux_l[layer_index + 1] = rho[j] * Math.Pow(Math.Pow(h_l[j] * u_l[j] ,2) + 0.50 * g * h_l[j] ,2);

                }

                // Add extra flux terms
                flux_r[3] = flux_r[3] + flux_transfer_r;
                flux_l[3] = flux_l[3] + flux_transfer_l;

                delta = flux_r.Select((x, i) => x - flux_l[i]).ToArray();

                // Momentum transfer and bathy terms
                delta[1] = delta[1] + momentum_transfer[0];
                delta[3] = delta[3] + momentum_transfer[1];

                // Wind forcing
                wind_speed = 0.50 * (w_l + w_r);
                tau = Wind_Drag(wind_speed) * rho_air * wind_speed;
                delta[1] = delta[1] - tau * wind_speed;


                // Solve system, solution is stored in delta
                A = eig_vec;
                double[] A_flattened = A.SelectMany(x => x).ToArray();
                var dgesv = new DGESV();
                dgesv.Run(4, 1, ref A_flattened, 0, 4, ref ipiv, 0, ref delta, 0, 4, ref info);
                if (!(info == 0)) 
                {
                    Console.WriteLine($"Location [i] = [{i}]");
                    Console.WriteLine($"Dry states, L= {dry_state_l[2]} R={dry_state_r[2]}");
                    Console.WriteLine($"h_l[2] = {h_l[2]} h_r[2] = {h_r[2]}");
                    Console.WriteLine($"Error solving R beta = delta, {info}");
                    Console.WriteLine($"Eigen-speeds:{s[i][0]}, {mwaves}");
                    Console.WriteLine($"Eigen-vectors:");
                    for (int j = 0; j < 4; j++)
                        Console.WriteLine($"[4d16.8] {eig_vec[j][0]}, {meqn}");
                }
                beta = delta;

                // Calculate waves
                 for (int k = 0; k < meqn; k++)
                    for (int mw = 0; mw < mwaves; mw++)
                        fwave[i][k][mw] = eig_vec[k][mw] * beta[mw];

            }

            for(int i=0; i < 5; i++)
            {
                // Calculate amdq and apdq
                for (int k = 0; k < meqn; k++)
                    for (int mw = 0; mw < mwaves; mw++)
                        if (s[i][k] > 0.0) apdq[i][k] += fwave[i][k][mw];
                        else amdq[i][k] += fwave[i][k][mw];
            }
            return 0;
        }

        private static void SetAll<T>(T[] refArr, T value)
        {
            for (int i = 0; i < refArr.Length; i++) refArr[i] = value;
        }

        private static (double[], double[][]) SingleLayerEigen(double[] h_l, double[] h_r, double[] u_l, double[] u_r, int mwaves = 4)
        {
            double[] s = new double[mwaves];
            double[][] eig_vec = new double[mwaves][];
            for (int i = 0; i < eig_vec.Length; i++) eig_vec[i] = new double[mwaves];

            s[0] = u_l[0] - Math.Sqrt(g * h_l[0]);
            s[1] = 0;
            s[2] = 0;
            s[3] = u_r[0] + Math.Sqrt(g * h_r[1]);

            SetAll(eig_vec[0], 1.0);
            eig_vec[1] = s;
            SetAll(eig_vec[2], 0.0);
            SetAll(eig_vec[3], 0.0);

            return (s, eig_vec);
        }

        private static (double[], double[][]) Linear_Eigen(double[] h_l, double[] h_r, double[] u_l, double[] u_r, int rare, int mwaves = 4)
        {

            double[] s = new double[mwaves];
            double[] alpha = new double[mwaves];
            double[][] eig_vec = new double[mwaves][];
            double r = Setprob.r;
            double g = Setprob.g;

            for (int i = 0; i < eig_vec.Length; i++) eig_vec[i] = new double[mwaves];

            
            double gamma_l, gamma_r;

            if (rare == 1)
            {
                gamma_l = h_l[1] / h_l[0];

                alpha[0] = 0.5 * (gamma_l - 1 + Math.Sqrt(Math.Pow(gamma_l - 1, 2) + 4 * r * gamma_l));
                alpha[1] = 0.5 * (gamma_l - 1 - Math.Sqrt(Math.Pow(gamma_l - 1, 2) + 4 * r * gamma_l));

                s[0] = -Math.Sqrt(g * h_l[0] * (1 + alpha[0]));
                s[1] = -Math.Sqrt(g * h_l[0] * (1 + alpha[2]));
                s[2] = u_l[1] + 2.0 * Math.Sqrt(g * h_l[1]);
                s[3] = u_r[0] + Math.Sqrt(g * h_r[0]);


                eig_vec[0][0] = 1.0;
                eig_vec[0][1] = 1.0;
                eig_vec[1][0] = s[0];
                eig_vec[1][1] = s[1];
                eig_vec[2][0] = alpha[0];
                eig_vec[2][1] = alpha[1];
                eig_vec[3][0] = s[0] * alpha[0];
                eig_vec[3][1] = s[1] * alpha[1];

            }
            else if (rare == 2) return default;
            else 
            {
                gamma_l = h_l[1] / h_l[0];
                gamma_r = h_r[1] / h_r[0];


                alpha[0] = 0.50 * (gamma_l - 1.0 + Math.Sqrt(Math.Pow(gamma_l - 1.0, 2) + 4.0 * r * gamma_l));
                alpha[1] = 0.50 * (gamma_l - 1.0 - Math.Sqrt(Math.Pow(gamma_l - 1.0, 2) + 4.0 * r * gamma_l));
                alpha[2] = 0.50 * (gamma_r - 1.0 - Math.Sqrt(Math.Pow(gamma_r - 1.0, 2) + 4.0 * r * gamma_r));
                alpha[3] = 0.50 * (gamma_r - 1.0 + Math.Sqrt(Math.Pow(gamma_r - 1.0, 2) + 4.0 * r * gamma_r));

                s[0] = -Math.Sqrt(g * h_l[0] * (1 + alpha[0]));
                s[1] = -Math.Sqrt(g * h_l[0] * (1 + alpha[1]));
                s[2] = -Math.Sqrt(g * h_r[0] * (1 + alpha[2]));
                s[3] = -Math.Sqrt(g * h_r[0] * (1 + alpha[3]));

                SetAll(eig_vec[0], 1.0);
                eig_vec[1] = s;
                eig_vec[2] = alpha;
                eig_vec[3] = s.Select((x, i) => x * alpha[i]).ToArray();
            }

            return (s, eig_vec);
        }

        private static (double[], double[][]) Exact_Eigen(double[] h_l, double[] h_r, double[] u_l, double[] u_r, int rare, int mwaves = 4)
        {

            double[] s = new double[mwaves];
            double[] real_evalues = new double[mwaves];
            double[] imag_evalues = new double[mwaves];
            double[] h_ave = new double[2];
            double[] u_ave = new double[2];
            double[] alpha = new double[mwaves];
            double[][] eig_vec = new double[mwaves][];
            double[][] A = new double[mwaves][];
            double r = Setprob.r;
            double g = Setprob.g;
            double[] empty = new double[mwaves * mwaves], work = new double[mwaves * mwaves];
            int info = default;
            int lwork = default;

            for (int i = 0; i < eig_vec.Length; i++) eig_vec[i] = new double[mwaves];
            for (int i = 0; i < A.Length; i++) A[i] = new double[mwaves];

            A[0] = new double[] { 0.0, 1.0, 0.0, 0.0 };
            A[1] = new double[] { -Math.Pow(u_ave[0],2) + g * h_ave[0], 2.0 * u_ave[0], g * h_ave[0], 0.0 };
            A[2] = new double[] { 0.0, 0.0, 0.0, 1.0};
            A[3] = new double[] { g * r * h_ave[1], 0.0, -Math.Pow(u_ave[1], 2) + g * h_ave[1], 2.0 * u_ave[1] };
            double[] A_flattened = A.SelectMany(x => x).ToArray();
            double[] eig_vec_flattened = eig_vec.SelectMany(x => x).ToArray();

            var dgeev = new DGEEV();
            //Call LAPACK eigen solver
            dgeev.Run("N", "V", 4, ref A_flattened, 0, 4, ref real_evalues, 0, ref imag_evalues, 0, ref empty, 0, 1, ref eig_vec_flattened, 0, 4, ref work, 0, lwork, ref info);

            info = (info < 0) ? -info : info;
            for (int i = 0; i < mwaves; i++) eig_vec[i] = eig_vec_flattened.Skip(4 * i).Take(4).ToArray();
            return (real_evalues, eig_vec);
        }

        private static (double[], double[][]) Velocity_Eigen(double[] h_l, double[] h_r, double[] u_l, double[] u_r, int rare, int mwaves = 4)
        {

            double[] s = new double[mwaves];
            double[] alpha = new double[mwaves];
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

    }
}
