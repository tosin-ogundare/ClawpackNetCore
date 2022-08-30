using GeoclawNetCore._1D;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ClassicClawPackNetCore._1D
{
    /// <summary>
    ///     Take one time step, updating q.
    /// </summary>
    public class Step1
    {

        public Step1(int num_waves, int num_ghost, int num_eqn, int num_aux, int mx, double dt, double dx, bool use_fwave)
        {
            this.num_waves = num_waves;
            this.num_ghost = num_ghost;
            this.mx = mx;
            this.dt = dt;
            this.dx = dx;
            this.num_eqn = num_eqn;
            this.use_fwave = use_fwave;
            q = new double[num_eqn][];
            aux = new double[num_aux][];
            f = new double[num_eqn][];
            s = new double[num_waves][];
            wave = new double[num_eqn][][];
            amdq = new double[num_eqn][];
            apdq = new double[num_eqn][];
            dtdx = new double[mx + 2*num_ghost];
            method = new int[7];
            mthlim = new int[num_waves];
            for (int i = 0; i < amdq.Length; i++) amdq[i] = new double[num_eqn];
            for (int i = 0; i < apdq.Length; i++) apdq[i] = new double[num_eqn];
            for (int i = 0; i < s.Length; i++) s[i] = new double[mx + 2 * num_ghost];
            for (int i = 0; i < f.Length; i++) f[i] = new double[mx + 2 * num_ghost];
            for (int i = 0; i < aux.Length; i++) aux[i] = new double[mx + 2 * num_ghost];
            for (int i = 0; i < q.Length; i++) q[i] = new double[mx + 2 * num_ghost];
            for (int i = 0; i < wave.Length; i++) wave[i] = new double[num_waves][];
            for (int i = 0; i < wave.Length; i++) for (int j = 0; j < wave[i].Length; j++) wave[i][j] = new double[mx + 2 * num_ghost];
        }

        int num_waves, num_ghost, mx, num_eqn;
        double dt, dx;
        public double[][] q;
        public double[][] aux;
        public double[][] f;
        public double[][] s;
        public double[][][] wave;
        public double[][] amdq;
        public double[][] apdq;
        public double[] dtdx;
        public int[] method;
        public int[] mthlim;
        public bool use_fwave, limit;

        ///<remarks>
        ///     method(1) = 1   ==>  Godunov method
        ///     method(1) = 2   ==>  Slope limiter method
        ///     mthlim(p)  controls what limiter is used in the pth family
        ///     amdq, apdq, wave, s, and f are used locally:
        ///     amdq(num_eqn, 1-num_ghost:mx+num_ghost) = left-going flux-differences
        ///     apdq(num_eqn, 1-num_ghost:mx+num_ghost) = right-going flux-differences
        ///        e.g. amdq(m,i) = m'th component of A^- \Delta q from i'th Riemann
        ///                         problem (between cells i-1 and i).
        ///     wave(num_eqn, num_waves, 1-num_ghost:mx+num_ghost) = waves from solution of
        ///                                           Riemann problems,
        ///            wave(m,mw,i) = mth component of jump in q across
        ///                           wave in family mw in Riemann problem between
        ///                           states i-1 and i.
        ///     s(num_waves, 1-num_ghost:mx+num_ghost) = wave speeds,
        ///            s(m,iw) = speed of wave in family mw in Riemann problem between
        ///                      states i-1 and i.
        ///     f(num_eqn, 1-num_ghost:mx+num_ghost) = correction fluxes for second order method
        ///            f(m,i) = mth component of flux at left edge of ith cell
        ///</remarks>
        public void Run()
        {
            // check if any limiters are used:
            limit = false;
            for (int mw = 0; mw < num_waves; mw++)
                if (mthlim[mw] > 0) limit = true;


            int index_capa = method[5];
            for (int i = 0; i < (mx + 2 * num_ghost); i++) 
            {
                if (index_capa > 0)
                {
                    if (aux[index_capa][i] <= 0.0)
                        throw new Exception("Error -- capa must be positive");

                    dtdx[i] = dt / (dx * aux[index_capa][i]);

                }
                else
                    dtdx[i] = dt / dx;
            }

            // solve Riemann problem at each interface
             var rpi = new Rp1(ql: q, qr: q, auxl: aux, auxr: aux, mx: mx, maxmx: mx );
            s = rpi.s;
            amdq = rpi.amdq;
            apdq = rpi.apdq;
            wave = rpi.fwave;
            // Modify q for Godunov update:
            // Note this may not correspond to a conservative flux-differencing
            // for equations not in conservation form.  It is conservative if
            // amdq + apdq = f(q(i)) - f(q(i-1)).

            for(int i = 1; i < mx + num_ghost + 1; i++)
            {
                // q(:, i - 1) is still in cache from last cycle of i loop, so
                // update it first
                for (int m = 0; m < num_eqn; m++) 
                { 
                    q[m][i - 1] = q[m][i - 1] - dtdx[i - 1] * amdq[m][i];
                    q[m][i] = q[m][i] - dtdx[i] * apdq[m][i];
                }
                 

            }

            // compute maximum wave speed:
            double cfl = 0.0;
            for (int i = 1; i <= mx + 1; i++)
                for (int mw = 0; mw < num_waves; mw++)
                    // if s>0 use dtdx(i) to compute CFL,
                    // if s<0 use dtdx(i-1) to compute CFL:
                    cfl = Math.Max(cfl, Math.Max(dtdx[i] * s[mw][i], -dtdx[i - 1] * s[mw][i]));

            if (method[1] == 1) return;

            // compute correction fluxes for second order q_{xx} terms:/
            // apply limiter to waves:
            var limiter = new Limiter(maxm: mx, num_eqn, num_waves, num_ghost, mx);
            if (limit) 
            {
                limiter.Run();
                wave = limiter.wave;
                s = limiter.s;
                mthlim = limiter.mthlim;
            }


            if (use_fwave == false)
            {
                for (int i = 1; i <= mx + 1; i++)
                {
                    for (int m = 0; m < num_eqn; m++)
                        f[m][i] = 0.0;

                    double dtdxave = 0.50 * (dtdx[i - 1] + dtdx[i]);
                    for (int mw = 0; mw < num_waves; mw++)
                        for (int m = 0; m < num_eqn; m++)
                            f[m][i] = f[m][i] + 0.50 * Math.Abs(s[mw][i]) * (1.0 - Math.Abs(s[mw][i]) * dtdxave) * wave[m][mw][i];
                } 
            }
            else
                for (int i = 1; i <= mx + 1; i++)
                {
                    for (int m = 0; m < num_eqn; m++)
                        f[m][i] = 0.0;

                    double dtdxave = 0.50 * (dtdx[i - 1] + dtdx[i]);
                    for (int mw = 0; mw < num_waves; mw++)
                        for (int m = 0; m < num_eqn; m++)
                            f[m][i] = f[m][i] + 0.50 * (0 > s[mw][i] ? -1.0 : 1.0) * (1.0 - Math.Abs(s[mw][i]) * dtdxave) * wave[m][mw][i];
                }

            // update q by differencing correction fluxes
            // (Note:  Godunov update has already been performed above)

            for (int i = 1; i <= mx + 1; i++)
                for (int m = 0; m < num_eqn; m++)
                    q[m][i] = q[m][i] - dtdx[i] * (f[m][i + 1] - f[m][i]);
        }
    }
}
