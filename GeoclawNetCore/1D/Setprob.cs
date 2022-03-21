namespace GeoclawNetCore._1D
{
    public static class Setprob
    {
        public static double[] rho = new double[2];
        public static double dry_tolerance;
        public static double PI;
        public static double g;
        public static double r;
        public static int eigen_method;
        public static double one_minus_r;
        public static double rho_air;
        public static double xlower;
        public static double ndim;
        public static double bulk;
        public static int inundation_method;
        public static bool entropy_fix;

        public static void Init(double[] rho, double rho_air, double dry_tolerance, int eigen_method, double xlower, double ndim, double bulk, int inundation_method, bool entropy_fix)
        {
            Setprob.rho = rho;
            Setprob.dry_tolerance = dry_tolerance;
            PI = 3.141592654;
            g = 9.8;
            r = rho[0] / rho[1];
            one_minus_r = 1.0 - r;
            Setprob.eigen_method = eigen_method;
            Setprob.rho_air = rho_air;
            Setprob.xlower = xlower;
            Setprob.ndim = ndim;
            Setprob.bulk = bulk;
            Setprob.inundation_method = inundation_method;
            Setprob.entropy_fix = entropy_fix;
    }

        public static double Wind_Drag(double wind_speed)
        {
            double wind_drag;

            if (wind_speed <= 11.0) wind_drag = 1.20;
            else if ((wind_speed > 11.0) && (wind_speed <= 25.0)) wind_drag = 0.490 + 0.0650 * wind_speed;
            else
                wind_drag = 0.49 + 0.0650 * 25.0;
            wind_drag = wind_drag * 10.0 - 3;

            return wind_drag;
        }
    }
}
