
namespace GeoclawNetCore._1D
{
    public static class datafile
    {
        public const double gravity = 9.8;
        public const int eigen_method = 4;
        public const int bulk = 1;
        public const int inundation_method = 2;
        public static double[] rho = new[] { 1025.0, 1028.0 };
        public const double rho_air = 1.15;
        public const double ambient_pressure = 1.15;
        public const int coordinate_system = 1;
        public const bool coriolis_forcing = true;
        public const double theta_0 = 45.0;
        public const bool friction_forcing = true;
        public const double dry_tolerance = 1e-3;
        public const double friction_depth = 1.0e6;
        public const double sea_level = 0.0;
        public const bool entropy_fix = false;
        public static double[] manning_coefficient = new double[] { 0.025 };

    }
}
