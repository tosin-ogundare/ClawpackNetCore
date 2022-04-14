using JetBrains.Annotations;

namespace ClassicClawPackNetCore._1D
{
    public class Philim
    {
        [Pure]
        public static double Run(double a, double b, int meth)
        {
            double r = b / a;
            double ret;

            switch (meth)
            {
                default: ret = r; break;
                case 1:
                    ret = Math.Max(0.0, Math.Min(1.0, r)); break;
                case 2: 
                    ret = Math.Max(0.0, Math.Max(Math.Min(1.0, 2.0 * r), Math.Min(2.0, r))); break;
                case 3: 
                    ret = (r + Math.Abs(r)) / (1.0 + Math.Abs(r)); break;
                case 4:
                    double c = (1.0 + r) / 2.0;
                    ret = Math.Max(0.0, Math.Min(Math.Min(c, 2.0), 2.0 * r)); break;
                case 5:
                    ret = r; break;

            };

            return ret;
        }
    }
}
