# ClawpackNetCore
C# Port of  Clawpack encapsulating finite volume methods for hyperbolic Partial Differential Equations (PDEs).

```
Setprob.Init(); Should be run to set up the default conditions.
```
The overloaded method should be run to modify the parameters
```
Init(double[] rho, double rho_air, double dry_tolerance, int eigen_method, double xlower, double ndim, double bulk, int inundation_method, bool entropy_fix)
```
The default values are listed in *datafile.cs*
