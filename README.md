Ambit Stochastics
=================

A collection of numerical codes to help perform simulations of ambit
fields.

We are working towards providing a library of useful C-codes.
Currently, the project is very much in its infancy as we experiment
and explore.

TO DO
=====

* Use `GSL_RND_SEED` to properly propagate seeds of the random number
  generators.
  
* The documentation is still somewhat incomplete.

Executables
===========

The directory `ambit-stochastics/software` contains a few executables
which are described below.


simulate-trawl-process
----------------------

To see a short summary of the use of the program, try

    ./simulate-trawl-process --help
    
With the `--cascade` option, the program simulates from the process
`X(t) = L(A + (0,t))` where `L` is a homogeneous Levy basis and `A` is
an ambit set given by

    A = {(x,t) | 0 <= t <= T and |x| <= g(t)},

where

    g(t) = ((1 - (t / T) ^ s) / (1 + (L t / T) ^ s)) ^ (1 / s).
    
Then the process `Y = exp(X)` displays scaling of correlators when
`T/L << lag << T`. The parameter `s` determines the smoothness of the
correlator when the lag is near `T` and `T/L`.

The Levy basis can be either normal (specified with the `--normal`
option) or generalised hyperbolic (specified with the
`--generalised-hyperbolic` option). The parameters of the
corresponding Levy seed are `lambda`, `alpha`, `beta`, `mu`, and
`delta` in the generalised hyperbolic case. In the normal case, the
mean and variance are given by `mu` and `delta`, respectively.

The time step of the simulation is specified using the `--resolution`
option.  The option `--samples` specifies the length of the
simulation, i.e., the number of samples to generate. The output is
stored as the dataset `/simulation` in the HDF5 file specified by
`--output`. The output file is overwritten if it already exists.

Suppose we want to generate one million samples from the process `X`
such that the process `Y = exp(X)` has the following properties.

* One million samples are generated.
* The output is stored in `$HOME/output.h5`.
* The decorrelation length `T` is 1.
* The cascade length `L` is 1000.
* The smoothness parameter `s` is 2.
* A normal Levy basis is used.
* The mean of `Y` is 1.
* The correlator of exponent `tau_11` of order (1,1) of `Y` is 0.1.
* The time step is 0.005.

First, we calculate that the parameters `mu` and `delta` of the Levy
basis must be

    delta = 1 / 2 * tau_11 * L / T =  50,
    mu    = -1 / 2 * delta         = -25.

Next, we run the simulation.

    ./simulate-trawl-process --samples=1000000 --output=$HOME/output.h5 \
        --cascade --decorrelation-length=1 --cascade-length=1000 \
        --smoothness=2 --normal --mu=-25 --delta=50 --resolution=0.005

The simulation takes a few seconds to complete.

Suppose now that we want to do a similar simulation but with a normal
inverse Gaussian Levy basis with shape parameters `alpha` = 2.25 and
`beta` = -1.75. Let `gamma(p) = sqrt(alpha^2 - (beta + p)^2)`. Then

    delta = 1 / 2 * L / T * tau_11 / (2 * gamma(1) - gamma(0) - gamma(2)) = 84.41,
    mu    = -delta * (gamma(0) - gamma(1))                                = 59.69

Next, we run the simulation.

    ./simulate-trawl-process --samples=1000000 --output=$HOME/output.h5 \
        --cascade --decorrelation-length=1 --cascade-length=1000 \
        --smoothness=2 --generalised-hyperbolic --lambda=-0.5 --alpha=2.25 \
        --beta=-1.75 --mu=59.69 --delta=84.41 --resolution=0.005

It also takes a few seconds to complete, though a bit longer than when
the normal Levy basis is used.

Other shapes of ambit sets are supported, but the documentation on how
to use this feature still has to be written.

Libraries
=========

generalised-inverse-gaussian 
----------------------------

Sample from the generalised inverse Gaussian distribution. The library
supports the boundary cases of the parameter region.

multivariate-generalised-hyperbolic
-----------------------------------

Sample from multivariate and univariate generalised hyperbolic
distributions. The library supports the boundary cases of the
parameter region.

multivariate-normal
-------------------

Sample from multivariate and univariate normal distributions

trawl-process
-------------

Generate time series from trawl processes.

utilities
---------

Various small functions used by the other libraries.

