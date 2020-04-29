# DUNE-PRISM Flux Fitting

## Overview

This package contains some useful scripts for making plots and simple analyses using the PRISM method.

## Dependencies

This package relies on a working installation of ROOT with Python support

[root.cern.ch](root.cern.ch)

Numpy and Matplotlib are also required.

Lastly, you will need some root files containing at lease one near detector flux (as a function of energy and off-axis angle) and one far detector flux (as a function of energy).  If you are unable to find these, please contact the author.

## Environment

The following environment variables *must* be set:

`ROOTSYS`: points to the install location of the CERN ROOT package.  Usually set by sourcing `/path/to/ROOT/installation/bin/thisroot.sh`.

These variables are not required, but highly recommended:

`PROB3ROOT`: points to the install location of the python-compatible Prob3 package (https://github.com/DanielMDouglas/Prob3).  Currently, there is no setup script, but this will come soon!  If this is unset, `oscProbs.py` will look for `./Prob3/` instead. 

`DP_FLUX_FILE`: points to the root file containing nominal and systematic fluxes.  If this is unset, `fluxes.py` will look for `../flux/syst/DUNE_Flux_OffAxis_Nov2017Review_syst_shifts_fine.root` instead.

## Components

### fluxes.py

This program defines a flux class and exposes a dictionary of nominal fluxes.  These fluxes can be accessed with the normal dictionary syntax, e.g. `ND_nominal[beamMode][flavor]`.  In order to avoid excessive loading times, the flux class does not load any data until the `load` method is called, which will return a `numpy` array (1-D for FD and 2-D for ND) containing the contents of the flux histograms.

Allowed beam modes are "nu" and "nubar" and allowed flavors are "nue", "numu", "nuebar", and "numubar".

### oscProbs.py

This program defines an oscillation probability class.  A new oscillation hypothesis is instantiated like `oscProb(fromFlavor, toFlavor, sin2theta23, dm32, deltaCP)`.  Other oscillation parameters can be added to the initializer trivially, but these tend to be the most useful ones, so this functionality has not been added yet.

Since calculating oscillation probability entails slightly more computation, the initializer will first check to see if these particular parameters/flavors have been calculated before (in which case they will exist inside `DP_OSC_PROB_DIR`) and, if not, will calculated them.  Like the flux object, this data will not actually by loaded until the `load` method is invoked.

### flux_fitter.py

This program contains the bulk of the "PRISM" methods, namely the `flux_fitter` object.  This object contains relevant data, including FD flux, oscillation hypothesis, ND flux, as well as parameters to the fit, such as energy range, out-of-region weights, and maximum off-axis distance.

Important methods of this object are:

#### `set_fit_region(energies = None, peaks = None)`

This method takes a list of length 2 of either energies or peaks.

If the keyword argument `energies` is supplied, the fit will simply use these energies as lower and upper bounds, respectively.

If the keyword argument `peaks` is supplied instead, the fitter will attempt to find the locations of these peaks and use those as the bounds instead.  Since the oscillation tends to be faster at lower energies, these peaks are counted from the right.  For example, `peaks = [1, 3]` will try to find the highest and third-highest energy peaks and use those.

#### `set_OOR(weights)`

This method sets the weights for the fit outside of the energy region described above.  The default value is `(0, 0)`, meaning that the value of the target flux outside of the region makes no difference to the fit at all.  Usually, however, a small weight on the low-energy side is preferred so that the combined flux does not become too large in the low-energy regions.  For this, about 5% is enough (e.g., `set_OOR((0.05, 0))`).

#### `set_oscHypothesis(Posc)`

Often, it is desireable to change the oscillation hypothesis without instantiating a new `flux_fitter` object.  This method takes a `numpy` array of oscillation probability (the output of a `oscProb.load()` call) and updates the `Posc`, `FD_oscillated`, and `target` data.

#### `set_maxOA(maxOA)`

This method is useful for studying the effects of differing ND hall size on the quality of fits.  Calling this method will update the `maxOA` and `ND` data.

#### `calc_coeffs(reg)`

This is perhaps the most interesting part of the flux fitting procedure.  The current version of the coefficient calculation uses [Tikhonov Regularization](https://en.wikipedia.org/wiki/Tikhonov_regularization) to find a solution for the underconstrained problem `dot(ND, c) = FD`.

The regularization factor, `reg`, quantifies the magnitude of the penalty for very rapidly varying coefficients.  A larger value for `reg` will produce a very smoothly varying (as a function of off-axis position) set of coefficients, while `reg = 0` will reduce to a simple matrix inversion: `c = np.dot(ND.I, FD)`, which will produce a very noisy set of coefficients.  For the usual FD and N normalization, a value of `reg` around 1.e-8 usually produces a sufficiently smooth set of coefficients.

## Contributing

Contributions are definitely welcome!  Please submit pull requests through Github or reach out to me at dougl215@msu.edu if you would like to be added to this project.