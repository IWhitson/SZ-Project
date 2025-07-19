# Contents of this folder

## Integrating signal over bandpasses

Since our bandpasses are quite wide, we have to average the frequency dependence of an underlying signal over the bandpass to predict the flux that we would measure with our instrument in the band.  For a real instrument, you would have a detailed (measured) bandpass shape; for AtLAST at the moment we are just assuming flat bandpasses within some range.

The methodology for averaging the bandpasses and converting between units is explained in [this Planck paper](https://arxiv.org/pdf/1303.5070) (see section 3.2.2).  It is implemented in `integrate_bp.py`.  I'm probably doing things the long way around in this code; in the past when dealing with Planck data we needed to convert between $y_\mathrm{SZ}$, $K_\mathrm{CMB}$ and MJy sr$^{-1}$ as explained in the paper; here we don't really need to do all of those conversions so we could simplify the code a bit but I have been too lazy to do so!

This code reads in the band parameters (upper and lower bounds) from `sensitivity_calculations.txt` which you'll also find in the folder.  It also contains resolutions and sensitivity estimates which will be handy later.

For the classical thermal SZ effect, everything is analytic and easy to calculate.  Things get more tricky for the relativistic SZ effect, which is still analytic but much more difficult to implement.  I use the [SZpack](https://www.jb.man.ac.uk/~jchluba/Science/SZpack/SZpack.html) implementation using [this python wrapper](https://toltec-astro.github.io/szpack_wrapper/).  I'm not sure how easy it would be to write a Julia wrapper so if this looks like it's in the too-hard basket just use the outputs of the code (also in the folder):

- `nu_eff.txt` contains effective frequencies for each band.  Since we have flat bandpasses these are just halfway inbetween the limits.  Mostly for plotting purposes.
- `KCMB2YSZ.txt` which contains factors to multiply by to convert from $K_\mathrm{CMB}$ to $y_\mathrm{tSZ}$ (the classical, thermal tSZ)
- `KCMB2YrSZ.txt` which contains the corresponding factors to convert to $y_\mathrm{tSZ}$ (relativistic thermal SZ) for a range of temperatures.  
- For computational efficiency, `YrSZ2KCMB_polyfits.txt` contains the coefficients of a polynomial fit to the rSZ conversion table.  You can see these plotted in `AtLAST_rSZ_signal.png`.  
- Finally, `KCMB2MJysr.txt` contains conversion factors to multiply $K_\mathrm{CMB}$ by to obtain MJy sr$^{-1}$, which are the units we want for a map of signal at a given frequency.

## Predicting SZ signal

`rSZ_spec.py` will predict SZ signal and errorbars for an AtLAST observation of a cluster, taking into account the sensitivities and resolutions of the different bands.  It normalizes to a requested SNR in a reference band (arbitrarily chosen to be band 8), so the $y$-parameter chosen doesn't matter.

This code uses the conversion factors in the files output in the previous step, so you won't need SZpack to run it.  It does also use SZpack just to make a spectrum for plotting purposes, but I've included a .npy file of the output you can use instead for a 10 keV cluster.  So you can just comment out the SZpack parts.