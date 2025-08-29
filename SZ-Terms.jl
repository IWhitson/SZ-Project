#Cell 1 Graph
# --- Load conversion tables ---
poly_pars = load_numeric_matrix("YrSZ2KCMB_polyfits.txt")
Tconv = load_numeric_matrix("KCMB2MJysr.txt")[:, 1]
yconv = load_numeric_matrix("KCMB2YSZ.txt")[:, 1]
	
# --- Load band definitions and error bars ---
bands = load_numeric_matrix("sensitivity_calculations.txt")
band_inds = bands[:, 1]
band_errs = bands[:, 6]
nband = length(band_inds)
y = 1e-4
Te = 10.0

#Cell 2 Chi^2
y_obs = 1.1e-4
observed = [1.0 / yconv[i] * Tconv[i] * y_obs for i in 1:length(band_inds)]
