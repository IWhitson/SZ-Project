using DelimitedFiles

function polynomial(params, x)
    npoly = length(params)
    y = zeros(length(x))
    for i in 1:npoly
        y .+= params[i] * x .^ (i - 1)
    end
    return y
end

poly_pars = readdlm("YrSZ2KCMB_polyfits.txt")
Tconv = readdlm("KCMB2MJysr.txt")
yconv = readdlm("KCMB2YSZ.txt")

# Constants
h_const = 6.62607015e-34 # Js
k_B = 1.38064852e-23 # J/K
mec2 = 511.0 # keV
c_light = 299792458.0 # m/s
Tcmb0 = 2.726
Dn_DI_conversion = 13.33914078 * Tcmb0^3

# AtLAST proposed bands
bands = readdlm("sensitivity_calculations.txt", ',', header=true)
band_inds = bands[:, 1]
nband = length(band_inds)
band_freq = zeros(nband, 2)
band_freq[:, 1] = bands[:, 2]
band_freq[:, 2] = bands[:, 3]

# Start with a relatively hot cluster
Te = 10.0 # keV

# Let's assume y=1e-4 - plot spectra. We are going to normalise everything by SNR in one band so this is arbitrary
y = 1e-4

# Predict signal in AtLAST bands
SZsig = zeros(nband)
nrSZsig = zeros(nband)
for (i, b) in enumerate(band_inds)
    SZsig[i] = polynomial(poly_pars[i, :], Te) * Tconv[i] * y
    nrSZsig[i] = 1.0 / yconv[i] * Tconv[i] * y
end

# Assume maximum signal (B8 - arbitrary choice) has a given SNR and give frequency-dependent errorbar to others based on instrumental simulations
ref = findfirst(==(8), band_inds)
instr_sim_results = readdlm("sensitivity_calculations.txt", ',', header=true)

# Error in mean SB in MJy/sr
refSNR = 50
err_ref = SZsig[ref] / refSNR

# Need to account for different resolutions in different bands
# Calculations give sensitivity for a point source in mJy/beam after a given integration time
# Want to ask the question: if I map the same area for the same amount of time at different frequencies, what will be the noise ratios?
sens = instr_sim_results[:, 4] * 1e-3 # flux error in Jy

# Also need to account for beam sizes. If I have flux = SB x Omega +/- D_flux then SB over one beam = flux/Omega +/- D_flux/Omega
# But then I'm taking a mean surface brightness over some area A so mean SB = sum(SB_i)/Nbeams +/- D_SB/sqrt(Nbeams) where Nbeams = A/Omega so mean SB error = D_flux/Omega / sqrt(A/Omega) = D_flux/sqrt(A*Omega)
# So mean SB err 1 = D_flux1/sqrt(A*Omega1)
# mean SB err 2 = D_flux2/sqrt(A*Omega2)
# mean SB err 1 / mean SB err 2 = D_flux1/D_flux2 / sqrt(Omega1/Omega2)
beam_sizes = instr_sim_results[:, 5]

err = err_ref * sens ./ sens[ref] ./ (beam_sizes ./ beam_sizes[ref]) # error in average SB in MJy/sr

# From here, can output these estimates and errors (maybe shuffle estimates by a random amount consistent with the error bars) and proceed with fitting
tSZsig = hcat(nrSZsig, err)
writedlm("Te_$(Te)_tSZ_signal.txt", tSZsig)
rSZsig = hcat(SZsig, err)
writedlm("Te_$(Te)_rSZ_signal.txt", rSZsig)
