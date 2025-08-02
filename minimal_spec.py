import numpy as np

def polynomial(params, x):
  npoly=len(params)
  y=np.zeros_like(x)
  for i in range(npoly):
    y+=params[i]*x**i
  return y

poly_pars=np.loadtxt('YrSZ2KCMB_polyfits.txt')
Tconv=np.loadtxt('KCMB2MJysr.txt')
yconv=np.loadtxt('KCMB2YSZ.txt')

# Constants
h_const=6.62607015e-34 # Js
k_B=1.38064852e-23 # J/K
mec2=511. # keV
c_light=299792458 # m/s
Tcmb0=2.726
Dn_DI_conversion=13.33914078*Tcmb0**3

# AtLAST proposed bands
bands=np.genfromtxt('sensitivity_calculations.txt', dtype=None, names=True)
band_inds=bands['Band']
nband=len(band_inds)
band_freq=np.zeros((nband,2))
band_freq[:,0]=bands['nu1']
band_freq[:,1]=bands['nu2']

# Start with a relatively hot cluster
Te=10. # keV

# Let's assume y=1e-4 - plot spectra.  We are going to normalise everything by SNR in one band so this is arbitrary
y=1e-4

# Predict signal in AtLAST bands
SZsig=np.zeros((nband))
nrSZsig=np.zeros_like(SZsig)
for i, b in enumerate(band_inds):
  SZsig[i]=polynomial(poly_pars[i,:], Te)*Tconv[i]*y
  nrSZsig[i]=1./yconv[i]*Tconv[i]*y

# Assume maximum signal (B8 - arbitrary choice) has a given SNR and give frequency-dependent errorbar to others based on instrumental simulations
ref=np.nonzero(band_inds==8)[0][0]
instr_sim_results=np.genfromtxt('sensitivity_calculations.txt', dtype=None, names=True)

# Error in mean SB in MJy/sr
refSNR=50
err_ref=SZsig[ref]/refSNR

# Need to account for different resolutions in different bands
# Calculations give sensitivity for a point source in mJy/beam after a given integration time
# Want to ask the question: if I map the same area for the same amount of time at different frequencies, what will be the noise ratios?
sens=instr_sim_results['sensitivity']*1e-3 # flux error in Jy

# Also need to account for beam sizes.  If I have flux = SB x Omega +/- D_flux then SB over one beam = flux/Omega +/- D_flux/Omega
# But then I'm taking a mean surface brightness over some area A so mean SB = sum(SB_i)/Nbeams +/- D_SB/sqrt(Nbeams) where Nbeams = A/Omega so mean SB error = D_flux/Omega / sqrt(A/Omega) = D_flux/sqrt(A*Omega)
# So mean SB err 1 = D_flux1/sqrt(A*Omega1)
# mean SB err 2 = D_flux2/sqrt(A*Omega2)
# mean SB err 1 / mean SB err 2 = D_flux1/D_flux2 / sqrt(Omega1/Omega2)
beam_sizes=instr_sim_results['beam']

err=err_ref * sens/sens[ref] / (beam_sizes/beam_sizes[ref]) # error in average SB in MJy/sr

# From here, can output these estimates and errors (maybe shuffle estimates by a random amount consistent with the error bars) and proceed with fitting
tSZsig = np.vstack((nrSZsig, err)).T
np.savetxt('Te_'+str(Te)+'_tSZ_signal.txt', tSZsig)
rSZsig = np.vstack((SZsig, err)).T
np.savetxt('Te_'+str(Te)+'_rSZ_signal.txt', rSZsig)
