import os, sys, glob
import numpy as np
#sys.path.append('/home/reinyv/software/SZpack.v1.1.1')
#import SZpack
from szpack_wrapper.sz import SZpack
from matplotlib import pyplot
from copy import deepcopy

def coth(x):
  return 1./np.tanh(x)

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
bands=np.genfromtxt('new_bands.txt', dtype=None, names=True)
band_inds=bands['Band']
nband=len(band_inds)
band_freq=np.zeros((nband,2))
band_freq[:,0]=bands['nu_i']
band_freq[:,1]=bands['nu_f']

band_x=h_const*band_freq*1e9/k_B/Tcmb0

# Start with a relatively hot cluster
Te=10. # keV

# Let's assume y=1e-4 - plot spectra
y=1e-4

nu=np.logspace(np.log10(25), np.log10(2000), 100)*1e9 # GHz
x=h_const*nu/k_B/Tcmb0

# Non-relativistic
Y0 = (2*h_const/c_light**2)*(k_B*Tcmb0/h_const)**3*x**4*np.exp(x)/(np.exp(x)-1)**2*(x*coth(x/2)-4)*y*1e26*1e-6 # MJy/sr

if os.path.exists('Te_'+str(Te)+'_Yrel.npy'):
  Yrel=np.load('Te_'+str(Te)+'_Yrel.npy')
else:
  # Relativistic
  Dtau=mec2*y/Te
  betac=0.0
  muc=1
  betao=0.0
  muo=0
  Yrel=np.zeros_like(nu)
  for i, xi in enumerate(x):
    Yrel[i]=SZpack.compute_combo(xi, Dtau, Te, betac, muc, betao, muo)

  Yrel*=Dn_DI_conversion*x**3 # MJy/sr
  np.save('Te_'+str(Te)+'_Yrel.npy', Yrel)

fig=pyplot.figure()
fig.subplots_adjust(top=0.95, right=0.95)
ax0=fig.add_subplot(3,1,(1,2))
ax1=fig.add_subplot(3,1,3,sharex=ax0)
ax0.plot(nu*1e-9, Y0, 'k', label=r'$T_{e}=0$')
ax0.plot(nu*1e-9, Yrel, 'k--', label=r'$T_{e}='+str(Te)+'$keV')
ax1.plot(nu*1e-9, Yrel-Y0, 'k--')
ax1.axhline(y=0., color='k')
ax0.legend(loc='best')

# Predict signal in AtLAST bands
global SZsig
SZsig=np.zeros((nband))
nrSZsig=np.zeros_like(SZsig)
for i, b in enumerate(band_inds):
  SZsig[i]=polynomial(poly_pars[i,:], Te)*Tconv[i]*y
  nrSZsig[i]=1./yconv[i]*Tconv[i]*y
  ax0.plot(band_freq[i,:], [SZsig[i], SZsig[i]], 'k')
  ax1.plot(band_freq[i,:], (SZsig[i]-nrSZsig[i])*np.ones((2)), 'k')

# Assume maximum signal (B8) has a given SNR and give frequency-dependent errorbar to others based on instrumental simulations
ref=np.nonzero(band_inds==8)[0][0]
instr_sim_results=np.genfromtxt('sensitivity_calculations.txt', dtype=None, names=True)

# Error in mean SB in MJy/sr
refSNR=50
err_ref=SZsig[ref]/refSNR

# Calculations give sensitivity for a point source in mJy/beam after a given integration time
# Want to ask the question: if I map the same area for the same amount of time at different frequencies, what will be the noise ratios?
sens=instr_sim_results['sensitivity']*1e-3 # flux error in Jy

# Also need to account for beam sizes.  If I have flux = SB x Omega +/- D_flux then SB over one beam = flux/Omega +/- D_flux/Omega
# But then I'm taking a mean surface brightness over some area A so mean SB = sum(SB_i)/Nbeams +/- D_SB/sqrt(Nbeams) where Nbeams = A/Omega so mean SB error = D_flux/Omega / sqrt(A/Omega) = D_flux/sqrt(A*Omega)
# So mean SB err 1 = D_flux1/sqrt(A*Omega1)
# mean SB err 2 = D_flux2/sqrt(A*Omega2)
# mean SB err 1 / mean SB err 2 = D_flux1/D_flux2 / sqrt(Omega1/Omega2)
beam_sizes=instr_sim_results['beam']

global err
err=err_ref * sens/sens[ref] / (beam_sizes/beam_sizes[ref]) # error in average SB in MJy/sr

nu_eff=np.loadtxt('nu_eff.txt')
ax0.errorbar(nu_eff, SZsig, yerr=err, fmt='k.', markersize=10)
ax1.errorbar(nu_eff, SZsig-nrSZsig, yerr=err, fmt='k.', markersize=10)
ax0.set_xscale('log')
ax1.set_xlabel('Freq / GHz')
ax0.set_ylabel(r'$I$ / (MJy sr$^{-1}$)')
ax1.set_ylabel(r'$\Delta I$ / (MJy sr$^{-1}$)')
ax0.tick_params(axis='x', direction='in', which='both')
pyplot.subplots_adjust(wspace=0, hspace=0)
pyplot.savefig('Te_'+str(Te)+'_sig_errs.png')
#pyplot.show()
pyplot.close('all')

# Adding a dust component consistent with Erler+ 2018, eq 17
def mod_BB(nu, beta, T, A):
  nu0=857e9 # Hz
  return A*(nu/nu0)**(beta+3)*(np.exp(h_const*nu0/k_B/T)-1)/(np.exp(h_const*nu/k_B/T)-1)

# Use Erler+ 2018 parameter fits from stacked Planck data
Adust_857=0.2*0.52 # MJy/sr, multiplying by fudge factor to account for redshift dependence in Erler stacked sample
Tdust=18.44 # K
beta=1.5

dust_sig=mod_BB(np.sum(band_freq, axis=1)/2*1e9, beta, Tdust, Adust_857) # MJy/sr

fig=pyplot.figure()
fig.subplots_adjust(top=0.95, right=0.95)
ax0=fig.add_subplot(1,1,1)
ax0.plot(nu*1e-9, Y0, 'k', label=r'$T_{e}=0$')
ax0.plot(nu*1e-9, Yrel, 'k--', label=r'$T_{e}='+str(Te)+'$keV')
ax0.plot(nu*1e-9, mod_BB(nu, beta, Tdust, Adust_857), 'r:', label='Dust')
ax0.plot(nu*1e-9, Y0+mod_BB(nu, beta, Tdust, Adust_857), 'r')
ax0.plot(nu*1e-9, Yrel+mod_BB(nu, beta, Tdust, Adust_857), 'r--')
ax0.legend(loc='best')

for i, b in enumerate(band_inds):
  ax0.plot(band_freq[i,:], [SZsig[i]+dust_sig[i], SZsig[i]+dust_sig[i]], 'r')

ax0.errorbar(nu_eff, SZsig+dust_sig, yerr=err, fmt='r.', markersize=10)
ax0.set_xscale('log')
ax0.set_xlabel('Freq / GHz')
ax0.set_ylabel(r'$I$ / (MJy sr$^{-1}$)')
pyplot.savefig('Te_'+str(Te)+'_dust_sig_errs.png')
#pyplot.show()
pyplot.close('all')
