import astropy.io.fits as pyfits
import numpy as np
from matplotlib import pyplot
from scipy.integrate import trapz
import sys
#sys.path.append('/home/reinyv/software/SZpack.v1.1.1')
#import SZpack
from szpack_wrapper.sz import SZpack

def coth(x):
  return 1./np.tanh(x)

def Planck_der(nu, T):
  x=h_const*nu/k_B/T
  return 2*h_const*nu**3/c_light**2/(np.exp(x)-1)**2*np.exp(x)*x/T

def dI_SZ(nu, T):
  x=h_const*nu/k_B/T
  return Planck_der(nu,T)*T*(x*(np.exp(x)+1)/(np.exp(x)-1)-4)

def dI_rSZ(nu, T, Te):
  x=h_const*nu/k_B/T
  Dtau=mec2/Te # y=1.
  dI=np.zeros_like(nu)
  for i, xi in enumerate(x):
    # Not sure why the float conversion is necessary but getting a type error without it
    dI[i]=SZpack.compute_combo(float(xi), Dtau, Te, betac, muc, betao, muo)*xi**3*Dn_DI_conversion/1e20
  return dI

# Methodology of integrating over bandpasses from Planck 2013 results. IX. HFI spectral response

def KCMB2YSZ(nu, T, bp):
  b_prime=Planck_der(nu, T)
  numer=b_prime*bp
  numer[np.isnan(numer)]=0.
  numer=trapz(numer, nu)
  denom=bp*dI_SZ(nu,T)
  denom[np.isnan(denom)]=0.
  denom=trapz(denom, nu)
  return numer/denom

def KCMB2YrSZ(nu, T, Te, bp):
  b_prime=Planck_der(nu, T)
  numer=b_prime*bp
  numer[np.isnan(numer)]=0.
  numer=trapz(numer, nu)
  denom=bp*dI_rSZ(nu,T,Te)
  denom[np.isnan(denom)]=0.
  denom=trapz(denom, nu)
  return numer/denom

def KCMB2MJysr(nu, nu_c, T, bp):
  b_prime=Planck_der(nu, T)
  numer=b_prime*bp
  numer[np.isnan(numer)]=0.
  numer=trapz(numer, nu)
  denom=bp*(nu_c/nu)
  denom=trapz(denom, nu)
  return numer/denom*1e20

h_const=6.62607015e-34 # Js
k_B=1.38064852e-23 # J/K
mec2=511. # keV
Tcmb0=2.726 # K
c_light=299792458. # m/s
Dn_DI_conversion=13.33914078*Tcmb0**3

# AtLAST proposed bands
bands=np.genfromtxt('sensitivity_calculations.txt', dtype=None, names=True)
band_inds=bands['Band']
nband=len(band_inds)
band_freq=np.zeros((nband,2))
band_freq[:,0]=bands['nu1']
band_freq[:,1]=bands['nu2']


betac=0.0
muc=1.
betao=0.0
muo=0.
eps=1e-4

nu_eff=np.zeros((nband), dtype=float)
yconv=np.zeros_like(nu_eff)
Tconv=np.zeros_like(nu_eff)
Te=np.linspace(0, 20., 1000) # grid of conversion factors for interpolation
yconv_rel=np.zeros((nband, len(Te)))
for i, b in enumerate(band_inds):
  freqGHz=np.linspace(band_freq[i,0], band_freq[i,1], 100)
  transmission=np.ones_like(freqGHz) # assume flat bandpass
  nu_eff[i]=trapz(freqGHz*transmission, freqGHz)/trapz(transmission, freqGHz) # effective frequency
  yconv[i]=KCMB2YSZ(freqGHz*1e9, Tcmb0, transmission) # classical tSZ conversion factor
  Tconv[i]=KCMB2MJysr(freqGHz*1e9, nu_eff[i]*1e9, Tcmb0, transmission)
  for j, Tej in enumerate(Te):
    if Tej==0:
      yconv_rel[i,j]=KCMB2YSZ(freqGHz*1e9, Tcmb0, transmission)
    else:
      yconv_rel[i,j]=KCMB2YrSZ(freqGHz*1e9, Tcmb0, Tej, transmission)


# Fit a polynomial to the rSZ conversion factors for each frequency
from lmfit import Minimizer, Parameters, report_fit
def polynomial_res(params, x, data):
  model=params['p0']+params['p1']*x+params['p2']*x**2+params['p3']*x**3+params['p4']*x**4+params['p5']*x**5
  return model - data

poly_pars=np.zeros((nband, 6))
for i,  b in enumerate(band_inds):
  params=Parameters()
  params.add('p0', value=1./yconv_rel[i,0], vary=True)
  params.add('p1', value=0., vary=True)
  params.add('p2', value=0., vary=False)
  params.add('p3', value=0., vary=False)
  params.add('p4', value=0., vary=False)
  params.add('p5', value=0., vary=False)
  minner=Minimizer(polynomial_res, params, fcn_args=(Te, 1./yconv_rel[i,:]))
  result1=minner.minimize()
  model1=1./yconv_rel[i,:]+result1.residual
  params['p2'].set(vary=True)
  result2=minner.minimize()
  model2=1./yconv_rel[i,:]+result2.residual
  params['p3'].set(vary=True)
  result3=minner.minimize()
  model3=1./yconv_rel[i,:]+result3.residual
  params['p4'].set(vary=True)
  result4=minner.minimize()
  model4=1./yconv_rel[i,:]+result4.residual
  params['p5'].set(vary=True)
  result5=minner.minimize()
  model5=1./yconv_rel[i,:]+result5.residual
  pyplot.figure()
  pyplot.plot(Te, 1./yconv_rel[i,:], '.')
  pyplot.plot(Te, model1)
  pyplot.plot(Te, model2)
  pyplot.plot(Te, model3)
  pyplot.plot(Te, model4)
  pyplot.plot(Te, model5)
  pyplot.title('B%d' % b)
  #pyplot.show()
  pyplot.close('all')
  # Check how many polynomial parameters we need, take as many as necessary
  if np.max(np.abs(result1.residual*yconv_rel[i,:]))<1e-4:
    result=result1
  elif np.max(np.abs(result2.residual*yconv_rel[i,:]))<1e-4:
    result=result2
  elif np.max(np.abs(result3.residual*yconv_rel[i,:]))<1e-4:
    result=result3
  elif np.max(np.abs(result4.residual*yconv_rel[i,:]))<1e-4:
    result=result4
  elif np.max(np.abs(result5.residual*yconv_rel[i,:]))<1e-4:
    result=result5
  else:
    print('More polynomial degrees required for B', b)
    result=result5
  poly_pars[i,0]=result.params['p0'].value
  poly_pars[i,1]=result.params['p1'].value
  poly_pars[i,2]=result.params['p2'].value
  poly_pars[i,3]=result.params['p3'].value
  poly_pars[i,4]=result.params['p4'].value
  poly_pars[i,5]=result.params['p5'].value

# Output tables
np.savetxt('KCMB2YSZ.txt', yconv, fmt='%.6e')
np.savetxt('KCMB2MJysr.txt', Tconv, fmt='%.6e')
np.savetxt('nu_eff.txt', nu_eff)

np.savetxt('KCMB2YrSZ.txt', yconv_rel, fmt='%.6e', header='%.2f '*len(Te) % tuple(Te))

# Output polynomial fit parameters
np.savetxt('YrSZ2KCMB_polyfits.txt', poly_pars, fmt='%15.6e', header='%15s'*6 % ('p0', 'p1', 'p2', 'p3', 'p4', 'p5'))


# Compare bands to rel/non-rel spectra
Te=15.
y=1e-4
Dtau=mec2*y/Te
betac=0.0
muc=1
betao=0.0
muo=0
nu=np.logspace(np.log10(60), np.log10(2000), 100)*1e9 # GHz
x=h_const*nu/k_B/Tcmb0
Yrel=np.zeros_like(nu)
for i, xi in enumerate(x):
  Yrel[i]=SZpack.compute_combo(xi, Dtau, Te, betac, muc, betao, muo)

Yrel*=Dn_DI_conversion*x**3 # MJy/sr

# Non-relativistic
Y0 = (2*h_const/c_light**2)*(k_B*Tcmb0/h_const)**3*x**4*np.exp(x)/(np.exp(x)-1)**2*(x*coth(x/2)-4)*y*1e26*1e-6 # MJy/sr

fig, ax1 = pyplot.subplots()
ax1.plot(nu*1e-9, Y0, 'k', label=r'$T_{e}=0$')
ax1.plot(nu*1e-9, Yrel, 'k--', label=r'$T_{e}=15$keV')
ax1.legend(loc='best')

ax2 = ax1.twinx()
ax2.set_prop_cycle('color',[pyplot.cm.viridis(i) for i in np.linspace(0, 1, nband)])
for i, b in enumerate(band_inds):
  freqGHz=np.linspace(band_freq[i,0], band_freq[i,1], 100)
  transmission=np.ones_like(freqGHz)
  ax2.plot(freqGHz, transmission)

ax1.set_xscale('log')
ax1.set_ylim(-0.12, 0.2)
ax1.set_xlabel(r'$\nu$ / GHz')
ax1.set_ylabel(r'$\Delta I_{SZ}$ / MJy sr$^{-1}$')
ax2.set_ylabel('Transmission')
pyplot.savefig('transmission.png')
#pyplot.show()
pyplot.close('all')

# Plot the rSZ conversion factors and polynomial fits
fig, ax=pyplot.subplots(int(np.ceil(nband/2)),2,sharex=True)
fig.set_size_inches(8.3, 11.7)
Te=np.linspace(0, 20., 1000)
for i, b in enumerate(band_inds):
  irow=int(i/2)
  icol=int(i % 2)
  ax[irow,icol].plot(Te, 1./yconv_rel[i,:]*Tconv[i], ',')
  params=poly_pars[i,:]
  model=params[0]+params[1]*Te+params[2]*Te**2+params[3]*Te**3+params[4]*Te**4+params[5]*Te**5
  ax[irow,icol].plot(Te, model*Tconv[i])
  ax[irow,icol].text(0.05, 0.9, 'B%d, %d - %d GHz' % (b, band_freq[i,0], band_freq[i,1]), transform=ax[irow,icol].transAxes)

fig.subplots_adjust(hspace=0.)
ax_big=fig.add_subplot(1,1,1) #axisbg='none')
# Turn off axis lines and ticks of the big subplot
ax_big.spines['top'].set_color('none')
ax_big.spines['bottom'].set_color('none')
ax_big.spines['left'].set_color('none')
ax_big.spines['right'].set_color('none')
ax_big.xaxis.set_ticks([])
ax_big.yaxis.set_ticks([])
ax_big.set_xlabel(r'$T_{\mathrm{e}}$ / keV', labelpad=20, fontsize=14)
ax_big.set_ylabel('SZ signal / (MJy/sr)', labelpad=50, fontsize=14)
ax_big.set_title('AtLAST rSZ\n', fontsize=14)

for i, b in enumerate(band_inds):
  irow=int(i/2)
  icol=int(i % 2)
  ax[irow,icol].set_zorder(10)

fig.savefig('AtLAST_rSZ_signal.png')
