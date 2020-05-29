from astropy.io import fits
from astropy.wcs import WCS
from astropy.utils.data import get_pkg_data_filename
from astropy.visualization import astropy_mpl_style
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from matplotlib.patches import Ellipse
from matplotlib import rcParams
import matplotlib.pyplot as plt
import numpy as np
import math
plt.style.use(astropy_mpl_style)
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman']

''' constant '''
#cosmo = FlatLambdaCDM(H0=67.8*u.km / u.s / u.Mpc, Om0=0.308)
cosmo = FlatLambdaCDM(H0=70*u.km / u.s / u.Mpc, Om0=0.3)
c = 3e5  #km/s
redshift = 0.0026
Mpc = 1e6*3.09*1e16*1e2  #cm
d_L = cosmo.luminosity_distance(redshift)  #Mpc
d_A = cosmo.angular_diameter_distance(redshift)  #Mpc
L_sun = 3.86*10**26 #W (J/s)
Z_sun = 0.02
erg = 10**(-7) #J
pix_ang = (0.2*math.pi)/(180*3600)
pix_area = (pix_ang*1000*d_A.value)**2  #kpc^2
Del_lamb = 1.25

#print(d_L.value)
#print(pix_area)

Ha_0 = 6563
Hb_0 = 4861
OIII4959_0 = 4959
OIII5007_0 = 5007
NII6548_0 = 6548
NII6583_0 = 6583
SII6717_0 = 6717
SII6731_0 = 6731

Ha_R = 2547
Hb_R = 1696
OIII4959_R = 1745
OIII5007_R = 1771
NII6548_R = 2541
NII6583_R = 2557
SII6717_R = 2623
SII6731_R = 2633



''' Define function '''
def emlines_read(number):
    flux = emlines_flux[number, :, :]
    flux_err = emlines_flux_err[number, :, :]
    center = emlines_center[number, :, :]
    sigma = emlines_sigma[number, :, :]
    ew = emlines_ew[number, :, :]
    return flux, flux_err, center, sigma, ew

def k_lamb(wave):
    Rv = 4.05
    lamb = wave/10000
    if 0.12<=lamb<0.63:
        k = 2.659*(-2.156+1.509/lamb-0.198/lamb**2+0.011/lamb**3)+Rv
    elif 0.63<=lamb<=2.20:
        k = 2.659*(-1.857+1.040/lamb)+Rv
    return k

def mask(A, A_err):
    mask = np.zeros(A.shape)
    [rows, cols] = snr.shape
    for i in range(rows):
        for j in range(cols):
            ratio = A[i,j]/A_err[i,j]
            if ratio > 3:
                mask[i,j] = 1
            else:
                mask[i,j] = 'nan'
    return mask

def ebvgas_mask(ebv_gas):
    ebvmask = np.zeros(ebv_gas.shape)
    [rows, cols] = ebvmask.shape
    for i in range(rows):
        for j in range(cols):
            ebvg = ebv_gas[i,j]
            if np.isnan(ebvg)==True:
                ebvmask[i,j] = 'nan'
            else:
                ebvmask[i,j] = 1
    return ebvmask

def moment(number, lamb_0, R):
    sigma_fsl = lamb_0/(2.355*R)
    print('FSL=', sigma_fsl)
    flux, flux_err, center, sigma, ew = emlines_read(number)
    Amask = mask(flux, flux_err)
    mom1_all = (c/lamb_0)*(center-np.ones(center.shape)*lamb_0)+ve
    mom1 = mom1_all*Amask
    sigma_tar = np.power(np.true_divide(np.power(sigma, 2), sigma_fsl**2), 0.5)
    mom2 = (c/lamb_0)*sigma_tar*Amask
    return mom1, mom2

def correct(number, lamb_0):
    u = k_lamb(lamb_0)
    print('k(lambda)=',u)
    flux, flux_err, center, sigma, ew = emlines_read(number)
    Bmask = mask(flux, flux_err)
    flux_corr = np.zeros(flux.shape)
    [rows, cols] = flux.shape
    for i in range(rows):
        for j in range(cols):
            ebvg = ebv_gas[i,j]
            if np.isnan(ebvg)==True:
                flux_corr[i,j] = 'nan'
            else:
                flux_corr[i,j] = flux[i,j]*10**(0.4*ebvg*u)
            flux_end = flux_corr*Bmask
    return flux_end

def color_excess():
    Hb_flux, Hb_flux_err, Hb_center, Hb_sigma, Hb_ew = emlines_read(8)
    Ha_flux, Ha_flux_err, Ha_center, Ha_sigma, Ha_ew = emlines_read(15)
    Ha_mask = mask(Ha_flux, Ha_flux_err)
    Hb_mask = mask(Hb_flux, Hb_flux_err)
    obs_ratio = np.true_divide(Ha_flux, Hb_flux)
    select_ratio = np.zeros(obs_ratio.shape)
    [rows, cols] = select_ratio.shape
    for i in range(rows):
        for j in range(cols):
            value = obs_ratio[i,j]
            if value>=2:
                select_ratio[i,j] = 1
            else:
                select_ratio[i,j] = 'nan'
    mask_ratio = obs_ratio*Ha_mask*Hb_mask*select_ratio
    ebv_gas = (2.5/1.27)*np.log10(np.true_divide(mask_ratio, 2.86))
    return ebv_gas

def cont_mask(cont_flux, cont_err):
    cmask = np.zeros(ebv.shape)
    cont_snr = np.zeros(ebv.shape)
    [rows, cols] = cmask.shape
    for i in range(rows):
        for j in range(cols):
            signal = np.sum(cont_flux[610:690, i, j])
            noise = np.sum(cont_err[610:690, i, j])
            cont_snr[i,j] = signal/noise
            if signal/noise>=50:
                cmask[i,j] = 1
            else:
                cmask[i,j] = 'nan'
    return cmask, cont_snr

def suf_brightness(flux, d_L, pix_area):
    factor = 1e-20*(4*math.pi*(d_L.value*Mpc)**2)/pix_area
    L = flux*factor
    return L
    

''' Input file '''
filename = get_pkg_data_filename('NGC3351-MAPS1.fits.gz')
hdu = fits.open(filename)
header = hdu[27].header
wcs = WCS(header).celestial

snr = hdu[1].data
nor = hdu[2].data
ve = hdu[3].data
ve_err = hdu[4].data
vd = hdu[5].data
weights = hdu[7].data
ebv = hdu[8].data
ebv_err = hdu[9].data
dof = hdu[10].data
chi2 = hdu[11].data
M = hdu[12].data
age_L = hdu[13].data
z_L = hdu[14].data
age_M = hdu[15].data
z_M = hdu[16].data
age_L2 = hdu[17].data
z_L2 = hdu[18].data
age_M2 = hdu[19].data
z_M2 = hdu[20].data
emlines_center = hdu[23].data
emlines_sigma = hdu[25].data
emlines_flux = hdu[27].data
emlines_flux_err = hdu[28].data
emlines_ew = hdu[29].data
lick_hd_a = hdu[34].data

''' Input OUTCUBE file '''
outcubename = get_pkg_data_filename('NGC3351-OUTCUBE.fits.gz')
spec = fits.open(outcubename)
cont_flux = spec[1].data
cont_err = spec[2].data

''' calculate continuum snr and mask ebv '''
cmask, cont_snr = cont_mask(cont_flux, cont_err)
# wave region is [5452, 5553]
mask_ebv = ebv*cmask

''' emission line '''
#Hb_flux, Hb_flux_err, Hb_center, Hb_sigma, Hb_ew = emlines_read(8)
#Ha_flux, Ha_flux_err, Ha_center, Ha_sigma, Ha_ew = emlines_read(15)

''' Calculate Balmer decrement and E(B-V) of ionized gas '''
#ebv_gas = color_excess()
#np.save('./picture0526/ebv_gas.npy', ebv_gas)

''' Stellar mass density '''
#M_dens = np.true_divide(M, pix_area)

''' Metallicity '''
#z_light = np.true_divide(z_L, Z_sun)
#z_mass = np.true_divide(z_M, Z_sun)
                         
''' Moments '''
ebv_gas = np.load('./picture0526/ebv_gas.npy')
ebvmask = ebvgas_mask(ebv_gas)
mom1, mom2 = moment(18, SII6731_0, SII6731_R)
mask_mom1 = mom1*ebvmask
mask_mom2 = mom2*ebvmask

''' Corrected emission line '''
#ebv_gas = np.load('./picture0526/ebv_gas.npy')
#mom0 = correct(18, SII6731_0)
#L_dens = suf_brightness(mom0, d_L, pix_area)  #erg/s/kpc^2
#L_dens_simple = ((4*math.pi*(1+redshift)**4)/pix_ang**2)*mom0*1e-20*(Mpc/1000)**2
#print(L_dens[150,152])
#print(L_dens_simple[150,152])



''' Plot '''
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection=wcs)
main = ax.imshow(mask_ebv, cmap=plt.get_cmap('jet'))
# map color: 'plasma', 'jet'
plt.grid(b=None)
plt.xlabel('RA', fontsize=14)
plt.ylabel('DEC', fontsize=14)
plt.tick_params(labelsize=14)
plt.title('Velocity',fontsize=16)
cb = plt.colorbar(main, label=r'$V_{SII6731}\ (km/s)$')
cb.ax.tick_params(labelsize=16)
#plt.savefig('./picture0526/SII6731_mom1.png')

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection=wcs)
#main = ax.imshow(mask_mom2, cmap=plt.get_cmap('jet'), vmax=180)
plt.grid(b=None)
plt.xlabel('RA',fontsize=14)
plt.ylabel('DEC',fontsize=14)
plt.tick_params(labelsize=14)
plt.title('Velocity dispersion',fontsize=16)
cb = plt.colorbar(main, label=r'$\sigma_{SII6731}\ (km/s)$')
cb.ax.tick_params(labelsize=16)
#plt.savefig('./picture0526/SII6731_mom2.png')


plt.show()

    