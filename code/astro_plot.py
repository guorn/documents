from astropy.io import fits
from astropy.wcs import WCS
from astropy.utils.data import get_pkg_data_filename
from astropy.visualization import astropy_mpl_style
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
import numpy as np
import math

# Define function
  # moment 0
def mom_0(data_cube, del_v):
    Nchan = data_cube.shape[1]
    Nx = data_cube.shape[2]
    Ny = data_cube.shape[3]
    mom0_data = np.zeros((Nx,Ny))
    for i in range(0,Nchan):
        mom0_data = mom0_data + data_cube[0,i]*abs(del_v)
    return mom0_data

# moment 1
def mom_1(data_cube, v0, del_v):
    Nchan = data_cube.shape[1]
    Nx = data_cube.shape[2]
    Ny = data_cube.shape[3]
    ''' calculate rms of a noise region for every channel '''
    noise = data_cube[0,:,250:390,10:150]
    sigma = np.zeros(data_cube.shape[1])
    for k in range(0,Nchan):
        sigma[k] = np.std(noise[k,:,:])
    level = 5*sigma
    ''' calculate mean velocity '''
    mom1_data = np.empty((Nx,Ny))
    mom1_data[:,:] = 'nan'
    for i in range(100,300):
        for j in range(100,300):
            A = 0
            B = 0
            for t in range(40,155):
                if data_cube[0,t,i,j] > level[t]:
                    v = v0 + t*del_v
                    A = A + v*data_cube[0,t,i,j]
                    B = B + data_cube[0,t,i,j]
                else:
                    A = A
                    B = B
            if B==0:
                mom1_data[i,j] = 'nan'
            else:
                mom1_data[i,j] = A/B
    return mom1_data, level

# moment 2
def mom_2(data_cube, v0, del_v):
    Nchan = data_cube.shape[1]
    Nx = data_cube.shape[2]
    Ny = data_cube.shape[3]
    mean_v, level = mom_1(data_cube, v0, del_v)
    mom2_data = np.empty((Nx,Ny))
    mom2_data[:,:] = 'nan'
    for i in range(100,300):
        for j in range(100,300):
            A = 0
            B = 0
            for t in range(40,155):
                if data_cube[0,t,i,j] > level[t]:
                    v_bar = mean_v[i,j]
                    v = v0 + t*del_v
                    A = A + data_cube[0,t,i,j]*(v-v_bar)**2
                    B = B + data_cube[0,t,i,j]
                else:
                    A = A
                    B = B
            if B==0:
                mom2_data[i,j] = 'nan'
            else:
                mom2_data[i,j] = math.sqrt(A/B)
    return mom2_data
    


# set the plot style
plt.style.use(astropy_mpl_style)

filename = get_pkg_data_filename('APEX_SB31_CO_6.fits')
hdu = fits.open(filename)

# look at file infomation
#print(hdu.info())

header = hdu[0].header
data = hdu[0].data
wcs = WCS(header).celestial
CRVAL3 = header['CRVAL3']/1000
CDELT3 = header['CDELT3']/1000
# set beam
BMAJ = header['BMAJ']
BMIN = header['BMIN']
BPA = header['BPA']
delta = abs(header['CDELT1'])
#print(BMAJ/delta,BMIN/delta)
beam = Ellipse(xy=(20,20), width=BMIN/delta, height=BMAJ/delta, angle=BPA, facecolor='None', edgecolor='black',alpha=1)

#look at header and data size
#print(header)
#print(data.shape[2])
# (1, 240, 400, 400)

# moments
mom0 = mom_0(data,CDELT3)
#mom1, sigma_5 = mom_1(data, CRVAL3, CDELT3)
#mom2 = mom_2(data, CRVAL3, CDELT3)
#print(mom1)

# Other marks
A = np.zeros(70)
for i in range(0,70):
    A[i] = i + 50
point_ra = [90, 95, 100, 105, 110, 115]
point_dec = []
for p in point_ra:
    point_dec.append(round(1.35*p-22))
print(point_ra)
print(point_dec)


# plot
fig = plt.figure()
ax = fig.add_subplot(111, projection=wcs)
#ax = fig.add_subplot(111)

ax.add_patch(beam)
#ax.imshow(data[0,100,:,:], cmap=plt.get_cmap('rainbow'))
main = ax.imshow(mom0[100:280,100:280], cmap=plt.get_cmap('jet'))
#plt.grid(b=None)

#ax.plot(A, 1.35*A-22, c='white', lw='1')
#ax.scatter(point_ra, point_dec, marker='o', c='white', s=20)

plt.xlabel('RA',fontsize=16)
plt.ylabel('DEC',fontsize=16)
plt.tick_params(labelsize=14)
plt.title('APEX_SB31 CO(2-1) moment 0')
plt.colorbar(main, label='(Jy/beam km/s)')


#plt.savefig('./moment0_mark.png')
plt.show()
