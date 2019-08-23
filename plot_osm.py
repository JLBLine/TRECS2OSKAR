'''Jack Line, CIRA, 23/08/19'''
from numpy import *
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.modeling.models import Gaussian2D

##Swaps between major/minor and standard deviation for gaussians 
factor = 2. * sqrt(2.*log(2.))

data = loadtxt('TRECS_sky.osm')

ras = data[:,0]
decs = data[:,1]
S150s = data[:,2]
majs = data[:,9]
mins = data[:,10]
PAs = data[:,11]

reso = 0.006

target_header = fits.Header()

nside = 4096

target_header['NAXIS'] = 2
target_header['NAXIS1'] = nside
target_header['NAXIS2'] = nside
target_header['CTYPE1'] = 'RA---SIN'
target_header['CRPIX1'] = nside / 2 + 1
target_header['CRVAL1'] = 0.2372726728
target_header['CDELT1'] = -reso
target_header['CUNIT1'] = 'deg '
target_header['CTYPE2'] = 'DEC--SIN'
target_header['CRPIX2'] = nside / 2 + 1
target_header['CRVAL2'] = -26.7932606719
target_header['CDELT2'] = reso
target_header['CUNIT2'] = 'deg '
target_header['EQUINOX']= 2000.

wcs = WCS(target_header)

x_locs,y_locs = wcs.all_world2pix(ras,decs,0)

edge_size = 50

good = where((1+ edge_size <= x_locs) & (x_locs < nside - edge_size- 1) & (1 + edge_size <= y_locs) & (y_locs < nside - edge_size -1))

sky = zeros((nside,nside))

def plot_gauss(x_loc,y_loc,S150,maj,min,PA):

    x_loc = int(round(x_loc))
    y_loc = int(round(y_loc))

    x_stddev = maj / (factor*reso*3600.0)
    y_stddev = min / (factor*reso*3600.0)

    gauss_func = Gaussian2D(amplitude=S150, x_mean=0, y_mean=0, x_stddev=x_stddev, y_stddev=y_stddev,theta=pi/2 + PA*(pi/180.0))


    x_mesh,y_mesh = meshgrid(linspace(-edge_size,edge_size,edge_size*2+1),linspace(-edge_size,edge_size,edge_size*2+1))

    gauss = gauss_func(x_mesh,y_mesh)

    gauss *=  S150 / gauss.sum()


    sky[-edge_size+y_loc:edge_size+1+y_loc,-edge_size+x_loc:edge_size+1+x_loc] += gauss


def plot_point(x_loc,y_loc,S150):
    x_loc = int(round(x_loc))
    y_loc = int(round(y_loc))

    sky[y_loc,x_loc] += S150

for ind in good[0]:
    if majs[ind] == 0.0:
        plot_point(x_locs[ind],y_locs[ind],S150s[ind])
    else:
        plot_gauss(x_locs[ind],y_locs[ind],S150s[ind],majs[ind],mins[ind],PAs[ind])

fits.writeto('TRECS_sky_raw.fits', sky, target_header, overwrite=True)

from scipy.signal import fftconvolve

BMAJ    =   0.0135435498617062
BMIN    =   0.00907275012433796
BPA     =      135.923459748969 * (pi/180.0)

x_stddev = BMAJ / (factor*reso)
y_stddev = BMIN / (factor*reso)

rest_gauss_func = Gaussian2D(amplitude=1, x_mean=0, y_mean=0, x_stddev=x_stddev, y_stddev=y_stddev,theta=pi/2 + BPA)

xrange = arange(-25,26)
yrange = arange(-25,26)

x_mesh, y_mesh = meshgrid(xrange,yrange)
rest_gauss_kern = rest_gauss_func(x_mesh,y_mesh)
rest_gauss_kern /= rest_gauss_kern.sum()

convolved_sky = fftconvolve(sky, rest_gauss_kern, 'same')

solid_beam = (pi*BMAJ*BMIN) / (4*log(2))
solid_pixel = reso**2

convert2pixel = solid_pixel/solid_beam

convolved_sky /= convert2pixel

fits.writeto('TRECS_sky_convolved.fits', convolved_sky, target_header, overwrite=True)
