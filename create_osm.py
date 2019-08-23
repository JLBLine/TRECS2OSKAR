'''Jack Line, CIRA, 23/08/19'''
from astropy.io import fits
from astropy.table import Table
import argparse
from numpy import *
from copy import deepcopy
from numpy import random
from astropy.modeling.models import Gaussian2D
import matplotlib.pyplot as plt
# from my_plotting_lib import *
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--TRECS_tables', help='A text file containing a list of reduced TRECS FITS files to convert e.g. --TRECS_tables=files_to_convert.txt ')
parser.add_argument('--osm_name', help='Name for output OSKAR sky model, e.g. --osm_name=TRECS_sky.osm')
parser.add_argument('--plot_AGN', help='Include to produce an example plot of the AGN included in the .osm',default=False,action='store_true')
parser.add_argument('--plot_SFG', help='Include to produce an example plot of the SFG included in the .osm',default=False,action='store_true')

args = parser.parse_args()

factor = 2. * sqrt(2.*log(2.))

VELC = 299792458.0
bmax = 65e+3
D2R = pi / 180.0
freq = 95e+6

tel_reso = (VELC / freq) / bmax

def create_AGNs(fitsfile,num_sources,outfile):
    '''Takes a reduced TRECS AGN FITS table and turns it into .osm entries'''
    ##Grab the data
    t = Table.read(fitsfile)
    ras = t['RA']
    decs = t['DEC']
    S150 = t['S150']
    S160 = t['S160']
    Rs = t['Rs']
    size = t['size']

    ##I needed consistent results when testing before, but you can change or
    ##remove the seed here
    random.seed(9083745)
    ##Make some random position angles for the orientation of the AGN on the sky
    pa = random.uniform(low=0,high=2*pi,size=len(ras))

    ##Use the TRECS Rs and size params to work out the distance of the hotspots
    ##from the core
    dist_core = 0.5 * Rs * size

    ##Set the major axis of gaussian for the lobe by measuring the distance of
    ##the hotspot to the edge of the source
    lobe_maj = 0.5 * size - dist_core

    ##Make some random ratios to scale the minor axis of the lobe gaussian,
    ##gives some slight variety in the shape of the lobes
    lobe_ratio = random.uniform(low=0.8,high=1.0,size=len(ras))
    lobe_min = lobe_maj*lobe_ratio

    ##Set up parameters for the hot spots. Make them 10 times smaller than the
    ##lobes. Scale found by trial and error looking at output plots
    hot_maj = size / 10
    hot_min = lobe_ratio*size / 10

    ##Set up hot-spot locations relative to core
    ##size is in arcsec, ras in degrees
    ras_l = -dist_core / 3600.0
    ras_r = dist_core / 3600.0

    ##rotate the hot spot locations about the core, using the random pa
    ##in line with dec so dec = 0

    dec_l, dec_r = 0, 0
    ras_1 = ras_l*cos(-(pa-pi/2)) - dec_l*sin(-(pa-pi/2)) + ras
    decs_1 = dec_l*cos(-(pa-pi/2)) - ras_l*sin(-(pa-pi/2)) + decs

    ras_2 = ras_r*cos(-(pa-pi/2)) - dec_r*sin(-(pa-pi/2)) + ras
    decs_2 = dec_r*cos(-(pa-pi/2)) - ras_r*sin(-(pa-pi/2)) + decs

    ##Divide up the flux of the source between the core, lobes, and hotspots
    ##Again, these were fiddled to look good on a plot; this could be
    ##done better with physical motivation / made to be random
    core_frac = 0.00005
    hot_frac = 0.01
    lobe_frac = (1.0 - core_frac - hot_frac*2) / 2.0

    ##Calculate the SI of the source
    SIs = log(S150 / S160) / log(150.0 / 160.0)

    if args.plot_AGN:
        ##Plot the 10 largest AGN
        large = where((size > 300))

        fig = plt.figure(figsize=(10,10))
        for plot,ind in enumerate(large[0][:9]):
            ax = fig.add_subplot(3,3,plot+1)

            num_pixels = 101
            plot_width = size[ind] / (2*3600.0)
            x_range = linspace(ras[ind] - plot_width,ras[ind] + plot_width,num_pixels)
            y_range = linspace(decs[ind] - plot_width,decs[ind] + plot_width,num_pixels)

            reso = x_range[1] - x_range[0]

            x_stddev = lobe_maj[ind] / (factor*3600.0)
            y_stddev = lobe_min[ind] / (factor*3600.0)

            lobe1_func = Gaussian2D(amplitude=S150[ind]*lobe_frac, x_mean=ras_1[ind], y_mean=decs_1[ind], x_stddev=x_stddev, y_stddev=y_stddev,theta=pi/2 + pa[ind])
            lobe2_func = Gaussian2D(amplitude=S150[ind]*lobe_frac, x_mean=ras_2[ind], y_mean=decs_2[ind], x_stddev=x_stddev, y_stddev=y_stddev,theta=pi/2 + pa[ind])

            x_stddev = hot_maj[ind] / (factor*3600.0)
            y_stddev = hot_min[ind] / (factor*3600.0)

            hot1_func = Gaussian2D(amplitude=S150[ind]*hot_frac, x_mean=ras_1[ind], y_mean=decs_1[ind], x_stddev=x_stddev, y_stddev=y_stddev,theta=pi/2 + pa[ind])
            hot2_func = Gaussian2D(amplitude=S150[ind]*hot_frac, x_mean=ras_2[ind], y_mean=decs_2[ind], x_stddev=x_stddev, y_stddev=y_stddev,theta=pi/2 + pa[ind])

            x_mesh, y_mesh = meshgrid(x_range,y_range)
            lobe1 = lobe1_func(x_mesh,y_mesh)
            lobe2 = lobe2_func(x_mesh,y_mesh)

            lobe1 *=  S150[ind]*lobe_frac / lobe1.sum()
            lobe2 *=  S150[ind]*lobe_frac / lobe2.sum()

            hot1 = hot1_func(x_mesh,y_mesh)
            hot2 = hot2_func(x_mesh,y_mesh)

            hot1 *=  S150[ind]*hot_frac / hot1.sum()
            hot2 *=  S150[ind]*hot_frac / hot2.sum()

            gal_array = hot1 + hot2 + lobe1 + lobe2
            gal_array[num_pixels // 2,num_pixels // 2] += S150[ind]*core_frac

            im = ax.imshow(gal_array,origin='lower',extent=[x_range[0],x_range[-1],y_range[0],y_range[-1]])

        plt.tight_layout()
        fig.savefig('example_AGNs_in_osm.png',bbox_inches='tight')

    ##Write the source into the .osm file
    for ind in arange(len(ras)):
        if size[ind] < 30:
            ##If less than 30 arcsec, just stick in a point source
            ##This was due to the resolution of my simulation - save some comutation with just a point source
            outfile.write('%.7f %.7f %.9f 0 0 0 150e+6 %.2f 0 0 0 0\n' %(ras[ind],decs[ind],S150[ind],SIs[ind]))
        else:

            ##Write out core - assume is a point source
            ##Scale the flux according the settings above
            outfile.write('%.7f %.7f %.9f 0 0 0 150e+6 %.2f 0 0 0 0\n' %(ras[ind],decs[ind],S150[ind]*core_frac,SIs[ind]))

            ##Then write out two lobes and hotspots either side of the code
            ##Scale the flux according the settings above
            outfile.write('%.7f %.7f %.9f 0 0 0 150e+6 %.2f 0 %.2f %.2f %.2f\n' %(ras_1[ind],decs_1[ind],S150[ind]*lobe_frac,SIs[ind],lobe_maj[ind],lobe_min[ind],pa[ind]*(180/pi)))
            outfile.write('%.7f %.7f %.9f 0 0 0 150e+6 %.2f 0 %.2f %.2f %.2f\n' %(ras_2[ind],decs_2[ind],S150[ind]*lobe_frac,SIs[ind],lobe_maj[ind],lobe_min[ind],pa[ind]*(180/pi)))
            outfile.write('%.7f %.7f %.9f 0 0 0 150e+6 %.2f 0 %.2f %.2f %.2f\n' %(ras_1[ind],decs_1[ind],S150[ind]*hot_frac,SIs[ind],hot_maj[ind],hot_min[ind],pa[ind]*(180/pi)))
            outfile.write('%.7f %.7f %.9f 0 0 0 150e+6 %.2f 0 %.2f %.2f %.2f\n' %(ras_2[ind],decs_2[ind],S150[ind]*hot_frac,SIs[ind],hot_maj[ind],hot_min[ind],pa[ind]*(180/pi)))
    #
    return num_sources + len(ras)

def create_SFGs(fitsfile,num_sources,outfile):
    '''Takes a reduced TRECS SFG FITS table and turns it into .osm entries'''
    ##Read in data
    t = Table.read(fitsfile)
    ras = t['RA']
    decs = t['DEC']
    S150 = t['S150']
    S160 = t['S160']
    ##Calc SI
    SIs = log(S150 / S160) / log(150.0 / 160.0)
    e1 = t['e1']
    e2 = t['e2']
    size = t['size']

    ##Give the sources random position angles
    random.seed(345673547)
    pa = random.uniform(low=0,high=2*pi,size=len(ras))

    ##Use TRECS e1,e2 values to make a gaussian maj/min axis
    ##No idea if this really makes sense to do
    maj_ax = size * abs(e1)
    min_ax = size * abs(e2)

    if args.plot_SFG:
        large = where((size > 100))
        fig = plt.figure(figsize=(10,10))
        for plot,ind in enumerate(large[0][:9]):
            ax = fig.add_subplot(3,3,plot+1)

            num_pixels = 101
            plot_width = size[ind] / (2*3600.0)
            x_range = linspace(ras[ind] - plot_width,ras[ind] + plot_width,num_pixels)
            y_range = linspace(decs[ind] - plot_width,decs[ind] + plot_width,num_pixels)
            x_mesh, y_mesh = meshgrid(x_range,y_range)

            reso = x_range[1] - x_range[0]

            x_stddev = maj_ax[ind] / (factor*3600.0)
            y_stddev = min_ax[ind] / (factor*3600.0)

            sfg_func = Gaussian2D(amplitude=S150[ind], x_mean=ras[ind], y_mean=decs[ind], x_stddev=x_stddev, y_stddev=y_stddev,theta=pi/2 + pa[ind])

            sfg = sfg_func(x_mesh,y_mesh)
            sfg *=  S150[ind]*sfg / sfg.sum()

            im = ax.imshow(sfg,origin='lower',extent=[x_range[0],x_range[-1],y_range[0],y_range[-1]])

        plt.tight_layout()
        fig.savefig('example_SFGs_in_osm.png',bbox_inches='tight')

    for ind in arange(len(ras)):

        if size[ind] < 30:
            ##If less than 30 arcsec, just stick in a point source
            ##Again, just
            outfile.write('%.7f %.7f %.9f 0 0 0 150e+6 %.2f 0 0 0 0\n' %(ras[ind],decs[ind],S150[ind],SIs[ind]))
        else:
            ##Then write a gaussian
            outfile.write('%.7f %.7f %.9f 0 0 0 150e+6 %.2f 0 %.2f %.2f %.2f\n' %(ras[ind],decs[ind],S150[ind],SIs[ind],maj_ax[ind],min_ax[ind],pa[ind]*(180/pi)))

    return num_sources + len(ras)

num_sources = 0
outfile = open(args.osm_name,'w+')

files = open(args.TRECS_tables,'r').read().split('\n')

for file in files:
    if 'AGN' in file:
        num_sources = create_AGNs(file,num_sources,outfile)
    elif 'SFG' in file:
        num_sources = create_SFGs(file,num_sources,outfile)

print('Sky model contains %d sources' %num_sources)
