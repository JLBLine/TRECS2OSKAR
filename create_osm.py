'''Jack Line, CIRA, 23/08/19'''
'''Modified by Anna Bonaldi SKAO 11/11/19'''

from astropy.io import fits
from astropy.table import Table
import argparse
import numpy as np
from copy import deepcopy
from numpy import random
from astropy.modeling.models import Gaussian2D
import matplotlib.pyplot as plt
# from my_plotting_lib import *
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--TRECS_tables', help='A text file containing a list of reduced TRECS FITS files to convert e.g. --TRECS_tables=files_to_convert.txt ')
parser.add_argument('--freq1', help='center frequency')
parser.add_argument('--freq2', help='nearby frequency for spectral index calculation')
parser.add_argument('--osm_name', help='Name for output OSKAR sky model, e.g. --osm_name=TRECS_sky.osm')
parser.add_argument('--nsources', help='For testing purposes, limit to a fixed number of sources in the .osm',default='-1')
parser.add_argument('--plot_AGN', help='Include to produce an example plot of the AGN included in the .osm',default=False,action='store_true')
parser.add_argument('--plot_SFG', help='Include to produce an example plot of the SFG included in the .osm',default=False,action='store_true')

#parser.add_argument('--nsources', help='For testing purposes, limit to a fixed number of sources in the .osm',default=False,action='store_true')
args = parser.parse_args()

factor = 2. * np.sqrt(2.*np.log(2.))   #conversion FWHM sigma
pi=np.pi
VELC = 299792458.0
bmax = 65e+3
D2R = pi / 180.0
#freq = 95e+6  ANNA: why is this not 150MHZ?

freqname1=args.freq1
freq1=float(freqname1)
freqname2=args.freq2
freq2=float(freqname2)

tel_reso = (VELC / freq1/1.e+6) / bmax
tel_reso_arcmin=tel_reso/D2R*3600.

def create_AGNs(t,num_sources,outfile):
    '''Takes a reduced TRECS AGN FITS table and turns it into .osm entries'''
    ##Grab the data
    #t = Table.read(fitsfile)
    ras = t['longitude'] 
    decs = t['latitude']
    Snu = t['I'+freqname1]/1000. # mJY-> Jy
    Sdnu = t['I'+freqname2]/1000. # mJY-> Jy
    Rs = t['Rs']
    size = t['size']
    pop=t['PopFlag']
    ##Calc Spectral indices
    SIs = np.log10(Snu / Sdnu) / np.log10(freq1 / freq2)

    if num_sources <0:
            num_sources=len(ras)
    ##I needed consistent results when testing before, but you can change or
    ##remove the seed here
    random.seed(9083745)
    ##Make some random position angles for the orientation of the AGN on the sky
    pa = random.uniform(low=0,high=2*pi,size=len(ras))

    ##Use the TRECS Rs and size params to work out the distance of the hotspots
    ##from the core

    dist_core = 0.5 * Rs * size   #used only for steep-spectrum AGNs 

    ##Set the major axis of gaussian for the lobe by measuring the distance of
    ##the hotspot to the edge of the source
    size_fwhm=size*factor/3.  #T-RECS size is LAS corresponding to 3 sigma of the Gaussian. 
    lobe_maj = 0.5 * size_fwhm #the size contains 2 lobes (used only for steep-spectrum)

    ##Make some random ratios to scale the minor axis of the lobe gaussian,
    ##gives some slight variety in the shape of the lobes
    lobe_ratio = random.uniform(low=0.3,high=1.0,size=len(ras))
    lobe_min = lobe_maj*lobe_ratio

    ##Set up parameters for the hot spots. Make them 10 times smaller than the
    ##lobes. Scale found by trial and error looking at output plots
    hot_maj = size / 10
    hot_min = lobe_ratio*size / 10

    ##Set up hot-spot locations relative to core
    ##size is in arcsec, ras in degrees
    ras_l_hot = -dist_core / 3600.0
    ras_r_hot = dist_core / 3600.0

    ras_l_lobe = -lobe_maj / 3600.0
    ras_r_lobe = lobe_maj / 3600.0

    
    ##rotate the hot spot locations about the core, using the random pa
    ##in line with dec so dec = 0

    dec_l_hot, dec_r_hot = 0, 0
    ras_1_hot = ras_l_hot*np.cos(-(pa-pi/2)) - dec_l_hot*np.sin(-(pa-pi/2)) + ras
    decs_1_hot = dec_l_hot*np.cos(-(pa-pi/2)) - ras_l_hot*np.sin(-(pa-pi/2)) + decs

    ras_2_hot = ras_r_hot*np.cos(-(pa-pi/2)) - dec_r_hot*np.sin(-(pa-pi/2)) + ras
    decs_2_hot = dec_r_hot*np.cos(-(pa-pi/2)) - ras_r_hot*np.sin(-(pa-pi/2)) + decs

    ##rotate the lobe locations about the core, using the random pa
    ##in line with dec so dec = 0
    dec_l_lobe, dec_r_lobe = 0, 0
    ras_1_lobe = ras_l_lobe*np.cos(-(pa-pi/2)) - dec_l_lobe*np.sin(-(pa-pi/2)) + ras
    decs_1_lobe = dec_l_lobe*np.cos(-(pa-pi/2)) - ras_l_lobe*np.sin(-(pa-pi/2)) + decs

    ras_2_lobe = ras_r_lobe*np.cos(-(pa-pi/2)) - dec_r_lobe*np.sin(-(pa-pi/2)) + ras
    decs_2_lobe = dec_r_lobe*np.cos(-(pa-pi/2)) - ras_r_lobe*np.sin(-(pa-pi/2)) + decs

    
    ##Divide up the flux of the source between the core, lobes, and hotspots
    ##Different for flat-spectrum and steep-spectrum as the steep-spectrum are seen perpendicular to the jet so much more diffuse emission
    sigma_flat=0.1
    mu_flat=0.75
    sigma_steep=0.005
    mu_steep=0.03
    core_frac = random.normal(0,1,len(ras)) #gauss2=randomn(seed,nrows)*0.1+0.75    ; IDL for core/total flux   for flat spectrum
    hot_frac=core_frac*0. 
    for ind in np.arange(len(ras)):
        if (pop[ind]==6): #steep-spectrum
            core_frac[ind]=core_frac[ind]*sigma_steep+mu_steep
        else: #flat_spectrum
            core_frac[ind]=core_frac[ind]*sigma_flat+mu_flat
       # hot_frac[ind] = random.uniform(low=0.01,high=1.-core_frac[ind]-0.1,size=len(ras))
        limit=1.-core_frac[ind]-0.8
        hot_frac[ind] = random.uniform(low=0.01,high=limit,size=1)


    lobe_frac = (1.0 - core_frac - hot_frac*2) / 2.0

 
    if args.plot_AGN:
        ##Find the order of the biggest galaxies, with largest first
        idx_order = np.argsort(size)[::-1]
        ##Plot the 9 largest galaxies
        fig = plt.figure(figsize=(10,10))
        for plot,ind in enumerate(idx_order[:9]):
            ax = fig.add_subplot(3,3,plot+1)

            num_pixels = 101
            plot_width = size[ind] / (1*3600.0)
            x_range = np.linspace(ras[ind] - plot_width,ras[ind] + plot_width,num_pixels)
            y_range = np.linspace(decs[ind] - plot_width,decs[ind] + plot_width,num_pixels)

            reso = x_range[1] - x_range[0]

            x_stddev = lobe_maj[ind] / (factor*3600.0)
            y_stddev = lobe_min[ind] / (factor*3600.0)

            lobe1_func = Gaussian2D(amplitude=Snu[ind]*lobe_frac[ind], x_mean=ras_1_lobe[ind], y_mean=decs_1_lobe[ind], x_stddev=x_stddev, y_stddev=y_stddev,theta=pi/2 + pa[ind])
            lobe2_func = Gaussian2D(amplitude=Snu[ind]*lobe_frac[ind], x_mean=ras_2_lobe[ind], y_mean=decs_2_lobe[ind], x_stddev=x_stddev, y_stddev=y_stddev,theta=pi/2 + pa[ind])

            x_stddev = hot_maj[ind] / (factor*3600.0)
            y_stddev = hot_min[ind] / (factor*3600.0)

            hot1_func = Gaussian2D(amplitude=Snu[ind]*hot_frac[ind], x_mean=ras_1_hot[ind], y_mean=decs_1_hot[ind], x_stddev=x_stddev, y_stddev=y_stddev,theta=pi/2 + pa[ind])
            hot2_func = Gaussian2D(amplitude=Snu[ind]*hot_frac[ind], x_mean=ras_2_hot[ind], y_mean=decs_2_hot[ind], x_stddev=x_stddev, y_stddev=y_stddev,theta=pi/2 + pa[ind])

            x_mesh, y_mesh = np.meshgrid(x_range,y_range)
            lobe1 = lobe1_func(x_mesh,y_mesh)
            lobe2 = lobe2_func(x_mesh,y_mesh)

            lobe1 *=  Snu[ind]*lobe_frac[ind] / lobe1.sum()
            lobe2 *=  Snu[ind]*lobe_frac[ind] / lobe2.sum()

            hot1 = hot1_func(x_mesh,y_mesh)
            hot2 = hot2_func(x_mesh,y_mesh)

            hot1 *=  Snu[ind]*hot_frac[ind] / hot1.sum()
            hot2 *=  Snu[ind]*hot_frac[ind] / hot2.sum()

            gal_array = hot1 + hot2 + lobe1 + lobe2
            gal_array[num_pixels // 2,num_pixels // 2] += Snu[ind]*core_frac[ind]

            im = ax.imshow(gal_array,origin='lower',extent=[x_range[0],x_range[-1],y_range[0],y_range[-1]],vmax=np.mean(gal_array)+2.*np.std(gal_array))
            
            if plot in [6,7,8]:
                ax.set_xlabel('RA (deg)')
            if plot in [0,3,6]:
                ax.set_ylabel('Dec (deg)')
        plt.tight_layout()
        fig.savefig('example_AGNs_in_osm.png',bbox_inches='tight')

    ##Write the source into the .osm file
    ##Go from brightest to dimmest
    idx_bright = np.argsort(Snu)[::-1]
    for ind in idx_bright[:num_sources]:
        if ((size[ind] < tel_reso_arcmin*60.) or (pop[ind]<6)):  #unresolved steep-spectrum sources or flat-spectrum (sources viewed along the jet axis)
            
            ##Point source for the core and a circular Gaussian for the lobe viewed side-on
             outfile.write('%.7f %.7f %.9f 0 0 0 150e+6 %.2f 0 0 0 0\n' %(ras[ind],decs[ind],Snu[ind]*core_frac[ind],SIs[ind]))
             outfile.write('%.7f %.7f %.9f 0 0 0 150e+6 %.2f 0 %.2f %.2f %.2f\n' %(ras[ind],decs[ind],Snu[ind]*(1.-core_frac[ind]),SIs[ind],size_fwhm[ind],size_fwhm[ind],pa[ind]*(180/pi)))  

        else:

            ##Write out core - assume is a point source
            ##Scale the flux according the settings above
            #point source for core, 2 gaussian lobes, two gaussian hot spots for the steep-spectrum resolved sources
            #print('AGN',size[ind],pa[ind]/D2R)
            outfile.write('%.7f %.7f %.9f 0 0 0 150e+6 %.2f 0 0 0 0\n' %(ras[ind],decs[ind],Snu[ind]*core_frac[ind],SIs[ind]))
            outfile.write('%.7f %.7f %.9f 0 0 0 150e+6 %.2f 0 %.2f %.2f %.2f\n' %(ras_1_lobe[ind],decs_1_lobe[ind],Snu[ind]*lobe_frac[ind],SIs[ind],lobe_maj[ind],lobe_min[ind],pa[ind]*(180/pi)))
            outfile.write('%.7f %.7f %.9f 0 0 0 150e+6 %.2f 0 %.2f %.2f %.2f\n' %(ras_2_lobe[ind],decs_2_lobe[ind],Snu[ind]*lobe_frac[ind],SIs[ind],lobe_maj[ind],lobe_min[ind],pa[ind]*(180/pi)))
            outfile.write('%.7f %.7f %.9f 0 0 0 150e+6 %.2f 0 %.2f %.2f %.2f\n' %(ras_1_hot[ind],decs_1_hot[ind],Snu[ind]*hot_frac[ind],SIs[ind],hot_maj[ind],hot_min[ind],pa[ind]*(180/pi)))
            outfile.write('%.7f %.7f %.9f 0 0 0 150e+6 %.2f 0 %.2f %.2f %.2f\n' %(ras_2_hot[ind],decs_2_hot[ind],Snu[ind]*hot_frac[ind],SIs[ind],hot_maj[ind],hot_min[ind],pa[ind]*(180/pi)))
    #
     #return num_sources + len(ras)
    return num_sources #+ num_sources

def create_SFGs(t,num_sources,outfile):
    '''Takes a reduced TRECS SFG FITS table and turns it into .osm entries'''
    ##Read in data
  
    
    #t = Table.read(fitsfile)
    ras = t['longitude'] 
    decs = t['latitude']
    Snu = t['I'+freqname1]/1000. #mJy -> Jy
    Sdnu = t['I'+freqname2]/1000. #mJy -> Jy
    e1 = t['e1']
    e2 = t['e2']
    size = t['size']
    size=size*np.sqrt(2.)  #conversion from exponential to gaussian FWHM. TODO: check if OSKAR wants sigma or FWHM
    ##Calc Spectral indices
    SIs = np.log10(Snu / Sdnu) / np.log10(freq1 / freq2)

    if num_sources <0:
            num_sources=len(ras)
    ##Calc PA and Bmaj/Bmin from ellipticity parameters e1,e2
    pa=size*0. 
    emod=np.sqrt(e1**2+e2**2)
    for ind in np.arange(len(ras)):
            theta=np.math.atan2(e2[ind],e1[ind])*180./pi
            pa[ind]=-theta/2. # change sign to get it to CW as AGNs. factor 2 is because rotations are of theta/2.
    q=1-emod #axis ratio
    maj_ax=np.sqrt(size**2/q) # major axis
    min_ax=q*maj_ax # minor axis
    
    if args.plot_SFG:
        ##Find the order of the biggest galaxies, with largest first
        idx_order = np.argsort(size)[::-1]
        ##Plot the 9 largest galaxies
        fig = plt.figure(figsize=(10,10))
        for plot,ind in enumerate(idx_order[:9]):
            ax = fig.add_subplot(3,3,plot+1)

            num_pixels = 101
            plot_width = size[ind] / (1.*3600.0)
            x_range = np.linspace(ras[ind] - plot_width,ras[ind] + plot_width,num_pixels)
            y_range = np.linspace(decs[ind] - plot_width,decs[ind] + plot_width,num_pixels)
            x_mesh, y_mesh = np.meshgrid(x_range,y_range)

            reso = x_range[1] - x_range[0]

            x_stddev = maj_ax[ind] / (factor*3600.0)
            y_stddev = min_ax[ind] / (factor*3600.0)

            sfg_func = Gaussian2D(amplitude=Snu[ind], x_mean=ras[ind], y_mean=decs[ind], x_stddev=x_stddev, y_stddev=y_stddev,theta=pi/2 + pa[ind]*(pi/180.0))  

            sfg = sfg_func(x_mesh,y_mesh)
            sfg *=  Snu[ind]*sfg / sfg.sum()

            im = ax.imshow(sfg,origin='lower',extent=[x_range[0],x_range[-1],y_range[0],y_range[-1]])
            
            if plot in [6,7,8]:
                ax.set_xlabel('RA (deg)')
            if plot in [0,3,6]:
                ax.set_ylabel('Dec (deg)')

        plt.tight_layout()
        fig.savefig('example_SFGs_in_osm.png',bbox_inches='tight')


    idx_bright = np.argsort(Snu)[::-1]
    for ind in idx_bright[:num_sources]:

        if size[ind] < 60.*tel_reso_arcmin:
            ##If less than 3 beams, just stick in a point source
             outfile.write('%.7f %.7f %.9f 0 0 0 150e+6 %.2f 0 0 0 0\n' %(ras[ind],decs[ind],Snu[ind],SIs[ind]))
        else:
            ##Then write a gaussian
            #print('sfg',size[ind],pa[ind])
            outfile.write('%.7f %.7f %.9f 0 0 0 150e+6 %.2f 0 %.2f %.2f %.2f\n' %(ras[ind],decs[ind],Snu[ind],SIs[ind],maj_ax[ind],min_ax[ind],pa[ind]))
    return num_sources #+ num_sources

num_sources = 0
outfile = open(args.osm_name,'w+')
print(args.nsources)
nsources=int(args.nsources)
if nsources >0:
    print('Limiting sky model to '+args.nsources+' sources')
 
    
with open(args.TRECS_tables) as f:files=f.readlines()

for file in files:
    print('file=',file)
    file=file.split('\n')
    t = Table.read(file[0])
    pop = t['PopFlag']

    if (max(pop) <=3.): #SFGs
        num_sources=nsources
        num_sources = create_SFGs(t,num_sources,outfile)
    else:               #AGNs
        num_sources=nsources
        num_sources = create_AGNs(t,num_sources,outfile)

    print('Sky model contains %d sources' %num_sources)
