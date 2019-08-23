'''Jack Line, CIRA, 23/08/19'''
from astropy.io import fits
from astropy.table import Table,Column
import argparse
from numpy import *
from copy import deepcopy

VELC = 299792458.0
D2R = pi / 180.0

def condon_fixbase(bmax=None,freq=None):
    '''Estimates the confusion noise at a given freq (Hz) for a maximum
    baseline length (m). Returns the value in milli Jy per beam
    Equation 27 in Condon et al 2012'''
    eightarcsec = (8. / 3600.0) * (pi / 180.)
    reso = (VELC / freq) / bmax
    nano2milli = 10e-3 ##Want mJy, orig equation returns in nJy
    noise = nano2milli * 1.2 * (freq / 3.02e+9)**-0.7 * (reso / eightarcsec) ** (10.0/3.0)
    return noise

parser = argparse.ArgumentParser()
parser.add_argument('--fitsfile', help='TRECS FITS table file to shift')
parser.add_argument('--max_baseline', help='Maximum baseline for confusion noise flux cutoff (m)', default=65e+3)
parser.add_argument('--flux_cut', help='Alternatively, instead of confusion noise, just set a hard flux cutoff (mJy)', default=False)


args = parser.parse_args()

bmax = float(args.max_baseline)

if args.flux_cut:
    flux_cut = float(args.flux_cut)
else:
    flux_cut = condon_fixbase(bmax=bmax,freq=150e+6)

t = Table.read(args.fitsfile)

##This should succeed if the file contains AGNs
try:
    ##Grab relevant data
    ras_i = t['longitude']
    decs_i = t['latitude']
    fluxes_150 = t['I150']
    fluxes_160 = t['I160']
    sizes = t['size']
    Rs = t['Rs']

    ##Shift everthing from dec 0 centre to dec -26.7 centre
    decs_f = array(decs_i) - 26.7

    ##Need to rescale the RAs depending on what declination they are now at
    ras_m = deepcopy(ras_i)
    go_neg = where(ras_m > 180.0)
    ras_m[go_neg] -= 360.0
    ras_f = (array(ras_m)) * (cos(array(decs_i)*D2R) / cos(array(decs_f)*D2R))
    ras_f[where(ras_f < 0)] += 360.0

    ##Find where fluxes are above the cut
    cut = where(fluxes_150 > flux_cut)

    ##Create new FITS table
    out = Table()

    t_ras = Column(data=ras_f[cut],name='RA',unit='deg')
    t_decs = Column(data=decs_f[cut],name='DEC',unit='deg')
    t_ras_old = Column(data=ras_i[cut],name='RA_orig',unit='deg')
    t_decs_old = Column(data=decs_i[cut],name='DEC_orig',unit='deg')

    t_150s = Column(data=fluxes_150[cut] / 1000.0,name='S150',unit='Jy')
    t_160s = Column(data=fluxes_160[cut] / 1000.0,name='S160',unit='Jy')

    t_sizes = Column(data=sizes[cut],name='size',unit='arcsec')
    t_RSs = Column(data=Rs[cut],name='Rs')

    out.add_columns([t_ras,t_decs,t_ras_old,t_decs_old,t_150s,t_160s,t_sizes,t_RSs])

##If not, it should work with SFG sources
except:
    ##Grab relevant data
    ras_i = t['longitude']
    decs_i = t['latitude']
    fluxes_150 = t['I150']
    fluxes_160 = t['I160']
    sizes = t['size']
    e1s = t['e1']
    e2s = t['e2']

    ##Shift everthing from dec 0 centre to dec -26.7 centre
    decs_f = array(decs_i) - 26.7

    ##Need to rescale the RAs depending on what declination they are now at
    ras_m = deepcopy(ras_i)
    go_neg = where(ras_m > 180.0)
    ras_m[go_neg] -= 360.0
    ras_f = (array(ras_m)) * (cos(array(decs_i)*D2R) / cos(array(decs_f)*D2R))
    ras_f[where(ras_f < 0)] += 360.0

    ##Find where fluxes are above the cut
    cut = where(fluxes_150 > flux_cut)

    ##Create new FITS table
    out = Table()

    t_ras = Column(data=ras_f[cut],name='RA',unit='deg')
    t_decs = Column(data=decs_f[cut],name='DEC',unit='deg')
    t_ras_old = Column(data=ras_i[cut],name='RA_orig',unit='deg')
    t_decs_old = Column(data=decs_i[cut],name='DEC_orig',unit='deg')

    t_150s = Column(data=fluxes_150[cut] / 1000.0,name='S150',unit='Jy')
    t_160s = Column(data=fluxes_160[cut] / 1000.0,name='S160',unit='Jy')

    t_sizes = Column(data=sizes[cut],name='size',unit='arcsec')
    t_e1s = Column(data=e1s[cut],name='e1')
    t_e2s = Column(data=e1s[cut],name='e2')

    out.add_columns([t_ras,t_decs,t_ras_old,t_decs_old,t_150s,t_160s,t_sizes,t_e1s,t_e2s])

out.write('%s_reduced.fits' %args.fitsfile.split('.')[0],overwrite=True)
