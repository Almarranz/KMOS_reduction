#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 18 17:13:17 2023

@author: alvaromartinez
"""

# Plot the spectrum from fits files for comparison

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys
import astropy.units as u
from astropy.utils.data import download_file
from astropy.io import fits  # We use fits to open the actual data file
from spectral_cube import SpectralCube
from matplotlib import rc
from matplotlib import rcParams
import glob
rcParams.update({'xtick.major.pad': '7.0'})
rcParams.update({'xtick.major.size': '7.5'})
rcParams.update({'xtick.major.width': '1.5'})
rcParams.update({'xtick.minor.pad': '7.0'})
rcParams.update({'xtick.minor.size': '3.5'})
rcParams.update({'xtick.minor.width': '1.0'})
rcParams.update({'ytick.major.pad': '7.0'})
rcParams.update({'ytick.major.size': '7.5'})
rcParams.update({'ytick.major.width': '1.5'})
rcParams.update({'ytick.minor.pad': '7.0'})
rcParams.update({'ytick.minor.size': '3.5'})
rcParams.update({'ytick.minor.width': '1.0'})
rcParams.update({'font.size': 20})
rcParams.update({'figure.figsize':(10,5)})
rcParams.update({

    "text.usetex": False,
    "font.family": "sans",
    "font.sans-serif": ["Palatino"]})
plt.rcParams["mathtext.fontset"] = 'dejavuserif'
rc('font',**{'family':'serif','serif':['Palatino']})
plt.rcParams.update({'figure.max_open_warning': 0})# 
# %%
pruebas = '/Users/amartinez/Desktop/PhD/KMOS/practice/'
spec_folder = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/regular_reduction/young_candidates/' 
ifu_sel = 6
d_spec, u_spec = 2.295,2.31
# d_spec, u_spec = 2.295,2.31 # Feautureless
with open(pruebas + 'yso_ifu%s.reg'%(ifu_sel),'w') as reg:
    reg.write('# Region file format: DS9 version 4.1'+"\n"+'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'+"\n"+'image'+'\n')
    reg.close
fig, ax = plt.subplots(1,1)

for st in range(13):
    
    name = glob.glob(spec_folder + 'ifu_%s/'%(ifu_sel) +'star*')[st]
   
    xy_in = name[name.find('star_')+5:name.find('star_')+5+7]
    f = fits.open(name)
    specdata = f[0].data 
    cab = f[0].header
    lam = np.array([cab['CRVAL1']+cab['CDELT1']*i for i in range(len(specdata))] )
    good = np.where((lam > d_spec) & (lam < u_spec))
    specdata, lam  = specdata[good], lam[good]
    ax.plot(lam,specdata + st*1*1e-16,  label ='%s,%s,SNR = %.2g'%(st,xy_in, np.mean(specdata)/np.std(specdata)))
    # ax.set_xlim(2.1,2.33)
    # ax.set_ylim(0,4e-17)
    print(np.mean(specdata))
    ax.legend(fontsize = 10)
    with open(pruebas + 'yso_ifu%s.reg'%(ifu_sel),'a') as reg:
        reg.write('# text(%s,%s) text={%s,SN = %.1f}\n'%(xy_in[0:3],xy_in[4:7],st,np.mean(specdata)/np.std(specdata)))
        reg.close
ax.set_title('IFU %s'%(ifu_sel))
sys.exit('61')
name = 'sky.fits'
print(name)
# fig, ax = plt.subplots(1,1)
f = fits.open(pruebas + name)
specdata_sky = f[0].data 
cab = f[0].header
lam = np.array([cab['CRVAL1']+cab['CDELT1']*i for i in range(len(specdata))] )
good = np.where((lam > d_spec) & (lam < u_spec))
specdata_sky, lam  = specdata_sky[good], lam[good]
ax.plot(lam,specdata_sky,  label ='sky', alpha = 0.5, ls ='dashed')
# ax.legend()
# # %% 
# fig, ax = plt.subplots(1,1)
# ax.scatter(np.arange(0,len(specdata)),specdata/specdata_sky)


# %%


































