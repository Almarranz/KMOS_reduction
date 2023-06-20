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

import astropy.units as u
from astropy.utils.data import download_file
from astropy.io import fits  # We use fits to open the actual data file
from spectral_cube import SpectralCube
from matplotlib import rc
from matplotlib import rcParams
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
pruebas = '/Users/alvaromartinez/Desktop/Phd/KMOS/pruebas/'
d_spec, u_spec = 2.10,2.4
fig, ax = plt.subplots(1,1)
for st in range(1,3):
    
    name = 'spec_%s.fits'%(st)
    print(name)
    f = fits.open(pruebas + name)
    specdata = f[0].data 
    
    cab = f[0].header
    lam = np.array([cab['CRVAL1']+cab['CDELT1']*i for i in range(len(specdata))] )
    good = np.where((lam > d_spec) & (lam < u_spec))
    specdata, lam  = specdata[good], lam[good]
    ax.plot(lam,specdata,  label ='%s')







































