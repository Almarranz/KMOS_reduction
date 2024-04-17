#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 13:26:13 2023

@author: amartinez
"""

# Fitting continuum channels to an algotithms 
import os
import numpy as np
import sys
import glob
from subprocess import call
import astropy.units as u
from astropy.utils.data import download_file
from astropy.io import fits  # We use fits to open the actual data file
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from matplotlib.colors import LogNorm
from matplotlib.tri import Triangulation
import regions
from matplotlib.patches import Rectangle
from matplotlib.patches import Circle
import IPython
import tkinter as tk
from tkinter import simpledialog
from astropy.stats import sigma_clip
from numpy import mean
from time import time
# %%plotting pa    metres
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
rcParams.update({'figure.figsize':(30,30)})
rcParams.update({
    "text.usetex": False,
    "font.family": "sans",
    "font.sans-serif": ["Palatino"]})
plt.rcParams["mathtext.fontset"] = 'dejavuserif'
rc('font',**{'family':'serif','serif':['Palatino']})
plt.rcParams.update({'figure.max_open_warning': 0})

# Enable automatic plotting mode
# IPython.get_ipython().run_line_magic('matplotlib', 'auto')
IPython.get_ipython().run_line_magic('matplotlib', 'inline')


reduction = 'ABC'
pruebas = '/Users/amartinez/Desktop/PhD/KMOS/practice/'
aling = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/ifu_alignment_%s/'%('ABC')
log_5 = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p105_%s/'%('ABC')
esorex_cube_5 = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p105_%s/COMBINE_SKY_TWEAK_mapping.fits'%(reduction)
esorex_ima_5= '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p105_%s/comoving_group_mosaic_K_Half1_COMBINED_IMAGE_mapping.fits'%(reduction)
log_5 = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p105_%s/'%('ABC')

# ifu_sel = 7
ifu_sel = 'all'
half_ifu = 2
if ifu_sel == 'all':    
    aligned_cube = esorex_cube_5
    lamb, ind, chunk = np.loadtxt(pruebas + 'continuum_all.txt'  ,unpack = True)

else:
    aligned_cube = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/%s_reduction/cubes/cube_ifu%s_half%s_mean.fits'%(reduction,ifu_sel,half_ifu)
    lamb, ind, chunk = np.loadtxt(pruebas + 'continuum_ifu%s_half%s.txt'%(ifu_sel, half_ifu),unpack = True)

cube = fits.open(aligned_cube)
flux = cube[1].data
cab = cube[1].header
all_lamb = np.array([cab['CRVAL3']+cab['CDELT3']*i for i in range(cab['NAXIS3'])])



cont_cube = np.empty(flux.shape)
ti = time()
for i in range(flux.shape[1]):  
    for j in range(flux.shape[2]):
        z = np.polyfit(lamb, flux[ind.astype(int),i,j], 4, rcond=None, full=False, w=None, cov=False)
        p = np.poly1d(z)
        cont_flux = p(all_lamb)
        cont_cube[:,i,j] = cont_flux
tf = time()
print('Fitting took %.0f seconds'%(tf-ti))
        
        
hdu_cont = fits.PrimaryHDU(data = cont_cube, header = cab )
cube_no_cont = flux -  cont_cube
hdu_no_cont = fits.PrimaryHDU(cube_no_cont, header = cab)

ti = time()
if ifu_sel == 'all':
    hdu_cont.writeto(pruebas + 'cont_cube_all.fits',overwrite= True)
    hdu_no_cont.writeto(pruebas + 'cube_no_cont_all.fits',overwrite= True)
else:
    hdu_cont.writeto(pruebas + 'cont_cube_ifu%s_half%s.fits'%(ifu_sel, half_ifu),overwrite= True)
    hdu_no_cont.writeto(pruebas + 'cube_no_cont_ifu%s_half%s.fits'%(ifu_sel, half_ifu),overwrite= True)
tf = time()
print('Save the fits file took %.0f seconds'%(tf-ti))

# %%

cont_ima = np.mean(cont_cube, axis = 0)
# %%
fig, ax = plt.subplots(1,1)
ax.imshow(cont_ima, origin = 'lower')
hdu_cont_imag = fits.PrimaryHDU(data = cont_ima, header = cab )
hdu_cont_imag.writeto(pruebas + 'Kmos_image.fits')















# l_a, l_b = 0,3000
# fig, ax = plt.subplots(1,1, figsize = (10,10))
# ax.plot(all_lamb, flux)
# ax.plot(lamb, flux[ind.astype(int)],'--',lw =2)
# HeI = 2.058
# COI =  2.29322
# COII = 2.32246
# Brg = 2.165
# He = 2.12
# HeII = 2.189
# # H2 = 2.12
# l_names = ['HeI', 'COI', 'COII','Br$\gamma$', 'He', 'HeII']
# lines = [HeI, COI, COII,Brg, He, HeII]
# for l in lines:
#     ax.axvline(l, color = 'grey', alpha = 0.5, ls = 'dashed')
# ax2 = ax.twiny()
# ax2.set_xticks(lines)
# tl = l_names
# ax2.set_xticklabels(tl)
# ax2.plot(all_lamb, flux, alpha = 0)
# # ax.set_xlim(0,2000)

# z = np.polyfit(lamb, flux[ind.astype(int)], 
#                   4, rcond=None, full=False, w=None, cov=False)
# p = np.poly1d(z)

# x_lamb = np.arange(min(all_lamb),max(all_lamb),
#                    (max(all_lamb)-min(all_lamb))/3072)
# ax.plot(x_lamb,p(x_lamb))
# # %%
# ima = fits.open(aligned_cube)
# mapa = WCS(ima[1].header).celestial

# h0 = ima[0].header
# h1 = ima[1].header

# cube_cont = ima[1].data-p(x_lamb)[:, np.newaxis, np.newaxis]
# fits.writeto(pruebas + 'cube_cont_star_%s_%s.fits'%(star[0],star[1]),cube_cont, overwrite = True)
# v_min =0.001e-16
# v_max = 0.05e-16
# # fig2, ax2 = plt.subplots(1,1,figsize = (10,10))
# # im = ax2.imshow(cube_cont, cmap='inferno',label='overlays',origin='lower',vmin=v_min,vmax = v_max,alpha =1)





