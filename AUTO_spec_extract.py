#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 26 11:30:04 2023

@author: amartinez
"""

# It re-extracts spetra from Kmos cubes using the pixels coordinates previously
# extracted by hand

# Extract the sky and the stars sky-subtracted
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
import shutil
from matplotlib.widgets import Slider
import re
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
# %%
# Enable automatic plotting mode
IPython.get_ipython().run_line_magic('matplotlib', 'auto')
# IPython.get_ipython().run_line_magic('matplotlib', 'inline')


reduction = 'ABC'
pruebas = '/Users/amartinez/Desktop/PhD/KMOS/practice/'
aling = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/ifu_alignment_%s/'%('ABC')
log_5 = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p105_%s/'%('ABC')
esorex_cube_5 = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p105_%s/COMBINE_SKY_TWEAK_mapping.fits'%(reduction)
esorex_ima_5= '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p105_%s/comoving_group_mosaic_K_Half1_COMBINED_IMAGE_mapping.fits'%(reduction)
log_5 = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p105_%s/'%('ABC')



# %%
half_ifu = 2
ifu_sel = 7
spec_folder = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/%s_reduction/cluster_spectra/ifu_%s/half_%s/'%(reduction, ifu_sel, half_ifu)
spec_young = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/%s_reduction/young_candidates/ifu_%s/half_%s/'%(reduction, ifu_sel, half_ifu)
aligned_cube = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/%s_reduction/cubes/cube_ifu%s_half%s_mean.fits'%(reduction,ifu_sel,half_ifu)


w1 = 2
w2 = 2

ima = fits.open(aligned_cube)
cube = ima[1].data
cube_plot = np.nanmedian(cube, axis = 0)

wide = 5.3#TODO FWHM of the standard star in the kmos pipeline
w = int(np.rint(wide/2))

fig, ax = plt.subplots(1,1,figsize =(10,10)) 
for i in range(len((glob.glob(spec_folder +'/spec*')))):
    name = os.path.basename(glob.glob(spec_folder +'/spec*')[i])
    coord = re.findall(r'\d+', name)
    x,y = np.array(coord[2]).astype(int), np.array(coord[3]).astype(int)
    print(x,y)
    # xx,yy = np.indices(cube_plot.shape)
    # distances = np.sqrt((xx - x)**2 + (yy - y)**2)
    # mask = np.where(distances <=w)    
    # spec = cube[:,mask[0],mask[1]]
    # spec_mean = np.mean(spec, axis = 1)
    # print(spec_mean.shape, np.nanmean(spec_mean))
    
    
    # mask_around = np.where((distances>w*2)&(distances<w*2 + w2))
    # circle = np.where(mask_around[1])
    # m = tuple(array[circle] for array in mask_around)
    # mask_rand = m
    # ring = cube[:, mask_rand[0], mask_rand[1]]
    # ring_mean = np.mean(ring, axis =1)
    # spec_mean_sky = spec_mean-ring_mean
    # print(np.nanmean(spec_mean_sky))
    
    ax.scatter(y,x, s = 100, color ='lime', marker = 'x')


ax.set_title('AUTO IFU %s, Half %s '%(ifu_sel, half_ifu))
v_min =np.nanmin(cube_plot)
v_max =np.nanmax(cube_plot)
im = ax.imshow(cube_plot, cmap='inferno',label='overlays',origin='lower',vmin=v_min,vmax = v_max/5,alpha =1)
fig.colorbar(im, ax=ax, orientation='vertical')
ax_vmin = plt.axes([0.2, 0.1, 0.3, 0.03], facecolor='lightgoldenrodyellow')
vmin_slider = Slider(ax_vmin, 'vmin', v_min, v_max, valinit=v_min)

# Add a slider for vmax
ax_vmax = plt.axes([0.2, 0.05, 0.3, 0.03], facecolor='lightgoldenrodyellow')
vmax_slider = Slider(ax_vmax, 'vmax', v_min, v_max, valinit=v_max)

def update(val):
    vmin = vmin_slider.val
    vmax = vmax_slider.val
    im.set_clim(vmin, vmax)
    plt.draw()

vmin_slider.on_changed(update)
vmax_slider.on_changed(update)

plt.draw()


























