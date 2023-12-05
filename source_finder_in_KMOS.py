#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 16:06:22 2023

@author: amartinez
"""

# Try to find sources in Kmos colapse image.
# Fitting continuum channels to an algotithms 
import sep

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
from matplotlib.widgets import Slider
import astroalign as aa
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
# %%

pruebas = '/Users/amartinez/Desktop/PhD/KMOS/practice/'
choped_ifus = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/ifu_alignment_ABC/'
fig, ax = plt.subplots(1,1,figsize =(10,10)) 
# im = fits.open(pruebas  + 'Kmos_image.fits')
ifu_sel = 17
half_ifu = 1
ima = fits.open(choped_ifus  + 'ifu%s_half%s_p105.fits'%(ifu_sel,half_ifu))
data = ima[1].data
cube_plot = data
v_min =np.nanmin(cube_plot)
v_max =np.nanmax(cube_plot)
im = ax.imshow(cube_plot, cmap='inferno',label='overlays',origin='lower',vmin=v_min,vmax = v_max,alpha =1)
fig.colorbar(im, ax=ax, orientation='vertical')

# Add a slider for vmin
ax_vmin = plt.axes([0.2, 0.1, 0.3, 0.03], facecolor='lightgoldenrodyellow')
vmin_slider = Slider(ax_vmin, 'vmin', v_min, v_max, valinit=v_min)

# Add a slider for vmax
ax_vmax = plt.axes([0.2, 0.05, 0.3, 0.03], facecolor='lightgoldenrodyellow')
vmax_slider = Slider(ax_vmax, 'vmax', v_min, v_max, valinit=v_max)


# %%
data = data.byteswap().newbyteorder()
bkg = sep.Background(data)

data_sub = data - bkg
# %%
obs = sep.extract(data_sub, 2, err=bkg.globalrms)

from matplotlib.patches import Ellipse

# plot background-subtracted image

# plot an ellipse for each object
# for i in range(len(objects)):
#     e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
#                 width=6*objects['a'][i],
#                 height=6*objects['b'][i],
#                 angle=objects['theta'][i] * 180. / np.pi)
#     e.set_facecolor('none')
#     e.set_edgecolor('lime')
#     ax.add_artist(e)
for i in range(len(obs)):
    ax.scatter(obs['x'][i],obs['y'][i], s= 200, marker = 'x', color ='lime')

def update(val):
    vmin = vmin_slider.val
    vmax = vmax_slider.val
    im.set_clim(vmin, vmax)
    plt.draw()

vmin_slider.on_changed(update)
vmax_slider.on_changed(update)

plt.draw()

# %%
coor = np.loadtxt(pruebas + 'text_coord6.txt')
x = obs['x']
y = obs['y']
# x = obs['y']
# y = obs['x']
xy = np.array([x,y]).T
wcs = WCS(ima[0].header)
coor_pix = wcs.wcs_world2pix(coor[:,0],coor[:,1],0)

fig, ax = plt.subplots(1,1)
ax.scatter(coor_pix[0],coor_pix[1],s=100)
ax.scatter(xy[:,0],xy[:,1],)
# sys.exit(137)
coor_pix = np.array([coor_pix[0],coor_pix[1]]).T
m, (source_list, target_list) = aa.find_transform(coor_pix, xy,max_control_points=500)
print('source_list',source_list)
print('target_list',target_list)

print("Translation: (x, y) = (%.2f, %.2f)"%(m.translation[0],m.translation[1]))
print("Rotation: %.3f degrees"%(m.rotation * 180.0 / np.pi))
print("Scale factor: %.4f"%(m.scale))

coor_pix_t = aa.matrix_transform(coor_pix, m.params)

ig, ax = plt.subplots(1,1)
ax.scatter(coor_pix_t[:,0],coor_pix_t[:,1],s=100)
ax.scatter(xy[:,0],xy[:,1])



