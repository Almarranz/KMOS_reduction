#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 12:26:51 2023

@author: amartinez
"""
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

# We are going to make a Br_g emision map
reduction = 'ABC'
pruebas = '/Users/amartinez/Desktop/PhD/KMOS/practice/'
aling = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/ifu_alignment_%s/'%('ABC')
log_5 = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p105_%s/'%('ABC')
esorex_cube_5 = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p105_%s/COMBINE_SKY_TWEAK_mapping.fits'%(reduction)

ima = fits.open(esorex_cube_5)
h0 = ima[0].header
h1 = ima[1].header
temp5 = np.zeros((ima[1].data.shape[1],ima[1].data.shape[2]))
# %%
data = ima[1].data
brg = np.mean(data[1286:1288,:,:], axis =0) 
brg_mean = np.nanmean((data - brg),axis = 0)
# %%
hdul = fits.HDUList()
hdul.append(fits.PrimaryHDU(header = h0))
hdul.append(fits.ImageHDU(brg_mean*-1,h1, name ='Brg emission'))
hdul.writeto(pruebas  + 'brg_mean.fits', overwrite = True)


# %%
fig, ax = plt.subplots(1, 1, figsize=(12, 4))
data = brg_mean*-1
# Plot the original data
v_min =4.68248e-20
v_max = 7.96049e-18
im = ax.imshow(data, cmap='viridis', origin='lower',
               norm = LogNorm(vmin=v_min,vmax = v_max))


cbar = fig.colorbar(im, orientation='vertical')
# %%

data = brg_mean*-1
# Plot the original data

# ax.contour(data, levels =niveles, color = 'k')  
# Replace NaN or infinite values with a fill value
fill_value = 0.0  # Replace with an appropriate fill value
data = np.nan_to_num(data, nan=fill_value, posinf=fill_value, neginf=fill_value)
# %%
fig, ax = plt.subplots(1, 1, figsize=(12, 4))
v_min =-0.18248e-17
v_max = 0.6049e-17
niveles = np.arange(-.0e-17,0.6e-17, 0.2e-17)
y, x = np.indices(data.shape)
x = x.flatten()
y = y.flatten()
z = data.flatten()
triangulation = Triangulation(x, y)
im = ax.imshow(data, cmap='Greys', origin='lower',vmin=v_min,vmax = v_max)
axcb =ax.tricontour(triangulation, data.ravel(), levels=niveles, colors='white',vmin=v_min, vmax=v_max)

# ax.tricontour(x,y,z, levels=niveles, colors='white')

fig.colorbar(axcb,ax =ax,shrink=0.7)
# cbar = fig.colorbar(im, orientation='vertical')

# %%
data = brg_mean*-1
import numpy as np
from scipy.interpolate import griddata

# Create example data (replace this with your actual data)
# Assuming your data is in the variable 'data'
# 'x' and 'y' should represent the coordinates where 'data' is defined
data= data[:,500:]
x, y = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[0]))

# Perform the interpolation
# Assuming 'data' contains your array with NaN values
# 'method' can be 'linear', 'cubic', or other interpolation methods
interpolated_data = griddata((x[~np.isnan(data)], y[~np.isnan(data)]), data[~np.isnan(data)],
                            (x, y), method='cubic')

# 'interpolated_data' now contains the interpolated values

# You can visualize the interpolated data using a heatmap or other methods
fig, ax = plt.subplots(1, 1, figsize=(12, 4))
# Plot the original data
v_min =4.68248e-20
v_max = 7.96049e-18
im = ax.imshow(data, cmap='viridis', origin='lower',
               norm = LogNorm(vmin=v_min,vmax = v_max))
cbar = fig.colorbar(im, orientation='vertical')

# %%
import numpy as np
import matplotlib.pyplot as plt

# Create example data (replace this with your actual data)
# Assuming your data is in the variable 'data'
  # Replace with your data

# Create a grid
x = np.arange(650)
y = np.arange(433)

# Plot the data as a heatmap
plt.figure(figsize=(8, 6))  # Adjust the figure size as needed
plt.imshow(data, extent=(x.min(), x.max(), y.min(), y.max()), origin='lower', cmap='viridis')
plt.colorbar(label='Value')
plt.title('Heatmap of Data')
plt.xlabel('X-Axis')
plt.ylabel('Y-Axis')
plt.show()

# %%
import regions


