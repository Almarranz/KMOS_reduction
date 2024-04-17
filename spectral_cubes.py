#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 12:19:06 2023

@author: amartinez
"""
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.utils import data
from spectral_cube import SpectralCube
from astropy.io import fits  
from spectral_cube import Projection
from spectral_cube import OneDSpectrum
import matplotlib.colors as mcolors
import matplotlib.colors as mcolors
from astropy.wcs import WCS

# %%

folder = '/Users/amartinez/Desktop/PhD/KMOS/practice/'
folder = '/Users/amartinez/Desktop/PhD/KMOS/cubes_sky_tweak/p105/'


# cube = SpectralCube.read(folder + 'example_cube.fits') 
cube = SpectralCube.read(folder + 'COMBINE_SKY_TWEAK_mapping.fits', hdu =1) 

print(cube)
# %%
hdul = fits.open(folder + 'COMBINE_SKY_TWEAK_mapping.fits')
# projection = Projection.from_hdu(hdul[2].data)
# projection = OneDSpectrum.from_hdu(hdul[1])
# %%
image = cube[1500]
# %%
slice_unmasked = cube.unmasked_data[1500,:,:]  
# %%
lambs = cube.spectral_axis  
# %%
lam, dec, ra = cube.world[:,20,21]  
# %%
# fig, ax = plt.subplots(1,1)
# ax.scatter(ra,dec, c= lam, norm = mcolors.LogNorm())

# %%
moment_0 = cube.moment(order=0)  
# %%
fig, ax = plt.subplots(1,1)
ax.plot(moment_0)
# %%
sigma_map = cube.linewidth_sigma()  

# %%

clus_obs = np.loadtxt()
hdul = fits.open(folder + 'COMBINE_SKY_TWEAK_mapping.fits')
map_ = WCS(hdul[1])
cor =[266.3864974, -28.9379862,0]
coor = map_.wcs_world2pix([cor],1)

# %%
fig, ax = plt.subplots(1,1)
# ax.scatter( 400,   400,zorder=3,s=50,color ='red')
plt.imshow(moment_0.value,vmin=-0.8e-20,vmax=0.8e-17,origin='lower',cmap = 'Greys')
ax.scatter( 400,   400,zorder=3,s=50,color ='red')
# plt.imshow(moment_0.value,norm = mcolors.LogNorm(vmin=33,vmax=56),origin='lower')



# %%













