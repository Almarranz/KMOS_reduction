#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 17:53:54 2023

@author: amartinez
"""

# Here we are going to extract spectra from the combines images.

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

data = '/Users/amartinez/Desktop/PhD/KMOS/practice/for_combination/'
pruebas = '/Users/amartinez/Desktop/PhD/KMOS/practice/'

hi_data = fits.open(data +'COMBINE_SKY_TWEAK_1.fits')
cube = SpectralCube.read(hi_data, hdu = 1)  # Initiate a SpectralCube
hi_data.close()  # Close the FITS file - we already read it in and don't need it anymore!


# plt.imshow(cube[1500,:,:].value, origin='lower', norm=LogNorm())

cube[1500, :, :].quicklook()
# cube[1500, 0:, 0:].quicklook()  
# %%
stars = [1,2]
stars_nB = [1,2]
CO_nB =[1]

tl=[2.16,2.295,2.325,2.352]
cns = [0,1,2.2,0.5,5]
color_ban = ['b','royalblue']
# color_noban = ['red','tomato']
color_noban = ['green','darkgreen']

fig, ax = plt.subplots(1,1, figsize = (10,10))
i =0
d_spec, u_spec = 2.10,2.4
for st in stars:
    
    name = 'ban%s_noCO.fits'%(st)
    f = fits.open(pruebas + name)
    specdata = f[0].data 
    
    cab = f[0].header
    lam = np.array([cab['CRVAL1']+cab['CDELT1']*i for i in range(len(specdata))] )
    good = np.where((lam > d_spec) & (lam < u_spec))
    specdata, lam  = specdata[good], lam[good]
    ax.plot(lam,specdata+1e-16*(cns[i]), color =color_ban[st-1], label ='Ban')
    i += 1


for st in stars_nB:
    name = 'noban%s_noCO.fits'%(st)
    f = fits.open(pruebas + name)
    specdata = f[0].data 
    cab = f[0].header
    lam = np.array([cab['CRVAL1']+cab['CDELT1']*i for i in range(len(specdata))] )
    good = np.where((lam > d_spec) & (lam < u_spec))
    specdata, lam  = specdata[good], lam[good]
    ax.plot(lam,specdata+1e-16*(cns[i]),color =color_noban[st-1],label ='No Ban')
    i += 1
for st in CO_nB:
    name = 'noban%s_CO.fits'%(st)
    f = fits.open(pruebas + name)
    specdata = f[0].data 
    cab = f[0].header
    lam = np.array([cab['CRVAL1']+cab['CDELT1']*i for i in range(len(specdata))] )
    good = np.where((lam > d_spec) & (lam < u_spec))
    specdata, lam  = specdata[good], lam[good]
    ax.plot(lam,specdata+1e-16*(cns[i]),color ='k', alpha = 0.5, label = 'CO line')
    i += 1

name = 'pachen_alpha.fits'
f = fits.open(pruebas + name)
specdata = f[0].data 
cab = f[0].header
lam = np.array([cab['CRVAL1']+cab['CDELT1']*i for i in range(len(specdata))] )
good = np.where((lam > d_spec) & (lam < u_spec))
specdata, lam  = specdata[good], lam[good]
ax.plot(lam,specdata+1e-16,color ='violet', alpha = 1, label = 'Ph alpha$')


ax.legend(fontsize = 10, loc =2)
ax.axvline(2.295, linewidth =3, alpha = 0.2, color ='green')
ax.axvline(2.16, linewidth =3, alpha = 0.2, color ='grey')
ax.axvline(2.325, linewidth =3, alpha = 0.2, color ='grey')
ax.axvline(2.352, linewidth =3, alpha = 0.2, color ='grey')

ax.set_xlim(d_spec,u_spec)
# ax.set_ylim(0,0.8e-15)
ticks = list(ax.get_xticks())
# del ticks[3]
ax.set_xticks(tl)
# tl=[2.16,2.295,2.325,2.352]
tl[1]='2.3\nCO(0,1)'
tl[2]='2.33\n\nCO(3,1)'
tl[3]='2.35\nCO(4,2)'

tl[0]='2.16\nBr$\gamma$'
ax.set_xticklabels(tl)

ax.set_xlabel('$\lambda[\mu m]$')
# ax.set_ylabel('$Flux\,[erg/(s \cdot cm^2 \cdot \AA)]$')
ax.set_yticks([])
ax.set_ylabel('Flux + constant')

ax.xaxis.set_label_coords(0.5, -0.02)

# %%
# We are going to try and combined differents parts from different spectra in a 
# single one

# Windows used in the esoReflex pipeline

A = [1.975,1.987,1.993,2.010,2.041,2.060]
B = [2.269,2.291,2.308,2.335,2.360,2.379]
C = [2.416,2.440,2.445,2.475] 

ABC = [A,B,C]
name = 'ban1_noCO.fits'
f = fits.open(pruebas + name)
specdata = f[0].data 

cab = f[0].header
lam = np.array([cab['CRVAL1']+cab['CDELT1']*j for j in range(len(specdata))] )
# %


diff_Aup = np.absolute(lam - A[-1])
A_up = np.argmin(diff_Aup)

lam_A = np.array([cab['CRVAL1']+cab['CDELT1']*j for j in range(A_up)] )
spec_A =specdata[0:A_up]




# %
good = np.where((lam > d_spec) & (lam < u_spec))
# specdata, lam  = specdata[good], lam[good]
fig, ax = plt.subplots(1,1, figsize=(10,10))
ax.plot(lam,specdata+1e-16, label ='Ban')
ax.plot(lam_A,spec_A+1e-16, label ='Ban')


    



















