#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 20:17:33 2023

@author: amartinez
"""

# Here we are going to analized the proper motions and positions of the stars
# we have spectra for
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
from scipy.spatial import distance
import re
import time
import random
from astropy.stats import sigma_clip
from matplotlib.ticker import FormatStrFormatter
from regions import Regions
from astropy.utils.data import get_pkg_data_filename
import pandas as pd
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
IPython.get_ipython().run_line_magic('matplotlib', 'auto')
# IPython.get_ipython().run_line_magic('matplotlib', 'inline')
# %%

ifu_sel = np.arange(1,25)
no_data = np.array([2,4,16,14,20,22,24,18])
ifu_sel = np.delete(ifu_sel,no_data-1)
half_ifu = [1,2]

reduction = 'ABC'
pruebas = '/Users/amartinez/Desktop/PhD/KMOS/practice/'
choped_ifus = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/ifu_alignment_ABC/'
gns_ls = '/Users/amartinez/Desktop/PhD/KMOS/GNS_lists/'
spec_folder = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/%s_reduction/cluster_spectra/ifu_%s/half_%s/'%(reduction, ifu_sel, half_ifu)



gns_young_all = np.empty((0,27))
gns_all= np.empty((0,25))
ifu_half = np.empty((0,2))
ifu_half_all = np.empty((0,2))
for ifu in ifu_sel:
    for half in half_ifu:
        spec_folder = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/%s_reduction/cluster_spectra/ifu_%s/half_%s/'%(reduction, ifu, half)
        spec_young = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/%s_reduction/young_candidates/ifu_%s/half_%s/'%(reduction, ifu, half)
        
        if os.path.isfile(spec_young + 'gns_lib_in_young_kmos_ifu%s_half%s.txt'%(ifu, half)):
            kmos_xy = np.loadtxt(spec_young + 'xy_coord_young_kmos_ifu%s_half%s.txt'%(ifu, half),)
            gns_young = np.loadtxt(spec_young + 'gns_lib_in_young_kmos_ifu%s_half%s.txt'%(ifu, half))
            if len(gns_young.shape) == 1 or len(kmos_xy.shape) ==1:
                gns_young = np.reshape(gns_young,(1,25))
                kmos_xy = np.reshape(kmos_xy,(1,2))
            # gns_young_all = np.append(gns_young_all,gns_young, axis = 0) 
            gns_young_all = np.append(gns_young_all, np.c_[gns_young,np.repeat(np.array([[ifu, half]]),len(gns_young),axis = 0)], axis = 0) 
            
            ifu_half = np.append(ifu_half,np.array([ifu, half]))
            print(ifu, half,kmos_xy)
        if os.path.isfile(spec_folder + 'gns_lib_in_kmos_ifu%s_half%s.txt'%(ifu, half)):
            # 0RA_gns 	1DE_gns 	2Jmag 	3Hmag 	4Ksmag 	5ra 	6Dec 	7x_c 	8y_c 	9mua 10dmua 	11mud 	12dmud 	13time 	14n1	15n2	16ID 	17mul 	18mub 	19dmul 	20dmub 	21m139 	
            gns_data = np.loadtxt(spec_folder + 'gns_lib_in_kmos_ifu%s_half%s.txt'%(ifu, half))
            gns_all = np.append(gns_all,gns_data,axis=0) 
            ifu_half_all = np.append(ifu_half_all, np.repeat(np.array([[ifu, half]]),len(gns_data),axis = 0), axis = 0)                           
        if  os.path.isfile(spec_folder + 'gns_lib_in_kmos_ifu%s_half%s.txt'%(ifu, half)) == False:
            print('No proper motions in ifu%s half%s'%(ifu, half))
gns_all = np.c_[gns_all, ifu_half_all]
# sys.exit()
# 0RA_gns 	1DE_gns 	2Jmag 	3Hmag 	4Ksmag 	5ra 	6Dec 	7x_c 	8y_c 	9mua 10dmua 	11mud 	12dmud 	13time 	14n1	15n2	16ID 	17mul 	18mub 	19dmul 	20dmub 	21m139 	22Separation
# %%
# This is a good set of colors
# ['#20B065', '#FD5EF9', '#41DE1E', '#E0F774', '#FBC379', '#AAF02C', '#68165F', '#EB103C', '#79452B', '#6C63E7', '#C16BA8', '#E3651E', '#855622', '#E2BFDE', '#E30EAC', '#FC2EC8', '#ED20E1', '#7E4A78', '#C08BEC', '#E6D0B5', '#AC52F5', '#B78ACE', '#6670CF', '#F7DBB5']

colorines = []
alp = 0.2#!!!
# bg = np.where(((gns_young_all[:,3]-gns_young_all[:,4])> 1.3)&
#               ((gns_young_all[:,3]-gns_young_all[:,4])<1.9))
bg = np.where(((gns_young_all[:,3]-gns_young_all[:,4])> 1.3))
              
gns_young_all = gns_young_all[bg]
gns_young_all = gns_young_all[gns_young_all[:,4].argsort()]


mag_lim = 14.5#!!!
bri = np.where(gns_young_all[:,4] < mag_lim)
young_good = gns_young_all[bri]
s_brig_ra = np.std(young_good[:,9])
s_brig_dec = np.std(young_good[:,11])
s_brig = np.std(np.sqrt(young_good[:,9]**2 + young_good[:,11]**2))
brig_color, s_brig_color = np.nanmean(young_good[:,3] - young_good[:,4]), np.nanstd(young_good[:,3] - young_good[:,4])

bg_all = np.where((gns_all[:,3]-gns_all[:,4])>1.3)
gns_all = gns_all[bg_all]


all_pm_sig=sigma_clip(np.sqrt(gns_all[:,9]**2 +gns_all[:,11]**2 ),
                      sigma=3,maxiters=10,cenfunc='median',
                      masked=True)
# all_pm_good=np.array(all_pm_sig)[all_pm_sig.mask==False]
gns_all = gns_all[all_pm_sig.mask==False]


s_all_ra = np.std(gns_all[:,9])
s_all_dec = np.std(gns_all[:,11])
all_color, s_all_color = np.nanmean(gns_all[:,3]-gns_all[:,4]), np.nanstd(gns_all[:,3]-gns_all[:,4])


for i in range(len(gns_young_all)):
    colorines.append("#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)]))
talla = [i*1 for i in gns_young_all[:,4]]
# sizs = np.array(np.where(np.array(sizs)< 15,2,0.2))
colorines = ['#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2','#7f7f7f', '#bcbd22', '#17becf', '#67AFAD', '#B018CC', '#D24FBD', '#A12322', '#4CC351', '#54DF4F', '#7389D2', '#898EE0', '#289C88', '#18EAA4', '#9ECC27', '#71A317', '#421256', '#A23C97', '#44302F']
# defa_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2','#7f7f7f', '#bcbd22', '#17becf']


fs = 10
fb = 40
sizs =[]
for i, t in enumerate(talla):
    if t <mag_lim:
        sizs.append(fb*talla[i]/(np.round(t) -10))
    else: 
        sizs.append(fs*talla[i]/(np.round(t) -10))
sizs = np.array(sizs)
# %%
# Br_gamma map
y_size = 100#!!!
Brg = fits.open(pruebas  +  'brg_emission.fits')
brg_d = Brg[0].data
brg_h = Brg[0].header
del brg_h['CRVAL3'], brg_h['CRPIX3'],brg_h['CDELT3'],brg_h['CUNIT3'],brg_h['CD3_1'],brg_h['CD3_2'],brg_h['CD3_3'],brg_h['CRDER3'],brg_h['CD1_3'],brg_h['CD2_3'],brg_h['CTYPE3 ']

wcs = WCS(brg_h)
# wcs = WCS(fits.getheader(pruebas + 'brg_emission.fits', ext=0)).celestial

fig, ax = plt.subplots(2,1, figsize = (5,10))

# r = Regions.read(pruebas + 'Brg_region_coord.reg',format='ds9')
r = Regions.read(pruebas + 'Brg_region.reg',format='ds9')

young_to_kpix = wcs.wcs_world2pix(gns_young_all[:,0],gns_young_all[:,1],1)
all_to_kpix = wcs.wcs_world2pix(gns_all[:,0],gns_all[:,1],1)

ax[0].scatter(gns_all[:,9],gns_all[:,11], alpha = alp)

# ax[0].scatter(gns_young_all[:,9][bri],gns_young_all[:,11][bri],c = np.array(defa_colors)[bri],s = 60, 
#               alpha = 1,lw =1 , edgecolor = 'k',marker = '^',)
ax[0].scatter(gns_young_all[:,9][bri],gns_young_all[:,11][bri],c = np.array(colorines)[bri],s = y_size, 
              alpha = 1,lw =1 , edgecolor = 'k',marker = '^',)
ax[0].scatter(gns_young_all[0,9],gns_young_all[0,11],facecolor = 'none',s = y_size, 
              lw =1 , edgecolor = 'k',alpha = 0,marker = '^',
              label = '$\sigma \mu_{ra} = %.2f$\n$\sigma \mu_{dec}$ = %.2f'%(s_brig_ra,s_brig_dec))

ax[0].axis('scaled')
ax[1].scatter(gns_all[:,3]-gns_all[:,4],gns_all[:,4], alpha = alp)
ax[1].scatter(gns_young_all[:,3][bri]-gns_young_all[:,4][bri],gns_young_all[:,4][bri],
              edgecolor = 'k',c =np.array(colorines)[bri],s = y_size + 40,marker = '^')
ls_young = []
# for t in range(len(gns_young_all[:,3])):
for t in range(len(bri[0])):
    # ax[1].text(gns_young_all[t,3]-gns_young_all[t,4]+0.3,gns_young_all[t,4],
    #               '[%.0f,%.0f] %.0f,%.0f'%(gns_young_all[t,-2],gns_young_all[t,-1],gns_young_all[t,-4],gns_young_all[t,-3]), 
    #               fontsize = 12,)
    if t>1 and t%2 >0:
        plus = -1
    else:
        plus = 0
    ax[1].text(gns_young_all[t,3]-gns_young_all[t,4]+0.3 + plus,gns_young_all[t,4],
                  'Y%s'%(t+1), 
                  fontsize = 12,)
    if gns_young_all[t,4] < mag_lim:
        print('---- %s'%(t%2))
        print(gns_young_all[t,-2],gns_young_all[t,-1],gns_young_all[t,-4],gns_young_all[t,-3],gns_young_all[t,4])
        ls_young.append((gns_young_all[t,-2],gns_young_all[t,-1],gns_young_all[t,-4],gns_young_all[t,-3],gns_young_all[t,4]))
np.savetxt(pruebas + 'ls_young.txt', np.array(ls_young), fmt = '%.0f '*4 + '%.2f ',header = 'ifu, half, x, y, Ks')
alp_lines = 0.5
ax[1].invert_yaxis()
ax[1].axvline(1.3, color = 'grey', ls = 'dashed', alpha = alp_lines)
# ax[1].axvline(1.9, color = 'grey', ls = 'dashed', alpha =alp_lines)
ax[1].axhline(mag_lim, color = 'grey', ls = 'dashed', alpha = alp_lines)

ax[1].axis('scaled')
ax[1].set_xlim(0.5,5)
# ax[1].axis('scaled')


ax[0].set_xlabel('$\mu_{RA}$ (mas/yr)', fontsize = 12,labelpad=-2)
ax[0].set_ylabel('$\mu_{Dec}$ (mas/yr)',labelpad=-10,fontsize = 12)
ax[1].set_ylabel('Ks',fontsize = 12)
ax[1].set_xlabel('H$-$Ks',fontsize = 12)

ax[0].scatter(gns_all[0,9],gns_all[0,11], alpha = 0, color= '#1f77b4',
              label = '$\sigma \mu_{ra} = %.2f$\n$\sigma \mu_{dec}$ = %.2f'%(s_all_ra,s_all_dec))
ax[1].scatter(gns_young[0,3]-gns_young[:,4],gns_young[:,4],alpha = 0,s =y_size,marker = '^',
              facecolor= 'none', edgecolor = 'k',
              label = '$\overline{H-Ks}$ = %.2f\n$\sigma$ = %.2f'%(brig_color,s_brig_color))
ax[1].scatter(gns_all[0,3]-gns_all[0,4],gns_all[0,4],alpha = 0,color =  '#1f77b4',
              label = '$\overline{H-Ks}$ = %.2f\n$\sigma$ = %.2f'%(all_color,s_all_color))

ax[0].invert_xaxis()
leg_siz = 8
leg0 = ax[0].legend(fontsize = leg_siz +2, loc =1) 
leg1 = ax[1].legend(fontsize = leg_siz) 

for lh in leg0.legend_handles: 
    lh.set_alpha(1)
for lh in leg1.legend_handles: 
    lh.set_alpha(1)
# ax[0].legend()
# ax[1].legend()
l_size = 10
ax[0].tick_params(length=4, width=0.5)
ax[1].tick_params(length=4, width=0.5)
ax[0].xaxis.set_tick_params(labelsize=l_size )
ax[0].yaxis.set_tick_params(labelsize=l_size )
ax[1].xaxis.set_tick_params(labelsize=l_size )
ax[1].yaxis.set_tick_params(labelsize=l_size )
plt.savefig(pruebas + 'pm_color_young.png', dpi =300, bbox_inches = 'tight')
# 

# sys.exit(236)
# fig.colorbar(im0, ax=ax[0], orientation='vertical')
# %%
fig, ax = plt.subplots(figsize =(8,8))
ax.scatter(all_to_kpix[0], all_to_kpix[1])
ax.axis('scaled')

# %%
lib_pruebas = '/Users/amartinez/Desktop/PhD/Libralato_data/pruebas/'
# Ra_cl, Dec_cl, mura_cl,mudec_cl, H, Ks
Ra_cl, Dec_cl, mura_cl,mudec_cl, H_cl, Ks_cl,ms_id= np.loadtxt(lib_pruebas + 'clus_14996_16_55.txt', unpack=True)
# lib_cl= np.loadtxt(lib_pruebas + 'clus_14996_16_55.txt', unpack=True)
fig, ax = plt.subplots(1,1)
ax.scatter(gns_young_all[:,0],gns_young_all[:,1],s = sizs, c= colorines)
ax.scatter(Ra_cl,Dec_cl, marker = 'x')
# %%


fig, ax = plt.subplots(1,2)
ax[0].scatter(gns_young_all[:,10], gns_young_all[:,12], s = sizs, c = colorines)
ax[0].set_xlabel('$\sigma \mu_{RA}$(mas/yr)')
ax[0].set_ylabel('$\sigma \mu_{Dec}$(mas/yr)')
ax[0].axis('scaled')
error = np.sqrt(gns_young_all[:,10]**2 + gns_young_all[:,12]**2)
ax[1].scatter(gns_young_all[:,4], error, s = sizs, c = colorines)
ax[1].set_ylabel('$\sigma \mu$(mas/yr)')
ax[1].set_xlabel('Ks')
ax[1].axis('scaled')


# %%
eso_pro = '/Users/amartinez/Desktop/PhD/ESO_proposals/KMOS'
esorex_ima_5= '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p105_%s/comoving_group_mosaic_K_Half1_COMBINED_IMAGE_mapping.fits'%(reduction)

# abc_one =

# Load your FITS file and create the WCS
# fits_file = get_pkg_data_filename('one_ABC.fits')
image_data = fits.getdata(pruebas + 'brg_emission.fits', ext=0)
image_data = np.squeeze(image_data) 
mapa = WCS(fits.getheader(pruebas + 'brg_emission.fits', ext=0)).celestial
fig, ax = plt.subplots(subplot_kw={'projection': mapa})  # Adjust the figsize as needed

fig = plt.figure()  # Adjust the figsize as needed
ax = plt.subplot(projection=mapa)
lon = ax.coords[0]
lat = ax.coords[1]
lon.set_ticks(spacing=5. * u.arcsec)
lat.set_ticks(spacing=5. * u.arcsec)
ax.tick_params(axis = 'y',which = 'both',labelright = True, labelleft = True)
ax.imshow(image_data, vmin=-0.8e-20, vmax=0.3e-17, origin='lower', cmap='Greys', label ='KMOS')

ax.grid()

lon.set_ticks_visible(True)
lon.set_ticklabel_visible(True)
# lat.set_ticklabel_visible(True)
# lat.set_ticklabel_visible(True)
# lat.set_ticks_visible(False)
lat.set_ticklabel(rotation='vertical')
ax.coords[0].set_axislabel('RA')
ax.coords[1].set_axislabel('Dec')
xticks = ax.get_xticks()
yticks = ax.get_yticks()
print("X-axis tick locations:", xticks)
print("Y-axis tick locations:", yticks)
sys.exit(343)
# ax.
# %%
# Create a latex table
df = pd.DataFrame({'ID':[],
                   'RA(º)':[],
                   'Dec(º)':[], 
                   '$\mu_{RA}$(mas/yr)':[],
                   '$\mu_{Dec}$(mas/yr)':[],
                   'H':[],
                   'Ks':[]})
ID_ls = ['Y%.0f'%(k+1) for k in range(len(bri[0]))]
df['ID'] = ID_ls
df['RA(º)'] = ['%.3f'%(gns_young_all[bri][i,0]) for i in range(len(bri[0]))]
df['Dec(º)'] = ['%.3f'%(gns_young_all[bri][i,1]) for i in range(len(bri[0]))]


# 0RA_gns 	1DE_gns 	2Jmag 	3Hmag 	4Ksmag 	5ra 	6Dec 	7x_c 	8y_c 	9mua 10dmua 	11mud 	12dmud 	13time 	14n1	15n2	16ID 	17mul 	18mub 	19dmul 	20dmub 	21m139 	22Separation
ls_mura = ['%.02f$\pm$%.02f'%(gns_young_all[bri][i,9],gns_young_all[bri][i,10]) for i in range(len(bri[0]))]
ls_mudec = ['%.02f$\pm$%.02f'%(gns_young_all[bri][i,11],gns_young_all[bri][i,12]) for i in range(len(bri[0]))]
ls_H = ['%.02f'%(gns_young_all[bri][i,3]) for i in range(len(bri[0]))]
ls_Ks = ['%.02f'%(gns_young_all[bri][i,4]) for i in range(len(bri[0]))]



df['$\mu_{RA}$(mas/yr)'] = ls_mura
df['$\mu_{Dec}$(mas/yr)'] = ls_mudec
df['H'] = ls_H
df['Ks'] = ls_Ks


print(df.to_latex(index = False, column_format = 'c'*len(df.columns), label = 'tab:ysos', 
                  caption = 'YSO parametres. Uncertainties in the propermotions as they appear in LIB21'))
# %%

# # Enable automatic plotting mode
# # IPython.get_ipython().run_line_magic('matplotlib', 'auto')
# IPython.get_ipython().run_line_magic('matplotlib', 'inline')

eso_pro = '/Users/amartinez/Desktop/PhD/ESO_proposals/KMOS'
esorex_ima_5= '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/p105_%s/comoving_group_mosaic_K_Half1_COMBINED_IMAGE_mapping.fits'%(reduction)

# abc_one =

# Load your FITS file and create the WCS
# fits_file = get_pkg_data_filename('one_ABC.fits')

brg_data = fits.getdata(pruebas + 'brg_emission.fits', ext=0)
im_data = fits.getdata(esorex_ima_5, ext=1)
brg_data = np.squeeze(brg_data) 
mapa = WCS(fits.getheader(pruebas + 'brg_emission.fits', ext=0)).celestial




fig, ax = plt.subplots(2,1,subplot_kw={'projection': mapa}, figsize = (12,12))  # Adjust the figsize as needed


# r = Regions.read(pruebas + 'Brg_region_coord.reg',format='ds9')
r = Regions.read(pruebas + 'Brg_region.reg',format='ds9')


# for i in range(11):
#     ax[0].plot(r[i].vertices.x,r[i].vertices.y, color = 'white',zorder = 1,alpha = 1,lw =3,)


lon = ax[0].coords[0]
lat = ax[0].coords[1]
lon1 = ax[1].coords[0]
lat1 = ax[1].coords[1]

lon.set_ticks(spacing=5. * u.arcsec)
lat.set_ticks(spacing=7. * u.arcsec)
lon1.set_ticks(spacing=5. * u.arcsec)
lat1.set_ticks(spacing=7. * u.arcsec)
ax[0].tick_params(axis = 'y',which = 'both',labelright = False, labelleft = True, fontsize = 20)
ax[1].tick_params(axis = 'y',which = 'both',labelright = True, labelleft = True)

ax[1].imshow(brg_data, vmin=-0.8e-20, vmax=0.3e-17, origin='lower', cmap='Greys', label ='KMOS')
ima1 = ax[0].imshow(im_data, vmin=-0.011e-16, vmax=0.136e-16, origin='lower', cmap='Greys', label ='KMOS')
ax[0].scatter(young_to_kpix[0][bri],young_to_kpix[1][bri], c = np.array(colorines)[bri],s = sizs[bri]*0.7,
              edgecolor = 'k')
ax[1].scatter(young_to_kpix[0][bri],young_to_kpix[1][bri],  c = np.array(colorines)[bri],s = sizs[bri]*0.7,
              edgecolor = 'k')

# ax[0].grid()
l_size = 20
lon.set_ticks_visible(True)
lon.set_ticklabel_visible(False)
# lat.set_ticklabel_visible(True)
# lat1.set_ticklabel_visible(False)
# lat1.set_ticks_visible(False)
lat.set_ticklabel(rotation='vertical', fontsize = l_size)
lat1.set_ticklabel(rotation='vertical',fontsize = l_size)
lon1.set_ticklabel(fontsize = l_size)

ax[0].coords[0].set_axislabel('RA',fontsize = l_size)
ax[0].coords[1].set_axislabel('Dec',fontsize = l_size)
ax[1].coords[0].set_axislabel('RA',fontsize = l_size)
ax[1].coords[1].set_axislabel('Dec ',fontsize = l_size)

xticks = ax[0].get_xticks()
yticks = ax[0].get_yticks()
print("X-axis tick locations:", xticks)
print("Y-axis tick locations:", yticks)

ax[0].tick_params(length=4, width=0.5)
ax[1].tick_params(length=4, width=0.5)

ax[0].xaxis.set_tick_params(labelsize=15)
ax[0].yaxis.set_tick_params(labelsize=15)
ax[1].xaxis.set_tick_params(labelsize=15)
ax[1].yaxis.set_tick_params(labelsize=15)

plt.savefig(pruebas + 'im_plus_brg.png', dpi =300, bbox_inches = 'tight')
# 
# =============================================================================
# # Add a slider for vmin
# v_min =np.nanmin(im_data)
# v_max =np.nanmax(im_data)
# ax_vmin = plt.axes([0.2, 0.1, 0.3, 0.03], facecolor='lightgoldenrodyellow')
# vmin_slider = Slider(ax_vmin, 'vmin', v_min, v_max, valinit=v_min)
# 
# # Add a slider for vmax
# ax_vmax = plt.axes([0.2, 0.05, 0.3, 0.03], facecolor='lightgoldenrodyellow')
# vmax_slider = Slider(ax_vmax, 'vmax', v_min, v_max, valinit=v_max)
# 
# def update(val):
#     vmin = vmin_slider.val
#     vmax = vmax_slider.val
#     ima1.set_clim(vmin, vmax)
#     plt.draw()
#     
# vmin_slider.on_changed(update)
# vmax_slider.on_changed(update)
# =============================================================================
# %%
# Find late typ stars with similar magnitudes of those of the YSO
ls_y = np.array(ls_young)
gns_all_but = gns_all
for y in ls_y:
    mag = y[4]
    # print(mag)
    equal_mag = np.where(gns_all[:,4] == mag)
    # print(equal_mag[0])
    gns_all_but = np.delete(gns_all_but, equal_mag[0], axis = 0) 
    # print(len(gns_all),len(gns_all_but))
    #rint(close_mag)

# %
late_type = []
for y in ls_y:
    mag = y[4]
    print('----')
    print(mag, y[0],y[1])
    close_mag = np.where(abs(gns_all_but[:,4]-mag) == min(abs(gns_all_but[:,4]-mag)))
    # print(min(abs(gns_all_but[:,4]-mag)))
    if min(abs(gns_all_but[:,4]-mag)) == 0:
        gns_all_but = np.delete(gns_all_but, close_mag[0], axis = 0) 
        close_mag = np.where(abs(gns_all_but[:,4]-mag) == min(abs(gns_all_but[:,4]-mag)))   
        # print('∫∫∫∫')
        # print(min(abs(gns_all_but[:,4]-mag)), close_mag)
        # print('∫∫∫∫')
    print(gns_all_but[close_mag][0][4],gns_all_but[close_mag][0][-2],gns_all_but[close_mag][0][-1])
    print('****')
    late_type.append((gns_all_but[close_mag][0][-2],gns_all_but[close_mag][0][-1],gns_all_but[close_mag][0][4]))
np.savetxt(pruebas + 'late_type_for_comparation.txt', np.array(late_type), fmt = 2*'%.0f ' + '%.8f',header = 'ifu, half, Ks')











