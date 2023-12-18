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
from scipy.spatial import distance
import re
from compar_lists import compare_lists
import time
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

ifu_sel =13#!!!
half_ifu =1#!!!
max_d = 5#!!!
min_d = 2#!!!
search_r = 0.0015#!!!


reduction = 'ABC'
pruebas = '/Users/amartinez/Desktop/PhD/KMOS/practice/'
choped_ifus = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/ifu_alignment_ABC/'
gns_ls = '/Users/amartinez/Desktop/PhD/KMOS/GNS_lists/'
spec_folder = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/%s_reduction/cluster_spectra/ifu_%s/half_%s/'%(reduction, ifu_sel, half_ifu)
spec_young = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/%s_reduction/young_candidates/ifu_%s/half_%s/'%(reduction, ifu_sel, half_ifu)
# sources = 'auto'#!!!
sources = 'hand'#!!!
np.savetxt(spec_folder + 'distance_for_aa_ifu%s_half%s.txt'%(ifu_sel,half_ifu), np.array([max_d,min_d]).T,fmt = '%.0f',header ='max distance,  min distance')

fig, ax = plt.subplots(1,1,figsize =(10,10)) 
# im = fits.open(pruebas  + 'Kmos_image.fits')

ima = fits.open(choped_ifus  + 'ifu%s_half%s_p105.fits'%(ifu_sel,half_ifu))
data = ima[1].data
cube_plot = data
v_min =np.nanmin(cube_plot)
v_max =np.nanmax(cube_plot)
im = ax.imshow(cube_plot, cmap='inferno',label='overlays',origin='lower',vmin=v_min,vmax = v_max,alpha =1)
fig.colorbar(im, ax=ax, orientation='vertical')
ax.set_title('IFU%s,half%s'%(ifu_sel, half_ifu))

# Add a slider for vmin
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


# %%
if sources =='auto':
    data = data.byteswap().newbyteorder()
    bkg = sep.Background(data)
    data_sub = data - bkg
    obs = sep.extract(data_sub, 3, err=bkg.globalrms)
    
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
    x = obs['x']
    y = obs['y']
# %%

if sources == 'hand':
    x = []
    y = []
    x_y = []
    y_y = []
    for i in range(len((glob.glob(spec_folder +'/spec*')))):
        name = os.path.basename(glob.glob(spec_folder +'/spec*')[i])
        # print(name)
        coord = re.findall(r'\d+', name)
        x0,y0 = np.array(coord[2]).astype(int), np.array(coord[3]).astype(int)
        x.append(y0)
        y.append(x0)
        print(x0,y0)
    for j in range(len((glob.glob(spec_young +'/spec*')))):
        name_y = os.path.basename(glob.glob(spec_young +'/spec*')[j])
        coord_y = re.findall(r'\d+', name_y)
        x0_y,y0_y = np.array(coord_y[2]).astype(int), np.array(coord_y[3]).astype(int)
        print('------')
        print(y0_y,x0_y)
        x_y.append(x0_y)
        y_y.append(y0_y)
        # print(name_y)
    
        
    x = np.array(x)
    y = np.array(y)
    x_y = np.array(x_y)
    y_y = np.array(y_y)
    ax.scatter(x,y, s = 100, color ='lime', marker = 'x')
    ax.scatter(x_y,y_y, s = 100, color ='lime', marker = 'o')
   
plt.draw()

# %%
wcs = WCS(ima[1].header)

x_mid, y_mid = ima[1].data.shape[1], ima[1].data.shape[0]
radec_m = wcs.wcs_pix2world(np.mean(x),np.mean(y),0)
lb_m = SkyCoord(ra = radec_m[0], dec = radec_m[1], unit = 'deg').galactic
np.savetxt(pruebas + 'middle_test.txt',np.array([radec_m[0],radec_m[1]]), fmt = '%.8f')
# 0RA_gns 	1DE_gns 	2Jmag 	3Hmag 	4Ksmag 	5ra 	6Dec 	7x_c 	8y_c 	9mua 	
#10dmua 	11mud 	12dmud 	13time 	14n1	15n2	16ID 	17mul 	18mub 	19dmul 	
#20dmub 	21m139 	22Separation
# gns = np.loadtxt(gns_ls + 'GNS_LIB_overKMOS_dmu5.txt')
# gns = np.loadtxt('/Users/amartinez/Desktop/PhD/Libralato_data/results/sec_All_match_GNS_and_WFC3IR_refined_galactic.txt')
gns = np.loadtxt('/Users/amartinez/Desktop/PhD/Libralato_data/results/relaxed_match_GNS_and_WFC3IR_refined_galactic.txt')

selec = np.where(np.sqrt((gns[:,0]-radec_m[0])**2 +(gns[:,1]-radec_m[1])**2)<search_r)


gns_=gns[selec[0]]
# gns_id = np.c_[gns_,np.arange(len(gns_))]

coor = np.array([gns_[:,0],gns_[:,1]]).T
np.savetxt(pruebas + 'around_ifu%s_half%s.txt'%(ifu_sel, half_ifu),coor)
xy = np.array([x,y]).T


coor_pix = wcs.wcs_world2pix(coor[:,0],coor[:,1],0)
gns_= np.c_[gns_,coor_pix[0],coor_pix[1]]
fig, ax = plt.subplots(1,1)
ax.set_title('IFU%s,half%s'%(ifu_sel, half_ifu))
ax.scatter(coor_pix[0],coor_pix[1],s=100, color = 'r', label = 'GNS')
ax.scatter(xy[:,0],xy[:,1],color = 'k', label = 'KMOS')

coor_pix = np.array([coor_pix[0],coor_pix[1]]).T
# %%

dist_mat=distance.cdist(coor_pix[:],xy[:], 'euclidean')

# close = np.where(dist_mat<max_d)[0]
# close_xy = np.where(dist_mat<max_d)[1]
close = np.where((dist_mat<max_d) &(dist_mat>min_d))[0]
close_xy = np.where((dist_mat<max_d)&(dist_mat>min_d))[1]


coor_pix_close = coor_pix[close]
coor_pix_close = coor_pix_close[np.unique(coor_pix_close,axis=0, return_index = True)[1]]
xy_close = xy[close_xy]
xy_close = xy_close[np.unique(xy_close,axis=0, return_index = True)[1]]
ax.scatter(coor_pix_close[:,0],coor_pix_close[:,1],s=200, color = 'g',alpha = 0.2)
ax.scatter(xy_close[:,0],xy_close[:,1],s=200, color = 'b',alpha = 0.2)
ax.set_xlabel('x pix')
ax.set_ylabel('y pix')
ax.legend()
# coor_pix_close = np.array(coor_pix_close, dtype="<f4")
# xy_close= np.array(xy_close, dtype="<f4")
coor_pix_close = coor_pix_close.astype(float)
xy_close= xy_close.astype(float)
save_lits = 0
m, (source_list, target_list) = aa.find_transform(coor_pix_close, xy_close,max_control_points=300)
save_lists = 1
# sys.exit(166)

# %%

# m, (source_list, target_list) = aa.find_transform(coor_pix, xy,max_control_points=300)

print('source_list',source_list)
print('target_list',target_list)

print("Translation: (x, y) = (%.2f, %.2f)"%(m.translation[0],m.translation[1]))
print("Rotation: %.3f degrees"%(m.rotation * 180.0 / np.pi))
print("Scale factor: %.4f"%(m.scale))

coor_pix_t = aa.matrix_transform(coor_pix, m.params)

fig, ax = plt.subplots(1,1)
ax.set_title('IFU%s,half%s'%(ifu_sel, half_ifu))
ax.scatter(coor_pix_t[:,0],coor_pix_t[:,1],s=100, color ='red', label ='aa_GNS', facecolor = 'none')
ax.scatter(xy[:,0],xy[:,1], color = 'k', label = 'KMOS')
ax.scatter(target_list[:,0], target_list[:,1], color = 'lime', s= 300, alpha = 0.5, label ='KMOS in GNS',zorder = 0)
ax.legend()
# %%
xy_y = np.array([x_y,y_y]).T
distancia = 2

lc1, lc2, l1 , l2 = compare_lists(coor_pix_t, xy,0,1,0,1,distancia)
lc1_y, lc2_y, l1_y, l2_y = compare_lists(coor_pix_t, xy_y, 0,1,0,1,distancia)

ax.scatter(lc2[:,0],lc2[:,1], 
            s= 300, marker = '^',alpha = 0.2, color = 'b', label ='Compare lsts\nmax_dis = %spix'%(distancia))    
ax.scatter(lc1[:,0],lc1[:,1], 
            s= 300, marker = '<',alpha = 0.2, color = 'b') 
# ax.scatter(lc1_y[:,0],lc1_y[:,1], 
#             s= 500, marker = 'o',alpha = 0.2, color = 'red')  
ax.scatter(xy_y[:,0],xy_y[:,1], 
            s= 500, marker = 'o',alpha = 0.2, color = 'red', label = 'Young = %s'%(len(xy_y)))  

# ax.scatter(xy[l2][:,0],xy[l2][:,1], 
#             s= 300, marker = '^',alpha = 0.2, color = 'b', label ='Compare lsts\nmax_dis = %spix'%(distancia))    
# ax.scatter(coor_pix_t[l1][:,0],coor_pix_t[l1][:,1], 
#             s= 300, marker = '<',alpha = 0.2, color = 'b')    
ax.set_xlabel('x pix')
ax.set_ylabel('y pix')
ax.legend()
# %%
fig, ax = plt.subplots(1,1)
ax.scatter(gns_[:,-2][l1], gns_[:,-1][l1])
ax.scatter(gns_[:,-2][l1_y], gns_[:,-1][l1_y],s = 100, color = 'red', alpha = 0.5)

ax.set_xlabel('x pix')
ax.set_ylabel('y pix')
fig, ax = plt.subplots(1,1)
ax.scatter(gns_[:,9][l1], gns_[:,11][l1])
ax.scatter(gns_[:,9][l1_y], gns_[:,11][l1_y],s = 100, color = 'red', alpha = 0.5)
ax.set_xlabel('mu_ra')
ax.set_ylabel('mu_dec')

youg_tosave = gns_[l1_y]
youg_tosave[:,-2:] = lc2_y

if save_lists == 1:
    gns_header = '0RA_gns 1DE_gns 	2Jmag 	3Hmag 	4Ksmag 	5ra 	6Dec 	7x_c 	8y_c 	9mua 10dmua 	11mud 	12dmud 	13time 	14n1	15n2	16ID 	17mul 	18mub 	19dmul 	20dmub 	21m139 	22Separation 23Xkmos 24Ykmos'
    np.savetxt(spec_folder + 'gns_lib_in_kmos_ifu%s_half%s.txt'%(ifu_sel, half_ifu),gns_[l1],header = gns_header, fmt = '%.8f')
    if len(xy_y) > 0:
        np.savetxt(spec_young + 'gns_lib_in_young_kmos_ifu%s_half%s.txt'%(ifu_sel, half_ifu),youg_tosave,header = gns_header, fmt = '%8f')
        np.savetxt(spec_young + 'xy_coord_young_kmos_ifu%s_half%s.txt'%(ifu_sel, half_ifu),xy_y, header = 'x, y',fmt ='%.8f')
# %%
# fig, ax = plt.subplots(1,1)
# test_g = np.loadtxt(spec_folder + 'gns_lib_in_kmos_ifu%s_half%s.txt'%(ifu_sel, half_ifu))
# test_p = np.loadtxt(spec_young + 'gns_lib_in_young_kmos_ifu%s_half%s.txt'%(ifu_sel, half_ifu))
# ax.scatter(test_g[:,9],test_g[:,11], s=200)
# ax.scatter(test_p[:,9],test_p[:,11], s=100)

    
    
    
    

