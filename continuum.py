
"""
Created on Tue Nov  7 09:57:30 2023

@author: amartinez
"""

# We are going to make a continuum map 
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
ifu_sel = 'all'
# ifu_sel = 7 
half_ifu = 2

aligned_cube = '/Users/amartinez/Desktop/PhD/KMOS/Kmos_iMac/%s_reduction/cubes/cube_ifu%s_half%s_mean.fits'%(reduction,ifu_sel,half_ifu)


if type(ifu_sel) == int:
    if os.path.isfile(pruebas  + 'sky_ifu%s_half%s.fits'%(ifu_sel,half_ifu)):
        print('There is sky spectrum for this region. Go On!')
    else:
        print('There is no sky computed for ifu = %s, half = %s'%(ifu_sel,half_ifu))
        sys.exit('DO IT!')
    dic_x = {}
    dic_y = {}
    for ifu in range(1,24):
        x_shifts = []
        y_shifts = []
        with open(log_5 + 'esorex.log', 'r') as f:
            for line in f:
                if 'IFU %s]'%(ifu) in line:
                    ind1 = line.find('in x:')
                    ind2 = line.find('in y:')
                    ind_x = ind1 + 6
                    ind_y = ind2 + 6
                    x_shifts.append(float(line[ind_x:ind_x+7]))
                    y_shifts.append(float(line[ind_y:ind_y+7]))
                    # print (line)
                    # break
                else:
                    continue
            x_shifts = np.array(x_shifts)
            y_shifts = np.array(y_shifts)
            if len(x_shifts)>0:
                dic_x['ifu%s'%(ifu)] = np.asarray(np.rint(217 - x_shifts),int)
                dic_y['ifu%s'%(ifu)] = np.asarray(np.rint(433 + y_shifts),int)
                # dic_x['ifu%s'%(ifu)] = np.asarray(np.rint(217 - x_shifts),int)
                # dic_y['ifu%s'%(ifu)] = np.asarray(np.rint(433 + y_shifts ),int)
    if half_ifu == 0:
        y_d = int(min(dic_y['ifu%s'%(ifu_sel)])-27)
        y_up = int(max(dic_y['ifu%s'%(ifu_sel)]))
        x_d = int(min(dic_x['ifu%s'%(ifu_sel)]))
        x_up = int(max(dic_x['ifu%s'%(ifu_sel)])+27)
    elif half_ifu == 1:
        y_d = int(min(dic_y['ifu%s'%(ifu_sel)])-27)
        y_up = np.unique(dic_y['ifu%s'%(ifu_sel)])[1]
        x_d = int(min(dic_x['ifu%s'%(ifu_sel)]))
        x_up = int(max(dic_x['ifu%s'%(ifu_sel)])+27)
    elif half_ifu == 2:   
        y_d = np.unique(dic_y['ifu%s'%(ifu_sel)])[1]
        y_up = np.unique(dic_y['ifu%s'%(ifu_sel)])[3]
        x_d = int(min(dic_x['ifu%s'%(ifu_sel)]))
        x_up = int(max(dic_x['ifu%s'%(ifu_sel)])+27)

# %%

# ima = fits.open(esorex_cube_5)

if ifu_sel == 'all':
    ima = fits.open(esorex_cube_5)
else:
    ima = fits.open(aligned_cube)
mapa = WCS(ima[1].header).celestial

h0 = ima[0].header
h1 = ima[1].header

# This is the header for the spectrum fits file
# h_esp = h1.copy()
model = '/Users/amartinez/Desktop/PhD/KMOS/practice/spec_header_model.fits'
h_esp = fits.open(model)[0].header
# h_esp['CRVAL1'] = h_esp['CRVAL3']
# h_esp['CDELT1'] = h_esp['CDELT3']
# h_esp['CTYPE1'] = 'WAVE'
# h_esp['CUNIT1'] = 'um'
# del h_esp['CRVAL3'], h_esp['CRVAL2'], h_esp['CDELT2'], h_esp['CDELT3']

ima_1 = fits.open(esorex_ima_5)
wcs = WCS(ima[1].header)
h0_1 = ima_1[0].header
h1_1 = ima_1[1].header

cube = ima[1].data

# cube = cube[:,y_d:y_up, x_d:x_up]
if ifu_sel == 'all':
    cube_im = ima_1[1].data
else:
    cube_im = ima_1[1].data[y_d:y_up, x_d:x_up]



fig, ax = plt.subplots(1,1,figsize =(10,10)) 
ax.set_title('Continuum')
#Uncomment these for the WCS projection 
# =============================================================================
# ax.axis("off")
# ax = plt.subplot(projection=wcs,slices=(1000, 'y', 'x'))
# =============================================================================



if ifu_sel == 'all':
    v_min =0.01e-16
    v_max = 0.1e-16
    im = ax.imshow(cube_im, cmap='inferno', origin='lower',vmin=v_min,vmax = v_max,alpha =1)
else:
    v_min =0.001e-16
    v_max = 0.05e-16
    cube_plot = np.nanmedian(cube, axis = 0)
    im = ax.imshow(cube_plot, cmap='inferno',label='overlays',origin='lower',vmin=v_min,vmax = v_max,alpha =1)
fig.colorbar(im, ax=ax, orientation='vertical')
# sys.exit(173)
# x, y = 62,43
# yy,xx = np.indices(cube_im.shape)
# distances = np.sqrt((xx - x)**2 + (yy - y)**2)
# mask = np.where(distances <=3)
# values = cube[:,mask[0],mask[1]]
# spec = np.mean(values, axis = 1)
# ax.scatter(mask[1],mask[0],color = 'r')

x = None
y = None
wide = 4#TODO
w = int(wide/2)
def extract_spec(cube,y,x):
    # yy,xx = np.indices(cube_im.shape)
    xx,yy = np.indices(cube_im.shape)
    distances = np.sqrt((xx - x)**2 + (yy - y)**2)
    mask = np.where(distances <=w)
    spec = cube[:,mask[0],mask[1]]
    spec_mean = np.mean(spec, axis = 1)
    # This subtract the sky that I don not know if it is correct
    # hdu_sky = fits.open(pruebas  + 'sky_ifu%s_half%s.fits'%(ifu_sel,half_ifu))
    # sky_data = hdu_sky[0].data
    # spec_mean = spec_mean - sky_data
    ax.scatter(mask[1],mask[0],color = 'green')
    
    # For normalizing, uncommnet
# =============================================================================
#     norm =  np.where((lam > d_norm) & (lam < u_norm))
#     filtered_data = sigma_clip(spec_mean[norm], sigma=1.5, maxiters=None,cenfunc='median', masked=False, copy=False)
#     spec_mean = spec_mean/np.mean(filtered_data)
# =============================================================================
    
    return spec_mean
    
    
   

# Function to update the plot with the clicked point
def update_plot(x,y):
    x = int(np.rint(x))
    y = int(np.rint(y))
    plt.scatter(x, y, color='blue', marker='x', s=100)  # Customize the marker style
    # square = Rectangle((x - wide/2-0.5, y - wide/2-0.5), wide, wide, linewidth=1, edgecolor='r', facecolor='none')
    circle = Circle((x, y),w, facecolor = 'none', edgecolor = 'green', linewidth=1)
    # ax.add_patch(square)
    ax.add_patch(circle)
    ax.text(x-w,y-w-1,'%.0f,%.0f'%(x,y), color = 'green',fontsize = 15)
    plt.draw()  # Redraw the plot with the updated point


cab = ima[1].header
lam = np.array([cab['CRVAL3']+cab['CDELT3']*i for i in range(cab['NAXIS3'])])



fig2, ax2 = plt.subplots(1,1,figsize =(20,10))
props = dict(boxstyle='round', facecolor='grey', alpha=0.5)
ax2.text(0.05, 0.95, 'D-click left: select interval limitis\nD-click righ: delete interval\n"g": save continuum', transform=ax.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)

HeI = 2.058
COI =  2.29322
COII = 2.32246
Brg = 2.165
He = 2.12
HeII = 2.189
# H2 = 2.12
l_names = ['HeI', 'COI', 'COII','Br$\gamma$', 'He', 'HeII']
lines = [HeI, COI, COII,Brg, He, HeII]
for l in lines:
    ax2.axvline(l, color = 'grey', alpha = 0.5, ls = 'dashed')
ax2b = ax2.twiny()
ax2b.set_xticks(lines)
tl = l_names
ax2b.set_xticklabels(tl)
d_norm, u_norm =  2.2,2.28
norm =  np.where((lam > d_norm) & (lam < u_norm))
def plot_spec(flux):
    
    factor = 5e-18
    # ax2.set_ylim(np.nanmean(flux)-sig, np.nanmean(flux) + sig)
    ax2.plot(lam,flux+ clicks*factor, label = '(%.0f,%.0f)'%(x,y))
    ax2.set_xlabel('$\lambda (\mu m)$')
    ax2.set_ylabel('flux')
    ax2.legend(fontsize=10)
    ax2b.plot(lam,flux + clicks*factor, alpha = 0)
    print(f'Clicks = {clicks}')

clicks = 0 
def onclick(event):
    global x, y, clicks  # Use global variables
    if event.xdata is not None and event.ydata is not None and event.dblclick:
        x = event.xdata
        y = event.ydata 
        x = int(np.rint(x))
        y = int(np.rint(y))
        flux = extract_spec(cube,int(np.rint(x)), int(np.rint(y)))
        update_plot(x,y)
        plot_spec(flux)
        clicks += 1
        sig = np.nanstd(flux)/2
        print(f'Clicked at x={x}, y={y}')
        
        return flux
chunk = 0     
no_lines = []   
continuum = []
continuum_inds = []
chunk_ind = []
def vertical(event):
    global x2, y2, chunk, no_lines  # Use global variables
    if event.dblclick and event.button ==1:
        no_lines.append(event.xdata)
        ax2.axvline(event.xdata, ls = 'dashed',color = 'r' )
        if len(no_lines) == 2:
            ind1 = np.where(abs(lam-no_lines[0])== min(abs(lam-no_lines[0])))
            ind2 = np.where(abs(lam-no_lines[1])== min(abs(lam-no_lines[1])))
            
            ax2.axvline(lam[ind1], ls = 'dashed',color = 'green' )
            ax2.axvline(lam[ind2], ls = 'dashed',color = 'green' )
            # print(extract_spec)
            continuum.append(lam[ind1[0][0]:ind2[0][0]+1])
            continuum_inds.append(np.arange(ind1[0][0], ind2[0][0]+1,1))
            chunk_ind.append(np.full((ind2[0][0]+1-ind1[0][0]),chunk))
            ax2.axvspan(lam[ind1][0], lam[ind2][0],color = 'green', alpha = 0.2)
            # ax2.axvspan(no_lines[0], no_lines[1], color = 'r', alpha = 0.2)
            no_lines = []
            chunk +=1
            print('chunk_ind', chunk_ind)
    if event.dblclick and event.button ==3:
            del continuum[-1]
            del continuum_inds[-1]
            del chunk_ind[-1]
            no_lines = []
    # if event.xdata is not None and event.ydata is not None:
    #     ax2.axvline(event.xdata)

def save_continuum(event):
    if event.key == 'g':
        cont = np.concatenate(continuum)
        cont_id = np.concatenate(continuum_inds)
        chunk_id = np.concatenate(chunk_ind)
        if ifu_sel == 'all':
            np.savetxt(pruebas + 'continuum_all.txt', np.c_[cont,cont_id,chunk_id], fmt = '%.8f %.0f %.0f')
        else:            
            np.savetxt(pruebas + 'continuum_ifu%s_half%s.txt'%(ifu_sel,half_ifu), np.c_[cont,cont_id,chunk_id], fmt = '%.8f %.0f %.0f')
        fig.canvas.mpl_disconnect(cid)
        fig2.canvas.mpl_disconnect(cidv)
        fig2.canvas.mpl_disconnect(save_continuum)
        plt.close()
        # plt.close()
    
cid = fig.canvas.mpl_connect('button_press_event',onclick)
cidv = fig2.canvas.mpl_connect('button_press_event',vertical)
cids = fig2.canvas.mpl_connect('key_press_event',save_continuum)

# %%
# fig, ax = plt.subplots()
# ax.plot(np.random.rand(10))

# def onclick(event):
#     if event.dblclick and event.button ==3:
#         print('yomamma')
#     print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
#           ('double' if event.dblclick else 'single', event.button,
#             event.x, event.y, event.xdata, event.ydata))

# cid = fig.canvas.mpl_connect('button_press_event', onclick)





