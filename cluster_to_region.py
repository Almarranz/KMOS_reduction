#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 17 12:06:39 2022

@author: amartinez
"""

# =============================================================================
# Creates a .reg file in order to check cluster members in DS9
# =============================================================================




import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import QTable
from matplotlib import rcParams
import os
import glob
import sys
import math
# pruebas = '/Users/amartinez/Desktop/PhD/HAWK/pm_gns1_gns2_off/pruebas/'
pruebas = '/Users/amartinez/Desktop/PhD/regions/KMOS/'

pm_folder ='/Users/amartinez/Desktop/PhD/HAWK/pm_gns1_gns2_off/'
# pm_folder ='/Users/amartinez/Desktop/PhD/HAWK/pm_gns1_gns2_wcs/'

field_one = 7
chip_one = 4
field_two = 7
chip_two =1 


# ra, dec, l, b, pml, pmb,J, H, Ks,x, y, Aks_mean, dAks_mean, radio("),cluster_ID
name = 'cluster0_0_0_knn_12_area_2.77_pm_color.txt' 
# cluster = '/Users/amartinez/Desktop/morralla/GNS1_f7c4_dmu1.4_at_GNS2_f7c1_2023-01-24 13:40:07.809861/cluster_num0_1_knn25_area1.40/cluster0_0_0_knn_25_area_1.40_pm_color.txt'
cluster = '/Users/amartinez/Desktop/PhD/KMOS/gns1_data/young_bg_stars.txt'

ra, dec = np.loadtxt(cluster,unpack=True, usecols=(0,1))
# pmra = pmra*-1
with open(pruebas+ 'young_bg_stars.reg', 'w') as f:
    f.write('# Region file format: DS9 version 4.1'+"\n"+'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'+"\n"+'fk5'+'\n')
    f.close

# with arrows    
# =============================================================================
# for i in range(len(ra)):
#     with open(pruebas+ 'GNSrel_reg_f%sc%s_clus%s_stars%s.reg'%(field_one,chip_one,name[7], len(ra)), 'a') as f:
#         
#         if pmdec[0]>0 and pmdec[0]>0:
#             f.write('\n'.join(('point(%s,%s) # point=x'%(float(ra[i]),float(dec[i])),'# vector(%s,%s,%s",%s)'%(float(ra[i]),float(dec[i]),np.sqrt(pmra[i]**2+pmdec[i]**2)*10,math.degrees(math.atan(pmdec[i]/pmra[i]))),'\n')))   
#             print('ssss')
#         elif pmra[0]<0 and pmdec[0]>0:
#             f.write('\n'.join(('point(%s,%s) # point=x'%(float(ra[i]),float(dec[i])),'# vector(%s,%s,%s",%s)'%(float(ra[i]),float(dec[i]),np.sqrt(pmra[i]**2+pmdec[i]**2)*10,180-math.degrees(math.atan(pmdec[i]/pmra[i]))),'\n')))   
#             print('ssss')
#         elif pmra[0]<0 and pmdec[0]<0:
#             f.write('\n'.join(('point(%s,%s) # point=x'%(float(ra[i]),float(dec[i])),'# vector(%s,%s,%s",%s)'%(float(ra[i]),float(dec[i]),np.sqrt(pmra[i]**2+pmdec[i]**2)*10,180 + math.degrees(math.atan(pmdec[i]/pmra[i]))),'\n')))   
#             print('ssss')
#         elif pmra[0]>0 and pmdec[0]<0:
#             print('ssss')
#             f.write('\n'.join(('point(%s,%s) # point=x'%(float(ra[i]),float(dec[i])),'# vector(%s,%s,%s",%s)'%(float(ra[i]),float(dec[i]),np.sqrt(pmra[i]**2+pmdec[i]**2)*10,(-1)*math.degrees(math.atan(pmdec[i]/pmra[i]))),'\n')))   
# f.close
# 
# =============================================================================

# without arrows
for i in range(len(ra)):
    with open(pruebas+ 'young_bg_stars.reg', 'a') as f:

        f.write('\n'.join(('point(%s,%s) # point=circle'%(float(ra[i]),float(dec[i])),'\n')))   
        print('ssss')    
    f.close

    
    
    
    