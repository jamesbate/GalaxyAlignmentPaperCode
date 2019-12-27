import numpy as np
from scipy.io import FortranFile
import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from collections import Counter
import time
import pandas as pd
import sys

import merger_code_4
import TidalField2
import Tidal_Plots_2
import z_avangle_plot


#general notes: snapshots must be in numerical order
#761 galaxy file must always be included


#redshifts_array = np.array([[3,15,131],[2,25,197],[1,47,343],[0.64,52,406],[0.42,62,570],[0.2,76,691],[0.06,-1,761]])
#masscut_array = [[9,9.75],[9.75,10.5],[10.5,11.25],[11.25,12],[12,12.75]]
#redshifts_array = np.array([[3,15,131],[2,25,197],[1,47,343],[0.64,52,406],[0.42,62,570],[0.31,68,638],[0.2,76,691],[0.12,84,733],[0.06,-1,761]])
list = []

for i in range(0,4):
    for j in range(0,4):
        list += ['4x4cutm' + str(i) + 'v' + str(j)]

#for i in range(0,3):
#    for j in range(0,3):
#        list += ['newcutm' + str(i) + 'v' + str(j)]

#for i in range(0,3):
#    list += ['mscutm' + str(i) + 'v0.6']

##################################################################
#you may assign the following variables:

#redshifts_array = np.array([[0.06,-1,761]])
#redshifts_array = np.array([[3,15,131],[2,25,197],[1,47,343],[0.31,68,638],[0.12,84,733],[0.06,-1,761]])#
#[0.64,52,406]

redshifts_array = np.array([[3,15,131],[2,25,197],[1,47,343],[0.5,59,519],[0.42,62,570],[0.31,68,638],[0.2,76,691],[0.12,84,733],[0.06,-1,761]])

#redshifts_array = np.array([[3,15,131],[1,47,343],[0.06,-1,761]])
smoothness = [2]
masscutbool = True
vdcutbool = True
#false ensures full population, make true with empty arrays for progen populaions

level_lock = True#level_lock index taken to be 19
masscut_array = [10.4,15]
#masscut_array = []
vdcut_array = [0,0.6]
#vdcut_array = []

bin_name = 'fullpop'

#redshiftplot_bins = np.array(['fullpop'])
redshiftplot_bins = np.array(list)

#redshiftplot_bins = np.array(['vsigrepromcsm'])

name_761_file = 'sample_full_761_allgal_wspin.dat'
#name_761_file = 'sample_full_761_lvlone.dat'

#redshiftplot_masscut_array = [[10.4,10.7],[10.7,11.06],[11.06,15]]
#redshiftplot_vdcut_array = [[0,0.6]]


#redshiftplot_masscut_array = [[0,9.5],[9.5,10.6],[10.6,15]]
#redshiftplot_vdcut_array = [[0,0.6],[0.6,0.8],[0.8,2]]
redshiftplot_masscut_array = [[0,9.4],[9.4,10],[10,10.6],[10.6,15]]
redshiftplot_vdcut_array = [[0,0.4],[0.4,0.7],[0.7,1],[1,2]]


#######################################################################

if masscutbool == False and vdcutbool == False:
    cut = False
else:
    cut = True

#List all the functions you want to use
#Note: the tidalfield2 function takes the most time and so shouldn't be run unnecessarily
#merger_code_4.merger_code_4(name_761_file,masscut_array = masscut_array,vdcut_array = vdcut_array,redshifts_array = redshifts_array,bin_name = bin_name,level_lock = level_lock)
#TidalField2.TidalField2(masscut_array,vdcut_array,redshifts_array[:,2].astype(int).tolist(),smoothness,name_761_file,cutname = bin_name,masscutbool = masscutbool,vdcutbool = vdcutbool,level_lock = level_lock,wspin = False)
#Tidal_Plots_2.Tidal_Plots_2(redshifts_array,smoothness,cut,name_761_file,cutname = bin_name,level_lock = level_lock)
#z_avangle_plot.z_avangle_plot(redshifts_array[:,[0,2]],smoothness,redshiftplot_bins,np.array([1,2,3]),True,'Full Population','major',True,masscutarray = redshiftplot_masscut_array,vsigcutarray=redshiftplot_vdcut_array,interpolater=False)
#z_avangle_plot.z_avangle_plot(redshifts_array[:,[0,2]],smoothness,redshiftplot_bins,np.array([4,5,6]),True,'Full Population','minor',True,masscutarray = redshiftplot_masscut_array,vsigcutarray=redshiftplot_vdcut_array,interpolater=False)
#z_avangle_plot.z_avangle_plot(redshifts_array[:,[0,2]],smoothness,redshiftplot_bins,np.array([7,8,9]),True,'Full Population','spin',True,masscutarray = redshiftplot_masscut_array,vsigcutarray=redshiftplot_vdcut_array,interpolater=False)

z_avangle_plot.z_avangle_plot(redshifts_array[:,[0,2]],smoothness,redshiftplot_bins,np.array([1,2,3]),False,'','major1',True,masscutarray = redshiftplot_masscut_array,vsigcutarray=redshiftplot_vdcut_array,interpolater=False)
#z_avangle_plot.z_avangle_plot(redshifts_array[:,[0,2]],smoothness,redshiftplot_bins,np.array([2]),True,'','major2',True,masscutarray = redshiftplot_masscut_array,vsigcutarray=redshiftplot_vdcut_array,interpolater=False)
#z_avangle_plot.z_avangle_plot(redshifts_array[:,[0,2]],smoothness,redshiftplot_bins,np.array([3]),True,'','major3',True,masscutarray = redshiftplot_masscut_array,vsigcutarray=redshiftplot_vdcut_array,interpolater=False)
