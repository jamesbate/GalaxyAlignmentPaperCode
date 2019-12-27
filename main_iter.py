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


#general notes: snapshots must be in numerical order
#761 galaxy file must always be included



#redshifts_array = np.array([[3,15,131],[2,25,197],[1,47,343],[0.64,52,406],[0.42,62,570],[0.2,76,691],[0.06,-1,761]])
#masscut_array = [[9,9.75],[9.75,10.5],[10.5,11.25],[11.25,12],[12,12.75]]


###################################################################################
#you may assign the following variables:
#redshifts_array = np.array([[3,15,131],[2,25,197],[1,47,343],[0.64,52,406],[0.42,62,570],[0.31,68,638],[0.2,76,691],[0.12,84,733],[0.06,-1,761]])
#redshifts_array = np.array([[3,15,131],[2,25,197],[1,47,343],[0.31,68,638],[0.12,84,733],[0.06,-1,761]])#
redshifts_array = np.array([[3,15,131],[2,25,197],[1,47,343],[0.5,59,519],[0.42,62,570],[0.31,68,638],[0.2,76,691],[0.12,84,733],[0.06,-1,761]])
#redshifts_array = np.array([[0.06,-1,761]])
smoothness = [2]
masscutbool = True
vdcutbool = True
level_lock = True

name_761_file = 'sample_full_761_allgal_wspin.dat'
####################################################################################


masscut_array = [float(sys.argv[1]),float(sys.argv[2])]
vdcut_array = [float(sys.argv[3]),float(sys.argv[4])]

bin_name = sys.argv[5]


if masscutbool == False and vdcutbool == False:
    cut = False
else:
    cut = True

#List all the functions you want to use
#Note: the tidalfield2 function takes the most time and so shouldn't be run unnecessarily


merger_code_4.merger_code_4(name_761_file,masscut_array = masscut_array,vdcut_array = vdcut_array,redshifts_array = redshifts_array,bin_name = bin_name,level_lock = level_lock)
TidalField2.TidalField2(masscut_array,vdcut_array,redshifts_array[:,2].astype(int).tolist(),smoothness,name_761_file,cutname = bin_name,masscutbool = masscutbool,vdcutbool = vdcutbool,level_lock = level_lock,wspin = False)
Tidal_Plots_2.Tidal_Plots_2(redshifts_array,smoothness,cut,name_761_file,cutname = bin_name,level_lock = level_lock)
#z_avangle_plot.z_avangle_plot(redshifts_array[:,[0,2]],smoothness,redshiftplot_bins,np.array([1,3]),True,'$\phi$ test','test',True)
