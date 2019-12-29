import numpy as np

import merger_code_4
import TidalField2
import Tidal_Plots_2
import z_avangle_plot
#preamble
###################################################################################

#general notes for use: redshift snapshots must be in numerical order
#761 galaxy file must always be included

###################################################################################
redshifts_array = np.array([[3,15,131],[2,25,197],[1,47,343],[0.5,59,519],[0.42,62,570],[0.31,68,638],[0.2,76,691],[0.12,84,733],[0.06,-1,761]])
#redshifts you would like to collect alignment data for in plots 
#form: [z_redshift,index,snapshotnumber]

smoothness = [2,4,8]
#iterate multipe plots on each figure for different smoothnesses
#of the Universe simulation grid size. 

masscutbool = True
vdcutbool = True
level_lock = True#level_lock index taken to be 19
#Boolean flags for certain functionalities
#masscutbool: choose whether or not to split population into mass subpopulations 
#vdcutbool: choose whether or not to split population into velocity disperion subpopulations 
#level_lock: choose whether to only include level 1 galaxies in z = 0 snapshot

name_761_file = 'sample_full_761_allgal_wspin.dat'
#name of z = 0 file from which progenitors are found

#parameters that can be assigned 
####################################################################################

masscut_array = [float(sys.argv[1]),float(sys.argv[2])]
vdcut_array = [float(sys.argv[3]),float(sys.argv[4])]
bin_name = sys.argv[5]

#This is the iterative version of the main file. It is designed to be iteratively called 
#from the main_iter_callfunction.py, where each iteration runs all the plots for a different
#subpopulation. The selection is made based on mass and/or velocity dispersion at z = 0.
#Here the particular suspopulation is read in from system arguments.


if masscutbool == False and vdcutbool == False:
    cut = False
else:
    cut = True

#####################################################################################
#List all the functions you want to use here
#Note: the tidalfield2 function takes the most time as it has to read though all data 
#in the full tidal file and so shouldn't be run unnecessarily


merger_code_4.merger_code_4(name_761_file,masscut_array = masscut_array,vdcut_array = vdcut_array,redshifts_array = redshifts_array,bin_name = bin_name,level_lock = level_lock)
TidalField2.TidalField2(masscut_array,vdcut_array,redshifts_array[:,2].astype(int).tolist(),smoothness,name_761_file,cutname = bin_name,masscutbool = masscutbool,vdcutbool = vdcutbool,level_lock = level_lock,wspin = False)
Tidal_Plots_2.Tidal_Plots_2(redshifts_array,smoothness,cut,name_761_file,cutname = bin_name,level_lock = level_lock)
z_avangle_plot.z_avangle_plot(redshifts_array[:,[0,2]],smoothness,redshiftplot_bins,np.array([1,3]),True,'$\phi$ test','test',True)
#example of set of functions to run. 