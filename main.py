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
list = []
#list of string name for each subpopulation cut

for i in range(0,4):
    for j in range(0,4):
        list += ['4x4cutm' + str(i) + 'v' + str(j)]

####################################################################################
redshifts_array = np.array([[3,15,131],[2,25,197],[1,47,343],[0.5,59,519],[0.42,62,570],[0.31,68,638],[0.2,76,691],[0.12,84,733],[0.06,-1,761]])
#redshifts you would like to collect alignment data for in plots 
#form: [z_redshift,index,snapshotnumber]


smoothness = [2,4,8]
#iterate multipe plots on each figure for different smoothnesses
#of the Universe simulation grid size. 

masscutbool = True
vdcutbool = True
#false ensures full population, make true with empty arrays for progen populations
level_lock = True#level_lock index taken to be 19
#Boolean flags for certain functionalities
#masscutbool: choose whether or not to split population into mass subpopulations 
#vdcutbool: choose whether or not to split population into velocity disperion subpopulations 
#level_lock: choose whether to only include level 1 galaxies in z = 0 snapshot

masscut_array = [[9,9.75],[9.75,10.5],[10.5,11.25],[11.25,12],[12,12.75]]
#masscut_array = []
vdcut_array = [[0,0.4],[0.4,0.7],[0.7,1],[1,2]]
#vdcut_array = []

#The galaxy population can be split into a subpopulation over a Velocity Dispersion/mass grid. 
#Plots are then done iteratively over this grid. 
#The grid is defined by masscut_array and vdcut_array

bin_name = 'fullpop'
#label specifier to use when saving files and data

#redshiftplot_bins = np.array(['fullpop'])
redshiftplot_bins = np.array(list)
#list of strings to label each grid bin. used for saving data and images. 


name_761_file = 'sample_full_761_allgal_wspin.dat'
#name of z = 0 file from which progenitors are found

redshiftplot_masscut_array = [[0,9.4],[9.4,10],[10,10.6],[10.6,15]]
redshiftplot_vdcut_array = [[0,0.4],[0.4,0.7],[0.7,1],[1,2]]
#A slight variation on the masscut_arry and vdcut_array, used in the plotting
#functions rather than the data extraction functions.

#parameters that you main asign
#################################################################################

if masscutbool == False and vdcutbool == False:
    cut = False
else:
    cut = True

#################################################################################
#List all the functions you want to use here
#Note: the tidalfield2 function takes the most time as it has to read though all data 
#in the full tidal file and so shouldn't be run unnecessarily


merger_code_4.merger_code_4(name_761_file,masscut_array = masscut_array,vdcut_array = vdcut_array,redshifts_array = redshifts_array,bin_name = bin_name,level_lock = level_lock)
TidalField2.TidalField2(masscut_array,vdcut_array,redshifts_array[:,2].astype(int).tolist(),smoothness,name_761_file,cutname = bin_name,masscutbool = masscutbool,vdcutbool = vdcutbool,level_lock = level_lock,wspin = False)
Tidal_Plots_2.Tidal_Plots_2(redshifts_array,smoothness,cut,name_761_file,cutname = bin_name,level_lock = level_lock)

z_avangle_plot.z_avangle_plot(redshifts_array[:,[0,2]],smoothness,redshiftplot_bins,np.array([1,2,3]),False,'','major1',True,masscutarray = redshiftplot_masscut_array,vsigcutarray=redshiftplot_vdcut_array,interpolater=False)
z_avangle_plot.z_avangle_plot(redshifts_array[:,[0,2]],smoothness,redshiftplot_bins,np.array([2]),True,'','major2',True,masscutarray = redshiftplot_masscut_array,vsigcutarray=redshiftplot_vdcut_array,interpolater=False)
z_avangle_plot.z_avangle_plot(redshifts_array[:,[0,2]],smoothness,redshiftplot_bins,np.array([3]),True,'','major3',True,masscutarray = redshiftplot_masscut_array,vsigcutarray=redshiftplot_vdcut_array,interpolater=False)
