import numpy as np
from scipy.io import FortranFile
import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

if not os.path.exists('Progen_Plots_2'):
    os.makedirs('Progen_Plots_2')
os.chdir('Progen_Plots_2')
#make directory for files and images

mass_cut = False
mass_split_3 = False
mass_split_2 = False
mass_split_1 = False

if mass_cut == True:
    file_string = 'masscut_'
    l_masscut = 10.4
    u_masscut = 20
elif mass_split_3 == True:
    file_string = 'masssplit3_'
    l_masscut = 11.06
    u_masscut = 20
elif mass_split_2 == True:
    file_string = 'masssplit2_'
    u_masscut = 11.06
    l_masscut = 10.7
elif mass_split_1 == True:
    file_string = 'masssplit1_'
    u_masscut = 10.7
    l_masscut = 10.4
else:
    file_string = ''
#does mass cut and saves new mass cut galaxy file

galaxy_761_file = open('/mnt/zfsusers/orie3568/Galaxy_Data/sample_full_761_allgal.dat','r')


merger_file = FortranFile('/mnt/zfsusers/orie3568/merger_tree_GM.dat','r')
#create all file objects

#merger_array = merger_file.read_record(dtype=np.float32)
#id_cat = merger_array.reshape((2,91,126361))
#id_cat[1,90,:] = np.arange(0,126361)
#id_cat[0,90,:] = 761
#read in merger file

#don't think I even need this merger file?

galaxy_761_data_list = galaxy_761_file.readlines()
galaxy_761_list_length = len(galaxy_761_data_list)
galaxy_761_data_array = np.zeros([galaxy_761_list_length,16])

for j,line in enumerate(galaxy_761_data_list,0):
    galaxy_761_data_array[j] = np.array([float(x) for x in line.split(" ")])
#read list into array

if mass_cut == True or mass_split_3 == True or mass_split_2 == True or mass_split_1 == True:
    mask = np.logical_and(galaxy_761_data_array[:,6] > l_masscut,galaxy_761_data_array[:,6] < u_masscut)
    galaxy_761_data_array = galaxy_761_data_array[mask]
    np.savetxt('/mnt/zfsusers/orie3568/Galaxy_Data/sample_full_'+file_string+'761.dat',galaxy_761_data_array, delimiter=" ")


mask = galaxy_761_data_array[:,4] < 0.55
galaxy_761_data_array = galaxy_761_data_array[mask]
#TEMPORARY

fig,ax = plt.subplots()
plt.subplots_adjust(left = 0.15,right = 0.85,top = 0.85, bottom = 0.15)

bin_size = 0.1

min = np.amin(galaxy_761_data_array[:,6])//bin_size*bin_size
max = np.amax(galaxy_761_data_array[:,6])//bin_size*bin_size + bin_size
n1, bins1, _ = plt.hist(galaxy_761_data_array[:,6], bins=np.arange(min, max,bin_size))
#plot all masses
plt.xlabel(r'$\log_{10}(M/{\rm M}_{\odot})$', fontsize = 17)
plt.axvline(x=10.4,color = 'k', linestyle = ':')
plt.xticks([9,10,11,12])
plt.yticks([1000,2000,3000,4000,5000,6000])
plt.ylabel('Galaxy Number',fontsize = 17)
ax.tick_params(axis='both', labelsize = 12, length = 5,width = 1)
plt.savefig('Mass Hist 761'+file_string+'.png')
plt.clf()
bin_size = 0.5
min = np.amin(galaxy_761_data_array[:,5])//bin_size*bin_size
max = np.amax(galaxy_761_data_array[:,5])//bin_size*bin_size + bin_size
n1, bins1, _ = plt.hist(galaxy_761_data_array[:,5], bins=np.arange(min,max, bin_size))
#plot all magnitudes
plt.savefig('Mag Hist 761'+file_string+'.png')
plt.clf()
bin_size = 0.05
min = np.amin(galaxy_761_data_array[:,4])//bin_size*bin_size
max = np.amax(galaxy_761_data_array[:,4])//bin_size*bin_size + bin_size
n1, bins1, _ = plt.hist(galaxy_761_data_array[:,4], bins=np.arange(min, max, bin_size))
#plot all masses
plt.savefig('Disp Hist 761'+file_string+'.png')
plt.clf()


bin_size = 0.1
min = np.amin(galaxy_761_data_array[:,7]/galaxy_761_data_array[:,9])//bin_size*bin_size
max = np.amax(galaxy_761_data_array[:,7]/galaxy_761_data_array[:,9])//bin_size*bin_size + bin_size
n1, bins1, _ = plt.hist(galaxy_761_data_array[:,7]/galaxy_761_data_array[:,9], bins=np.arange(min, max, bin_size))
#plot all masses
plt.savefig('c_over_a_761_'+file_string+'.png')
plt.clf()
bin_size = 0.1
min = np.amin(galaxy_761_data_array[:,8]/galaxy_761_data_array[:,9])//bin_size*bin_size
max = np.amax(galaxy_761_data_array[:,8]/galaxy_761_data_array[:,9])//bin_size*bin_size + bin_size
n1, bins1, _ = plt.hist(galaxy_761_data_array[:,8]/galaxy_761_data_array[:,9], bins=np.arange(min, max, bin_size))
#plot all masses
plt.savefig('b_over_a_761_'+file_string+'.png')
plt.clf()
