import numpy as np
from scipy.io import FortranFile
import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from collections import Counter

#change how images are saved for masscut

def hist_plot(galaxy,data,progen_data,bin_size,xlab):
    #normalise
    #descriptive labels

    min = np.amin(data)//bin_size*bin_size
    max = np.amax(data)//bin_size*bin_size
    if np.amax(data)%bin_size != 0:
        max += bin_size

    if all(np.array(-1*data) == np.ones(data.size)):
        print 'Histogram Plot Error:\tInvalid data for '+xlab+str(i)
        return np.array([0])
    #if all data -1, call error

    n1, bins1, _ = plt.hist(data, bins=np.arange(min, max,bin_size))

    if n1[n1.size-1] == 0:
        n1 = np.delete(n1,n1.size-1)
    #sometimes last element gives 0
    #this removes this 0

    if xlab == 'Mass':
        plt.xlabel(r'$log(M/M_{sun})$')
    else:
        plt.xlabel(r'$'+xlab+'$')
    plt.ylabel(r'Number of Galaxies')
    plt.title(r'Histrogram Plot of '+xlab)
    #plt.grid(True)

    if xlab == 'b/a':
        plt.savefig('b_over_a'+'_hist_'+file_string+str(galaxy.snapshot)+'.png')
    elif xlab == 'c/a':
        plt.savefig('c_over_a'+'_hist_'+file_string+str(galaxy.snapshot)+'.png')
    else:
        plt.savefig(xlab+'_hist_'+file_string+str(galaxy.snapshot)+'.png')


    plt.clf()
    return np.array(n1)

def ratio_plot(galaxy,n1,data,progen_data,bin_size,lab):

    min = np.amin(data)//bin_size*bin_size
    max = np.amax(data)//bin_size*bin_size
    if np.amax(data)%bin_size != 0:
        max += bin_size

    if n1.size == 1 and n1 == [0]:
        print 'Ratio Plot Error:\tInvalid data for '+lab+str(i)
        return [np.array([0]),np.array([0])]
    #ensures error if [0]

    n2, bins2,bin_num = plt.hist(progen_data, bins=np.arange(min, max,bin_size),color = 'r')

    if lab == 'b/a':
        plt.savefig('b_over_a'+'_hist_progen_plot_'+file_string+str(galaxy.snapshot)+'.png')
    elif lab == 'c/a':
        plt.savefig('c_over_a'+'_hist_progen_plot_'+file_string+str(galaxy.snapshot)+'.png')
    else:
        plt.savefig(str(lab)+'_hist_progen_plot_'+file_string+str(galaxy.snapshot)+'.png')

    n2 = n2[0:n1.size]
    bins2 = bins2[0:n1.size+1]

    yerror = np.sqrt(n2*(1/n1)**2 + n1*(n2/n1**2)**2)

    plt.clf()

    plt.xlabel(lab)
    plt.ylabel('Galaxy Number Ratio')
    plt.title(r'Histrogram Plot of '+lab)
    plt.grid(True)

    bin_centers = 0.5*(bins2[:-1] + bins2[1:])
    #bins for both plots should be the same

    if bin_centers.any() != True or n2.any() != True or n1.any() != True:
        print 'Ratio Plot Error:\tInvalid data for '+lab+str(i)
        return [np.array([0]),np.array([0])]
    else:
        plt.plot(bin_centers,n2/n1)
        plt.errorbar(bin_centers,n2/n1,yerr = yerror,fmt='bo',ms = 4, capthick=1,capsize = 8,elinewidth = 1)
        if lab == 'b/a':
            plt.savefig('b_over_a'+'_ratio_plot_'+file_string+str(galaxy.snapshot)+'.png')
        elif lab == 'c/a':
            plt.savefig('c_over_a'+'_ratio_plot_'+file_string+str(galaxy.snapshot)+'.png')
        else:
            plt.savefig(str(lab)+'_ratio_plot_'+file_string+str(galaxy.snapshot)+'.png')
        plt.clf()
        return np.array([n2,bin_centers])
        #plot mass ratios
#add errors

def joint_ratio_plot(lab):
    if lab == 'Mass':
        n = 0
    elif lab == 'Magnitude':
        n = 1
    elif lab == 'Velocity Dispersion':
        n = 2
    elif lab == 'c/a':
        n = 3
    elif lab == 'b/a':
        n = 4

    a = []
    colours = ['m','c','b','g','r']

    fig,ax = plt.subplots()
    plt.subplots_adjust(left = 0.15,right = 0.85,top = 0.85, bottom = 0.15)

    for i in range(0,snapshot_num):
        if galaxy_z[i].n2[n].size <= 1:
            continue
        if lab == 'Velocity Dispersion':
            if galaxy_z[i].bin_centers[n][-1] > 1.75:
                temp = plt.plot(galaxy_z[i].bin_centers[n][:-1],galaxy_z[i].n2[:][n][:-1]/galaxy_z[i].n1[:][n][:-1],str(colours[i])+'--',label = 'z = '+ str(redshifts[i]))
                #plt.scatter(galaxy_z[i].bin_centers[n][:-1],galaxy_z[i].n2[:][n][:-1]/galaxy_z[i].n1[:][n][:-1],s = 6)
                n2 = galaxy_z[i].n2[:][n][:-1]
                n1 = galaxy_z[i].n1[:][n][:-1]
                yerror = np.sqrt(n2*(1/n1)**2 + n1*(n2/n1**2)**2)
                plt.errorbar(galaxy_z[i].bin_centers[n][:-1],galaxy_z[i].n2[:][n][:-1]/galaxy_z[i].n1[:][n][:-1],fmt = str(colours[i])+'o',yerr = yerror,ms = 4, capthick=1,capsize = 8,elinewidth = 1)
            else:
                temp = plt.plot(galaxy_z[i].bin_centers[n],galaxy_z[i].n2[:][n]/galaxy_z[i].n1[:][n],str(colours[i])+'--',label = 'z = '+ str(redshifts[i]))
                #plt.scatter(galaxy_z[i].bin_centers[n],galaxy_z[i].n2[:][n]/galaxy_z[i].n1[:][n],s = 6)
                n2 = galaxy_z[i].n2[:][n]
                n1 = galaxy_z[i].n1[:][n]
                yerror = np.sqrt(n2*(1/n1)**2 + n1*(n2/n1**2)**2)
                plt.errorbar(galaxy_z[i].bin_centers[n][:],galaxy_z[i].n2[:][n][:]/galaxy_z[i].n1[:][n][:],fmt = str(colours[i])+'o',yerr = yerror,ms = 4, capthick=1,capsize = 8,elinewidth = 1)
            #TEMPORARY
        else:
            temp = plt.plot(galaxy_z[i].bin_centers[n],galaxy_z[i].n2[:][n]/galaxy_z[i].n1[:][n],str(colours[i])+'--',label = 'z = '+ str(redshifts[i]))
            #plt.scatter(galaxy_z[i].bin_centers[n],galaxy_z[i].n2[:][n]/galaxy_z[i].n1[:][n],s = 6)
            n2 = galaxy_z[i].n2[:][n]
            n1 = galaxy_z[i].n1[:][n]
            yerror = np.sqrt(n2*(1/n1)**2 + n1*(n2/n1**2)**2)
            plt.errorbar(galaxy_z[i].bin_centers[n][:],galaxy_z[i].n2[:][n][:]/galaxy_z[i].n1[:][n][:],fmt = str(colours[i])+'o',yerr = yerror,ms = 4, capthick=1,capsize = 8,elinewidth = 1)

        a.append(temp[0])

    ax.tick_params(axis='both', labelsize = 12 , length = 5,width = 1)

    plt.xlabel(lab)
    if lab == 'Mass':
        plt.xlabel(r'$\log_{10}(M/{\rm M}_{\odot})$',fontsize = 17)
        plt.axvline(x=10.4,color = 'k', linestyle = ':')
    if lab == 'Velocity Dispersion':
        plt.xlabel(r'$V_{\theta}/{\sigma}$',fontsize = 17)

    plt.ylabel(r'$\frac{Number of Progenitors}{Number of Galaxies}$',fontsize = 17)
    #plt.title(r'Histrogram Ratio Plot of '+lab)
    #plt.grid(True)
    plt.legend(handles = [i for i in a])

    if lab == 'c/a':
        plt.savefig('ratio_plot_joint_'+file_string+'_c_over_a_'+'.png')
    elif lab == 'b/a':
        plt.savefig('ratio_plot_joint_'+file_string+'_b_over_a_'+'.png')
    else:
        plt.savefig('ratio_plot_joint_'+file_string+str(lab)+'.png')
    #can you save with space?
    plt.clf()
    return

mass_cut = True
#set to true if mass cut desired

level_lock = True
#set to true if only level 1 galaxies wanted

mass_split_1 = False
#10.4 to 10.7
mass_split_2 = False
#10.7 to 11.06
mass_split_3 =  False
#11.06 and above


if mass_cut == True:
    file_string = 'masscut_'
    cut_lower = 10.4
    cut_upper = 15
elif mass_split_1 == True:
    file_string = 'masssplit1_'
    cut_lower = 10.4
    cut_upper = 10.7
elif mass_split_2 == True:
    file_string = 'masssplit2_'
    cut_lower = 10.7
    cut_upper = 11.6
elif mass_split_3 == True:
    file_string = 'masssplit3_'
    cut_lower = 11.06
    cut_upper = 15
else:
    file_string = ''
if level_lock == True:
    lvl = 1
    file_string += 'lvlone_'


if not os.path.exists('/mnt/zfsusers/orie3568/Progen_Plots_4_r'):
    os.makedirs('/mnt/zfsusers/orie3568/Progen_Plots_4_r')
os.chdir('/mnt/zfsusers/orie3568/Progen_Plots_4_r')
#make directory for files and images

class galaxy:
    def __init__(self,label,snapshot):
        self.snapshot = snapshot
        self.label = label
        return

snapshots = [[15,131],[25,197],[47,343],[59,519],[84,733]]
redshifts = [3,2,1,0.5,0.12]
#snapshots = [[15,131],[25,197],[47,343],[52,406],[59,519],[62,570],[68,638],[76,691],[84,733]]
#redshifts = [3,2,1,0.64,0.5,0.42,0.31,0.2,0.12]
#redshifts = [10]
#snapshots = [[52,406]]
#specify number of snapshots required

snapshot_num = len(snapshots)

galaxy_761 = galaxy(90,761)
galaxy_z_file = []
galaxy_z =[]
file_write=[]
#initialisers

for i in range(0,snapshot_num):
    galaxy_z.append(galaxy(snapshots[i][0],snapshots[i][1]))
    galaxy_z_file.append(open('/mnt/zfsusers/orie3568/Galaxy_Data/sample_full_'+str(snapshots[i][1])+'.dat','r'))

    file_write.append(open('progen_'+file_string+'761_'+str(snapshots[i][1])+'.dat','w'))#relative to current directory
galaxy_761_file = open('/mnt/zfsusers/orie3568/Galaxy_Data/sample_full_761_wlevel.dat','r')
merger_file = FortranFile('/mnt/zfsusers/orie3568/merger_tree_GM.dat','r')
#create all file objects

galaxy_761.data_list = galaxy_761_file.readlines()
galaxy_761.list_length = len(galaxy_761.data_list)
galaxy_761.data_array = np.zeros([galaxy_761.list_length,17])
#read in 761 file and prime array

for i,line in enumerate(galaxy_761.data_list,0):
    galaxy_761.data_array[i] = np.array([float(x) for x in line.split(" ")])
#read in galaxy files

merger_array = merger_file.read_record(dtype=np.float32)
id_cat = merger_array.reshape((2,91,126361))
id_cat[1,90,:] = np.arange(0,126361)
id_cat[0,90,:] = 761
#read in merger file


for i in range(0,snapshot_num):

    print 'Plotting Snapshot ' + str(galaxy_z[i].snapshot)

    galaxy_z[i].n1 = []
    galaxy_z[i].n2 = []
    galaxy_z[i].n3 = []
    galaxy_z[i].bin_size = np.zeros(5)
    galaxy_z[i].bin_centers = []
    #initialisers

    galaxy_z[i].data_list = galaxy_z_file[i].readlines()
    galaxy_z[i].list_length = len(galaxy_z[i].data_list)
    galaxy_z[i].data_array = np.zeros([galaxy_z[i].list_length,16])
    #read in current galaxy file and prime array

    galaxy_z[i].progen_masses = np.zeros(galaxy_761.list_length+1)
    galaxy_z[i].progen_mag = np.zeros(galaxy_761.list_length+1)
    galaxy_z[i].progen_disp = np.zeros(galaxy_761.list_length+1)
    galaxy_z[i].progen_ax = np.zeros([galaxy_761.list_length+1,3])
    #prime desired plots

    for j,line in enumerate(galaxy_z[i].data_list,0):
        galaxy_z[i].data_array[j] = np.array([float(x) for x in line.split(" ")])
    #read galaxy list into array

    for j,index in enumerate(galaxy_761.data_array[:,0],0):
        index = int(index)
        progen_index = int(id_cat[1,galaxy_z[i].label,index-1])
        if (progen_index == -1) or (int(id_cat[0,galaxy_z[i].label,index-1]) == 0):
            continue
        #no progenitor in this snap shot
        line, = np.where(galaxy_z[i].data_array[:,0] == int(progen_index+1))
        if line.tolist() == []:
            continue
        #progenitor not included in file

        if mass_cut == True or mass_split_1 == True or mass_split_2 == True or mass_split_3 == True:
            if not cut_lower < galaxy_761.data_array[j,6] < cut_upper:
                continue
        if level_lock == True:
            if galaxy_761.data_array[j,16] != lvl:
                continue
        #its fine to not fill this plot array because zeros are shaved off afterwards

        galaxy_z[i].progen_masses[j] = galaxy_z[i].data_array[line.tolist(),6]
        galaxy_z[i].progen_mag[j] = galaxy_z[i].data_array[line.tolist(),5]
        galaxy_z[i].progen_disp[j] = galaxy_z[i].data_array[line.tolist(),4]
        for z in range(0,3):
            galaxy_z[i].progen_ax[j][z] = galaxy_z[i].data_array[line.tolist(),z+7]#7,8,9

        file_write[i].write('{} {}\n'.format(str(index),str(progen_index+1)))
    #create progen file and fill plot arrays by going through galaxies
    #in galaxy file and checking for progen in progen file

    galaxy_z[i].masses_plot = np.zeros(np.count_nonzero(galaxy_z[i].progen_masses))
    galaxy_z[i].mag_plot = np.zeros(np.count_nonzero(galaxy_z[i].progen_mag))
    galaxy_z[i].disp_plot = np.zeros(np.count_nonzero(galaxy_z[i].progen_disp))
    galaxy_z[i].ax_plot = np.zeros([np.count_nonzero(galaxy_z[i].progen_ax[:,z]),3])
    galaxy_z[i].masses_plot = [x for x in galaxy_z[i].progen_masses if x != 0]
    galaxy_z[i].mag_plot = [x for x in galaxy_z[i].progen_mag if x != 0]
    galaxy_z[i].disp_plot = [x for x in galaxy_z[i].progen_disp if x != 0]
    for z in range(0,3):
        galaxy_z[i].ax_plot[:,z] = [x for x in galaxy_z[i].progen_ax[:,z] if x != 0]
    #remove all zeros



    galaxy_z[i].bin_size[0] = 0.5
    galaxy_z[i].n1.append(hist_plot(galaxy_z[i],galaxy_z[i].data_array[:,6],galaxy_z[i].masses_plot,galaxy_z[i].bin_size[0],'Mass'))

    temp = ratio_plot(galaxy_z[i],galaxy_z[i].n1[0],galaxy_z[i].data_array[:,6],galaxy_z[i].masses_plot,galaxy_z[i].bin_size[0],'Mass')
    galaxy_z[i].n2.append(temp[0])
    galaxy_z[i].bin_centers.append(temp[1])

    galaxy_z[i].bin_size[1] = 0.5
    galaxy_z[i].n1.append(hist_plot(galaxy_z[i],galaxy_z[i].data_array[:,5],galaxy_z[i].mag_plot,galaxy_z[i].bin_size[1],'Magnitude'))

    temp = ratio_plot(galaxy_z[i],galaxy_z[i].n1[1],galaxy_z[i].data_array[:,5],galaxy_z[i].mag_plot,galaxy_z[i].bin_size[1],'Magnitude')
    galaxy_z[i].n2.append(temp[0])##!!!!!
    galaxy_z[i].bin_centers.append(temp[1])


    galaxy_z[i].bin_size[2] = 0.25
    galaxy_z[i].n1.append(hist_plot(galaxy_z[i],galaxy_z[i].data_array[:,4],galaxy_z[i].disp_plot,galaxy_z[i].bin_size[2],'Velocity Dispersion'))

    temp = ratio_plot(galaxy_z[i],galaxy_z[i].n1[2],galaxy_z[i].data_array[:,4],galaxy_z[i].disp_plot,galaxy_z[i].bin_size[2],'Velocity Dispersion')
    galaxy_z[i].n2.append(temp[0])
    galaxy_z[i].bin_centers.append(temp[1])

    #axis ratio

    galaxy_z[i].bin_size[3] = 0.1
    galaxy_z[i].n1.append(hist_plot(galaxy_z[i],galaxy_z[i].data_array[:,7]/galaxy_z[i].data_array[:,9],galaxy_z[i].ax_plot[:,0]/galaxy_z[i].ax_plot[:,2],galaxy_z[i].bin_size[3],'c/a'))

    temp = ratio_plot(galaxy_z[i],galaxy_z[i].n1[3],galaxy_z[i].data_array[:,7]/galaxy_z[i].data_array[:,9],galaxy_z[i].ax_plot[:,0]/galaxy_z[i].ax_plot[:,2],galaxy_z[i].bin_size[3],'c/a')
    galaxy_z[i].n2.append(temp[0])
    galaxy_z[i].bin_centers.append(temp[1])

    galaxy_z[i].bin_size[4] = 0.1
    galaxy_z[i].n1.append(hist_plot(galaxy_z[i],galaxy_z[i].data_array[:,8]/galaxy_z[i].data_array[:,9],galaxy_z[i].ax_plot[:,1]/galaxy_z[i].ax_plot[:,2],galaxy_z[i].bin_size[4],'b/a'))

    temp = ratio_plot(galaxy_z[i],galaxy_z[i].n1[4],galaxy_z[i].data_array[:,8]/galaxy_z[i].data_array[:,9],galaxy_z[i].ax_plot[:,1]/galaxy_z[i].ax_plot[:,2],galaxy_z[i].bin_size[4],'b/a')
    galaxy_z[i].n2.append(temp[0])
    galaxy_z[i].bin_centers.append(temp[1])

    #do all of the plots

joint_ratio_plot('Mass')
joint_ratio_plot('Magnitude')
joint_ratio_plot('Velocity Dispersion')

joint_ratio_plot('c/a')
joint_ratio_plot('b/a')
#need to make plots look prettier

galaxy_761_file.close()
for i in range(0,snapshot_num/2):
    galaxy_z_file[i].close()
    file_write[i].close()
merger_file.close()
#close all relavent files
