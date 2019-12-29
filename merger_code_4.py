import numpy as np
from scipy.io import FortranFile
import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
#preamble


def merger_code_4(name_761_file,masscut_array,vdcut_array,redshifts_array,bin_name = '',level_lock=True):

    """
    This function reads in the 761 galaxy file of all the galaxies and the galaxy files at each snapshot. It then applies desired cuts on the 761 population,
    and uses the mergerfile to relate the index of a galaxy at 761 to its progenitors at different
    specified redshifts, producing 'linker files'. ratio plots of the progenitor population to the full population at each redshift are made for various galactic properties.

    The inputs are an array of the desired massbins, [lowercut,uppercut], the bin name (for file saving),the desired redshifts
    [[redshift,mergertreeindex,,snapshot],...], and
    the level_lock (i.e. whether to only select galaxies which are level 1).

    The outputs are progen linker files which list 761 galaxy indices to their progenitors and
    each z, as well as plots of galaxy progenitors at each redshift.

    Note: to not use cuts, use empty lists for cut arrays
    """

    Master_dir = os.getcwd()

    def hist_plot(galaxy,data,progen_data,bin_size,xlab):
        """
        Plots and saves histogram for given data
        """

        #PLEASE NOT THAT THIS PLOTS THE FUNCTION FOR ALL ORIGINAL THE DATA, IT DOES NOT APPLY CUTS
        #progen data included as we are likely to extend this function, although it is currently 
        #not used

        #data bins
        min = np.amin(data)//bin_size*bin_size
        max = np.amax(data)//bin_size*bin_size
        if np.amax(data)%bin_size != 0:
            max += bin_size

        if all(np.array(-1*data) == np.ones(data.size)):
            print 'Histogram Plot Error:\tInvalid data for '+xlab+str(i)
            return np.array([0])
        #if all data -1, raise error

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
        #different labels for different data

        if xlab == 'b/a':
            plt.savefig('b_over_a'+'_hist_'+file_string+str(galaxy.snapshot)+'nocut.png')
            plt.savefig('b_over_a'+'_hist_'+file_string+str(galaxy.snapshot)+'nocut.eps',format='eps', dpi=1000)
        elif xlab == 'c/a':
            plt.savefig('c_over_a'+'_hist_'+file_string+str(galaxy.snapshot)+'nocut.png')
            plt.savefig('c_over_a'+'_hist_'+file_string+str(galaxy.snapshot)+'nocut.eps',format='eps', dpi=1000)
        else:
            plt.savefig(xlab+'_hist_'+file_string+str(galaxy.snapshot)+'nocut.png')
            plt.savefig(xlab+'_hist_'+file_string+str(galaxy.snapshot)+'nocut.eps',format='eps', dpi=1000)
        #saving for all the different types of data

        plt.clf()
        return np.array(n1)

    def ratio_plot(galaxy,n1,data,progen_data,bin_size,lab):
        """
        Plots and saves ratios of the data and progenitor data for 
        each bin.
        """
        min = np.amin(data)//bin_size*bin_size
        max = np.amax(data)//bin_size*bin_size
        if np.amax(data)%bin_size != 0:
            max += bin_size
        #data bin sizes 

        if n1.size == 1 and n1 == [0]:
            print 'Ratio Plot Error:\tInvalid data for '+lab+str(i)
            return [np.array([0]),np.array([0])]
        #ensures error if [0]

        n2, bins2,bin_num = plt.hist(progen_data, bins=np.arange(min, max,bin_size),color = 'r')
        #histogram plot


        if lab == 'b/a':
            plt.savefig('b_over_a'+'_hist_progen_plot_'+file_string+str(galaxy.snapshot)+'.png')

            plt.savefig('b_over_a'+'_hist_progen_plot_'+file_string+str(galaxy.snapshot)+'.eps',format='eps', dpi=1000)
        elif lab == 'c/a':
            plt.savefig('c_over_a'+'_hist_progen_plot_'+file_string+str(galaxy.snapshot)+'.png')
            plt.savefig('c_over_a'+'_hist_progen_plot_'+file_string+str(galaxy.snapshot)+'.eps',format='eps', dpi=1000)
        else:
            plt.savefig(str(lab)+'_hist_progen_plot_'+file_string+str(galaxy.snapshot)+'.png')
            plt.savefig(str(lab)+'_hist_progen_plot_'+file_string+str(galaxy.snapshot)+'.eps',format='eps', dpi=1000)

        #save data

        n2 = n2[0:n1.size]
        bins2 = bins2[0:n1.size+1]
        yerror = np.sqrt(n2*(1/n1)**2 + n1*(n2/n1**2)**2)
        #errorbar values

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
                plt.savefig('b_over_a'+'_ratio_plot_'+file_string+str(galaxy.snapshot)+'.eps',format='eps', dpi=1000)
            elif lab == 'c/a':
                plt.savefig('c_over_a'+'_ratio_plot_'+file_string+str(galaxy.snapshot)+'.png')
                plt.savefig('c_over_a'+'_ratio_plot_'+file_string+str(galaxy.snapshot)+'.eps',format='eps', dpi=1000)
            else:
                plt.savefig(str(lab)+'_ratio_plot_'+file_string+str(galaxy.snapshot)+'.png')
                plt.savefig(str(lab)+'_ratio_plot_'+file_string+str(galaxy.snapshot)+'.eps',format='eps', dpi=1000)
            plt.clf()
            return np.array([n2,bin_centers])
            #plot mass ratios

    def joint_ratio_plot(lab):
        """
        Applies ratio plots to all the data saved by the above two functions simultaneously
        for each redshift
        """
        uncut = False
        #apply no cuts to the data

        if uncut == True:
            uncut_string = 'uncut'
        else:
            uncut_string = ''


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
        #associate data label to data index

        a = []
        #for legend

        fig,ax = plt.subplots()
        plt.subplots_adjust(left = 0.15,right = 0.85,top = 0.85, bottom = 0.15)
        if snapshots[-1][1] == 761:
            snapshot_loop = snapshot_num -1
        else:
            snapshot_loop = snapshot_num
        #ratio plot needs to act slightly different if z = 0 is included

        for i in range(0,snapshot_loop):
            if galaxy_z[i].n2[n].size <= 1:
                continue
            if lab == 'Velocity Dispersion':
                if galaxy_z[i].bin_centers[n][-1] > 1.75:
                    temp = plt.plot(galaxy_z[i].bin_centers[n][:-1],galaxy_z[i].n2[:][n][:-1]/galaxy_z[i].n1[:][n][:-1],colour_dict[str(galaxy_z[i].snapshot)]+'--',label = 'z = '+ str(redshifts[i]))
                    plt.scatter(galaxy_z[i].bin_centers[n][:-1],galaxy_z[i].n2[:][n][:-1]/galaxy_z[i].n1[:][n][:-1],c = colour_dict[str(galaxy_z[i].snapshot)],s = 6)
                else:
                    temp = plt.plot(galaxy_z[i].bin_centers[n],galaxy_z[i].n2[:][n]/galaxy_z[i].n1[:][n],colour_dict[str(galaxy_z[i].snapshot)]+'--',label = 'z = '+ str(redshifts[i]))
                    plt.scatter(galaxy_z[i].bin_centers[n],galaxy_z[i].n2[:][n]/galaxy_z[i].n1[:][n],c = colour_dict[str(galaxy_z[i].snapshot)],s = 6)
                               
            else:
                temp = plt.plot(galaxy_z[i].bin_centers[n],galaxy_z[i].n2[:][n]/galaxy_z[i].n1[:][n],colour_dict[str(galaxy_z[i].snapshot)]+'--',label = 'z = '+ str(redshifts[i]))
                plt.scatter(galaxy_z[i].bin_centers[n],galaxy_z[i].n2[:][n]/galaxy_z[i].n1[:][n],c = colour_dict[str(galaxy_z[i].snapshot)],s = 6)
               
            a.append(temp[0])

        ax.tick_params(axis='both', labelsize = 12 , length = 5,width = 1)

        plt.xlabel(lab)
        if lab == 'Mass':
            plt.xlabel(r'$\log_{10}(M/{\rm M}_{\odot})$',fontsize = 17)
            plt.axvline(x=10.4,color = 'k', linestyle = ':')
        if lab == 'Velocity Dispersion':
            plt.xlabel(r'$V_{\theta}/{\sigma}$',fontsize = 17)

        plt.ylabel(r'$\frac{\rm Number\ of\ main\ progenitors}{\rm Total\ number\ of\ galaxies}$',fontsize = 17)
        plt.legend(handles = [i for i in a],frameon = False)

        if lab == 'c/a':
            plt.savefig('ratio_plot_joint_'+file_string+uncut_string+'_c_over_a_'+'.png')
            plt.savefig('ratio_plot_joint_'+file_string+uncut_string+'_c_over_a_'+'.eps',format='eps', dpi=1000)
        elif lab == 'b/a':
            plt.savefig('ratio_plot_joint_'+file_string+uncut_string+'_b_over_a_'+'.png')
            plt.savefig('ratio_plot_joint_'+file_string+uncut_string+'_b_over_a_'+'.eps',format='eps', dpi=1000)
        else:
            plt.savefig('ratio_plot_joint_'+file_string+uncut_string+str(lab)+'.png')
            plt.savefig('ratio_plot_joint_'+file_string+uncut_string+str(lab)+'.eps',format='eps', dpi=1000)
        plt.clf()
        return

    #These are the function which plot all the galaxy data
#################################################################################################################


    if masscut_array == []:
        mcut = False
    else:
        mcut = True

    if vdcut_array == []:
        vdcut = False
    else:
        vdcut = True
    #create mroe boolean flags to identify different cases

    masscolumns = 6
    vdcolumns = 4
    #index where mass ad vd can be found

    level_lock_index = 20
    column_num = level_lock_index
    #I HAVE ASSUMED HERE THAT THE LEVEL IS ALWAYS THE LAST INDEX

    file_string = bin_name + '_'

    if mcut == True:
        mcut_upper = masscut_array[1]
        mcut_lower = masscut_array[0]
    if level_lock == True:
        lvl = 1
        file_string += 'lvlone_'

    if vdcut == True:
        vdcut_upper = vdcut_array[1]
        vdcut_lower = vdcut_array[0]

    if not os.path.exists(Master_dir + '/Progen_Plots_new/'+bin_name):
        os.makedirs(Master_dir + '/Progen_Plots_new/'+bin_name)
    os.chdir(Master_dir + '/Progen_Plots_new/'+bin_name)
    #make directory for files and images

    class galaxy:
        def __init__(self,label,snapshot):
            self.snapshot = snapshot
            self.label = label
            return
    #each galaxy file has a class associated with it
    colour_dict = {
    "761" : "C0",
    "733" : "C8",
    "691" : "C2",
    "638" : "C3",
    "570" : "C4",
    "519" : "C5",
    "406" : "C6",
    "343" : "C7",
    "197" : "C1",
    "131" : "C9"
    }#standardise colours 


    redshifts = redshifts_array[:,0].tolist()
    snapshots = redshifts_array[:,1:].astype(int).tolist()
    #desired snapshots entered as argument 

    snapshot_num = len(snapshots)

    galaxy_761 = galaxy(90,761)
    galaxy_z_file = []
    galaxy_z =[]
    file_write=[]
    progen_uncut_file_write = []
    #initialisers

    for i in range(0,snapshot_num):
        if snapshots[i][1] == 761:
            galaxy_z.append(galaxy(snapshots[i][0],snapshots[i][1]))
            galaxy_z_file.append(open(Master_dir + '/Galaxy_Data_new/'+name_761_file,'r'))
        else:
            galaxy_z.append(galaxy(snapshots[i][0],snapshots[i][1]))
            #galaxy_z_file.append(open(Master_dir + '/Galaxy_Data_new/sample_full_'+str(snapshots[i][1])+'_allgal_wspin.dat','r'))#CHANGED TO NEW FILE
            galaxy_z_file.append(open(Master_dir + '/Galaxy_Data_new/sample_full_'+str(snapshots[i][1])+'_allgal_wspin.dat','r'))#CHANGED TO NEW FILE
            file_write.append(open('progen_'+file_string+'761_'+str(snapshots[i][1])+'.dat','w'))#relative to current directory
            progen_uncut_file_write.append(open('sample_progen_'+file_string+'_uncut_'+str(snapshots[i][1])+'.dat','w'))#relative to current directory
    galaxy_761_file = open(Master_dir + '/Galaxy_Data_new/'+name_761_file,'r')
    #changed from wlevel!!!
    merger_file = FortranFile('/mnt/zfsusers/orie3568/merger_tree_GM.dat','r')
    #create all file objects

    galaxy_761.data_list = galaxy_761_file.readlines()
    galaxy_761.list_length = len(galaxy_761.data_list)


    galaxy_761.data_array = np.zeros([galaxy_761.list_length,column_num])
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
        #Plots and save all galaxy data


        print 'Plotting Snapshot ' + str(galaxy_z[i].snapshot)

        galaxy_z[i].n1 = []
        galaxy_z[i].n2 = []
        galaxy_z[i].n3 = []#USED FOR UNCUT PLOTS
        galaxy_z[i].bin_size = np.zeros(5)
        galaxy_z[i].bin_centers = []
        #initialisers

        galaxy_z[i].data_list = galaxy_z_file[i].readlines()
        galaxy_z[i].list_length = len(galaxy_z[i].data_list)
        galaxy_z[i].data_array = np.zeros([galaxy_z[i].list_length,column_num])
        #read in current galaxy file and prime array

        galaxy_z[i].progen_masses = np.zeros(galaxy_761.list_length+1)
        galaxy_z[i].progen_mag = np.zeros(galaxy_761.list_length+1)
        galaxy_z[i].progen_disp = np.zeros(galaxy_761.list_length+1)
        galaxy_z[i].progen_ax = np.zeros([galaxy_761.list_length,3])
        #prime desired plots

        galaxy_z[i].progen_masses_uncut = np.zeros(galaxy_761.list_length+1)
        galaxy_z[i].progen_mag_uncut = np.zeros(galaxy_761.list_length+1)
        galaxy_z[i].progen_disp_uncut = np.zeros(galaxy_761.list_length+1)
        galaxy_z[i].progen_ax_uncut = np.zeros([galaxy_761.list_length,3])
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


            #I DECIDED TO ADD UNCUT PROGEN PLOTS BY DEFAULT

            galaxy_z[i].progen_masses_uncut[j] = galaxy_z[i].data_array[line.tolist(),6]
            galaxy_z[i].progen_mag_uncut[j] = galaxy_z[i].data_array[line.tolist(),5]
            galaxy_z[i].progen_disp_uncut[j] = galaxy_z[i].data_array[line.tolist(),4]
            for z in range(0,3):
                galaxy_z[i].progen_ax_uncut[j][z] = galaxy_z[i].data_array[line.tolist(),z+7]#7,8,9

            if galaxy_z[i].snapshot != 761:
                progen_uncut_file_write[i].write(line)


            if mcut == True:
                if not mcut_lower < galaxy_761.data_array[j,masscolumns] < mcut_upper:
                    continue
            if vdcut == True:
                if not vdcut_lower < galaxy_761.data_array[j,vdcolumns] < vdcut_upper:
                    continue
            if level_lock == True:
                if galaxy_761.data_array[j,level_lock_index - 1] != lvl:
                    continue
            #applying all required cuts
            #its fine to not fill this plot array because zeros are shaved off afterwards

            galaxy_z[i].progen_masses[j] = galaxy_z[i].data_array[line.tolist(),6]
            galaxy_z[i].progen_mag[j] = galaxy_z[i].data_array[line.tolist(),5]
            galaxy_z[i].progen_disp[j] = galaxy_z[i].data_array[line.tolist(),4]
            for z in range(0,3):
                galaxy_z[i].progen_ax[j][z] = galaxy_z[i].data_array[line.tolist(),z+7]#7,8,9

            if galaxy_z[i].snapshot != 761:
                file_write[i].write('{} {}\n'.format(str(index),str(progen_index+1)))

            #create progen file and fill plot arrays by going through galaxies
            #in galaxy file and checking for progen in progen file

        galaxy_z[i].masses_plot = np.zeros(np.count_nonzero(galaxy_z[i].progen_masses))
        galaxy_z[i].mag_plot = np.zeros(np.count_nonzero(galaxy_z[i].progen_mag))
        galaxy_z[i].disp_plot = np.zeros(np.count_nonzero(galaxy_z[i].progen_disp))
        galaxy_z[i].masses_plot = np.array([x for x in galaxy_z[i].progen_masses if x != 0])
        galaxy_z[i].mag_plot = np.array([x for x in galaxy_z[i].progen_mag if x != 0])#CHANGED TO ARRAY
        galaxy_z[i].disp_plot = np.array([x for x in galaxy_z[i].progen_disp if x != 0])
        galaxy_z[i].ax_plot = np.array([np.zeros(np.count_nonzero(galaxy_z[i].progen_ax[:,0])),np.zeros(np.count_nonzero(galaxy_z[i].progen_ax[:,1])),np.zeros(np.count_nonzero(galaxy_z[i].progen_ax[:,2]))])
        for z in range(0,3):
            #galaxy_z[i].ax_plot[z] = np.zeros(np.count_nonzero(galaxy_z[i].progen_ax[:,z]))
            galaxy_z[i].ax_plot[z] = [x for x in galaxy_z[i].progen_ax[:,z] if x != 0]
        #remove all zeros


        galaxy_z[i].masses_plot_uncut = np.zeros(np.count_nonzero(galaxy_z[i].progen_masses_uncut))
        galaxy_z[i].mag_plot_uncut = np.zeros(np.count_nonzero(galaxy_z[i].progen_mag_uncut))
        galaxy_z[i].disp_plot_uncut = np.zeros(np.count_nonzero(galaxy_z[i].progen_disp_uncut))
        galaxy_z[i].masses_plot_uncut = np.array([x for x in galaxy_z[i].progen_masses_uncut if x != 0])
        galaxy_z[i].mag_plot_uncut = np.array([x for x in galaxy_z[i].progen_mag_uncut if x != 0])#CHANGED TO ARRAY
        galaxy_z[i].disp_plot_uncut = np.array([x for x in galaxy_z[i].progen_disp_uncut if x != 0])
        galaxy_z[i].ax_plot_uncut = np.array([np.zeros(np.count_nonzero(galaxy_z[i].progen_ax_uncut[:,0])),np.zeros(np.count_nonzero(galaxy_z[i].progen_ax_uncut[:,1])),np.zeros(np.count_nonzero(galaxy_z[i].progen_ax_uncut[:,2]))])
        for z in range(0,3):
            #galaxy_z[i].ax_plot[z] = np.zeros(np.count_nonzero(galaxy_z[i].progen_ax[:,z]))
            galaxy_z[i].ax_plot_uncut[z] = [x for x in galaxy_z[i].progen_ax_uncut[:,z] if x != 0]
        #remove all zeros

        galaxy_z[i].bin_size[0] = 0.5
        galaxy_z[i].n1.append(hist_plot(galaxy_z[i],galaxy_z[i].data_array[:,6],galaxy_z[i].masses_plot,galaxy_z[i].bin_size[0],'Mass'))

        temp = ratio_plot(galaxy_z[i],galaxy_z[i].n1[0],galaxy_z[i].data_array[:,6],galaxy_z[i].masses_plot,galaxy_z[i].bin_size[0],'Mass')
        galaxy_z[i].n2.append(temp[0])
        galaxy_z[i].bin_centers.append(temp[1])
        temp2 = ratio_plot(galaxy_z[i],galaxy_z[i].n1[0],galaxy_z[i].data_array[:,6],galaxy_z[i].masses_plot_uncut,galaxy_z[i].bin_size[0],'Massuncut')
        galaxy_z[i].n3.append(temp2[0])
        #I WILL ASSUME THE BINCENTERS ARE THE SAME FOR BOTH

        galaxy_z[i].bin_size[1] = 0.5
        galaxy_z[i].n1.append(hist_plot(galaxy_z[i],galaxy_z[i].data_array[:,5],galaxy_z[i].mag_plot,galaxy_z[i].bin_size[1],'Magnitude'))

        temp = ratio_plot(galaxy_z[i],galaxy_z[i].n1[1],galaxy_z[i].data_array[:,5],galaxy_z[i].mag_plot,galaxy_z[i].bin_size[1],'Magnitude')
        galaxy_z[i].n2.append(temp[0])
        galaxy_z[i].bin_centers.append(temp[1])
        temp2 = ratio_plot(galaxy_z[i],galaxy_z[i].n1[1],galaxy_z[i].data_array[:,5],galaxy_z[i].mag_plot_uncut,galaxy_z[i].bin_size[1],'Magnitudeuncut')
        galaxy_z[i].n3.append(temp2[0])

        galaxy_z[i].bin_size[2] = 0.25
        galaxy_z[i].n1.append(hist_plot(galaxy_z[i],galaxy_z[i].data_array[:,4],galaxy_z[i].disp_plot,galaxy_z[i].bin_size[2],'Velocity Dispersion'))

        temp = ratio_plot(galaxy_z[i],galaxy_z[i].n1[2],galaxy_z[i].data_array[:,4],galaxy_z[i].disp_plot_uncut,galaxy_z[i].bin_size[2],'Velocity Dispersion')
        galaxy_z[i].n2.append(temp[0])
        galaxy_z[i].bin_centers.append(temp[1])
        temp2 = ratio_plot(galaxy_z[i],galaxy_z[i].n1[2],galaxy_z[i].data_array[:,4],galaxy_z[i].disp_plot,galaxy_z[i].bin_size[2],'Velocity Dispersionuncut')
        galaxy_z[i].n3.append(temp2[0])

        #axis ratio
        galaxy_z[i].bin_size[3] = 0.1
        galaxy_z[i].n1.append(hist_plot(galaxy_z[i],galaxy_z[i].data_array[:,7]/galaxy_z[i].data_array[:,9],galaxy_z[i].ax_plot[0]/galaxy_z[i].ax_plot[2],galaxy_z[i].bin_size[3],'c/a'))

        temp = ratio_plot(galaxy_z[i],galaxy_z[i].n1[3],galaxy_z[i].data_array[:,7]/galaxy_z[i].data_array[:,9],galaxy_z[i].ax_plot[0]/galaxy_z[i].ax_plot[2],galaxy_z[i].bin_size[3],'c/a')
        galaxy_z[i].n2.append(temp[0])
        galaxy_z[i].bin_centers.append(temp[1])


        galaxy_z[i].bin_size[4] = 0.1
        galaxy_z[i].n1.append(hist_plot(galaxy_z[i],galaxy_z[i].data_array[:,8]/galaxy_z[i].data_array[:,9],galaxy_z[i].ax_plot[1]/galaxy_z[i].ax_plot[2],galaxy_z[i].bin_size[4],'b/a'))

        temp = ratio_plot(galaxy_z[i],galaxy_z[i].n1[4],galaxy_z[i].data_array[:,8]/galaxy_z[i].data_array[:,9],galaxy_z[i].ax_plot[1]/galaxy_z[i].ax_plot[2],galaxy_z[i].bin_size[4],'b/a')
        galaxy_z[i].n2.append(temp[0])
        galaxy_z[i].bin_centers.append(temp[1])

        #do all of the plots

    joint_ratio_plot('Mass')
    joint_ratio_plot('Magnitude')
    joint_ratio_plot('Velocity Dispersion')
    #final joint plots with all acquired data

    galaxy_761_file.close()
    for i in range(0,snapshot_num/2):
        galaxy_z_file[i].close()
        file_write[i].close()
        progen_uncut_file_write[i].close()
    merger_file.close()
    #close all relavent files

    os.chdir(Master_dir)
    #move back to oringal directory

    return
