import numpy as np
import os
import matplotlib as mpl
mpl.use('Agg')
from scipy import stats
import matplotlib.pyplot as plt
from collections import Counter
from matplotlib import rc
#preamble

def Tidal_Plots_2(snapshots_array,smoothness_array,cut,name_761_file,cutname,level_lock = True):
    """
    This function uses the ToPlot file to plot information on the alignment of galaxies with the tidal field. It reads in the
    galaxy files form the Galaxy_Data_new dorectory produced by the TidalField2 function, and along with the ToPlot Tidal data is
    produces alignments plots against z.

    The inputs are an array of the desired mass bins, an array of the snapshots adn redshifts [redshift,tree index,snapshot], an array of
    the desired smoothness vales, a cut name for file reading and sharing, as well as the boolean flags cut, for whether a mass cut
    has been made, level_lock for whether only level 1 galaxies are accepted, and wspin which is irrelevant now.

    The outputs are alignemnt plots in the Tidal_Images_new_%cutname% directory

    Note: if no cut, looks for full sample file in Galaxy_Data_new

    """
    Master_dir = os.getcwd()
    Galaxy_Data_dir = Master_dir + '/Galaxy_Data_new'
    Tidal_Data_dir = Master_dir + '/Tidal_Data_new'

    level_lock_index = 20

    ToPlotBase = Tidal_Data_dir+'/'+cutname+'/Tidal_Fields_Data_ToPlot_progen_' + cutname
    if cut == True:
        GalaxyBase ='/sample_full_progen_' + cutname
    else:
        GalaxyBase ='/sample_full_'

    if level_lock == True:
        GalaxyBase += '_lvlone'
        ToPlotBase += '_lvlone'

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
    }#standardise colour scheme


    if not os.path.exists(Master_dir + '/Tidal_Image_Plots/Tidal_Images_new_'+cutname):
        os.makedirs(Master_dir + '/Tidal_Image_Plots/Tidal_Images_new_'+cutname)
    os.chdir(Master_dir + '/Tidal_Image_Plots/Tidal_Images_new_'+cutname)
    #make directory for images

    class galaxy:
        def __init__(self,snapshot):
            self.snapshot = snapshot
            return

    class tide:
        def __init__(self,snapshot):
            self.isempty = False
            self.snapshot = snapshot
            return

    def weight_finder(c,b,a):
        #apply ellipticity weighting to see how much it changes plots
        return 1 - c/a

    def bin_plot_angle_prob(angles,galaxy,snap,z,string):
        """
        Plots pdf distribution of galaxy angles with the tidal field
        """
        if len(angles) == 0 or galaxy.array.size == 0:
            print('No tidal data for ' + str(galaxy.snapshot) +' binplot')
            return
        plt.clf()
        galaxy_num = angles[0].size
        #assume each angle columns is the same size
        n_bins = 20

        fig_0 = plt.figure('Full')

        colours = ['m','c','b','g','r'][0:smooth_num]#if desired, use colours dict

        b = []
        #legend

        for i,angle in enumerate(angles,0):
            plt.figure('temp')
            probability, bins, patches = plt.hist(angle, n_bins,normed = True)
            plt.clf()
            plt.figure('Full')
            bincenters = 0.5*(bins[1:]+bins[:-1])
            bin_width = bins[1] - bins[0]
            #Use histogram to get binned data
            yerror = np.ones(n_bins)*np.sqrt((probability*galaxy_num*bin_width))/(galaxy_num*bin_width)

            temp = plt.plot(bincenters,probability,colours[i]+'--',label = 'Actual Probabilities Smoothness '+str(smoothness[i]) )
            b.append(temp[0])
            plt.errorbar(bincenters,probability,yerr = yerror,fmt = colours[i]+'o',ms = 2,ecolor = colours[i], capthick=1,capsize = 4,elinewidth = 1)
        #plots

        x = np.linspace(0,np.pi/2,20)
        a, = plt.plot(x,np.sin(x),'k:',label = 'Expected if Randomly Distributed')
        #sinusoidal basline
        plt.xlabel(string)
        plt.ylabel('Probability Density')
        plt.title('Histogram Plot of Angular Distribution Redshift '+str(z))
        plt.grid(True)
        plt.xticks(np.linspace(0,1.6,9),['0','','0.4','','0.8','','1.2','','1.6'])

        plt.legend(handles = [a]+[k for k in b])
        plt.savefig(string+'_'+str(snap) + '_sm'+str(smoothness)+'_angleplot.png')
        plt.savefig(string+'_'+str(snap) + '_sm'+str(smoothness)+'_angleplot.eps',format='eps', dpi=1000)
        plt.clf()
        plt.close()
        return

    def scatter_plot_inertia_avangle(mass,angles,snap,z,string):
        """
        scatter plot of mass against average angle of the bin 
        (Ended up using not using this plot)
        """
        if len(angles) == 0 or len(mass)== 0:
            print('No tidal data for massplot')
            return

        plt.clf()
        n_bins = 12
        colours = ['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9'][0:smooth_num]
        #if we were to use plots, make sure to use the global colour dict 

        a = []
        #legend

        fig = plt.figure(figsize = (8,6))

        for i,angle in enumerate(angles,0):
            avangle,edges,bin_num = stats.binned_statistic(mass, angle, statistic = 'mean', bins=n_bins)
            n = Counter(bin_num)
            yerror = np.ones(n_bins)
            for j in range(0,n_bins):
                yerror[j] *= np.std(angle[j])/np.sqrt(n[j+1])
            bincenters = 0.5*(edges[1:]+edges[:-1])
            plt.errorbar(bincenters,avangle,yerr = yerror,fmt=colours[i]+'o',ms = 4, ecolor=colours[i], capthick=1,capsize = 8,elinewidth = 1)
            temp =  plt.plot(bincenters,avangle,colours[i]+'--',label = 'Measured PDF Smoothness '+str(smoothness[i]))
            a.append(temp[0])
        #plot with errors

        plt.grid()
        plt.xticks([9,10,11,12,13])
        fig.suptitle('Binned Plot of Average Angle against Mass Redshift '+str(z), fontsize=18)
        plt.xlabel('Mass', fontsize=10)
        plt.ylabel('Average angle ' + string, fontsize=10)

        x = np.linspace(edges[0],edges[len(edges)-1],20)
        b, = plt.plot(x,np.ones(len(x)),'k:',label = 'Expected if Randomly Distributed')
        #random baseline
        plt.legend(handles = [l for l in a]+[b])
        plt.savefig(string+'_'+str(snap) + '_sm'+str(smoothness)+'-massangleplot.png')
        plt.savefig(string+'_'+str(snap) + '_sm'+str(smoothness)+'-massangleplot.eps',format='eps', dpi=1000)
        plt.clf()
        plt.close()
        return

    def bin_plot_angle_prob_diff(angles,galaxy,snap,z,string):
        """
        Pdf of galaxy, plotting difference from mean 
        """

        if len(angles) == 0 or galaxy.array.size == 0:
            print('No tidal data for ' + str(galaxy.snapshot) +' binplotdiff')
            return
        #close figures
        plt.clf()
        galaxy_num = angles[0].size
        n_bins = 7
        #n_bins = 20

        fig_0 = plt.figure('Full')

        colours = ['m','c','b','g','r'][0:smooth_num]
        colours = ['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9'][0:smooth_num]



        b = []

        for i,angle in enumerate(angles,0):
            plt.figure('temp')
            probability, bins, patches = plt.hist(np.cos(angle), n_bins,normed = True)
            plt.clf()
            plt.figure('Full')
            bincenters = 0.5*(bins[1:]+bins[:-1])
            bin_width = bins[1] - bins[0]
            #Use histogram to get binned data
            yerror = np.ones(n_bins)*np.sqrt((probability*galaxy_num*bin_width))/(galaxy_num*bin_width)


            #temp = plt.plot(bincenters,probability,'c--',label = 'Actual Probabilities Smoothness '+str(smoothness[i]) )
            #b.append(temp[0])
            #plt.errorbar(bincenters,probability,yerr = yerror,fmt='co',ms = 2, ecolor='c', capthick=1,capsize = 4,elinewidth = 1)
            temp = plt.plot(bincenters,probability-1,colours[i]+'--',label = 'Actual Probabilities Smoothness '+str(smoothness[i]) )
            b.append(temp[0])
            plt.errorbar(bincenters,probability-1,yerr = yerror,fmt = colours[i]+'o',ms = 2,ecolor = colours[i], capthick=1,capsize = 4,elinewidth = 1)

        x = np.linspace(0,1,100)
        a, = plt.plot(x,np.zeros(100),'k:',label = 'Expected if Randomly Distributed Redshift '+str(z))
        #sinusoidal basline
        plt.xlabel('cos('+string+')')
        plt.ylabel('Probability Density')
        plt.title('Histogram Plot of Angular Distribution')
        plt.grid(True)


        plt.legend(handles = [a]+[k for k in b])
        plt.savefig(string+'_'+str(snap) + '_sm'+str(smoothness)+'_diffangleplot.png')
        plt.savefig(string+'_'+str(snap) + '_sm'+str(smoothness)+'_diffangleplot.eps',format='eps', dpi=1000)
        plt.clf()
        plt.close()
        return

    #The above plots are designed to plot single angles, with different
    #smoothnesses if desired

    def bin_plot_angle_prob_joint(angles,snap,z,galaxies,string):
        #or galaxies.array.size == 0 -> could be adapted if needed (need to look at all galaxies in list)
        if len(angles) == 0:
            print('No tidal data for ' + snap + ' binplotjoint')
            return
        plt.clf()

        if string == 'theta_1_minor':
            lat_string = '\\theta^{\\rm minor}_{1}'
        elif string == 'theta_2_minor':
            lat_string = '\\theta^{\\rm minor}_{2}'
        elif string == 'theta_3_minor':
            lat_string = '\\theta^{\\rm minor}_{3}'
        elif string == 'theta_1_major':
            lat_string = '\\theta^{\\rm{major}}_{1}'
        elif string == 'theta_2_major':
            lat_string = '\\theta^{\\rm major}_{2}'
        elif string == 'theta_3_major':
            lat_string = '\\theta^{\\rm major}_{3}'
        elif string == 'theta_1_spin':
            lat_string = '\\theta^{\\rm spin}_{1}'
        elif string == 'theta_2_spin':
            lat_string = '\\theta^{\\rm spin}_{2}'
        elif string == 'theta_3_spin':
            lat_string = '\\theta^{\\rm spin}_{3}'
        else:
            print 'string unmatched'
            exit()


        galaxy_num = []
        for n in range(0,len(angles)):
            galaxy_num.append(len(angles[n]))

        n_bins = 20


        fig_0 = plt.figure('Full')

        b = []


        for i,angle in enumerate(angles,0):
            plt.figure('temp')
            probability, bins, patches = plt.hist(angles[i], n_bins,normed = True)
            plt.clf()
            fig_0 = plt.figure('Full')
            bincenters = 0.5*(bins[1:]+bins[:-1])
            bin_width = bins[1] - bins[0]
            #Use histogram to get binned data
            yerror = np.ones(n_bins)*np.sqrt((probability*galaxy_num[i]*bin_width))/(galaxy_num[i]*bin_width)

            temp = plt.plot(bincenters,probability,colour_dict[str(snap[i])]+'--',label = 'Redshift ' + str(z[i]) )

            b.append(temp[0])
            plt.errorbar(bincenters,probability,yerr = yerror,fmt = colour_dict[str(snap[i])]+'o',ms = 2,ecolor = colour_dict[str(snap[i])], capthick=1,capsize = 4,elinewidth = 1)

        x = np.linspace(0,np.pi/2,20)
        a, = plt.plot(x,np.sin(x),'k:',label = 'Random')
        #sinusoidal basline
        plt.xlabel(r'$'+lat_string+'$')
        plt.ylabel('Probability Density')
        plt.title('Joint Histogram Plot of Angular Distribution')
        plt.grid(True)
        plt.xticks(np.linspace(0,1.6,9),['0','','0.4','','0.8','','1.2','','1.6'])

        plt.legend(handles = [a]+[k for k in b])
        plt.savefig(string+'_angleplot_joint.png')
        plt.savefig(string+'_angleplot_joint.eps',format='eps', dpi=1000)
        plt.clf()
        plt.close()
        return

    def bin_plot_angle_prob_diff_joint(angles,snap,z,galaxies,string):
        # or galaxies.array.size == 0
        if len(angles) == 0:
            print('No tidal data for ' + snap +' binplotdiffjoint')
            return
        #close figures

        n_bins = 10

        if string == 'theta_1_minor':
            lat_string = '\\theta^{minor}_{1}'
        elif string == 'theta_2_minor':
            lat_string = '\\theta^{minor}_{2}'
        elif string == 'theta_3_minor':
            lat_string = '\\theta^{minor}_{3}'
        elif string == 'theta_1_major':
            lat_string = '\\theta^{major}_{1}'
        elif string == 'theta_2_major':
            lat_string = '\\theta^{major}_{2}'
        elif string == 'theta_3_major':
            lat_string = '\\theta^{major}_{3}'
        elif string == 'theta_1_spin':
            lat_string = '\\theta^{spin}_{1}'
        elif string == 'theta_2_spin':
            lat_string = '\\theta^{spin}_{2}'
        elif string == 'theta_3_spin':
            lat_string = '\\theta^{spin}_{3}'
        else:
            print 'string unmatched'
            exit()



        plt.clf()
        galaxy_num = []
        for n in range(0,len(angles)):
            galaxy_num.append(len(angles[0]))

        b = []
        fig_0,ax = plt.subplots(num = "Full")

        for i,angle in enumerate(angles,0):
            plt.figure('temp')
            probability, bins, patches = plt.hist(np.cos(angles[i]), n_bins,normed = True)
            plt.clf()
            plt.figure('Full')

            plt.subplots_adjust(left = 0.2,right = 0.8,top = 0.85, bottom = 0.15)#FIDDLE

            bincenters = 0.5*(bins[1:]+bins[:-1])
            bin_width = bins[1] - bins[0]
            #Use histogram to get binned data
            yerror = np.ones(n_bins)*np.sqrt((probability*galaxy_num[i]*bin_width))/(galaxy_num[i]*bin_width)


            temp = plt.plot(bincenters,probability,colour_dict[str(snap[i])]+'--',label = 'z =  ' + str(z[i]))
            b.append(temp[0])
            plt.errorbar(bincenters,probability,yerr = yerror,fmt = colour_dict[str(snap[i])]+'o',ms = 2,ecolor = colour_dict[str(snap[i])], capthick=1,capsize = 4,elinewidth = 1)

        x = np.linspace(0,1,100)
        a, = plt.plot(x,np.ones(100),'k:',label = 'Random')
        #sinusoidal basline
        #ax = plt.axes()
        ax.tick_params(axis='both', labelsize = 14 , length = 5,width = 1)

        plt.xlabel(r'$\cos({\rm '+lat_string+'})$',fontsize = 17)
        plt.ylabel(r'$1 + \xi$',fontsize = 17)

        #plt.title('Joint Histogram Plot of Angular Distribution')


        #plt.title('Full Population',fontsize = 17)
        plt.title('High Mass Elliptical Progenitors',fontsize = 17)

        # Shrink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        # Put a legend to the right of the current axis

        ax.legend(handles = [a]+[k for k in b],loc='center left', bbox_to_anchor=(1, 0.8))

        plt.savefig(string+'_diffangleplot_joint.png')
        plt.savefig(string+'_diffangleplot_joint.eps',format='eps', dpi=1000)
        plt.clf()
        plt.close()
        return

    def scatter_plot_inertia_avangle_joint(masses,angles,string):
        if len(angles) == 0 or len(masses) == 0:
            print('No tidal data for binplotmassjoint')
            return
        plt.clf()
        n_bins = 12

        a = []


        fig = plt.figure(figsize = (8,6))


        for i,(angle,mass) in enumerate(zip(angles,masses),0):
            avangle,edges,bin_num = stats.binned_statistic(mass, angle, statistic = 'mean', bins=n_bins)
            n = Counter(bin_num)
            yerror = np.ones(n_bins)*np.sqrt(np.std(angle))
            for j in range(0,n_bins):
                yerror[j] /= np.sqrt(n[j+1])
            bincenters = 0.5*(edges[1:]+edges[:-1])
            #plt.errorbar(bincenters,avangle,yerr = yerror,fmt=colours[i]+'o',ms = 4, ecolor=colours[i], capthick=1,capsize = 8,elinewidth = 1)
            plt.scatter(bincenters,avangle,c = colour_dict[str(snapshots[i])])
            temp =  plt.plot(bincenters,avangle,colour_dict[str(snapshots[i])]+'--',label = 'Snapshot '+str(snapshots[i]))
            a.append(temp[0])

        plt.grid()
        plt.xticks([9,10,11,12,13])
        fig.suptitle('Binned Plot of Average Angle against Mass ', fontsize=18)
        plt.xlabel('Mass', fontsize=10)
        plt.ylabel('Average angle ' + string, fontsize=10)

        x = np.linspace(edges[0],edges[len(edges)-1],20)
        b, = plt.plot(x,np.ones(len(x)),'k:',label = 'Expected if Randomly Distributed')
        #random baseline
        plt.legend(handles = [l for l in a]+[b])
        plt.savefig(string + '_massangleplot_joint.png')
        plt.savefig(string + '_massangleplot_joint.eps',format='eps', dpi=1000)
        plt.clf()
        plt.close()
        return

    #These joint plots are designed to plot angles of different snapshots simultaneously


    snapshots = snapshots_array[:,2].astype(int).tolist()
    redshifts = snapshots_array[:,0].astype(float).tolist()

    smoothness = smoothness_array

    snap_num = len(snapshots)
    smooth_num = len(smoothness)
    #read in desired snapshots and smoothnesses

    galaxy_files = np.zeros(snap_num,dtype = object)
    galaxy_objects = np.zeros(snap_num,dtype = object)
    tidal_files = np.zeros([snap_num,smooth_num],dtype = object)
    tidal_objects = np.zeros([snap_num,smooth_num],dtype = object)
    #initialisers

    for n,i in enumerate(snapshots,0):
        if cut == True:
            galaxy_files[n] = open(Galaxy_Data_dir+GalaxyBase+'_'+str(i)+'.dat','r')
        else:
            if i == 761:
                galaxy_files[n] = open(Galaxy_Data_dir+'/'+name_761_file,'r')
            else:
                galaxy_files[n] = open(Galaxy_Data_dir+GalaxyBase+str(i)+'_allgal_wspin.dat','r')

        galaxy_objects[n] = galaxy(i)
        for m,j in enumerate(smoothness,0):
            tidal_files[n][m] = open(ToPlotBase+'_smooth'+str(j)+'_'+str(i)+'.dat','r')
            tidal_objects[n][m] = tide(i)
    #open all required files



    for i in range(0,snap_num):
        galaxy_objects[i].list = galaxy_files[i].readlines()
        #read in galaxy file

        column_number = level_lock_index

        galaxy_objects[i].array = np.zeros([len(galaxy_objects[i].list),column_number])
        #prime galaxy array

        for j,line in enumerate(galaxy_objects[i].list,0):
            line_array = [float(x) for x in line.split(" ")]

            galaxy_objects[i].array[j] = np.array(line_array)
        #read in galaxy files


        for j in range(0,smooth_num):
            tidal_objects[i][j].list = tidal_files[i][j].readlines()
            #read in tidal data

            if tidal_objects[i][j].list == []:
                print str(tidal_objects[i][j].snapshot) + ' is empty...'
                tidal_objects[i][j].isempty = True
                continue
            #continue loop if tidal file empty

            file_length = 10
            tidal_objects[i][j].array = np.zeros([len(tidal_objects[i][j].list),file_length])
            #prime tidal array

            for k,line in enumerate(tidal_objects[i][j].list,0):
                tidal_objects[i][j].array[k] = np.array([float(x) for x in line.split(" ")])
            #read in tidal file

            args = np.argsort(tidal_objects[i][j].array[:,0],kind = 'quicksort')
            tidal_objects[i][j].array = tidal_objects[i][j].array[args]
            #sort tidal data in terms of index
        #reads in tidal file and sorts in terms of index

        for j in range(0,smooth_num):
            if tidal_objects[i][j].isempty == False:
                tidal_objects[i][j].theta_minor = np.zeros([len(tidal_objects[i][j].list),3])
                tidal_objects[i][j].theta_major = np.zeros([len(tidal_objects[i][j].list),3])
                tidal_objects[i][j].theta_spin = np.zeros([len(tidal_objects[i][j].list),3])
                tidal_objects[i][j].theta_minor = np.transpose(tidal_objects[i][j].array[:,1:4])
                tidal_objects[i][j].theta_major = np.transpose(tidal_objects[i][j].array[:,4:7])
                tidal_objects[i][j].theta_spin = np.transpose(tidal_objects[i][j].array[:,7:10])
            else:
                tidal_objects[i][j].theta_minor = np.array([])
                tidal_objects[i][j].theta_major = np.array([])
                tidal_objects[i][j].theta_minor = np.array([])
                tidal_objects[i][j].theta_major = np.array([])
                tidal_objects[i][j].theta_spin = np.array([])
        #prime minor major axis angle arrays

        bin_plot_angle_prob([x.theta_minor[0] for x in tidal_objects[i] if x.isempty == False],galaxy_objects[i],snapshots[i],redshifts[i],'theta_1_major')#theta 1 minor
        bin_plot_angle_prob([x.theta_minor[1] for x in tidal_objects[i] if x.isempty == False],galaxy_objects[i],snapshots[i],redshifts[i],'theta_2_major')#theta 2 minor
        bin_plot_angle_prob([x.theta_minor[2] for x in tidal_objects[i] if x.isempty == False],galaxy_objects[i],snapshots[i],redshifts[i],'theta_3_major')#theta 3 minor
        bin_plot_angle_prob([x.theta_major[0] for x in tidal_objects[i] if x.isempty == False],galaxy_objects[i],snapshots[i],redshifts[i],'theta_1_minor')#theta 1 major
        bin_plot_angle_prob([x.theta_major[1] for x in tidal_objects[i] if x.isempty == False],galaxy_objects[i],snapshots[i],redshifts[i],'theta_2_minor')#theta 2 major
        bin_plot_angle_prob([x.theta_major[2] for x in tidal_objects[i] if x.isempty == False],galaxy_objects[i],snapshots[i],redshifts[i],'theta_3_minor')#theta 3 major
        bin_plot_angle_prob([x.theta_spin[0] for x in tidal_objects[i] if x.isempty == False],galaxy_objects[i],snapshots[i],redshifts[i],'theta_1_spin')#theta 1 minor
        bin_plot_angle_prob([x.theta_spin[1] for x in tidal_objects[i] if x.isempty == False],galaxy_objects[i],snapshots[i],redshifts[i],'theta_2_spin')#theta 2 minor
        bin_plot_angle_prob([x.theta_spin[2] for x in tidal_objects[i] if x.isempty == False],galaxy_objects[i],snapshots[i],redshifts[i],'theta_3_spin')#theta 3 minor

        scatter_plot_inertia_avangle(galaxy_objects[i].array[:,6],[x.theta_minor[0] for x in tidal_objects[i] if x.isempty == False],snapshots[i],redshifts[i],'theta_1_major')
        scatter_plot_inertia_avangle(galaxy_objects[i].array[:,6],[x.theta_minor[1] for x in tidal_objects[i] if x.isempty == False],snapshots[i],redshifts[i],'theta_2_major')
        scatter_plot_inertia_avangle(galaxy_objects[i].array[:,6],[x.theta_minor[2] for x in tidal_objects[i] if x.isempty == False],snapshots[i],redshifts[i],'theta_3_major')
        scatter_plot_inertia_avangle(galaxy_objects[i].array[:,6],[x.theta_major[0] for x in tidal_objects[i] if x.isempty == False],snapshots[i],redshifts[i],'theta_1_minor')
        scatter_plot_inertia_avangle(galaxy_objects[i].array[:,6],[x.theta_major[1] for x in tidal_objects[i] if x.isempty == False],snapshots[i],redshifts[i],'theta_2_minor')
        scatter_plot_inertia_avangle(galaxy_objects[i].array[:,6],[x.theta_major[2] for x in tidal_objects[i] if x.isempty == False],snapshots[i],redshifts[i],'theta_3_minor')
        scatter_plot_inertia_avangle(galaxy_objects[i].array[:,6],[x.theta_spin[0] for x in tidal_objects[i] if x.isempty == False],snapshots[i],redshifts[i],'theta_1_spin')
        scatter_plot_inertia_avangle(galaxy_objects[i].array[:,6],[x.theta_spin[1] for x in tidal_objects[i] if x.isempty == False],snapshots[i],redshifts[i],'theta_2_spin')
        scatter_plot_inertia_avangle(galaxy_objects[i].array[:,6],[x.theta_spin[2] for x in tidal_objects[i] if x.isempty == False],snapshots[i],redshifts[i],'theta_3_spin')

        bin_plot_angle_prob_diff([x.theta_minor[0] for x in tidal_objects[i] if x.isempty == False],galaxy_objects[i],snapshots[i],redshifts[i],'theta_1_major')#theta 1 minor
        bin_plot_angle_prob_diff([x.theta_minor[1] for x in tidal_objects[i] if x.isempty == False],galaxy_objects[i],snapshots[i],redshifts[i],'theta_2_major')#theta 2 minor
        bin_plot_angle_prob_diff([x.theta_minor[2] for x in tidal_objects[i] if x.isempty == False],galaxy_objects[i],snapshots[i],redshifts[i],'theta_3_major')#theta 3 minor
        bin_plot_angle_prob_diff([x.theta_major[0] for x in tidal_objects[i] if x.isempty == False],galaxy_objects[i],snapshots[i],redshifts[i],'theta_1_minor')#theta 1 minor
        bin_plot_angle_prob_diff([x.theta_major[1] for x in tidal_objects[i] if x.isempty == False],galaxy_objects[i],snapshots[i],redshifts[i],'theta_2_minor')#theta 2 minor
        bin_plot_angle_prob_diff([x.theta_major[2] for x in tidal_objects[i] if x.isempty == False],galaxy_objects[i],snapshots[i],redshifts[i],'theta_3_minor')#theta 3 minor
        bin_plot_angle_prob_diff([x.theta_spin[0] for x in tidal_objects[i] if x.isempty == False],galaxy_objects[i],snapshots[i],redshifts[i],'theta_1_spin')#theta 1 minor
        bin_plot_angle_prob_diff([x.theta_spin[1] for x in tidal_objects[i] if x.isempty == False],galaxy_objects[i],snapshots[i],redshifts[i],'theta_2_spin')#theta 2 minor
        bin_plot_angle_prob_diff([x.theta_spin[2] for x in tidal_objects[i] if x.isempty == False],galaxy_objects[i],snapshots[i],redshifts[i],'theta_3_spin')#theta 3 minor
        #plots all smoothnesses

    if snap_num > 1:
        for z in range(0,3):
            pass

            bin_plot_angle_prob_joint([x.theta_minor[z] for x in np.transpose([tidal_objects[y] for y in range(0,5)])[0]],[snapshots[k] for k in range(0,5)],[redshifts[k] for k in range(0,5)],[galaxy_objects[i] for i in range(0,5)],'theta_'+str(z+1)+'_major')
            bin_plot_angle_prob_diff_joint([x.theta_minor[z] for x in np.transpose([tidal_objects[y] for y in [0,1,2,5,8]])[0]],[snapshots[k] for k in [0,1,2,5,8]],[redshifts[k] for k in [0,1,2,5,8]],[galaxy_objects[i] for i in [0,1,2,5,8]],'theta_'+str(z+1)+'_major')

            bin_plot_angle_prob_joint([x.theta_spin[z] for x in np.transpose([tidal_objects[y] for y in range(0,5)])[0]],[snapshots[k] for k in range(0,5)],[redshifts[k] for k in range(0,5)],[galaxy_objects[i] for i in range(0,5)],'theta_'+str(z+1)+'_spin')
            bin_plot_angle_prob_diff_joint([x.theta_spin[z] for x in np.transpose([tidal_objects[y] for y in [0,1,2,5,7]])[0]],[snapshots[k] for k in [0,1,2,5,7]],[redshifts[k] for k in [0,1,2,5,7]],[galaxy_objects[i] for i in [0,1,2,5,7]],'theta_'+str(z+1)+'_spin')

            bin_plot_angle_prob_joint([x.theta_major[z] for x in np.transpose([tidal_objects[y] for y in range(0,5)])[0]],[snapshots[k] for k in range(0,5)],[redshifts[k] for k in range(0,5)],[galaxy_objects[i] for i in range(0,5)],'theta_'+str(z+1)+'_minor')
            bin_plot_angle_prob_diff_joint([x.theta_major[z] for x in np.transpose([tidal_objects[y] for y in [0,1,2,5,8]])[0]],[snapshots[k] for k in [0,1,2,5,8]],[redshifts[k] for k in [0,1,2,5,8]],[galaxy_objects[i] for i in [0,1,2,5,8]],'theta_'+str(z+1)+'_minor')


    #plots desired snapshots

    os.chdir(Master_dir)
    return
