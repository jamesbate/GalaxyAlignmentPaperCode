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

def z_avangle_plot(redshifts_array,smoothings_array,cut_name_array,index_plot_array,joint,title,file_string,lvllock,masscutarray,vsigcutarray,interpolater):

    Master_dir = os.getcwd()

    if not os.path.exists(Master_dir + '/Tidal_Data_new'):
        print 'Error' + Master_dir + '/Tidal_Data_new' + ' does not exist'
        exit()
    os.chdir(Master_dir + '/Tidal_Data_new')

    #define tidal data class
    class tide:
        def __init__(self,snap):
                self.snapshot = snap
                return

    def linecrosschecker(points):
        line_y = 1
        point1 = points[0]
        point2 = points[1]
        if point1[0] < point2[0]:
            xl = point1[0]
            yl = point1[1]
            xu = point2[0]
            yu = point2[1]
        if point1[0] > point2[0]:
            xu = point1[0]
            yu = point1[1]
            xl = point2[0]
            yl = point2[1]
        if point1[0] == point2[0]:
            return None
        def interp(y):
            grad = (xu-xl)/(yu-yl)
            if y < yu and y > yl:
                return xl + grad*(y-yl)
            else:
                return None

        return interp(line_y)

    if interpolater == True:
        interpoint = np.zeros([5,5])

    #initialise subplots
    num = cut_name_array.size
    l = np.sqrt(num)

    if l != int(l):
        b_1 = int(l)
        b_2 = int(l) + 1

        v_1 = b_1 * b_1
        v_2 = b_1 * b_2
        v_3 = b_2 * b_2

        if l < v_2:
            x = int(l)
            e = int(l) + 1
        elif l < v_3:
            x = int(l) + 1
            e = int(l) + 1
    else:
        x = l
        e = l

    colours = ['m','c','b','g','r','m','c','b','g','r','m','c','b','g','r','m','c','b','g','r','m','c','b','g','r','m','c','b','g','r']
    ticks = ['-','--',':','-.','-','--',':','-','--',':','-','--',':']
    labels = {1:r'$\theta_1^{\rm major}$',2:r'$\theta_2^{\rm major}$',3:r'$\theta_3^{\rm major}$',4:r'$\theta_1^{\rm minor}$',5:r'$\theta_2^{\rm minor}$',6:r'$\theta_3^{\rm minor}$',7:r'$\theta_1^{\rm spin}$',8:r'$\theta_2^{\rm spin}$',9:r'$\theta_3^{\rm spin}$'}


    if joint == False:
        mpl.rcParams.update({'font.size': 12})
        fig, axes = plt.subplots(int(x),int(e),constrained_layout=False,figsize=(9,9))#boolean logic,
        #fig.suptitle('Alignment Plot over 3x3 Grid',fontsize = 16)
        fig.text(0.5, 0.02, r'Z', ha='center',size=12)
        fig.text(0.02, 0.5, r'$\langle \theta \rangle $ (radians)', va='center', rotation='vertical',size=12)

        for u in range(int(x)):
            for r in range(int(e)):
                if num == 2:
                    ax = axes[r]
                else:
                    ax = axes[u,r]

                plt.subplots_adjust(left = 0.075,right = 0.92,top = 0.90, bottom = 0.075)
                ax.tick_params(pad=1,length=1)
                ax.set_ylim(0.8, 1.2)

                cutname = cut_name_array[u*int(x) + r]

                #-1 for INDEX

                #load in to plot files
                tidal_objects =  np.zeros([redshifts_array[:,0].size,len(smoothings_array)],dtype = object)
                To_Plot_Base = cutname+'/Tidal_Fields_Data_ToPlot_progen_' + cutname + '_'
                if lvllock == True:
                    To_Plot_Base += 'lvlone_'

                for p,i in enumerate(redshifts_array[:,1],0):
                    for m,j in enumerate(smoothings_array,0):
                        tidal_objects[p][m] = tide(i)
                        tidal_objects[p][m].file_array = np.array(np.loadtxt(To_Plot_Base + 'smooth' + str(j) + '_' + str(int(i)) + '.dat',delimiter = ' '))
                        if tidal_objects[p][m].file_array.size == 0:
                            tidal_objects[p][m].isempty = True
                        else:
                            tidal_objects[p][m].isempty = False


                for q,i in enumerate(index_plot_array,0):

                    colour = colours[q]
                    b = []
                    for c,sm in enumerate(smoothings_array,0):
                        tick = ticks[c]
                        for y in tidal_objects[:,c]:
                            if len(y.file_array.shape) == 1:
                                y.file_array = np.array([y.file_array])
                        #reduce redshift array
                        mask = [not y.isempty for y in tidal_objects[:,c]]
                        redshifts_array_reduced = redshifts_array[:,0][mask]
                        temp = ax.plot(redshifts_array_reduced ,[np.mean(y.file_array[:,i]) for y in tidal_objects[:,c] if y.isempty == False],str(colour)+tick,label = labels[i] + 'sm ' + str(sm))#can remove smoothness
                        ax.scatter(redshifts_array_reduced,[np.mean(y.file_array[:,i]) for y in tidal_objects[:,c] if y.isempty == False],c = str(colour),marker = 'o',s = 3)
                        yerror = np.ones(redshifts_array_reduced.size)*[np.std(y.file_array[:,i])/np.sqrt(y.file_array[:,i].size) for y in tidal_objects[:,c] if y.isempty == False]
                        ax.errorbar(redshifts_array_reduced,[np.mean(y.file_array[:,i]) for y in tidal_objects[:,c] if y.isempty == False],ecolor = str(colour),fmt = colour+'o',yerr = yerror,ms = 3, capthick=1,capsize = 4,elinewidth = 1)
                        b.append(temp[0])
                        ax.plot(np.linspace(0,3,3),np.ones(3),'k--')
                        ax.set_title(r'$\log_{10}(M/{\rm M}_{\odot})  ' + str(masscutarray[u]) + 'V_{\\theta}/\\sigma ' + str(vsigcutarray[r]) + '$',fontsize = 7,pad = 2)

                        if interpolater == True:
                            points = np.transpose([redshifts_array_reduced,[np.mean(y.file_array[:,i]) for y in tidal_objects[:,c] if y.isempty == False]])
                            if points[:,0].size > 1:
                                points1 = points[1:]
                                points2 = points[:-1]
                                z = zip(points1.tolist(),points2.tolist())
                                z = [list(a) for a in z]
                                d = map(linecrosschecker,z)
                                crosses = d[d != None]
                                if isinstance(crosses,float):
                                    ax.scatter(crosses,1,s = 30,c = 'r')
                                    interpoint[u,r] = crosses
                                elif crosses != None:
                                    y_0 = np.ones(len(crosses))
                                    ax.scatter(crosses,y_0,s = 30, c = 'r')
                                    interpoint[u,r] = crosses[0]#take first one FOR NOW
                                else:
                                    interpoint[u,r] = None



    #plt.xlabel('redshift',fontsize = 10)
    #plt.ylabel('mean angle',fontsize = 10)
    plt.savefig(Master_dir + '/Tidal_Image_Plots/mean_angleplot_'+str(999)+'_'+file_string+'4x4.png', dpi=1000)
    plt.savefig(Master_dir + '/Tidal_Image_Plots/mean_angleplot_'+str(999)+'_'+file_string+'4x4.eps',format='eps', dpi=1000)

    if interpolater== True:
        masscut_array = np.array([[9,9.75],[9.75,10.5],[10.5,11.25],[11.25,12],[12,12.75]])
        vdcut_array = np.array([[0,0.5],[0.5,0.75],[0.75,1],[1,1.5],[1.5,2.0]])
        plt.figure('new')
        xx,yy = np.meshgrid(masscut_array[:,0],vdcut_array[:,0])
        plt.contourf(xx,yy,interpoint)
        plt.savefig(Master_dir +'/Tidal_Image_Plots/cock.png')

    if joint == True:
        fig_0,ax = plt.subplots(num = "Full")
        #plt.subplots_adjust(left = 0.15,right = 0.75,top = 0.85, bottom = 0.15)
        #if you want legend outside
        b = []
        for n,cutname in enumerate(cut_name_array,1):
            tick = ticks[n-1]
            tidal_objects =  np.zeros([redshifts_array[:,0].size,len(smoothings_array)],dtype = object)
            To_Plot_Base = cutname + '/Tidal_Fields_Data_ToPlot_progen_' + cutname + '_'
            if lvllock == True:
                To_Plot_Base += 'lvlone_'

            for p,i in enumerate(redshifts_array[:,1],0):
                for m,j in enumerate(smoothings_array,0):
                    tidal_objects[p][m] = tide(i)
                    tidal_objects[p][m].file_array = np.array(np.loadtxt(To_Plot_Base +'smooth' + str(j) + '_' + str(int(i)) + '.dat',delimiter = ' '))

            #TEMPORARY
            #smoothings_array_label = [0.4,0.8,1.6,3.2]
            for q,i in enumerate(index_plot_array,0):
                for c,sm in enumerate(smoothings_array,0):
                    if len(cut_name_array) == 1:
                        tick = ticks[c]
                    #colour = colours[q]
                    colour = colours[i-1]
                    lab = labels[i]
                    if len(smoothings_array) > 1:
                        lab += ' sm ' + str(sm)

                    #temp = plt.plot(redshifts_array[:,0] ,[np.mean(y.file_array[:,i]) for y in tidal_objects[:,c]],str(colour)+tick,label = lab)
                    temp = plt.plot(redshifts_array[:,0] ,[np.mean(y.file_array[:,i]) for y in tidal_objects[:,c]],str(colour)+tick,label = r'$\log_{10}(M/{\rm M}_{\odot})  '+str(masscutarray[n-1])+'$')
                    plt.plot(redshifts_array[:,0],[np.mean(y.file_array[:,i]) for y in tidal_objects[:,c]],str(colour)+'o')
                    yerror = np.ones(redshifts_array[:,0].size)*[np.std(y.file_array[:,i])/np.sqrt(y.file_array[:,i].size) for y in tidal_objects[:,c]]
                    plt.errorbar(redshifts_array[:,0],[np.mean(y.file_array[:,i]) for y in tidal_objects[:,c]],ecolor = str(colour),fmt = colour+'o',yerr = yerror,ms = 4, capthick=1,capsize = 8,elinewidth = 1)
                    b.append(temp[0])
                    plt.plot(np.linspace(0,3,3),np.ones(3),'k--')



        #plt.legend(handles = [k for k in b],frameon = False,loc = "center right",bbox_to_anchor=(1.4,0.5))
        plt.legend(handles = [k for k in b])
        plt.xlabel('z', fontsize = 15)
        plt.ylabel(r'$\langle \theta \rangle $ (radians)', fontsize = 15)
        ax.tick_params(labelsize=14)
        #ax.set_ylim([0.85,1.15])
        ax.set_ylim([0.8,1.2])
        #print('$\langle \theta^{\rm '+ str(labels[index_plot_array[0]]),'} \rangle $ (radians))
        #plt.ylabel(r'$\langle \theta^{\rm '+str(labels[index_plot_array[0]][1:-1])+'} \rangle $ (radians)', fontsize = 17)
        plt.title(r''+title, fontsize = 17)

    plt.savefig(Master_dir + '/Tidal_Image_Plots/mean_angleplot_'+cutname+'_'+file_string+'joint.png',dpi=1000)
    plt.savefig(Master_dir + '/Tidal_Image_Plots/mean_angleplot_'+cutname+'_'+file_string+'joint.eps',format='eps', dpi=1000)
    plt.clf
    plt.close()
    os.chdir(Master_dir)
    return
