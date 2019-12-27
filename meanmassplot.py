import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os

if not os.path.exists('MeanPlots'):
    os.makedirs('MeanPlots')
os.chdir('MeanPlots')

def scatter_plot_inertia_avangle_joint(masses,redshifts,string):
    plt.clf()
    fig_0,ax = plt.subplots(num = "Full")
    plt.subplots_adjust(left = 0.2,right = 0.8,top = 0.85, bottom = 0.15)


    means = np.array([np.mean(masses[i]) for i in range(0,len(redshifts))])
    plt.plot(redshifts,means)

    std_joint = np.array([np.std(masses[i]) for i in range(0,len(redshifts))])

    yerror = std_joint/np.sqrt([masses[i].size for i in range(0,len(redshifts))])
    #error for frations?

    plt.errorbar(redshifts,means,yerr = yerror,ms = 4,marker = 'o', ecolor = 'b',fmt = 'b',capthick=1,capsize = 8,elinewidth = 1)

    plt.xlabel('z', fontsize=10)

    ax.tick_params(axis='both', labelsize = 14 , length = 5,width = 1)

    plt.xlabel(r'$z$',fontsize = 17)
    if string == 'mass':
        plt.ylabel(r'$\langle M \rangle$',fontsize = 17)
    if string == 'caratio':
        plt.ylabel(r'$\langle \frac{b}{a} \rangle$',fontsize = 17)
    if string == 'baratio':
        plt.ylabel(r'$\langle \frac{c}{a} \rangle$',fontsize = 17)

    plt.savefig(string + '_meansplot.png')
    return

snapshots = [131,197,343,519,570,638,691,733]
redshifts = [3,2,1,0.5,0.42,0.31,0.2,0.12]

masses = np.array([np.loadtxt('/mnt/zfsusers/orie3568/Galaxy_Data/sample_full_progen_masscut_'+str(i)+'.dat',usecols = 6) for i in snapshots])
c_a = np.array([np.loadtxt('/mnt/zfsusers/orie3568/Galaxy_Data/sample_full_progen_masscut_'+str(i)+'.dat',usecols = 7)/np.loadtxt('/mnt/zfsusers/orie3568/Galaxy_Data/sample_full_progen_masscut_'+str(i)+'.dat',usecols = 9) for i in snapshots])
b_a = np.array([np.loadtxt('/mnt/zfsusers/orie3568/Galaxy_Data/sample_full_progen_masscut_'+str(i)+'.dat',usecols = 7)/np.loadtxt('/mnt/zfsusers/orie3568/Galaxy_Data/sample_full_progen_masscut_'+str(i)+'.dat',usecols = 8) for i in snapshots])
#b_a = np.array([np.loadtxt('Galaxy_Data/sample_full_progen_masscut_'+str(i)+'.dat',usecols = 6) for i in snapshots])
#c_a = np.array([np.loadtxt('Galaxy_Data/sample_full_progen_masscut_'+str(i)+'.dat',usecols = 6) for i in snapshots])

scatter_plot_inertia_avangle_joint(masses,redshifts,'mass')
scatter_plot_inertia_avangle_joint(c_a,redshifts,'caratio')
scatter_plot_inertia_avangle_joint(b_a,redshifts,'baratio')
