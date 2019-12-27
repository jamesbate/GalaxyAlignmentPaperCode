import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

masscut = False
progen = False
mass_split_1 =False
mass_split_2 = False
mass_split_3 = False
allmasssplit = True
allmasssplitjoint = False
level_lock = True

if masscut == True:
    file_string = '_progen_masscut'
elif progen == True:
    file_string = '_progen'
elif mass_split_1 == True:
    file_string = '_progen_masssplit1'
elif mass_split_2 == True:
    file_string = '_progen_masssplit2'
elif mass_split_3 == True:
    file_string = '_progen_masssplit3'
elif allmasssplit == True:
    file_string_1 = '_progen_masssplit1'
    file_string_2 = '_progen_masssplit2'
    file_string_3 = '_progen_masssplit3'
    file_string = 'joint_masssplit'
elif allmasssplitjoint == True:
    #file_string_1 = '_progen_masssplit1'
    #file_string_2 = '_progen_masssplit2'
    #file_string_3 = '_progen_masssplit3'
    #file_string = 'joint_masssplit'
    file_string_1 = ''
    file_string_2 = '_progen_masscut'
    file_string_3 = '_progen_masssplit3'
    file_string = 'cust_masssplit'
else:
    file_string = ''
if level_lock == True:
    file_string += '_lvlone'
    #if this is relavent
    file_string_1 += '_lvlone'
    file_string_2 += '_lvlone'
    file_string_3 += '_lvlone'


class galaxy:
    def __init__(self,snap):
            self.snapshot = snap
            return

class tide:
    def __init__(self,snap):
            self.snapshot = snap
            return


snapshots = [[3,131],[2,197],[1,343],[0.5,519],[0.42,570],[0.31,638],[0.2,691],[0.12,733]]
#add 761 full

if allmasssplitjoint == True or allmasssplit == True:
    tidal_objects = np.zeros([len(snapshots),3],dtype = object)
else:
    tidal_objects = np.zeros(len(snapshots),dtype = object)

for i,snap in enumerate(snapshots,0):

    if allmasssplitjoint == True or allmasssplit == True:
            tidal_objects[i][0] = tide(snap[1])
            tidal_objects[i][0].file_array = np.array(np.loadtxt('/mnt/zfsusers/orie3568/Tidal_Data/Tidal_Fields_Data_ToPlot'+file_string_1+'_smooth2_'+str(snap[1])+'.dat',delimiter = ' '))
            tidal_objects[i][1] = tide(snap[1])
            tidal_objects[i][1].file_array = np.array(np.loadtxt('/mnt/zfsusers/orie3568/Tidal_Data/Tidal_Fields_Data_ToPlot'+file_string_2+'_smooth2_'+str(snap[1])+'.dat',delimiter = ' '))
            tidal_objects[i][2] = tide(snap[1])
            tidal_objects[i][2].file_array = np.array(np.loadtxt('/mnt/zfsusers/orie3568/Tidal_Data/Tidal_Fields_Data_ToPlot'+file_string_3+'_smooth2_'+str(snap[1])+'.dat',delimiter = ' '))
    else:
        tidal_objects[i] = tide(snap[1])
        tidal_objects[i].file_array = np.array(np.loadtxt('/mnt/zfsusers/orie3568/Tidal_Data/Tidal_Fields_Data_ToPlot'+file_string+'_smooth2_'+str(snap[1])+'.dat',delimiter = ' '))

colours = ['m','c','b','g','r','m','c','b','g','r','m','c','b','g','r']
b = []
fig_0,ax = plt.subplots(num = "Full")
plt.subplots_adjust(left = 0.15,right = 0.85,top = 0.85, bottom = 0.15)
for i in range(1,4):
#for i in range(4,7):
#for i in range(1,2):

    colour = colours[i]
    if allmasssplit == True:
        plt.clf()
        plt.figure("Full")
        b = []
    if allmasssplit == True or allmasssplitjoint == True:
        temp1 = plt.plot([x[0] for x in snapshots],[np.mean(y.file_array[:,i]) for y in tidal_objects[:,0]],str(colour)+'-',label =  r'$10.4<\log(M/M_{sun})<10.7$')
        #temp1 = plt.plot([x[0] for x in snapshots],[np.mean(y.file_array[:,i]) for y in tidal_objects[:,0]],str(colour)+'-',label = r'$\theta_'+str(i)+'^{minor}$ full population')
        plt.plot([x[0] for x in snapshots],[np.mean(y.file_array[:,i]) for y in tidal_objects[:,0]],str(colour)+'o')
        yerror1 = np.ones(len(snapshots))*[np.std(y.file_array[:,i])/np.sqrt(y.file_array[:,i].size) for y in tidal_objects[:,0]]
        plt.errorbar([x[0] for x in snapshots],[np.mean(y.file_array[:,i]) for y in tidal_objects[:,0]],ecolor = str(colour),fmt = colour+'o',yerr = yerror1,ms = 4, capthick=1,capsize = 8,elinewidth = 1)

        temp2 = plt.plot([x[0] for x in snapshots],[np.mean(y.file_array[:,i]) for y in tidal_objects[:,1]],str(colour)+'--',label = r'$10.7<\log(M/M_{sun})<11.06$')
        #temp2 = plt.plot([x[0] for x in snapshots],[np.mean(y.file_array[:,i]) for y in tidal_objects[:,1]],str(colour)+'--',label = r'$\theta_'+str(i)+'^{minor} \log(M/M_{sun}) > 10.4$')
        plt.plot([x[0] for x in snapshots],[np.mean(y.file_array[:,i]) for y in tidal_objects[:,1]],str(colour)+'o')
        yerror2 = np.ones(len(snapshots))*[np.std(y.file_array[:,i])/np.sqrt(y.file_array[:,i].size) for y in tidal_objects[:,1]]
        plt.errorbar([x[0] for x in snapshots],[np.mean(y.file_array[:,i]) for y in tidal_objects[:,1]],ecolor = str(colour),fmt = colour+'o',yerr = yerror2,ms = 4, capthick=1,capsize = 8,elinewidth = 1)

        temp3 = plt.plot([x[0] for x in snapshots],[np.mean(y.file_array[:,i]) for y in tidal_objects[:,2]],str(colour)+':',label = r'$11.06<\log(M/M_{sun})$')
        #temp3 = plt.plot([x[0] for x in snapshots],[np.mean(y.file_array[:,i]) for y in tidal_objects[:,2]],str(colour)+':',label = r'$\theta_'+str(i)+'^{minor} \log(M/M_{sun}) > 11.06$')
        plt.plot([x[0] for x in snapshots],[np.mean(y.file_array[:,i]) for y in tidal_objects[:,2]],str(colour)+'o')
        yerror3 = np.ones(len(snapshots))*[np.std(y.file_array[:,i])/np.sqrt(y.file_array[:,i].size) for y in tidal_objects[:,2]]
        plt.errorbar([x[0] for x in snapshots],[np.mean(y.file_array[:,i]) for y in tidal_objects[:,2]],ecolor = str(colour),fmt = colour+'o',yerr = yerror3,ms = 4, capthick=1,capsize = 8,elinewidth = 1)


        b.append(temp1[0])
        b.append(temp2[0])
        b.append(temp3[0])

    else:
        temp = plt.plot([x[0] for x in snapshots],[np.mean(y.file_array[:,i]) for y in tidal_objects],colour,label = r'$\theta_'+str(i)+'^{major}$')
        plt.plot([x[0] for x in snapshots],[np.mean(y.file_array[:,i]) for y in tidal_objects],str(colour)+'o')
        yerror = np.ones(len(snapshots))*[np.std(y.file_array[:,i])/np.sqrt(y.file_array[:,i].size) for y in tidal_objects]
        plt.errorbar([x[0] for x in snapshots],[np.mean(y.file_array[:,i]) for y in tidal_objects],ecolor = str(colour),fmt = colour+'o',yerr = yerror,ms = 4, capthick=1,capsize = 8,elinewidth = 1)
        b.append(temp[0])


    plt.legend(handles = [k for k in b],frameon = False)
    plt.xlabel('z', fontsize = 17)

    if allmasssplit == True:
        plt.ylabel(r'$\langle \theta_'+str(i)+r'^{major} \rangle $ (radians)', fontsize = 17)

    else:
        plt.ylabel(r'$\langle \theta \rangle (radians)$', fontsize = 17)

    if masscut == True:
        plt.title(r'High mass elliptical progenitors')
    if progen == True:
        plt.title('Progenitors')
    if mass_split_1 == True:
        plt.title(r'Progenitors of Galaxies with $10.4<\log(M/M_{sun})<10.7$ at z = 0')
    if mass_split_2 == True:
        plt.title(r'Progenitors of Galaxies with $10.7<\log(M/M_{sun})<11.06$ at z = 0')
    if mass_split_3 == True:
        plt.title(r'Progenitors of Galaxies with $11.06<\log(M/M_{sun})$ at z = 0')
    if allmasssplit == True or allmasssplitjoint == True:
        #plt.title('mass split progenitor population')
        #plt.title('Average Angle Against Redshift')
        pass

    #plt.title('Full Population')
    #TEMPORARY

    plt.plot(np.linspace(0,3,3),np.ones(3),'k--')

    #ax = plt.axes()
    ax.tick_params(axis='both', labelsize = 12 , length = 5,width = 1)

    if allmasssplit == True:
        plt.savefig('/mnt/zfsusers/orie3568/RedshiftPlots/mean_angleplot'+file_string+str(i)+'.png')
        plt.savefig('/mnt/zfsusers/orie3568/RedshiftPlots/mean_angleplot'+file_string+str(i)+'.eps',format='eps', dpi=1000)

if allmasssplit == False:
    plt.savefig('/mnt/zfsusers/orie3568/RedshiftPlots/mean_angleplot_alljoint_'+file_string+'.png')
    plt.savefig('/mnt/zfsusers/orie3568/RedshiftPlots/mean_angleplot_alljoint_'+file_string+'.eps',format='eps', dpi=1000)
