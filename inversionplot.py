import numpy as np
import os

os.chdir('/mnt/zfsusers/orie3568/Tidal_Data_new')

file_write = open('inversions.dat','w')
snapshots_array = [131,197,343,406,570,638,691,761]
redshifts_array = [3,2,1,0.64,0.42,0.31,0.2,0.06]
inversions = np.zeros(5)


def linear_interpolate(coords,y):

    mask = coords[1] < y
    #< > doesn't matter, only change matters

    flag = coords[1][0]#syntax

    for a,i in enumerate(coords,0):
        flag_comp = i[1]#syntax
        if flag != flag_comp:
            ind_p = a
            ind_n = a-1

    grad = (coords[ind_p][1] - coords[ind_n][1])/(coords[ind_p][0] - coords[ind_n][0])

    y_m = coords[ind_n][1] + ((coords[ind_p][0] - coords[ind_n][0])/2)*grad

    return y_m


for m in range(0,5):
    for v in range(0,5):

        filename_base = 'Tidal_Fields_Data_ToPlot_progen_cutm' + str(m) + 'v' + str(v) + '_smooth2_'

        angle_averages = np.array([])


        for snap in snapshots_array:
            filename = filename_base + str(snap) + '.dat'

            file = open(filename,'r')
            data_list = file.readlines()

            angles = np.zeros([1,len(data_list)])

            for k,line in enumerate(data_list,0):
                angles[k] = np.array([float(x) for x in line.split(" ")])[4]

            angle_averages += np.mean(angles)

        inversions[m][v] = linear_interpolate(np.array([redshifts_array,angle_averages.tolist()]),0)

        file_write.write('{} {} {}\n'.format(m,v,inversions[m][v]))



file_write.close()

            #find averages and add to averages array
            #make linear interpolation
            #plot
