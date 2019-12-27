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

def TidalField2(masscut_array,vdcut_array,snapshots_array,smoothness_array,name_761_file,cutname = '',masscutbool = False,vdcutbool = False,level_lock = True,wspin = False):
    """
    For each desired redshift, this function reads the relavent galaxy files from 'Galaxy_Data' directory, finds the galaxy indices
    of progenotors using the 'linker file' progduced by the merger_code_4 function, and then using the merger tree data produces 'ToPlot' files of the
    tidal information in the region of each galaxy. This is saved in the Tidal_Plots_new directory.

    The inputs are an array of the desired mass bins, [lowercut,uppercut], a snapshot array
    providing the names of the snapshots desired, smoothness array for which smoothness
    levels are desired, cutname for file saving/reading, and
    the level_lock (i.e. whether to only select galaxies which are level 1). (wspin is
    largely irrelevant).

    The outputs are ToPlot dat files in the Tidal_Data_new directory, ready to be plotted by the Tidal_Plots_2 function, and
    galaxy data files which includes data of the progenitor galaxies in the Galaxy_Data_new directory

    """


    #packages

    Master_dir = os.getcwd()

    level_lock_index = 20
    column_num = level_lock_index

    Galaxy_Data_dir = Master_dir + '/Galaxy_Data_new'
    Tidal_Data_dir = Master_dir + '/Tidal_Data_new/' + cutname
    if not os.path.exists(Tidal_Data_dir):
        os.makedirs(Tidal_Data_dir)
    os.chdir(Tidal_Data_dir)
    #make directory for files and images

    #still need to do mass cut for 761

    class field:
        def __init__(self,snap,smooth):
            self.snapshot = snap
            self.smooth = smooth
            return

    class galaxy:
        def __init__(self,snap):
            self.snapshot = snap
            return

    #def reducer_761(galaxy,lvlone,cut_lower,cut_upper,string,columns):
    def reducer_761(galaxy,lvlone,string,*argv):
        #ARGS MUST BE LISTS OF [[LOWERCUT,UPPERCUT,COLUMN]]

        sets = len(argv[0])


        for [cut_lower,cut_upper,columns] in argv[0]:
            #cut_lower = argv[2*i+0]
            #cut_upper = argv[2*i+1]
            #columns = argv[2*i+2]
            #This applies masscut to teh 761 galaxy files, this should match the
            #cut done by the merger code function

            #mass_columns = 6

            level_index = level_lock_index - 1 #ADD EXTRA COLUMN, CHECK THIS
            #mass_cut = 10.4
            if lvlone == True:
                #mask = np.logical_and(galaxy.data_array[:,columns] > cut_lower,galaxy.data_array[:,level_index] == 1,galaxy.data_array[:,columns] < cut_upper)
                mask = np.logical_and(galaxy.data_array[:,columns] > cut_lower,galaxy.data_array[:,columns] < cut_upper)
                galaxy.data_array = galaxy.data_array[mask]
                galaxy.data_array = galaxy.data_array[galaxy.data_array[:,level_index] == 1]

            else:
                mask = np.logical_and(galaxy.data_array[:,columns] > cut_lower,galaxy.data_array[:,columns] < cut_upper)
                galaxy.data_array = galaxy.data_array[mask]






            galaxy.file_length = galaxy.data_array[:,0].size


        if galaxy.data_array.size == 0:
            print(str(galaxy.snapshot) + 'empty' )
        else:
            np.savetxt(Galaxy_Data_dir+'/sample_full_'+string+'_761.dat',galaxy.data_array[:,0:level_lock_index], delimiter=" ")#16

        return galaxy.data_array

    def progen_cutter(galaxy_object,progen_file):
        #This applies cut to data in memory to match cut done by merger code


        progen_file_list = progen_file.readlines()

        progen_file_array = np.zeros([len(progen_file_list),2])

        for j,line in enumerate(progen_file_list,0):
            progen_file_array[j] = np.array([float(x) for x in line.split(" ")])



        s = pd.Series(galaxy_object.data_array[:,0])
        mask = s.isin(progen_file_array[:,1])



        galaxy_object.data_array = galaxy_object.data_array[mask]


        galaxy_object.file_length = galaxy_object.data_array[:,0].size
        #new array smaller in size


        args = np.argsort(np.transpose(galaxy_object.data_array)[0],kind = 'quicksort')
        #order of indices which sorts array in terms of index
        galaxy_object.data_array = galaxy_object.data_array[args]


        if masscutbool == True or vdcutbool == True:
            np.savetxt(Galaxy_Data_dir+'/sample_full_'+string+'_'+str(galaxy_object.snapshot)+'.dat',galaxy_object.data_array[:,0:level_lock_index], delimiter=" ")


        return galaxy_object.data_array

    def line_finder(pos):
        pos += 50/0.704
        pixel_length = 100/(0.704*512)
        pos_floor = pos//pixel_length

        nearest_pixel = pos_floor[2]
        #data file goes through z first

        nearest_pixel += 512*pos_floor[1]
        #data file then goes through y

        nearest_pixel += 512*512*pos_floor[0]
        #data file then goes through x

        return int(nearest_pixel)
    #returns pixel location of galaxy

    def find_angle(tidal_eigenvector,minor):
        x = np.arccos(np.dot(minor,tidal_eigenvector)/((np.linalg.norm(minor)*np.linalg.norm(tidal_eigenvector))))
        if x <= np.pi/2:
            return x
        elif x > np.pi/2 and x <= np.pi:
            return np.pi - x
        #should be normalised
    #returns angle between two vectors < 90

    #tidal_list = [131,197,343,519,570,638,691,733,761]
    #tidal_list = [406]
    #tidal_list = [302]
    #make sure to put in ascending order
    #smoothness = [2]
    #alter the smoothness and snapshots you want to read in

    tidal_list = snapshots_array
    smoothness = smoothness_array

    if masscutbool == True:
        masscut_lower = masscut_array[0]
        masscut_upper = masscut_array[1]
    if vdcutbool == True:
        vdcut_lower = vdcut_array[0]
        vdcut_upper = vdcut_array[1]

    #MUST BE IN NUMERICAL ORDER

    masscolumns = 6
    vdcolumns = 4

    #progen_cut = False
    #progen_mass_cut = False
    #progen_mass_split_1 = False
    #progen_mass_split_2 = False
    #progen_mass_split_3 = False
    #level_lock =False
    #wspin = True


    #progen_cut = sys.argv[1]
    #progen_mass_cut = sys.argv[2]
    #progen_mass_split_1 = sys.argv[3]
    #progen_mass_split_2 = sys.argv[4]
    #progen_mass_split_3 = sys.argv[5]


    #make true if you require either one of these cuts. Make sure required files
    #are already saved in the right dirctory. Do not set both to true

    tidal_num = len(tidal_list)
    smooth_num = len(smoothness)

    tidal_files = np.zeros([tidal_num,smooth_num],dtype = object)
    tidal_objects = np.zeros([tidal_num,smooth_num],dtype = object)
    galaxy_files = np.zeros(tidal_num,dtype = object)
    #if cut == True: DIDNT KNOW WHY THIS WAS HERE
    if 761 in tidal_list:
        galaxy_progen_files = np.zeros(tidal_num - 1,dtype = object)
    else:
        galaxy_progen_files = np.zeros(tidal_num ,dtype = object)
        #there isn't a progen file for 761


    galaxy_objects = np.zeros(tidal_num,dtype = object)
    file_write_toplot = np.zeros([tidal_num,smooth_num],dtype = object)
    file_write_other = np.zeros([tidal_num,smooth_num],dtype = object)
    #initialisers

    #if masscutbool == True or vdcutbool == True: DONT KNOW WHY I HAD THIS CONDITION
    string = 'progen_' + cutname
    #else:
    #        string = ''
    if level_lock == True:
        string += '_lvlone'


    for i,x in enumerate(tidal_list,0):
        if x != 761:
            #if cut == True:  I DIDNT KNOW WHY I PUT THIS THERE
            if level_lock == True:
                galaxy_progen_files[i] = open(Master_dir + '/Progen_Plots_new/'+cutname+'/'+string+'_761_'+str(x)+'.dat','r')
            else:
                galaxy_progen_files[i] = open(Master_dir + '/Progen_Plots_new/'+cutname+'/'+string+'_761_'+str(x)+'.dat','r')
        #progenitors

        if x == 761:
            galaxy_files[i] = open(Master_dir + '/Galaxy_Data_new/' + name_761_file,'r')
        else:
            galaxy_files[i] = open(Master_dir + '/Galaxy_Data_new/sample_full_'+str(x)+'_allgal_wspin.dat','r')


        galaxy_objects[i] = galaxy(x)
        for j,y in enumerate(smoothness,0):

            tidal_files[i][j] = open('/mnt/extraspace/chisari/whenia/tidaleigenvectors-hzAGN-smooth'+str(y)+'-snap'+str(x)+'.dat','r')
            tidal_objects[i][j] = field(x,y)
            file_write_toplot[i][j] = open(Tidal_Data_dir+'/Tidal_Fields_Data_ToPlot_'+string+'_smooth'+str(y)+'_'+str(x)+'.dat','w')
            file_write_other[i][j] = open(Tidal_Data_dir+'/Tidal_Fields_Data_Other_'+string+'_smooth'+str(y)+'_'+str(x)+'.dat','w')
            #need more file writes
    #create field and galaxy and file objects from list
    #makae sure files are in the right directories, alter directories as required

    start = time.time()

    print('Extracting Tidal Field Data...')

    for i in range(0,tidal_num):
        print 'Extracting Snapshot ' + str(galaxy_objects[i].snapshot)
        galaxy_objects[i].data_list = galaxy_files[i].readlines()
        galaxy_files[i].close()
        galaxy_objects[i].file_length = len(galaxy_objects[i].data_list)

        #if galaxy_objects[i].snapshot == 761:
        #    if level_lock == True:
        #        columns = 17
        #    else:
        #        columns = 16#added
        #elif wspin == True:
        #    columns = 17
        #else:
        #    columns = 16

        columns = level_lock_index


        galaxy_objects[i].data_array = np.zeros([galaxy_objects[i].file_length,columns + 1])
        #read in galaxy file to list, initialise array
        #16 columns in dat file, plus one for pixels

        for counter,line in enumerate(galaxy_objects[i].data_list,0):
            current_line = np.array([float(x) for x in line.split(" ")])
            pos = current_line[1:4]
            nearest_pixel = line_finder(pos)
            galaxy_objects[i].data_array[counter][0:columns] = current_line
            galaxy_objects[i].data_array[counter][columns] = nearest_pixel
        #create array of split up galaxy data and nearest pixel



        #THIS NEXT SET OF CONDITIONAL STATEMENTS IS DEFINITELY NOT OPTIMAL
        if galaxy_objects[i].snapshot == 761:
            if masscutbool == True and vdcutbool == False:
                if level_lock == True:
                    galaxy_objects[i].data_array = reducer_761(galaxy_objects[i],True,string,[[masscut_lower,masscut_upper,masscolumns]])#NEED TO AD COLUMNS
                else:
                    galaxy_objects[i].data_array = reducer_761(galaxy_objects[i],False,string,[[masscut_lower,masscut_upper,masscolumns]])

            if vdcutbool == True and masscutbool == False:
                if level_lock == True:
                    galaxy_objects[i].data_array = reducer_761(galaxy_objects[i],True,string,[[vdcut_lower,vdcut_upper,vdcolumns]])#NEED TO AD COLUMNS
                else:
                    galaxy_objects[i].data_array = reducer_761(galaxy_objects[i],False,string,[[vdcut_lower,vdcut_upper,vdcolumns]])

            if vdcutbool == True and masscutbool == True:
                if level_lock == True:
                    galaxy_objects[i].data_array = reducer_761(galaxy_objects[i],True,string,[[masscut_lower,masscut_upper,masscolumns],[vdcut_lower,vdcut_upper,vdcolumns]])#NEED TO AD COLUMNS
                else:
                    galaxy_objects[i].data_array = reducer_761(galaxy_objects[i],False,string,[[masscut_lower,masscut_upper,masscolumns],[vdcut_lower,vdcut_upper,vdcolumns]])#set to true?

        if galaxy_objects[i].snapshot != 761:
            if masscutbool == True or vdcutbool == True:
                galaxy_objects[i].data_array = progen_cutter(galaxy_objects[i],galaxy_progen_files[i])





        #adjusts array if particular cut desired

        args = np.argsort(np.transpose(galaxy_objects[i].data_array)[columns],kind = 'quicksort')
        #order of indices which sorts array in terms of pixel
        galaxy_objects[i].data_array_sort = galaxy_objects[i].data_array[args]
        #sorts in terms of galaxy pixel necessary for the following for loop
        #which loops through tidal file until the next galaxy is found from this list



        for j in range(0,smooth_num):
            line_counter = -1
            flag = False
            #important to initalise this way for following for loop

            print 'Extracting Smoothness '+ str(smoothness[j])


            for k,galaxy_index in enumerate(np.transpose(galaxy_objects[i].data_array_sort)[0],0):
                #rememeber galaxy_index are floating point, and start at 1.0


                if flag == False:
                    line_counter += 1
                #flag being true corresponds to this galaxy being in the same pixel as the
                #last one
                while(line_counter != int(galaxy_objects[i].data_array_sort[k][columns])):
                    #read until next pixel with galaxy in
                    tidal_files[i][j].readline()
                    line_counter += 1
                    if line_counter > 512**3:
                        print 'Overrun'
                        exit()
                #loop round until the next pixel corresponds to galaxy

                if flag == False:
                    line = tidal_files[i][j].readline()
                #read in the next line

                flag = False

                if k < galaxy_objects[i].file_length - 1:
                    if int(galaxy_objects[i].data_array_sort[k][columns]) == int(galaxy_objects[i].data_array_sort[k+1][columns]):
                        flag = True
                #so long as this is not the last galaxy, test whether the next
                #galaxy is in the same pixel as this one, if so set flag to true

                pixel_central = np.array([float(x) for x in line.split(" ")])
                #zeroth order approximation, use whatever pixel the galaxy is
                #in to approximate the tidal data

                v1 = np.array(pixel_central[3:6],dtype = float)
                v2 = np.array(pixel_central[6:9],dtype = float)
                v3 = np.array(pixel_central[9:12],dtype = float)
                #zeroth order estimate of tidal field eigenvctors

                pos = galaxy_objects[i].data_array_sort[k][1:4]
                #load in position of galaxy
                minor = galaxy_objects[i].data_array_sort[k][10:13]
                major = galaxy_objects[i].data_array_sort[k][13:16]
                spin = galaxy_objects[i].data_array_sort[k][16:19]
                #load minor and major axis of galaxy tensor

                #REMEMBER THAT I MEAN A MISTAKE HERE, MAJOR MINOR REVERSED

                theta_1_minor = find_angle(v1,minor)
                theta_1_major = find_angle(v1,major)
                theta_2_minor = find_angle(v2,minor)
                theta_2_major = find_angle(v2,major)
                theta_3_minor = find_angle(v3,minor)
                theta_3_major = find_angle(v3,major)

                theta_1_spin = find_angle(v1,spin)
                theta_2_spin = find_angle(v2,spin)
                theta_3_spin = find_angle(v3,spin)
                #returns angle <90

                file_write_other[i][j].write('{} {} {} {} {} {} {} {} {} {}\n'.format(galaxy_index,v1[0],v1[1],v1[2],v2[0],v2[1],v2[2],v3[0],v3[1],v3[2]))
                file_write_toplot[i][j].write('{} {} {} {} {} {} {} {} {} {}\n'.format(galaxy_index,theta_1_minor,theta_2_minor,theta_3_minor,theta_1_major,theta_2_major,theta_3_major,theta_1_spin,theta_2_spin,theta_3_spin))
                #write data to ToPlot and Other file names
            #loops through all galaxies in galaxy file, and passes their pixel
            #position into the next while loop
        #for given snapshot and smoothness finds tidal data for each galaxy and writes
        #it to a file

    for i in range(0,tidal_num - 1):
        galaxy_progen_files[i].close()
        galaxy_files[i].close()
        for j in range(0,smooth_num):
            tidal_files[i][j].close()

            file_write_toplot[i][j].close()
            file_write_other[i][j].close()
    #close all remaining files


    end = time.time()



    print('Time Elapsed:\t', end-start)

    os.chdir(Master_dir)

    return
