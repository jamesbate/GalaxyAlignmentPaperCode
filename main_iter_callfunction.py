import os

#masscut_array = [[9,9.75],[9.75,10.5],[10.5,11.25],[11.25,12],[12,12.75]]
#vdcut_array = [[0,0.5],[0.5,0.75],[0.75,1],[1,1.5],[1.5,2.0]]

#masscut_array = [[0,9.4],[9.4,10],[10,10.6],[10.6,15]]
#vdcut_array = [[0,0.4],[0.4,0.7],[0.7,1],[1,2]]

masscut_array = [[10.4,10.7],[10.7,11.06],[11.06,15]]

#masscut_array = [[10.5,11.25]]
#vdcut_array = [[0.75,1]]
#for i in range(4):
#    for j in range(4):
#        string_command = 'addqueue -q berg -m 8 -n 1x1 /usr/bin/python2.7 main_iter.py ' + str(masscut_array[i][0])+ ' ' + str(masscut_array[i][1])+ ' ' + str(vdcut_array[j][0])+ ' ' + str(vdcut_array[j][1])+ ' 4x4cutm' + str(i) + 'v' + str(j)
#        os.system(string_command)

for i in range(3):
    string_command = 'addqueue -q berg -m 8 -n 1x1 /usr/bin/python2.7 main_iter.py ' + str(masscut_array[i][0])+ ' ' + str(masscut_array[i][1])+ ' ' + str(0)+ ' ' + str(0.6)+ ' mscutm' + str(i) + 'v0.6'
    os.system(string_command)
