#This is an example provided by Garreth on how to read the merger file
import numpy as np
from scipy.io import FortranFile

#Read the merger tree file from Garreth
f = FortranFile('merger_tree_GM.dat','r')
data = f.read_record(dtype=np.float32)
#Do some re-arranging to create id_cat, a matrix of these dimensions:
#i=2: i=0 is the label of the z=0 snapshot; i=1 is the index of the z=0 galaxy.
#j=91: number of snapshots for which progenitors of z=0 galaxies are available.
#k=126361: total number of galaxies in snapshot 761 at z=0
id_cat = data.reshape((2,91,126361))
f.close()

#And make some additions/replacements
#which populate i=0,1 as defined above:
id_cat[1,90,:] = np.arange(0,126361)
id_cat[0,90,:] = 761

#Some clarifications
#The second index specifies the snapshot and the third index specifies the id of the main descendant at snapshot 761. The first index is just to switch between the list of snapshot numbers (index 0) and galaxy id (index 1). Please note that I started the galaxy ids from 0 (instead of from 1 as they are given in Yohan's catalogs)

#An example:
print "Snapshots of main progenitors for galaxy 24:"
print np.int32(id_cat[0,:,23])
#gives the all of the snapshots of the main progenitors the galaxy with galaxy id = 23 (official id = 24) at snapshot 761

#Another example:
print "Galaxy ids of main progenitors for galaxy 24:"
print np.int32(id_cat[1,:,23])
#gives the galaxy id of the main progenitor for each of these snapshots (again official ids will be +1).

#Important notes
#If the value of the snapshot is 0, then the galaxy has no progenitor at that snapshot, probably because they are not massive enough.
#You should also disregard snapshots where the value of the index is -1 for a similar reason.

exit()
