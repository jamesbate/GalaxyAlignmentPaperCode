import numpy as np

mass = np.loadtxt('/mnt/zfsusers/orie3568/Galaxy_Data/sample_full_761.dat',usecols = 6)

mass = np.sort(mass,kind = 'quicksort')

mass = mass[mass>10.4]

n = mass.size

print mass[n//3]
print mass[2*n//3]
