import os
import sys
nfile = 100
filename1 = '100events_D21_Hybrid'
filename2 = '_numEvent100.root'
for x in xrange(1, nfile+1):
    filename=filename1+str(x)+filename2
    os.system('rm -rf '+filename)
