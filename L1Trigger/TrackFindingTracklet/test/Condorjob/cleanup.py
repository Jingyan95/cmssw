import os
import sys
nfile = 10
filename1 = '1000events_D21_Hybrid'
filename2 = '_numEvent1000.root'
for x in xrange(1, nfile+1):
    filename=filename1+str(x)+filename2
    os.system('rm -rf '+filename)
