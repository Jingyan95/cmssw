import os
import sys
from itertools import islice
submit = 'universe = vanilla\n' ## writing .sub file
submit += 'arguments = "$(argument)"\n'
submit += 'output = submit01.out\n'
submit += 'error = submit01.err\n'
submit += 'log = submit01.log\n'
submit += '+JobFlavour = "tomorrow"\n' 
submit += 'queue\n'
submitName = 'submit01.sub'
sub1 = open(submitName,'w')
sub1.write(submit+'\n')
sub1.close() ##finish writing .sub file
nfile = 1 ##number of sample file processed per each condor job
filename = 'D49_TTbar_PU200.txt' ## location of the txt file used to help locate samples
with open(filename,'r') as f:
    counter1 = 1 ## the nth job
    while True:
        lines = list(islice(f, nfile))
        if not lines:
            break
        counter2 = 1 ## the nth file in one job
        for line in lines:
            if counter2 == 1:
               input='root://cms-xrd-global.cern.ch//'+line.rstrip() ## Remove  'root:..'+  if you don't have grid certificate. 
            else:
               input =input+',root://cms-xrd-global.cern.ch//'+line.rstrip() ## Remove  ',root:..'+  if you don't have grid certificate.
            counter2+=1
        create = '#!/bin/bash\n' ##writng .sh file
        create += 'export CMSSW_PROJECT_SRC=/afs/cern.ch/user/j/jingyan/CMSSW_11_3_0_pre3/src\n' 
        create += 'cd $CMSSW_PROJECT_SRC\n' ## go to the src directly 
        create += 'eval `scramv1 runtime -sh`\n' ## this is effectively "cmsenv"
        create += 'export X509_USER_PROXY=/afs/cern.ch/user/j/jingyan/x509up_u122075\n' ## exporting the proxy
        create += 'cd /afs/cern.ch/user/j/jingyan/CMSSW_11_3_0_pre3/src/L1Trigger/TrackFindingTracklet/test\n'## go the directly where you have the producer and python config file
        create += 'cmsRun L1TrackNtupleMaker_cfg.py inputFiles='+input+' outputFile='+str(100*nfile)+'events_D49_Hybrid'+str(counter1)+'.root  maxEvents='+str(10000*nfile)+'\n'
        createName = 'submit'+str(counter1)+'.sh'
        sub2 = open(createName,'w')
        sub2.write(create+'\n')
        sub2.close() ## finish writing .sh file
        counter1+=1
        os.system('chmod 755 '+createName) ## make .sh file executable
        os.system('condor_submit '+ submitName+' executable='+createName) ## submit the job using condor_submit command
