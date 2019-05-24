#!/bin/bash
export CMSSW_PROJECT_SRC=/afs/cern.ch/user/j/jingyan/CMS_NEW/CMSSW_10_5_0_pre1/src
cd $CMSSW_PROJECT_SRC
eval `scramv1 runtime -sh`
export X509_USER_PROXY=/afs/cern.ch/user/j/jingyan/x509up_u122075
cd /afs/cern.ch/user/j/jingyan/CMS_NEW/CMSSW_10_5_0_pre1/src/L1Trigger/TrackFindingTracklet/test/Condorjob
cmsRun L1TrackNtupleMakerTestParsing_cfg.py inputFiles=root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/5CAC6671-5050-A04E-B495-F9D6791B4E48.root outputFile=100events_D21_Hybrid68.root  maxEvents=100

