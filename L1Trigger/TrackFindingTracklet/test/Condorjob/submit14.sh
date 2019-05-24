#!/bin/bash
export CMSSW_PROJECT_SRC=/afs/cern.ch/user/j/jingyan/CMS_NEW/CMSSW_10_5_0_pre1/src
cd $CMSSW_PROJECT_SRC
eval `scramv1 runtime -sh`
export X509_USER_PROXY=/afs/cern.ch/user/j/jingyan/x509up_u122075
cd /afs/cern.ch/user/j/jingyan/CMS_NEW/CMSSW_10_5_0_pre1/src/L1Trigger/TrackFindingTracklet/test/Condorjob
cmsRun L1TrackNtupleMakerTestParsing_cfg.py inputFiles=root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0_pre2/RelValNuGun/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/A30CE693-1AF3-154E-9ADA-B7DE22C26B58.root,root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0_pre2/RelValNuGun/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/A2160289-E780-4044-86BD-0CC2DFCDBF9A.root,root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0_pre2/RelValNuGun/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/9FD2652A-F6F2-F34A-A528-EC45F7927E28.root,root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0_pre2/RelValNuGun/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/9D7A3398-293F-4643-A590-7185FBF31912.root,root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0_pre2/RelValNuGun/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/9D078A8A-815B-8147-BF51-54E185413DE2.root,root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0_pre2/RelValNuGun/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/9A91C278-E4B6-4E4D-8D49-D3E0F3D2E246.root,root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0_pre2/RelValNuGun/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/98C35E02-9129-9042-A624-DB2CFE53E394.root,root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0_pre2/RelValNuGun/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/97DB3A19-6D18-5C49-8905-0A0E93747ADC.root,root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0_pre2/RelValNuGun/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/9400A38C-05BF-8248-B34D-2B5E48A03779.root,root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0_pre2/RelValNuGun/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/9087C285-0186-BF4F-895F-E6B272A1CA9A.root outputFile=1000events_D21_Hybrid14.root  maxEvents=1000

