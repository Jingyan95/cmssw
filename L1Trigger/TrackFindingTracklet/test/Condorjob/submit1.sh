#!/bin/bash
export CMSSW_PROJECT_SRC=/afs/cern.ch/user/j/jingyan/CMS_NEW/CMSSW_10_5_0_pre1/src
cd $CMSSW_PROJECT_SRC
eval `scramv1 runtime -sh`
export X509_USER_PROXY=/afs/cern.ch/user/j/jingyan/x509up_u122075
cd /afs/cern.ch/user/j/jingyan/CMS_NEW/CMSSW_10_5_0_pre1/src/L1Trigger/TrackFindingTracklet/test/Condorjob
cmsRun L1TrackNtupleMakerTestParsing_cfg.py inputFiles=root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/FFF48AB4-E5E6-3842-8A5B-20E2B7E497BC.root,root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/FFB04F2D-0C40-2040-94CD-587D1F5AC4A7.root,root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/FEDE3D03-2B54-BF4C-8D63-0C8C4DBAC6F0.root,root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/FE522E6A-09E8-6B4B-B111-C5F6322B71D1.root,root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/FD5A1B6D-C998-C149-8277-B278047E0011.root,root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/FBBA43FA-9AAB-D349-BE35-18F74BEBA975.root,root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/F7A323B1-9683-9846-9139-358A93E6E5AD.root,root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/F70A7F1D-EDF4-7A47-8A13-0710A5964D03.root,root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/F6FCF82F-C6A7-F744-99C2-CBCEA15FF124.root,root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/EED79AA9-1E2A-7A4B-8A74-507CBA91E250.root outputFile=1000events_D21_Hybrid1.root  maxEvents=1000

