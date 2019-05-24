#!/bin/bash
export CMSSW_PROJECT_SRC=/afs/cern.ch/user/j/jingyan/CMS_NEW/CMSSW_10_5_0_pre1/src
cd $CMSSW_PROJECT_SRC
eval `scramv1 runtime -sh`
export X509_USER_PROXY=/afs/cern.ch/user/j/jingyan/x509up_u122075
cd /afs/cern.ch/user/j/jingyan/CMS_NEW/CMSSW_10_5_0_pre1/src/L1Trigger/TrackFindingTracklet/test/Condorjob
cmsRun L1TrackNtupleMakerTestParsing_cfg.py inputFiles=root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0_pre2/RelValNuGun/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/C2A92EA2-DA6E-CB4A-8975-B0151AFCC8B6.root,root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0_pre2/RelValNuGun/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/C2898E3F-F469-C545-8CB5-9E047CF52390.root,root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0_pre2/RelValNuGun/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/C022A8EE-0873-8449-A8E2-1EB5915CCB8E.root,root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0_pre2/RelValNuGun/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/BFFFC3CD-A623-CD48-845D-DA9A479568F3.root,root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0_pre2/RelValNuGun/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/BFF15079-AB46-9540-B4EC-00624D1A9E54.root,root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0_pre2/RelValNuGun/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/BFABA9A6-67A4-5943-99CA-B4F7579AECD8.root,root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0_pre2/RelValNuGun/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/BF8ABFE4-8D5D-0948-8EF5-F806FE55DECC.root,root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0_pre2/RelValNuGun/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/BD0F7BB4-EC86-BD4F-A84A-528A0A9F866C.root,root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0_pre2/RelValNuGun/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/BB3102C3-43F8-8844-86BD-EEE8FBCBFCB8.root,root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0_pre2/RelValNuGun/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/B9B666D4-91D8-C342-927C-7D426E3C05EE.root outputFile=1000events_D21_Hybrid12.root  maxEvents=1000

