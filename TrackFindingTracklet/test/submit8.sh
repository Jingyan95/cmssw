#!/bin/bash
export CMSSW_PROJECT_SRC=/afs/cern.ch/user/j/jingyan/CMSSW_10_4_0/src
cd $CMSSW_PROJECT_SRC
eval `scramv1 runtime -sh`
export X509_USER_PROXY=/afs/cern.ch/user/j/jingyan/x509up_u122075
cd /afs/cern.ch/user/j/jingyan/CMSSW_10_4_0/src/L1Trigger/TrackFindingTracklet/test
cmsRun L1TrackNtupleMakerTestParsing_cfg.py inputFiles=root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/49E94AB3-4EDE-6B40-BD87-4FA47FDA1BCD.root,root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/4725E564-4DEA-FB4F-B9FE-834999AFDC4D.root,root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/45C4B37D-7CD2-0541-B399-181CB2C7633A.root,root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/456BF325-22AC-8D4E-A21C-8EB53E89AE14.root,root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/445584F5-57DA-6F4D-814C-B4A1A9C17981.root,root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/42A97B14-003B-1A4A-8D86-215DD230F088.root,root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/3FD7A367-96FD-D249-A362-E0DBD4D03E27.root,root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/37B90FE4-88D4-294D-B2EC-8CF002337A65.root,root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/379D9981-C7E6-CA48-8722-E5E022A35D15.root,root://cms-xrd-global.cern.ch///store/relval/CMSSW_10_4_0/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_103X_upgrade2023_realistic_v2_2023D21PU200-v1/20000/345B6106-BD75-8944-A9E7-3C6893C0317B.root outputFile=1000events_D21_Hybrid8.root  maxEvents=1000

