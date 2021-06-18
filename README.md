## This is a simple script Jack sets up as an example to run L1Tracking code through Condor

## I Setup 

```sh
cmsrel CMSSW_11_3_0_pre3
cd CMSSW_11_3_0_pre3/src/
cmsenv
git cms-checkout-topic -u Jingyan95:L1TK-dev-11_3_0_pre3_Condor
scram b -j 8
cd L1Trigger/TrackFindingTracklet/test/
```

## II Directory 

Modify the following two lines to set up your working space
<br />https://github.com/Jingyan95/cmssw/blob/L1TK-dev-11_3_0_pre3_Condor/L1Trigger/TrackFindingTracklet/test/testing01.py#L31
<br />https://github.com/Jingyan95/cmssw/blob/L1TK-dev-11_3_0_pre3_Condor/L1Trigger/TrackFindingTracklet/test/testing01.py#L35

## III Grid Proxy 

Modify the following line to export your proxy
<br />https://github.com/Jingyan95/cmssw/blob/L1TK-dev-11_3_0_pre3_Condor/L1Trigger/TrackFindingTracklet/test/testing01.py#L34
<br />I usually copy the proxy to my home directory
```sh
voms-proxy-init --voms cms -valid 192:00
cp /tmp/x509up_u122075 ~/x509up_u122075
```
## IV Make & Submit Condor Jobs

```sh
python testing01.py
```
