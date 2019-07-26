# L1 Tracking  

The L1 tracking can either be compiled standalone (outside CMSSW) or within CMSSW.

## To compile & run within CMSSW

```sh
cmsrel CMSSW_10_6_0
cd CMSSW_10_6_0/src
cmsenv

git init
git clone (-b myBranch) ssh://git@gitlab.cern.ch:7999/cms-tracker-phase2-backend-development/BE_software/L1Tracking.git L1Trigger

Optional: to update your branch with changes made to the master, since your branch was created: git pull --rebase origin master)

scram b -j 8
cd L1Trigger/TrackFindingTracklet/test/ 

cmsRun L1TrackNtupleMaker_cfg.py 

By default, the above runs on D21 samples. To run on D41 geometry, need to both change the flag in L1TrackNtupleMaker_cfg.py as well as change flag "geomTDR" to false in TrackFindingTracklet/interface/Constants.h

By default, it runs "Hybrid" L1 tracking emulation. Edit parameter USEHYBRID is Constants.h to instead use "Tracklet" emulation. Or edit parameter L1TRKALGO in L1TrackNtupleMaker_cfg.py to use either "TMTT" emulation, or simple floating point simulations of the Hybrid or Tracklet L1 tracking.
(N.B. To run the simple floating point Hybrid, you need to checkout additional code, as explained in
https://gitlab.cern.ch/cms-tracker-phase2-backend-development/BE_software/hybridfloat/blob/master/README.md).

```

## To compile & run standalone (currently not supported for Hybrid)

```sh
git clone https://gitlab.cern.ch/cms-tracker-phase2-backend-development/BE_software/L1Tracking.git
cd L1Tracking/TrackFindingTracklet/test/
make 
```

then to run over a file of single muon events for tilted barrel, do:

```sh
./fpga evlist_MuPt10_PU0_D4geom.txt 100 0
```

or to run and make selection on truth tracks (for efficiency / resolution plots), do: 
```sh
./fpga evlist_MuPt10_PU0_D4geom.txt 100 1
```

## Configuration parameters

Are hard-coded in TrackFindingTracklet/interface/Constants.h & TrackFindingTMTT/src/Settings.cc

## PLOTS 

### EFFICIENCY / RESOLUTION 

If running inside CMSSW, a ROOT TTree is created from the output TTracks & truth, with typical file name TTbar_PU200_hybrid.root . To make efficiency & resolution plots, and print out performance summary:

```sh
cd TrackFindingTracklet/test/
mkdir TrkPlots
root
.x L1TrackNtuplePlot.C("TTbar_PU200_hybrid")
```

(You can make similar plots for the tracks available before the track fit or before duplicate removal by disabling these algos.
To do this, edit Constants.h, setting fakefit=true and/or RemovalType="").

If running stand-alone, to make efficiency and resolution plots (that compare the integer based emulation to the floating point algorithm), set "writeResEff=true" in Constants.h to write a .txt file with the info,
and process it using macros in TrackFindingTracklet/test/PlotMacros/ . Warning: the efficiency isn't defined in the standard way.

```sh
root -l trackres.cc
root -l trackeff.cc
```

### DETAILED PERFORMANCE PLOTS (via hist utility)

Anders Ryd talk https://indico.cern.ch/event/781634/sessions/305061/#20190314 explains how to make plots & access truth/assoc
info at any point in code. Enabled via "bookHisto=true" in Constants.h, and booked in class "HistImp". But truth matching only partially available.

### DETAILED PERFORMANCE PLOTS (via txt files)

To generate performance plots you need to enable the relevant output (by editing cfg param named below in Constants.h), and after run the root script (from TrackFindingTracklet/test/PlotMacros/) to generate the plots.

```sh
Constants.h            root script
==================================================================
writeResEff                .L trackres/eff.cc++           trackres/eff()
writeStubsLayer            .L stubslayer.cc++             stubslayer()
writeStubsLayerperSector   .L stubslayerpersector.cc++    stubslayerpersector()
writeVMOccupancy           .L vmstubs.cc++                vmstubs()
writeTE                    .L trackletengine.cc++         trackletengine()
writeAllStubs              .L allstubs.cc++               allstubs()
writeTrackletCalculator    .L trackletcalculator.cc++     trackletcalculator()
writeMatchCalculator       .L matchcalculator.cc++        matchcalculator()
writeTrackletPars          .L trackletpars.cc++           trackletpars() (*)
writeNeighborProj          .L neighborproj.cc++           neighborproj()
writeAllProjections        .L allprojections.cc++         allprojections()
writeVMProjections         .L vmprojections.cc++          vmprojections()
writeTrackProj             .L trackproj.cc++              trackproj()
writeME                    .L matchengine.cc++            matchengine()
writeProjectionTransceiver .L projectiontransceiver.cc++  projectiontransceiver()
writeMatchTransceiver      .L matchtransceiver.cc++       matchtransceiver()
writeNMatches              .L nmatches.cc++               nmatches()
writez0andrinv             .L z0_and_rinv.cc++            z0_and_rinv()
```

(*) Needs some fixing/plots not meaningful

```sh
root -l stubs.cc
root -l stub_layer.cc
root -l stubpairs.cc
root -l trackletcands.cc
root -l trackletlayers.cc
root -l neighborproj.cc
root -l vmprojections.cc
root -l vmmatches.cc
```

## OTHER (stand-alone only?)

To turn on/off writing the files that dump the memory content
of the FPGA memories change the 'writememfiles' variable in the
Constants.h file.

To clean up all the output files that were produced do:

```sh
make clean
```

To generate the fitpattern.txt file do:

```sh
sort hitpattern.txt | uniq | sort -r -n -k 3 > fitpatter.txt
```

## PRODUCING ROOT Tree (historic option only?)

To produce a ROOT-Tree with the output of the emulation:

   1) Search Constants.h for "USEROOT" (at the top) and uncomment "#define USEROOT"
		
   2) Search Makefile.inc for "ROOT-Tree" and uncomment the loading of the FPGAEvent_cxx.so

   3) Building the FPGAEvent classes in ROOT. In a ROOT session:
   	gROOT->ProcessLine(".L FPGAEvent.cxx+");

   4) Build the fpga emulation code with:
       make fpga
		 
   5) Run the code as normal.

The produced file, "myTest.root", will be a ROOT tree with the class structure defined in FPGAEvent.h.  
It includes all the found tracks, mc particle information, stub information.


## Mixing PU events with signal (NOT WORKING)

program mixPU mixes events and dumps them to stdout.
fpga.cc now has an option of picking the input file from stdin (you just need to use stdin as the input file name)
This way one can pipe the huge files from the mixing straight into the fpga.cc without too much of a hassle.

To run on 3 muon events mixing each with two PU events do

```sh
make mixPU
./mixPU evlist_muminus_2_10_20000.txt 3 evlist_minbias_140PU_100.txt 2 | ./fpga stdin 6 1
```
note that in this example 6 is 2*3, so you run on all 6 events produced by the mixer.

