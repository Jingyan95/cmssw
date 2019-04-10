# L1 Tracking 

The L1 tracking can either be compiled standalone (outside CMSSW) or within CMSSW.

## To compile & run within CMSSW

```sh
cmsrel CMSSW_10_4_0
cd CMSSW_10_4_0/src
cmsenv

git init
git clone https://gitlab.cern.ch/cms-tracker-phase2-backend-development/BE_software/L1Tracking.git L1Trigger
scramv1 b -j 8
cd L1Trigger/TrackFindingTracklet/test/ 

cmsRun L1TrackNtupleMaker_cfg.py 
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

## PLOTS 

### EFFICIENCY / RESOLUTION 

If running inside CMSSW, a ROOT TTree is created from the output TTracks & truth, with typical file name TTbar_PU200_hybrid.root . To make efficiency & resolution plots, and print out performance summary:

```sh
cd TrackFindingTracklet/test/
mkdir TrkPlots
root
.x L1TrackNtuplePlot.C("TTbar_PU200_hybrid")
```

If running stand-alone, to make efficiency and resolution plots (that compare the integer based emulation to the floating point algorithm), set "writeResEff=true" in FPGAConstants.h to write a .txt file with the info,
and process it using macros in TrackFindingTracklet/test/PlotMacros/ . Warning: the efficiency isn't defined in the standard way.

```sh
root -l trackres.cc
root -l trackeff.cc
```

### DETAILED PERFORMANCE PLOTS 

To generate performance plots you need to enable the relevant output (by editing cfg param named below in FPGAConstants.hh), and after run the root script (from TrackFindingTracklet/test/PlotMacros/) to generate the plots.

FPGAConstants.hh            root script
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


(*) Needs some fixing/plots not meaningful

root -l stubs.cc
root -l stub_layer.cc
root -l stubpairs.cc
root -l trackletcands.cc
root -l trackletlayers.cc
root -l neighborproj.cc
root -l vmprojections.cc
root -l vmmatches.cc

## OTHER (stand-alone only?)

To turn on/off writing the files that dump the memory content
of the FPGA memories change the 'writememfiles' variable in the
FPGAConstants.hh file.

To clean up all the output files that were produced do:

make clean

To generate the fitpattern.txt file do:
sort hitpattern.txt | uniq | sort -r -n -k 3 > fitpatter.txt


## PRODUCING ROOT Tree (historic option only?)

To produce a ROOT-Tree with the output of the emulation:

   1) Search FPGAConstants.hh for "USEROOT" (at the top) and uncomment "#define USEROOT"
		
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

make mixPU
./mixPU evlist_muminus_2_10_20000.txt 3 evlist_minbias_140PU_100.txt 2 | ./fpga stdin 6 1

note that in this example 6 is 2*3, so you run on all 6 events produced by the mixer.

