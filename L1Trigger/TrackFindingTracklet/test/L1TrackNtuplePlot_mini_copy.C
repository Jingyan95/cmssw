// ----------------------------------------------------------------------------------------------------------------
// Basic example script for making tracking performance plots using the ntuples produced by L1TrackNtupleMaker.cc
// By Louise Skinnari, June 2013  
// ----------------------------------------------------------------------------------------------------------------

#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"
#include <TError.h>


#include <iostream>
#include <string>
#include <vector>

using namespace std;

void SetPlotStyle();
void mySmallText(Double_t x,Double_t y,Color_t color,char *text);

double getIntervalContainingFractionOfEntries( TH1* histogram, double interval );
void makeResidualIntervalPlot( TString type, TString dir, TString variable, TH1F* h_68, TH1F* h_90, TH1F* h_99, double minY, double maxY );


// ----------------------------------------------------------------------------------------------------------------
// Main script
// ----------------------------------------------------------------------------------------------------------------


void L1TrackNtuplePlot_mini(TString type, int TP_select_pdgid=0, int TP_select_eventid=0,
			    float TP_minPt=2.0, float TP_maxPt=100.0, float TP_maxEta=2.4) {

  // type:              this is the input file you want to process (minus ".root" extension)
  // TP_select_pdgid:   if non-zero, only select TPs with a given PDG ID
  // TP_select_eventid: if zero, only look at TPs from primary interaction, else, include TPs from pileup
  // TP_minPt:          only look at TPs with pt > X GeV
  // TP_maxPt:          only look at TPs with pt < X GeV
  // TP_maxEta:         only look at TPs with |eta| < X
 
  gROOT->SetBatch();
  gErrorIgnoreLevel = kWarning;
  
  SetPlotStyle();
    
  TFile* fout = new TFile("output_"+type+".root","recreate"); 
  
  
  // ----------------------------------------------------------------------------------------------------------------
  // define input options

  // baseline cuts for efficiency and rate plots ==> configure as appropriate
  int L1Tk_minNstub = 4;  
  float L1Tk_maxChi2 = 999999;  
  float L1Tk_maxChi2dof = 999999.;  
  
  //some counters for integrated efficiencies
  int n_all_eta2p5 = 0;
  int n_all_eta1p75 = 0;
  int n_all_eta1p0 = 0;
  int n_match_eta2p5 = 0;
  int n_match_eta1p75 = 0;
  int n_match_eta1p0 = 0;

  // counters for total track rates 
  int ntrk_pt2 = 0;
  int ntrk_pt3 = 0;
  int ntrk_pt10 = 0;
  int ntp_pt2 = 0;
  int ntp_pt3 = 0;
  int ntp_pt10 = 0;


  // ----------------------------------------------------------------------------------------------------------------
  // read ntuples
  TChain* tree = new TChain("L1TrackNtuple/eventTree");
  tree->Add(type+".root");
  
  if (tree->GetEntries() == 0) {
    cout << "File doesn't exist or is empty, returning..." << endl;
    return;
  }
  

  // ----------------------------------------------------------------------------------------------------------------
  // define leafs & branches

  // tracking particles
  vector<float>* tp_pt;
  vector<float>* tp_eta;
  vector<float>* tp_phi;
  vector<float>* tp_dxy;
  vector<float>* tp_z0;
  vector<float>* tp_d0;
  vector<int>*   tp_pdgid;
  vector<int>*   tp_nmatch;
  vector<int>*   tp_nstub;
  vector<int>*   tp_eventid;
  
  // *L1 track* properties, for tracking particles matched to a L1 track
  vector<float>* matchtrk_pt;
  vector<float>* matchtrk_eta;
  vector<float>* matchtrk_phi;
  vector<float>* matchtrk_d0;
  vector<float>* matchtrk_z0;
  vector<float>* matchtrk_chi2;
  vector<int>*   matchtrk_nstub;
  //vector<int>*   matchtrk_seed;

  // all L1 tracks
  vector<float>* trk_pt;
  vector<float>* trk_eta;
  vector<float>* trk_phi;
  vector<float>* trk_chi2;
  vector<int>*   trk_nstub;
  //vector<int>*   trk_seed;
  vector<int>*   trk_fake;
  vector<int>*   trk_genuine;
  vector<int>*   trk_loose;

  TBranch* b_tp_pt;
  TBranch* b_tp_eta;
  TBranch* b_tp_phi;
  TBranch* b_tp_dxy;
  TBranch* b_tp_z0;
  TBranch* b_tp_d0;
  TBranch* b_tp_pdgid;
  TBranch* b_tp_nmatch;
  TBranch* b_tp_nstub;
  TBranch* b_tp_eventid;

  TBranch* b_matchtrk_pt;
  TBranch* b_matchtrk_eta;
  TBranch* b_matchtrk_phi;
  TBranch* b_matchtrk_d0;
  TBranch* b_matchtrk_z0;
  TBranch* b_matchtrk_chi2; 
  TBranch* b_matchtrk_nstub;
  //TBranch* b_matchtrk_seed;

  TBranch* b_trk_pt; 
  TBranch* b_trk_eta; 
  TBranch* b_trk_phi; 
  TBranch* b_trk_chi2;
  TBranch* b_trk_chi2_fake;
  TBranch* b_trk_nstub; 
  //TBranch* b_trk_seed; 
  TBranch* b_trk_fake; 
  TBranch* b_trk_genuine; 
  TBranch* b_trk_loose; 

  tp_pt  = 0;
  tp_eta = 0;
  tp_phi = 0;
  tp_dxy = 0;
  tp_z0  = 0;
  tp_d0  = 0;
  tp_pdgid = 0;
  tp_nmatch = 0;
  tp_nstub = 0;
  tp_eventid = 0;

  matchtrk_pt  = 0;
  matchtrk_eta = 0;
  matchtrk_phi = 0;
  matchtrk_d0  = 0;
  matchtrk_z0  = 0;
  matchtrk_chi2  = 0; 
  matchtrk_nstub = 0;
  //matchtrk_seed = 0;

  trk_pt = 0; 
  trk_eta = 0; 
  trk_phi = 0; 
  trk_chi2 = 0; 
  trk_nstub = 0; 
  //trk_seed = 0; 
  trk_fake = 0; 
  trk_genuine = 0; 
  trk_loose = 0; 


  tree->SetBranchAddress("tp_pt",     &tp_pt,     &b_tp_pt);
  tree->SetBranchAddress("tp_eta",    &tp_eta,    &b_tp_eta);
  tree->SetBranchAddress("tp_phi",    &tp_phi,    &b_tp_phi);
  tree->SetBranchAddress("tp_dxy",    &tp_dxy,    &b_tp_dxy);
  tree->SetBranchAddress("tp_z0",     &tp_z0,     &b_tp_z0);
  tree->SetBranchAddress("tp_d0",     &tp_d0,     &b_tp_d0);
  tree->SetBranchAddress("tp_pdgid",  &tp_pdgid,  &b_tp_pdgid);
  tree->SetBranchAddress("tp_nmatch", &tp_nmatch, &b_tp_nmatch);
  tree->SetBranchAddress("tp_nstub",      &tp_nstub,      &b_tp_nstub);
  tree->SetBranchAddress("tp_eventid",    &tp_eventid,    &b_tp_eventid);

  tree->SetBranchAddress("matchtrk_pt",    &matchtrk_pt,    &b_matchtrk_pt);
  tree->SetBranchAddress("matchtrk_eta",   &matchtrk_eta,   &b_matchtrk_eta);
  tree->SetBranchAddress("matchtrk_phi",   &matchtrk_phi,   &b_matchtrk_phi);
  tree->SetBranchAddress("matchtrk_d0",    &matchtrk_d0,    &b_matchtrk_d0);
  tree->SetBranchAddress("matchtrk_z0",    &matchtrk_z0,    &b_matchtrk_z0);
  tree->SetBranchAddress("matchtrk_chi2",  &matchtrk_chi2,  &b_matchtrk_chi2);
  tree->SetBranchAddress("matchtrk_nstub", &matchtrk_nstub, &b_matchtrk_nstub);
  //tree->SetBranchAddress("matchtrk_seed",  &matchtrk_seed,  &b_matchtrk_seed);
  
  tree->SetBranchAddress("trk_pt",   &trk_pt,   &b_trk_pt);
  tree->SetBranchAddress("trk_eta",  &trk_eta,  &b_trk_eta);
  tree->SetBranchAddress("trk_phi",  &trk_phi,  &b_trk_phi);
  tree->SetBranchAddress("trk_chi2", &trk_chi2, &b_trk_chi2);
  tree->SetBranchAddress("trk_nstub",   &trk_nstub,   &b_trk_nstub);
  //tree->SetBranchAddress("trk_seed",    &trk_seed,    &b_trk_seed);
  tree->SetBranchAddress("trk_fake",    &trk_fake,    &b_trk_fake);
  tree->SetBranchAddress("trk_genuine", &trk_genuine, &b_trk_genuine);
  tree->SetBranchAddress("trk_loose",   &trk_loose,   &b_trk_loose);
    
  // Make TTree object
  TTree *t_fake = new TTree("t_fake","tree for fake track variables");
  TTree *t_real = new TTree("t_real","tree for real track variables");
  // variable
  float pt_fake;
  float eta_fake;
  float chi2_fake;
  int   nstub_fake;

  float pt_real;
  float eta_real;
  float chi2_real;
  int   nstub_real;

  // Tree branch
    
  t_fake->Branch("pt", &pt_fake, "Transverse Momentum/F");
  t_fake->Branch("eta", &eta_fake, "Pseudorapidity/F");
  t_fake->Branch("nstub", &nstub_fake, "Number of Stubs/I");
  t_fake->Branch("chi2", &chi2_fake, "Chi squared/F");
    
  t_real->Branch("pt", &pt_real, "Transverse Momentum/F");
  t_real->Branch("eta", &eta_real, "Pseudorapidity/F");
  t_real->Branch("nstub", &nstub_real, "Number of Stubs/I");
  t_real->Branch("chi2", &chi2_real, "Chi squared/F");
  

  // ----------------------------------------------------------------------------------------------------------------
  // histograms
  // ----------------------------------------------------------------------------------------------------------------

  /////////////////////////////////////////////////
  // NOTATION:                                   //
  // 'C' - Central eta range, |eta|<0.8          //
  // 'I' - Intermediate eta range, 0.8<|eta|<1.6 //
  // 'F' - Forward eta range, |eta|>1.6          //
  //                                             //
  // 'L' - Low pt range,  pt<8 GeV               //
  // 'H' - High pt range, pt>8 GeV               //
  /////////////////////////////////////////////////



  // ----------------------------------------------------------------------------------------------------------------
  // efficiencies

  TH1F* h_tp_pt    = new TH1F("tp_pt",   ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100,  0,   100.0);
  TH1F* h_tp_pt_L  = new TH1F("tp_pt_L", ";Tracking particle p_{T} [GeV]; Tracking particles / 0.1 GeV",  80,  0,     8.0);
  TH1F* h_tp_pt_H  = new TH1F("tp_pt_H", ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV",  92,  8.0, 100.0);
  TH1F* h_tp_eta   = new TH1F("tp_eta",  ";Tracking particle #eta; Tracking particles / 0.1",             50, -2.5,   2.5);
  TH1F* h_tp_eta_L = new TH1F("tp_eta_L",";Tracking particle #eta; Tracking particles / 0.1",             50, -2.5,   2.5);
  TH1F* h_tp_eta_H = new TH1F("tp_eta_H",";Tracking particle #eta; Tracking particles / 0.1",             50, -2.5,   2.5);

  TH1F* h_match_tp_pt    = new TH1F("match_tp_pt",   ";Tracking particle p_{T} [GeV]; Tracking particles / 1.0 GeV", 100,  0,   100.0);
  TH1F* h_match_tp_pt_L  = new TH1F("match_tp_pt_L", ";Tracking particle p_{T} [GeV]; Tracking particles / 0.1 GeV",  80,  0,     8.0);
  TH1F* h_match_tp_pt_H  = new TH1F("match_tp_pt_H", ";Tracking particle p_{T} [GeV]; Tracking particles / 0.1 GeV",  92,  8.0, 100.0);
  TH1F* h_match_tp_eta   = new TH1F("match_tp_eta",  ";Tracking particle #eta; Tracking particles / 0.1",             50, -2.5,   2.5);
  TH1F* h_match_tp_eta_L = new TH1F("match_tp_eta_L",";Tracking particle #eta; Tracking particles / 0.1",             50, -2.5,   2.5);
  TH1F* h_match_tp_eta_H = new TH1F("match_tp_eta_H",";Tracking particle #eta; Tracking particles / 0.1",             50, -2.5,   2.5);


  // ----------------------------------------------------------------------------------------------------------------
  // resolution vs eta histograms

  const int nETARANGE = 24;
  TString etarange[nETARANGE] = {"0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0",
				 "1.1","1.2","1.3","1.4","1.5","1.6","1.7","1.8","1.9","2.0",
				 "2.1","2.2","2.3","2.4"};

  TH1F* h_absResVsEta_eta[nETARANGE];
  TH1F* h_absResVsEta_z0[nETARANGE];
  TH1F* h_absResVsEta_phi[nETARANGE];
  TH1F* h_absResVsEta_ptRel[nETARANGE];

  // ----------------------------------------------
  // for new versions of resolution vs pt/eta plots
  unsigned int nBinsPtRelRes = 1000;
  double maxPtRelRes = 10.;

  unsigned int nBinsEtaRes = 500;
  double maxEtaRes = 0.1;

  unsigned int nBinsPhiRes = 500;
  double maxPhiRes = 0.2;

  unsigned int nBinsZ0Res = 100;
  double maxZ0Res = 4.0;
  // ----------------------------------------------

  for (int i=0; i<nETARANGE; i++) {
    h_absResVsEta_eta[i]   = new TH1F("absResVsEta_eta_"+etarange[i],  ";#eta residual (L1 - sim) [GeV]; L1 tracks / 0.1",     nBinsEtaRes,   0, maxEtaRes);
    h_absResVsEta_z0[i]    = new TH1F("absResVsEta_z0_"+etarange[i],   ";|z_{0} residual (L1 - sim)| [cm]; L1 tracks / 0.01",  nBinsZ0Res,    0, maxZ0Res);
    h_absResVsEta_phi[i]   = new TH1F("absResVsEta_phi_"+etarange[i],  ";#phi residual (L1 - sim) [GeV]; L1 tracks / 0.1",     nBinsPhiRes,   0, maxPhiRes);    
    h_absResVsEta_ptRel[i] = new TH1F("absResVsEta_ptRel_"+etarange[i],";p_{T} residual (L1 - sim) / p_{T}; L1 tracks / 0.02", nBinsPtRelRes, 0, maxPtRelRes);
  }

  // ----------------------------------------------------------------------------------------------------------------
  // total track rates 

  TH1F* h_trk_all_vspt        = new TH1F("trk_all_vspt",       ";Track p_{T} [GeV]; ",50,0,25);
  TH1F* h_trk_genuine_vspt    = new TH1F("trk_genuine_vspt",   ";Track p_{T} [GeV]; ",50,0,25);
  TH1F* h_tp_vspt             = new TH1F("tp_vspt",            ";TP p_{T} [GeV]; ",   50,0,25);

  // ----------------------------------------------------------------------------------------------------------------
  // number of tracks per event 

  // all tracks
  TH1F* h_ntrk_pt2  = new TH1F("ntrk_pt2",  ";# tracks (p_{T} > 2 GeV) / event; Events",  400, 0, 400.0);
  TH1F* h_ntrk_pt3  = new TH1F("ntrk_pt3",  ";# tracks (p_{T} > 3 GeV) / event; Events",  300, 0, 300.0);
  TH1F* h_ntrk_pt10 = new TH1F("ntrk_pt10", ";# tracks (p_{T} > 10 GeV) / event; Events", 100, 0, 100.0);
    
  // ----------------------------------------------------------------------------------------------------------------
  // Practice: fake/pri/sec
  TH1F* h_trk_eta_fake   = new TH1F("trk_eta_fake",       ";Fake Tracks #eta; ",50,-2.5,2.5);
  TH1F* h_trk_eta_pri    = new TH1F("trk_eta_pri",   ";Primary interaction #eta; ",50,-2.5,2.5);
  TH1F* h_trk_eta_sec    = new TH1F("trk_eta_sec",            ";Pileup #eta; ",   50,-2.5,2.5);
  TH1F* h_trk_chi2       = new TH1F("trk_chi2",            ";All tracks #chi^{2}; ",   50,0,150);
  TH1F* h_matchtrk_chi2  = new TH1F("matchtrk_chi2",            ";Matched tracks #chi^{2}; ",   50,0,150);
  TH1F* h_trk_nstub      = new TH1F("trk_nstub",            ";number of stubs (all tracks); ",   50,0,10);
  TH1F* h_tp_nstub       = new TH1F("tp_nstub",            ";number of stubs (tracking particles); ",   50,0,15);
  TH1F* h_matchtrk_nstub = new TH1F("matchtrk_nstub",            ";number of stubs (matched tracks); ",   50,0,10);
    
  // ----------------------------------------------------------------------------------------------------------------
  //   Practice 2: fake/pri/sec
  TH1F* h_trk_chi2_fake   = new TH1F("trk_chi2_fake",       ";Fake Tracks #chi^{2}; ",50,0,150);
  TH1F* h_trk_chi2_pri    = new TH1F("trk_chi2_pri",   ";Primary interaction #chi^{2}; ",50,0,150);
  TH1F* h_trk_chi2_sec    = new TH1F("trk_chi2_sec",            ";Pileup #chi^{2}; ",   50,0,150);
  TH1F* h_trk_nstub_fake  = new TH1F("trk_nstub_fake",       ";Fake Tracks nstub; ",50,0,10);
  TH1F* h_trk_nstub_pri   = new TH1F("trk_nstub_pri",   ";Primary interaction nstub; ",50,0,10);
  TH1F* h_trk_nstub_sec   = new TH1F("trk_nstub_sec",            ";Pileup nstub; ",   50,0,10);
  TH1F* h_trk_eta_H_fake  = new TH1F("trk_eta_H_fake",       ";Fake Tracks #eta; ",50,-2.5,2.5);
  TH1F* h_trk_eta_H_pri   = new TH1F("trk_eta_H_pri",   ";Primary interaction #eta; ",50,-2.5,2.5);
  TH1F* h_trk_eta_H_sec   = new TH1F("trk_eta_H_sec",            ";Pileup #eta; ",   50,-2.5,2.5);
  TH1F* h_trk_eta_L_fake  = new TH1F("trk_eta_L_fake",       ";Fake Tracks #eta; ",50,-2.5,2.5);
  TH1F* h_trk_eta_L_pri   = new TH1F("trk_eta_L_pri",   ";Primary interaction #eta; ",50,-2.5,2.5);
  TH1F* h_trk_eta_L_sec   = new TH1F("trk_eta_L_sec",            ";Pileup #eta; ",   50,-2.5,2.5);
    

  // ----------------------------------------------------------------------------------------------------------------
  //        * * * * *     S T A R T   O F   A C T U A L   R U N N I N G   O N   E V E N T S     * * * * *
  // ----------------------------------------------------------------------------------------------------------------
  
  int nevt = tree->GetEntries();
  cout << "number of events = " << nevt << endl;


  // ----------------------------------------------------------------------------------------------------------------
  // event loop
  for (int i=0; i<nevt; i++) {

    tree->GetEntry(i,0);
  
  
    // ----------------------------------------------------------------------------------------------------------------
    // loop over all L1 tracks in the event

    // counters for number of tracks / event with pt > 2/3/10 GeV for filling histograms
    int ntrkevt_pt2 = 0;
    int ntrkevt_pt3 = 0;
    int ntrkevt_pt10 = 0;

    for (int it=0; it<(int)trk_pt->size(); it++) {

      // basic kinematic selection criteria 
      if (fabs(trk_eta->at(it)) > TP_maxEta) continue;
      //only machine learning
      if (trk_pt->at(it) > TP_maxPt) continue;
      if (trk_pt->at(it) < TP_minPt) continue;

      if (trk_pt->at(it) > 2.0) {
	ntrk_pt2++;
	ntrkevt_pt2++;
	h_trk_all_vspt->Fill(trk_pt->at(it));
    h_trk_chi2->Fill(trk_chi2->at(it));
    h_trk_nstub->Fill(trk_nstub->at(it));
	if (trk_genuine->at(it) == 1) {
	  h_trk_genuine_vspt->Fill(trk_pt->at(it));
	}
          if (trk_fake->at(it)==0) {h_trk_eta_fake->Fill(trk_eta->at(it));
      h_trk_chi2_fake->Fill(trk_chi2->at(it));

              h_trk_nstub_fake->Fill(trk_nstub->at(it));
          }
          if (trk_fake->at(it)==1) {h_trk_eta_pri->Fill(trk_eta->at(it));
      h_trk_chi2_pri->Fill(trk_chi2->at(it));
              h_trk_nstub_pri->Fill(trk_nstub->at(it));
          }
          if (trk_fake->at(it)==2) {h_trk_eta_sec->Fill(trk_eta->at(it));
      h_trk_chi2_sec->Fill(trk_chi2->at(it));
              h_trk_nstub_sec->Fill(trk_nstub->at(it));
          }
          
          if (trk_fake->at(it)==0) {eta_fake=(trk_eta->at(it));
              chi2_fake=trk_chi2->at(it);
              nstub_fake=trk_nstub->at(it);
              pt_fake=trk_pt->at(it);
              t_fake->Fill();
          }
          if (trk_fake->at(it)>0) {eta_real=(trk_eta->at(it));
              chi2_real=trk_chi2->at(it);
              nstub_real=trk_nstub->at(it);
              pt_real=trk_pt->at(it);
              t_real->Fill();
          }
          
          if (trk_pt->at(it)>=5.0) {
          if (trk_fake->at(it)==0) h_trk_eta_H_fake->Fill(trk_eta->at(it));
          else if (trk_fake->at(it)==1) h_trk_eta_H_pri->Fill(trk_eta->at(it));
          else h_trk_eta_H_sec->Fill(trk_eta->at(it));}
          if (trk_pt->at(it)<=5.0) {
          if (trk_fake->at(it)==0) h_trk_eta_L_fake->Fill(trk_eta->at(it));
          else if (trk_fake->at(it)==1) h_trk_eta_L_pri->Fill(trk_eta->at(it));
          else h_trk_eta_L_sec->Fill(trk_eta->at(it));}
      
      }
      if (trk_pt->at(it) > 3.0) {
	ntrk_pt3++;
      }
      if (trk_pt->at(it) > 10.0) {
	ntrk_pt10++;
      }

    }

    h_ntrk_pt2->Fill(ntrkevt_pt2);
    h_ntrk_pt3->Fill(ntrkevt_pt3);
    h_ntrk_pt10->Fill(ntrkevt_pt10);


    // ----------------------------------------------------------------------------------------------------------------
    // loop over all tracking particle in the event

    for (int it=0; it<(int)tp_pt->size(); it++) {
      
      // cut on PDG ID at plot stage?
      if (TP_select_pdgid != 0) {
	if (abs(tp_pdgid->at(it)) != abs(TP_select_pdgid)) continue;
      }

      // kinematic cuts
      if (tp_dxy->at(it) > 1) continue;
      if (tp_pt->at(it) < 0.2) continue;
      if (tp_pt->at(it) > TP_maxPt) continue;
      if (fabs(tp_eta->at(it)) > TP_maxEta) continue;

      
      // total track rates
      if (tp_pt->at(it) > TP_minPt) {
	if (tp_pt->at(it) > 2.0) {
	  ntp_pt2++;
	  h_tp_vspt->Fill(tp_pt->at(it));
	}
	if (tp_pt->at(it) > 3.0) ntp_pt3++;
	if (tp_pt->at(it) > 10.0) ntp_pt10++;
      }

      // cut on event ID (eventid=0 means the TP is from the primary interaction, so *not* selecting only eventid=0 means including stuff from pileup)
      if (TP_select_eventid == 0 && tp_eventid->at(it) != 0) continue;


      h_tp_pt->Fill(tp_pt->at(it));
      if (tp_pt->at(it) < 8.0) h_tp_pt_L->Fill(tp_pt->at(it));
      else h_tp_pt_H->Fill(tp_pt->at(it));
      
      if (tp_pt->at(it) > TP_minPt) {
		
	if (fabs(tp_eta->at(it)) < 1.0) n_all_eta1p0++;
	else if (fabs(tp_eta->at(it)) < 1.75) n_all_eta1p75++;
	else n_all_eta2p5++;
	
	h_tp_eta->Fill(tp_eta->at(it));
    h_tp_nstub->Fill(tp_nstub->at(it));
	if (tp_pt->at(it) < 8.0) h_tp_eta_L->Fill(tp_eta->at(it));
	else h_tp_eta_H->Fill(tp_eta->at(it));
	
      }
      
      
      // ----------------------------------------------------------------------------------------------------------------
      // was the tracking particle matched to a L1 track?
      if (tp_nmatch->at(it) < 1) continue;
      
      
      // use only tracks with min X stubs
      if (matchtrk_nstub->at(it) < L1Tk_minNstub) continue;

      // cut on chi2?
      if (matchtrk_chi2->at(it) > L1Tk_maxChi2) continue;
      int ndof = 2*matchtrk_nstub->at(it)-4; //number of degrees of freedom
      if (matchtrk_chi2->at(it)/ndof > L1Tk_maxChi2dof) continue;
      
      
      // ----------------------------------------------------------------------------------------------------------------
      // fill matched track histograms
      
      h_matchtrk_chi2->Fill(matchtrk_chi2->at(it));
      h_matchtrk_nstub->Fill(matchtrk_nstub->at(it));
      h_match_tp_pt->Fill(tp_pt->at(it));
      if (tp_pt->at(it) < 8.0) h_match_tp_pt_L->Fill(tp_pt->at(it));
      else h_match_tp_pt_H->Fill(tp_pt->at(it));
      
      // only look at tracking particles with pt > TP_minPt
      if (tp_pt->at(it) < TP_minPt) continue;

     
      if (fabs(tp_eta->at(it)) < 1.0) n_match_eta1p0++;
      else if (fabs(tp_eta->at(it)) < 1.75) n_match_eta1p75++;
      else n_match_eta2p5++;
     
      h_match_tp_eta->Fill(tp_eta->at(it));
      
      if (tp_pt->at(it) < 8.0) h_match_tp_eta_L->Fill(tp_eta->at(it));
      else h_match_tp_eta_H->Fill(tp_eta->at(it));
      

      // ----------------------------------------------------------------------------------------------------------------
      // fill resolution vs eta histograms

      for (int im=0; im<nETARANGE; im++) {
	if ( (fabs(tp_eta->at(it)) > (float)im*0.1) && (fabs(tp_eta->at(it)) < (float)im*0.1+0.1) ) {
	  h_absResVsEta_ptRel[im]->Fill( fabs( (matchtrk_pt->at(it) - tp_pt->at(it)) )/tp_pt->at(it));
	  h_absResVsEta_eta[im]  ->Fill( fabs( matchtrk_eta->at(it) - tp_eta->at(it) ) );
	  h_absResVsEta_phi[im]  ->Fill( fabs( matchtrk_phi->at(it) - tp_phi->at(it) ) );
	  h_absResVsEta_z0[im]   ->Fill( fabs( matchtrk_z0->at(it)  - tp_z0->at(it) ) );	 
	}
      }

      
    } // end of matched track loop
    
  } // end of event loop
  // ----------------------------------------------------------------------------------------------------------------
  

  // ----------------------------------------------------------------------------------------------------------------
  // 2D plots  
  // ----------------------------------------------------------------------------------------------------------------

  // resolution vs. eta histograms (68 / 90 / 99% residuals)
  TH1F* h2_resVsEta_eta_68  = new TH1F("resVsEta_eta_68",   ";Tracking particle |#eta|; #eta resolution", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_eta_90  = new TH1F("resVsEta_eta_90",   ";Tracking particle |#eta|; #eta resolution", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_eta_99  = new TH1F("resVsEta_eta_99",   ";Tracking particle |#eta|; #e ta resolution", nETARANGE,0,2.4);

  TH1F* h2_resVsEta_z0_68   = new TH1F("resVsEta_z0_68",   ";Tracking particle |#eta|; z_{0} resolution [cm]", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_z0_90   = new TH1F("resVsEta_z0_90",   ";Tracking particle |#eta|; z_{0} resolution [cm]", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_z0_99   = new TH1F("resVsEta_z0_99",   ";Tracking particle |#eta|; z_{0} resolution [cm]", nETARANGE,0,2.4);

  TH1F* h2_resVsEta_phi_68   = new TH1F("resVsEta_phi_68",   ";Tracking particle |#eta|; #phi resolution [rad]", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_phi_90   = new TH1F("resVsEta_phi_90",   ";Tracking particle |#eta|; #phi resolution [rad]", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_phi_99   = new TH1F("resVsEta_phi_99",   ";Tracking particle |#eta|; #phi resolution [rad]", nETARANGE,0,2.4);

  TH1F* h2_resVsEta_ptRel_68   = new TH1F("resVsEta_ptRel_68",   ";Tracking particle |#eta|; p_{T} resolution / p_{T}", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_ptRel_90   = new TH1F("resVsEta_ptRel_90",   ";Tracking particle |#eta|; p_{T} resolution / p_{T}", nETARANGE,0,2.4);
  TH1F* h2_resVsEta_ptRel_99   = new TH1F("resVsEta_ptRel_99",   ";Tracking particle |#eta|; p_{T} resolution / p_{T}", nETARANGE,0,2.4);


  for (int i=0; i<nETARANGE; i++) {

    h2_resVsEta_eta_68->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_eta[i], 0.68 ));
    h2_resVsEta_eta_90->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_eta[i], 0.90 ));
    h2_resVsEta_eta_99->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_eta[i], 0.99 ));

    h2_resVsEta_z0_68->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_z0[i], 0.68 ));
    h2_resVsEta_z0_90->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_z0[i], 0.90 ));
    h2_resVsEta_z0_99->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_z0[i], 0.99 ));

    h2_resVsEta_phi_68->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_phi[i], 0.68 ));
    h2_resVsEta_phi_90->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_phi[i], 0.90 ));
    h2_resVsEta_phi_99->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_phi[i], 0.99 ));

    h2_resVsEta_ptRel_68->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_ptRel[i], 0.68 ));
    h2_resVsEta_ptRel_90->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_ptRel[i], 0.90 ));
    h2_resVsEta_ptRel_99->SetBinContent(i+1, getIntervalContainingFractionOfEntries( h_absResVsEta_ptRel[i], 0.99 ));
    
  }


  // -------------------------------------------------------------------------------------------
  // output file for histograms
  // -------------------------------------------------------------------------------------------
 
  if (TP_select_pdgid != 0) {
    char pdgidtxt[500];
    sprintf(pdgidtxt,"_pdgid%i",TP_select_pdgid);
    type = type+pdgidtxt;
  }

  if (TP_select_eventid != 0) type = type+"_wpu";

  if (TP_minPt > 2.0) {
    char pttxt[500];
    sprintf(pttxt,"_pt%.0f",TP_minPt);
    type = type+pttxt;
  }

  
  // -------------------------------------------------------------------------------------------
  // draw and save plots
  // -------------------------------------------------------------------------------------------

  char ctxt[500];
  TCanvas c;

  TString DIR = "TrkPlots/";

  // plots overlaying 68, 90, 99% confidence levels]

  float max_eta_ptRel = 0.2;
  float max_z0 = 2.0;
  float max_phi = 0.01;
  float max_eta = 0.03;

  if (type.Contains("El")) {
    max_eta_ptRel = 1.0;
    max_phi = 0.1;
  }


  // ----------------------------------------------------------------------------------------------------------
  // resolution vs eta plots
  // ----------------------------------------------------------------------------------------------------------

  // makeResidualIntervalPlot will save the individual plots to the root file  
  makeResidualIntervalPlot( type, DIR, "resVsEta_eta", h2_resVsEta_eta_68, h2_resVsEta_eta_90, h2_resVsEta_eta_99, 0, max_eta );
  makeResidualIntervalPlot( type, DIR, "resVsEta_z0", h2_resVsEta_z0_68, h2_resVsEta_z0_90, h2_resVsEta_z0_99, 0, max_z0 );
  makeResidualIntervalPlot( type, DIR, "resVsEta_phi", h2_resVsEta_phi_68, h2_resVsEta_phi_90, h2_resVsEta_phi_99, 0, max_phi );
  makeResidualIntervalPlot( type, DIR, "resVsEta_ptRel", h2_resVsEta_ptRel_68, h2_resVsEta_ptRel_90, h2_resVsEta_ptRel_99, 0, max_eta_ptRel );

  
  // ----------------------------------------------------------------------------------------------------------------
  // efficiency plots  
  // ----------------------------------------------------------------------------------------------------------------

  // rebin pt/phi plots
  h_tp_pt->Rebin(4);
  h_match_tp_pt->Rebin(4);
  h_tp_pt_L->Rebin(2);
  h_match_tp_pt_L->Rebin(2);
  h_tp_pt_H->Rebin(2);
  h_match_tp_pt_H->Rebin(2);

  h_tp_eta->Rebin(2);
  h_match_tp_eta->Rebin(2);
  h_tp_eta_L->Rebin(2);
  h_match_tp_eta_L->Rebin(2);
  h_tp_eta_H->Rebin(2);
  h_match_tp_eta_H->Rebin(2);

  // calculate the effeciency
  h_match_tp_pt->Sumw2();
  h_tp_pt->Sumw2();
  TH1F* h_eff_pt = (TH1F*) h_match_tp_pt->Clone();
  h_eff_pt->SetName("eff_pt");
  h_eff_pt->GetYaxis()->SetTitle("Efficiency");
  h_eff_pt->Divide(h_match_tp_pt, h_tp_pt, 1.0, 1.0, "B");

  h_match_tp_pt_L->Sumw2();
  h_tp_pt_L->Sumw2();
  TH1F* h_eff_pt_L = (TH1F*) h_match_tp_pt_L->Clone();
  h_eff_pt_L->SetName("eff_pt_L");
  h_eff_pt_L->GetYaxis()->SetTitle("Efficiency");
  h_eff_pt_L->Divide(h_match_tp_pt_L, h_tp_pt_L, 1.0, 1.0, "B");

  h_match_tp_pt_H->Sumw2();
  h_tp_pt_H->Sumw2();
  TH1F* h_eff_pt_H = (TH1F*) h_match_tp_pt_H->Clone();
  h_eff_pt_H->SetName("eff_pt_H");
  h_eff_pt_H->GetYaxis()->SetTitle("Efficiency");
  h_eff_pt_H->Divide(h_match_tp_pt_H, h_tp_pt_H, 1.0, 1.0, "B");

  h_match_tp_eta->Sumw2();
  h_tp_eta->Sumw2();
  TH1F* h_eff_eta = (TH1F*) h_match_tp_eta->Clone();
  h_eff_eta->SetName("eff_eta");
  h_eff_eta->GetYaxis()->SetTitle("Efficiency");
  h_eff_eta->Divide(h_match_tp_eta, h_tp_eta, 1.0, 1.0, "B");
  
  h_match_tp_eta_L->Sumw2();
  h_tp_eta_L->Sumw2();
  TH1F* h_eff_eta_L = (TH1F*) h_match_tp_eta_L->Clone();
  h_eff_eta_L->SetName("eff_eta_L");
  h_eff_eta_L->GetYaxis()->SetTitle("Efficiency");
  h_eff_eta_L->Divide(h_match_tp_eta_L, h_tp_eta_L, 1.0, 1.0, "B");

  h_match_tp_eta_H->Sumw2();
  h_tp_eta_H->Sumw2();
  TH1F* h_eff_eta_H = (TH1F*) h_match_tp_eta_H->Clone();
  h_eff_eta_H->SetName("eff_eta_H");
  h_eff_eta_H->GetYaxis()->SetTitle("Efficiency");
  h_eff_eta_H->Divide(h_match_tp_eta_H, h_tp_eta_H, 1.0, 1.0, "B");


  // set the axis range
  h_eff_pt  ->SetAxisRange(0,1.1,"Y");
  h_eff_pt_L->SetAxisRange(0,1.1,"Y");
  h_eff_pt_H->SetAxisRange(0,1.1,"Y");
  h_eff_eta ->SetAxisRange(0,1.1,"Y");
  h_eff_eta_L ->SetAxisRange(0,1.1,"Y");
  h_eff_eta_H ->SetAxisRange(0,1.1,"Y");

  gPad->SetGridx();
  gPad->SetGridy();

//  // draw and save plots
//  h_eff_pt->Draw();
//  h_eff_pt->Write();
//  c.SaveAs(DIR+type+"_eff_pt.pdf");
//
//  if (type.Contains("Mu")) {
//    h_eff_pt->GetYaxis()->SetRangeUser(0.8,1.01); // zoomed-in plot
//    c.SaveAs(DIR+type+"_eff_pt_zoom.pdf");
//  }
//
//  h_eff_pt_L->Draw();
//  h_eff_pt_L->Write();
//  sprintf(ctxt,"p_{T} < 8 GeV");
//  mySmallText(0.45,0.5,1,ctxt);
//  c.SaveAs(DIR+type+"_eff_pt_L.pdf");
//
//  h_eff_pt_H->Draw();
//  h_eff_pt_H->Write();
//  sprintf(ctxt,"p_{T} > 8 GeV");
//  mySmallText(0.45,0.5,1,ctxt);
//  c.SaveAs(DIR+type+"_eff_pt_H.pdf");
//
//  h_eff_eta->Draw();
//  h_eff_eta->Write();
//  c.SaveAs(DIR+type+"_eff_eta.pdf");
//
//  if (type.Contains("Mu")) {
//    h_eff_eta->GetYaxis()->SetRangeUser(0.8,1.01); // zoomed-in plot
//    c.SaveAs(DIR+type+"_eff_eta_zoom.pdf");
//  }
//
//  h_eff_eta_L->Draw();
//  h_eff_eta_L->Write();
//  sprintf(ctxt,"p_{T} < 8 GeV");
//  mySmallText(0.45,0.5,1,ctxt);
//  c.SaveAs(DIR+type+"_eff_eta_L.pdf");
//
//  h_eff_eta_H->Draw();
//  h_eff_eta_H->Write();
//  sprintf(ctxt,"p_{T} > 8 GeV");
//  mySmallText(0.45,0.5,1,ctxt);
//  c.SaveAs(DIR+type+"_eff_eta_H.pdf");
//
//  h_trk_eta_fake->Draw();
//  h_trk_eta_fake->Write();
//  c.SaveAs(DIR+type+"_eta_fake.pdf");
//
//  h_trk_eta_pri->Draw();
//  h_trk_eta_pri->Write();
//  c.SaveAs(DIR+type+"_eta_pri.pdf");
//
//  h_trk_eta_sec->Draw();
//  h_trk_eta_sec->Write();
//  c.SaveAs(DIR+type+"_eta_sec.pdf");
//
//  h_tp_eta_H->Draw("hist");
//  h_tp_eta_H->Write();
//  sprintf(ctxt,"p_{T} > 8 GeV");
//  mySmallText(0.45,0.5,1,ctxt);
//  c.SaveAs(DIR+type+"_tp_eta_H.pdf");
//
//  h_tp_eta_L->Draw("hist");
//  h_tp_eta_L->Write();
//  sprintf(ctxt,"p_{T} < 8 GeV");
//  mySmallText(0.45,0.5,1,ctxt);
//  c.SaveAs(DIR+type+"_tp_eta_L.pdf");
//
//  h_trk_chi2->Draw("hist");
//  h_trk_chi2->Write();
//  c.SaveAs(DIR+type+"_trk_chi2.pdf");
//
//  h_matchtrk_chi2->Draw("hist");
//  h_matchtrk_chi2->Write();
//  c.SaveAs(DIR+type+"_matchtrk_chi2.pdf");
//
//  h_trk_nstub->Draw("hist");
//  h_trk_nstub->Write();
//  c.SaveAs(DIR+type+"_trk_nstub.pdf");
//
//  h_tp_nstub->Draw("hist");
//  h_tp_nstub->Write();
//  c.SaveAs(DIR+type+"_tp_nstub.pdf");
//
//  h_matchtrk_nstub->Draw("hist");
//  h_matchtrk_nstub->Write();
//  c.SaveAs(DIR+type+"_matchtrk_nstub.pdf");
//
////  TCanvas *c1 = new TCanvas("c1","#chi^{2}",500, 600,500,800);
////  c1->Divide(1,3);
////  c1->cd(1);
//  h_trk_chi2_fake->Draw("hist");
//  h_trk_chi2_fake->Write();
//  c.SaveAs(DIR+type+"_trk_chi2_fake.png");
////  c1->cd(2);
//  h_trk_chi2_pri->Draw("hist");
//  h_trk_chi2_pri->Write();
//  c.SaveAs(DIR+type+"_trk_chi2_pri.png");
////  c1->cd(3);
//  h_trk_chi2_sec->Draw("hist");
//  h_trk_chi2_sec->Write();
//  c.SaveAs(DIR+type+"_trk_chi2_sec.png");
//
//    h_trk_nstub_fake->Draw("hist");
//    h_trk_nstub_fake->Write();
//    c.SaveAs(DIR+type+"_trk_nstub_fake.png");
//    //  c1->cd(2);
//    h_trk_nstub_pri->Draw("hist");
//    h_trk_nstub_pri->Write();
//    c.SaveAs(DIR+type+"_trk_nstub_pri.png");
//    //  c1->cd(3);
//    h_trk_nstub_sec->Draw("hist");
//    h_trk_nstub_sec->Write();
//    c.SaveAs(DIR+type+"_trk_nstub_sec.png");
//
//    h_trk_eta_H_fake->Draw("hist");
//    h_trk_eta_H_fake->Write();
//    sprintf(ctxt,"p_{T} > 5 GeV");
//    mySmallText(0.45,0.5,1,ctxt);
//    c.SaveAs(DIR+type+"_trk_eta_H_fake.png");
//    //  c1->cd(2);
//    h_trk_eta_H_pri->Draw("hist");
//    h_trk_eta_H_pri->Write();
//    sprintf(ctxt,"p_{T} > 5 GeV");
//    mySmallText(0.45,0.5,1,ctxt);
//    c.SaveAs(DIR+type+"_trk_eta_H_pri.png");
//    //  c1->cd(3);
//    h_trk_eta_H_sec->Draw("hist");
//    h_trk_eta_H_sec->Write();
//    sprintf(ctxt,"p_{T} > 5 GeV");
//    mySmallText(0.45,0.5,1,ctxt);
//    c.SaveAs(DIR+type+"_trk_eta_H_sec.png");
//
//    h_trk_eta_L_fake->Draw("hist");
//    h_trk_eta_L_fake->Write();
//    sprintf(ctxt,"p_{T} < 5 GeV");
//    mySmallText(0.45,0.5,1,ctxt);
//    c.SaveAs(DIR+type+"_trk_eta_L_fake.png");
//    //  c1->cd(2);
//    h_trk_eta_L_pri->Draw("hist");
//    h_trk_eta_L_pri->Write();
//    sprintf(ctxt,"p_{T} < 5 GeV");
//    mySmallText(0.45,0.5,1,ctxt);
//    c.SaveAs(DIR+type+"_trk_eta_L_pri.png");
//    //  c1->cd(3);
//    h_trk_eta_L_sec->Draw("hist");
//    h_trk_eta_L_sec->Write();
//    sprintf(ctxt,"p_{T} < 5 GeV");
//    mySmallText(0.45,0.5,1,ctxt);
//    c.SaveAs(DIR+type+"_trk_eta_L_sec.png");
    
//  TCanvas *c2 = new TCanvas("c2","nstub",900,700);
//  TCanvas *c3 = new TCanvas("c3","#eta_H",900,700);
//  TCanvas *c4 = new TCanvas("c4","#eta_L",900,700);


  gPad->SetGridx(0);
  gPad->SetGridy(0);

  

  // ---------------------------------------------------------------------------------------------------------
  // track rate plots
  // ---------------------------------------------------------------------------------------------------------

  h_trk_all_vspt->Sumw2();
  h_trk_genuine_vspt->Sumw2();
  h_tp_vspt->Sumw2();

  h_trk_all_vspt->Scale(1.0/nevt);
  h_trk_genuine_vspt->Scale(1.0/nevt);
  h_tp_vspt->Scale(1.0/nevt);

  h_tp_vspt->GetYaxis()->SetTitle("Tracks / event");
  h_tp_vspt->GetXaxis()->SetTitle("Track p_{T} [GeV]");
  h_tp_vspt->SetLineColor(4);
  h_tp_vspt->SetLineStyle(2);

  float max = h_tp_vspt->GetMaximum();
  if (h_trk_all_vspt->GetMaximum() > max) max = h_trk_all_vspt->GetMaximum();
  h_tp_vspt->SetAxisRange(0.001,max*1.5,"Y");  

  h_tp_vspt->Draw("hist");
  h_trk_all_vspt->Draw("same,hist");
  h_tp_vspt->Draw("same,hist");

//  h_trk_all_vspt->Write();
//  h_trk_genuine_vspt->Write();
//  h_tp_vspt->Write();

  char txt[500];
  sprintf(txt,"# tracks/event = %.1f",h_trk_all_vspt->GetSum());
  mySmallText(0.5,0.85,1,txt);
  char txt3[500];
  sprintf(txt3,"# TPs(stubs in #geq 4 layers)/");
  char txt2[500];
  sprintf(txt2,"event = %.1f",h_tp_vspt->GetSum());
  mySmallText(0.5,0.79,4,txt3);
  mySmallText(0.5,0.74,4,txt2);

  c.SaveAs(DIR+type+"_trackrate_vspt.pdf");

  gPad->SetLogy();
  c.SaveAs(DIR+type+"_trackrate_vspt_log.pdf");
  gPad->SetLogy(0);



  // nbr tracks per event 
//  h_ntrk_pt2->Write();
//  h_ntrk_pt3->Write();
//  h_ntrk_pt10->Write();
    
   // TMVA!!!!
    t_fake->Write();
    t_real->Write();
    
    
  
    

  // close output root ntuple file
  fout->Close();



  // ---------------------------------------------------------------------------------------------------------
  //some printouts
  // ---------------------------------------------------------------------------------------------------------

  float k = (float)n_match_eta1p0;
  float N = (float)n_all_eta1p0;
  if (fabs(N)>0) cout << endl << "efficiency for |eta| < 1.0 = " << k/N*100.0 << " +- " << 1.0/N*sqrt(k*(1.0 - k/N))*100.0 << endl;
  k = (float)n_match_eta1p75;
  N = (float)n_all_eta1p75;
  if (fabs(N)>0) cout << "efficiency for 1.0 < |eta| < 1.75 = " << k/N*100.0 << " +- " << 1.0/N*sqrt(k*(1.0 - k/N))*100.0 << endl;
  k = (float)n_match_eta2p5;
  N = (float)n_all_eta2p5;
  if (fabs(N)>0) cout << "efficiency for 1.75 < |eta| < 2.5 = " << k/N*100.0 << " +- " << 1.0/N*sqrt(k*(1.0 - k/N))*100.0 << endl;
  N = (float) n_all_eta1p0 + n_all_eta1p75 + n_all_eta2p5;
  k = (float) n_match_eta1p0 + n_match_eta1p75 + n_match_eta2p5;
  if (fabs(N)>0) cout << "combined efficiency for |eta| < 2.5 = " << k/N*100.0 << " +- " << 1.0/N*sqrt(k*(1.0 - k/N))*100.0 << endl << endl;

  // track rates
  cout << "# TP/event (pt > 2.0) = " << (float)ntp_pt2/nevt << endl;
  cout << "# TP/event (pt > 3.0) = " << (float)ntp_pt3/nevt << endl;
  cout << "# TP/event (pt > 10.0) = " << (float)ntp_pt10/nevt << endl;

  cout << "# tracks/event (pt > 2.0) = " << (float)ntrk_pt2/nevt << endl;
  cout << "# tracks/event (pt > 3.0) = " << (float)ntrk_pt3/nevt << endl;
  cout << "# tracks/event (pt > 10.0) = " << (float)ntrk_pt10/nevt << endl;

}


void SetPlotStyle() {

  // from ATLAS plot style macro

  // use plain black on white colors
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetHistLineColor(1);

  gStyle->SetPalette(1);

  // set the paper & margin sizes
  gStyle->SetPaperSize(20,26);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  gStyle->SetTitleXOffset(1.4);
  gStyle->SetTitleYOffset(1.4);

  // use large fonts
  gStyle->SetTextFont(42);
  gStyle->SetTextSize(0.05);
  gStyle->SetLabelFont(42,"x");
  gStyle->SetTitleFont(42,"x");
  gStyle->SetLabelFont(42,"y");
  gStyle->SetTitleFont(42,"y");
  gStyle->SetLabelFont(42,"z");
  gStyle->SetTitleFont(42,"z");
  gStyle->SetLabelSize(0.05,"x");
  gStyle->SetTitleSize(0.05,"x");
  gStyle->SetLabelSize(0.05,"y");
  gStyle->SetTitleSize(0.05,"y");
  gStyle->SetLabelSize(0.05,"z");
  gStyle->SetTitleSize(0.05,"z");

  // use bold lines and markers
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetHistLineWidth(2.);
  gStyle->SetLineStyleString(2,"[12 12]");

  // get rid of error bar caps
  gStyle->SetEndErrorSize(0.);

  // do not display any of the standard histogram decorations
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

}


void mySmallText(Double_t x,Double_t y,Color_t color,char *text) {
  Double_t tsize=0.044;
  TLatex l;
  l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x,y,text);
}

double getIntervalContainingFractionOfEntries( TH1* absResidualHistogram, double quantileToCalculate ) {
  
  double totalIntegral = absResidualHistogram->Integral( 0, absResidualHistogram->GetNbinsX() + 1 );

  // Check that the interval is not somewhere in the overflow bin
  double maxAllowedEntriesInOverflow = totalIntegral * ( 1 - quantileToCalculate );
  double nEntriesInOverflow = absResidualHistogram->GetBinContent( absResidualHistogram->GetNbinsX() + 1 );
  if ( nEntriesInOverflow > maxAllowedEntriesInOverflow ) {
        // cout << "WARNING : Cannot compute range corresponding to interval, as it is in the overflow bin" << endl;
        return absResidualHistogram->GetXaxis()->GetXmax()  * 1.2;
  }

  // Calculate quantile for given interval
  double interval[1];
  double quantile[1] = { quantileToCalculate };
  absResidualHistogram->GetQuantiles( 1, interval, quantile);

  return interval[0];
}

void makeResidualIntervalPlot( TString type, TString dir, TString variable, TH1F* h_68, TH1F* h_90, TH1F* h_99, double minY, double maxY ) {

  TCanvas c;

  h_68->SetMinimum( minY );
  h_90->SetMinimum( minY );
  h_99->SetMinimum( minY );

  h_68->SetMaximum( maxY );
  h_90->SetMaximum( maxY );
  h_99->SetMaximum( maxY );

  h_68->SetMarkerStyle(20);
  h_90->SetMarkerStyle(26);
  h_99->SetMarkerStyle(24);

//  h_68->Draw("P");
//  h_68->Write();
//  h_90->Draw("P same");
//  h_90->Write();
//  h_99->Draw("P same");
//  h_99->Write();

  TLegend* l = new TLegend(0.65,0.65,0.85,0.85);
  l->SetFillStyle(0);
  l->SetBorderSize(0);
  l->SetTextSize(0.04);
  l->AddEntry(h_99,"99%","p");
  l->AddEntry(h_90,"90%","p");
  l->AddEntry(h_68,"68%","p");
  l->SetTextFont(42);
  l->Draw();  

  //c.SaveAs(dir+type+"_"+variable+"_interval.png");
  c.SaveAs(dir+type+"_"+variable+"_interval.pdf");

  delete l;
}


