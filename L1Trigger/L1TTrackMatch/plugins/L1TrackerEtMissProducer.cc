// Original Author:  Emmanuelle Perez,40 1-A28,+41227671915,
//         Created:  Tue Nov 12 17:03:19 CET 2013
//Modified by Emily MacDonald, 30 Nov 2018

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkEtMissParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkEtMissParticleFwd.h"
#include "DataFormats/L1TVertex/interface/Vertex.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkPrimaryVertex.h" //new
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTTrackAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
// detector geometry
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
// TMVA
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

using namespace l1t;
using namespace TMVA;

class L1TrackerEtMissProducer : public edm::EDProducer {
public:

  typedef TTTrack< Ref_Phase2TrackerDigi_ >  L1TTTrackType;
  typedef std::vector< L1TTTrackType > L1TTTrackCollectionType;

  explicit L1TrackerEtMissProducer(const edm::ParameterSet&);
  ~L1TrackerEtMissProducer();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  bool BDTMethod;

  float L1Tk_maxZ0;	    // in cm
  float L1Tk_maxDeltaZ;	    // in cm
  float L1Tk_maxEta;
  float L1Tk_maxChi2dof;
  float L1Tk_maxBendchi2;
  float L1Tk_minPt;	    // in GeV
  int L1Tk_minNStubs;
  int L1Tk_minNStubsPS;  // minimum number of stubs in PS modules
  float L1Tk_maxPt;	    // in GeV
  int L1Tk_HighPtTracks; // saturate or truncate

  float TP_minPt;
  float TP_maxEta;
  float TP_maxZ0;
  float TP_maxDxy;
  int TP_minNStubs;
  int TP_minNLayers;

 // const edm::EDGetTokenT< VertexCollection > pvToken;
  const edm::EDGetTokenT< L1TkPrimaryVertexCollection > pvToken;
  const edm::EDGetTokenT<std::vector<TTTrack< Ref_Phase2TrackerDigi_ > > > trackToken;
  const edm::EDGetTokenT< std::vector < TrackingParticle > > trackingParticleToken;
  const edm::EDGetTokenT< TTStubAssociationMap < Ref_Phase2TrackerDigi_ > > ttStubMCTruthToken;
  const edm::EDGetTokenT< TTClusterAssociationMap< Ref_Phase2TrackerDigi_ > > ttClusterMCTruthToken;
  const edm::EDGetTokenT< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > > ttTrackMCTruthToken;
  const edm::EDGetTokenT< edmNew::DetSetVector < TTStub < Ref_Phase2TrackerDigi_ > > > ttStubToken;
};

///////////////
//constructor//
///////////////
L1TrackerEtMissProducer::L1TrackerEtMissProducer(const edm::ParameterSet& iConfig) :
pvToken(consumes<L1TkPrimaryVertexCollection>(iConfig.getParameter<edm::InputTag>("L1VertexInputTag"))),
trackToken(consumes< std::vector<TTTrack< Ref_Phase2TrackerDigi_> > > (iConfig.getParameter<edm::InputTag>("L1TrackInputTag"))),
trackingParticleToken(consumes< std::vector < TrackingParticle > >(iConfig.getParameter<edm::InputTag>("TrackingParticleInputTag"))),
ttStubMCTruthToken(consumes< TTStubAssociationMap < Ref_Phase2TrackerDigi_ > > (iConfig.getParameter<edm::InputTag>("MCTruthStubInputTag"))),
ttClusterMCTruthToken(consumes< TTClusterAssociationMap< Ref_Phase2TrackerDigi_ > >(iConfig.getParameter<edm::InputTag>("MCTruthClusterInputTag"))),
ttTrackMCTruthToken(consumes< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > > (iConfig.getParameter<edm::InputTag>("MCTruthTrackInputTag"))),
ttStubToken(consumes< edmNew::DetSetVector < TTStub < Ref_Phase2TrackerDigi_ > > > (iConfig.getParameter<edm::InputTag>("L1StubInputTag")))
{
  
  BDTMethod = (bool)iConfig.getParameter<bool>("BDTMethod");
   
  L1Tk_maxZ0 = (float)iConfig.getParameter<double>("L1Tk_maxZ0");
  L1Tk_maxDeltaZ = (float)iConfig.getParameter<double>("L1Tk_maxDeltaZ");
  L1Tk_maxChi2dof = (float)iConfig.getParameter<double>("L1Tk_maxChi2dof");
  L1Tk_maxBendchi2 = (float)iConfig.getParameter<double>("L1Tk_maxBendchi2");
  L1Tk_minPt = (float)iConfig.getParameter<double>("L1Tk_minPt");
  L1Tk_minNStubs = iConfig.getParameter<int>("L1Tk_minNStubs");
  L1Tk_minNStubsPS = iConfig.getParameter<int>("L1Tk_minNStubsPS");
  L1Tk_maxPt = (float)iConfig.getParameter<double>("L1Tk_maxPt");
  L1Tk_maxEta = (float)iConfig.getParameter<double>("L1Tk_maxEta");
  L1Tk_HighPtTracks = iConfig.getParameter<int>("L1Tk_HighPtTracks");
    
  produces<L1TkEtMissParticleCollection>("trkMET");
    
  TP_minPt = (float)iConfig.getParameter<double>("TP_minPt");
  TP_maxEta = (float)iConfig.getParameter<double>("TP_maxEta");
  TP_maxZ0 = (float)iConfig.getParameter<double>("TP_maxZ0");
  TP_maxDxy = (float)iConfig.getParameter<double>("TP_maxDxy");
  TP_minNStubs = iConfig.getParameter<int>("TP_minNStubs");
  TP_minNLayers = iConfig.getParameter<int>("TP_minNLayers");
    
  produces<L1TkEtMissParticleCollection>("trueMET");
}

//////////////
//destructor//
//////////////
L1TrackerEtMissProducer::~L1TrackerEtMissProducer() {
}

////////////
//producer//
////////////
void L1TrackerEtMissProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  
  // Tracker Topology
  edm::ESHandle<TrackerTopology> tTopoHandle_;
  iSetup.get<TrackerTopologyRcd>().get(tTopoHandle_);
  const TrackerTopology* tTopo = tTopoHandle_.product();

  std::unique_ptr<L1TkEtMissParticleCollection> METCollection(new L1TkEtMissParticleCollection);
  std::unique_ptr<L1TkEtMissParticleCollection> tp_METCollection(new L1TkEtMissParticleCollection);

  edm::Handle<L1TkPrimaryVertexCollection> L1VertexHandle;
  iEvent.getByToken(pvToken,L1VertexHandle);

  edm::Handle<L1TTTrackCollectionType> L1TTTrackHandle;
  iEvent.getByToken(trackToken, L1TTTrackHandle);
  L1TTTrackCollectionType::const_iterator trackIter;
    
  edm::Handle< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > > TTStubHandle;
  iEvent.getByToken(ttStubToken, TTStubHandle);
    
  // Truth Association Maps
  edm::Handle< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > > MCTruthTTTrackHandle;
  iEvent.getByToken(ttTrackMCTruthToken, MCTruthTTTrackHandle);
    
  edm::Handle< TTClusterAssociationMap< Ref_Phase2TrackerDigi_ > > MCTruthTTClusterHandle;
  iEvent.getByToken(ttClusterMCTruthToken, MCTruthTTClusterHandle);
    
  edm::Handle< TTStubAssociationMap< Ref_Phase2TrackerDigi_ > > MCTruthTTStubHandle;
  iEvent.getByToken(ttStubMCTruthToken, MCTruthTTStubHandle);
    
  // Tracking Particles
  edm::Handle< std::vector < TrackingParticle > > trackingParticleHandle;
  iEvent.getByToken(trackingParticleToken, trackingParticleHandle);
  std::vector< TrackingParticle >::const_iterator iterTP;
  

  if( !L1TTTrackHandle.isValid() ) {
    LogError("L1TrackerEtMissProducer")<< "\nWarning: L1TTTrackCollection not found in the event. Exit"<< std::endl;
    return;
  }

  if( !L1VertexHandle.isValid() ) {
    LogError("L1TrackerEtMissProducer")<< "\nWarning: VertexCollection not found in the event. Exit"<< std::endl;
    return;
  }

  if( !TTStubHandle.isValid() ) {
      LogError("L1TrackerEtMissProducer")<< "\nWarning: TTStubCollection not found in the event. Exit"<< std::endl;
      return;
  }
    
  if( !MCTruthTTTrackHandle.isValid() ) {
      LogError("L1TrackerEtMissProducer")<< "\nWarning: MCTruthTTTrackCollection not found in the event. Exit"<< std::endl;
      return;
  }
    
  if( !MCTruthTTClusterHandle.isValid() ) {
      LogError("L1TrackerEtMissProducer")<< "\nWarning: MCTruthTTClusterCollection not found in the event. Exit"<< std::endl;
      return;
  }
    
  if( !MCTruthTTStubHandle.isValid() ) {
      LogError("L1TrackerEtMissProducer")<< "\nWarning: MCTruthTTStubCollection not found in the event. Exit"<< std::endl;
      return;
  }
    
  if( !trackingParticleHandle.isValid() ) {
      LogError("L1TrackerEtMissProducer")<< "\nWarning: trackingParticleCollection not found in the event. Exit"<< std::endl;
      return;
  }
    
 // if( !L1TTTrackHandle.isValid() ) {
   // LogError("L1TrackerEtMissProducer")<< "\nWarning: L1TTTrackCollection not found in the event. Exit"<< std::endl;
   // return;
 // }


  float sumPx = 0;
  float sumPy = 0;
  float etTot = 0;
  double sumPx_PU = 0;
  double sumPy_PU = 0;
  double etTot_PU = 0;

 // float zVTX = L1VertexHandle->begin()->z0();
  float zVTX = L1VertexHandle->begin()->getZvertex();
    
    float pt;
    float eta;
    float nstub;
    float chi2dof;
    float bendchi2;
    float z0;
    float DeltaZ;
    float seed;
    
    float BDT_response;
    float BDT_response1;
    
    // TMVA Library+
    TMVA::Tools::Instance();
    
    TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
    reader->AddVariable( "trk_pt := trk_pt",&pt);
    reader->AddVariable( "trk_eta := trk_eta", &eta);
    reader->AddVariable( "trk_nstub := trk_nstub", &nstub);
    reader->AddVariable( "trk_chi2 := trk_chi2", &chi2dof);
    reader->AddVariable( "trk_bend_chi2 := trk_bend_chi2", &bendchi2);
    reader->AddVariable( "trk_z0 := trk_z0",&z0);
    reader->AddVariable( "trk_DeltaZ := trk_DeltaZ",&DeltaZ);
    reader->AddVariable( "trk_seed := trk_seed",&seed);
    reader->BookMVA( "BDT method", "/afs/cern.ch/user/j/jingyan/CMS_NEW/CMSSW_10_5_0_pre1/src/L1Trigger/TrackFindingTracklet/test/TMVAClassification_BDT850_D21_hybrid_4_1.weights.xml");
    
    TMVA::Reader *reader1 = new TMVA::Reader( "!Color:!Silent" );
    reader1->AddVariable( "trk_pt := trk_pt",&pt);
    reader1->AddVariable( "trk_eta := trk_eta", &eta);
    reader1->AddVariable( "trk_nstub := trk_nstub", &nstub);
    reader1->AddVariable( "trk_chi2 := trk_chi2", &chi2dof);
    reader1->AddVariable( "trk_bend_chi2 := trk_bend_chi2", &bendchi2);
    reader1->AddVariable( "trk_z0 := trk_z0",&z0);
    reader1->AddVariable( "trk_DeltaZ := trk_DeltaZ",&DeltaZ);
    reader1->BookMVA( "BDT method1", "/afs/cern.ch/user/j/jingyan/CMS_NEW/CMSSW_10_5_0_pre1/src/L1Trigger/TrackFindingTracklet/test/TMVAClassification_BDT850_D21_hybrid_4_2.weights.xml");

  for (trackIter = L1TTTrackHandle->begin(); trackIter != L1TTTrackHandle->end(); ++trackIter) {
    pt = trackIter->getMomentum().perp();
    float phi = trackIter->getMomentum().phi();
    eta = trackIter->getMomentum().eta();
    std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > >  theStubs = trackIter -> getStubRefs() ;
    int nstubs = (int) theStubs.size();
    nstub = nstubs;

    float chi2 = trackIter->getChi2();
    chi2dof = chi2 / (2*nstubs-4);
    bendchi2 = trackIter->getStubPtConsistency();
    z0  = trackIter->getPOCA().z();
    DeltaZ = fabs(z0-zVTX);
    int seeds = (int) trackIter->getWedge();
    seed = seeds;

    if (pt < L1Tk_minPt) continue;
    if (fabs(eta) > L1Tk_maxEta) continue;
    if (fabs(z0) > L1Tk_maxZ0) continue;
    if (nstubs < L1Tk_minNStubs) continue;
      
    int nPS = 0.;     // number of stubs in PS modules
    // loop over the stubs
    for (unsigned int istub=0; istub<(unsigned int)theStubs.size(); istub++) {
          DetId detId( theStubs.at(istub)->getDetId() );
          if (detId.det() == DetId::Detector::Tracker) {
              if ( (detId.subdetId() == StripSubdetector::TOB && tTopo->tobLayer(detId) <= 3) || (detId.subdetId() == StripSubdetector::TID && tTopo->tidRing(detId) <= 9) ) nPS++;
          }
      }
      
    if (nPS < L1Tk_minNStubsPS) continue;
      
    if ( L1Tk_maxPt > 0 && pt > L1Tk_maxPt)  {
       if (L1Tk_HighPtTracks == 0)  continue;    // ignore these very high PT tracks: truncate
       if (L1Tk_HighPtTracks == 1)  pt = L1Tk_maxPt; // saturate
    }
    
    if(BDTMethod){
      BDT_response = reader->EvaluateMVA("BDT method");
      if (BDT_response<0.0991) continue;
      BDT_response1=1;
      //BDT_response1 = reader1->EvaluateMVA("BDT method1");
      if ( BDT_response1>0) {
         sumPx += pt*cos(phi);
         sumPy += pt*sin(phi);
         etTot += pt ;
        }
      else {    // PU sums
           sumPx_PU += pt*cos(phi);
           sumPy_PU += pt*sin(phi);
           etTot_PU += pt ;
        }
      }
    else{
    if (chi2dof > L1Tk_maxChi2dof) continue;
    if (bendchi2 > L1Tk_maxBendchi2) continue;
//    BDT_response = reader->EvaluateMVA("BDT method");
//    if (BDT_response<0.0648) continue;
    // construct deltaZ cut to be based on track eta
    if      ( fabs(eta)>=0   &&  fabs(eta)<0.7)  L1Tk_maxDeltaZ = 0.4;
    else if ( fabs(eta)>=0.7 &&  fabs(eta)<1.0)  L1Tk_maxDeltaZ = 0.6;
    else if ( fabs(eta)>=1.0 &&  fabs(eta)<1.2)  L1Tk_maxDeltaZ = 0.76;
    else if ( fabs(eta)>=1.2 &&  fabs(eta)<1.6)  L1Tk_maxDeltaZ = 1.0;
    else if ( fabs(eta)>=1.6 &&  fabs(eta)<2.0)  L1Tk_maxDeltaZ = 1.7;
    else if ( fabs(eta)>=2.0 &&  fabs(eta)<=2.4) L1Tk_maxDeltaZ = 2.2;

    if ( DeltaZ <= L1Tk_maxDeltaZ) {
      sumPx += pt*cos(phi);
      sumPy += pt*sin(phi);
      etTot += pt ;
    }
    else {    // PU sums
      sumPx_PU += pt*cos(phi);
      sumPy_PU += pt*sin(phi);
      etTot_PU += pt ;
    }
   }
  } // end loop over tracks
    
  float tp_sumPx = 0;
  float tp_sumPy = 0;
  float tp_etTot = 0;
  double tp_sumPx_PU = 0;
  double tp_sumPy_PU = 0;
  double tp_etTot_PU = 0;
    
  int this_tp = 0;
  
  float tmp_tp_pt;
  float tmp_tp_eta;
  float tmp_tp_nstub;
  float tmp_tp_chi2dof;
  float tmp_tp_bendchi2;
  float tmp_tp_z0;
  float tmp_tp_DeltaZ;
  float tmp_tp_seed;
    
    float BDT_response2;
    float BDT_response3;
    
    TMVA::Reader *reader2 = new TMVA::Reader( "!Color:!Silent" );
    reader2->AddVariable( "trk_pt := trk_pt",&tmp_tp_pt);
    reader2->AddVariable( "trk_eta := trk_eta", &tmp_tp_eta);
    reader2->AddVariable( "trk_nstub := trk_nstub", &tmp_tp_nstub);
    reader2->AddVariable( "trk_chi2 := trk_chi2", &tmp_tp_chi2dof);
    reader2->AddVariable( "trk_bend_chi2 := trk_bend_chi2", &tmp_tp_bendchi2);
    reader2->AddVariable( "trk_z0 := trk_z0",&tmp_tp_z0);
    reader2->AddVariable( "trk_DeltaZ := trk_DeltaZ",&tmp_tp_DeltaZ);
    reader2->AddVariable( "trk_seed := trk_seed",&tmp_tp_seed);
    reader2->BookMVA( "BDT method2", "/afs/cern.ch/user/j/jingyan/CMS_NEW/CMSSW_10_5_0_pre1/src/L1Trigger/TrackFindingTracklet/test/TMVAClassification_BDT850_D21_hybrid_4_1.weights.xml");
    
    TMVA::Reader *reader3 = new TMVA::Reader( "!Color:!Silent" );
    reader3->AddVariable( "trk_pt := trk_pt",&tmp_tp_pt);
    reader3->AddVariable( "trk_eta := trk_eta", &tmp_tp_eta);
    reader3->AddVariable( "trk_nstub := trk_nstub", &tmp_tp_nstub);
    reader3->AddVariable( "trk_chi2 := trk_chi2", &tmp_tp_chi2dof);
    reader3->AddVariable( "trk_bend_chi2 := trk_bend_chi2", &tmp_tp_bendchi2);
    reader3->AddVariable( "trk_z0 := trk_z0",&tmp_tp_z0);
    reader3->AddVariable( "trk_DeltaZ := trk_DeltaZ",&tmp_tp_DeltaZ);
    reader3->BookMVA( "BDT method3", "/afs/cern.ch/user/j/jingyan/CMS_NEW/CMSSW_10_5_0_pre1/src/L1Trigger/TrackFindingTracklet/test/TMVAClassification_BDT850_D21_hybrid_4_2.weights.xml");
    
    for (iterTP = trackingParticleHandle->begin(); iterTP != trackingParticleHandle->end(); ++iterTP) {
        edm::Ptr< TrackingParticle > tp_ptr(trackingParticleHandle, this_tp);
        this_tp++;

        bool tmp_tp_signal = (iterTP->eventId().event() == 0);

        tmp_tp_pt  = iterTP->pt();
        tmp_tp_eta = iterTP->eta();
        float tmp_tp_phi = iterTP->phi();
        float tmp_tp_vz  = iterTP->vz();
        float tmp_tp_vx  = iterTP->vx();
        float tmp_tp_vy  = iterTP->vy();

        if (tmp_tp_pt < TP_minPt) continue;
        if (fabs(tmp_tp_eta) > TP_maxEta) continue;

        float tmp_tp_t = tan(2.0*atan(1.0)-2.0*atan(exp(-tmp_tp_eta)));
        float delx = -tmp_tp_vx;
        float dely = -tmp_tp_vy;
        float A = 0.01*0.5696;
        float Kmagnitude = A / tmp_tp_pt;
        float tmp_tp_charge = iterTP->charge();

        if (tmp_tp_charge==0) continue;

        float K = Kmagnitude * tmp_tp_charge;
        float tmp_tp_x0p = delx - (1./(2. * K)*sin(tmp_tp_phi));
        float tmp_tp_y0p = dely + (1./(2. * K)*cos(tmp_tp_phi));
        float tmp_tp_rp = sqrt(tmp_tp_x0p*tmp_tp_x0p + tmp_tp_y0p*tmp_tp_y0p);
        float tmp_tp_d0 = tmp_tp_charge*tmp_tp_rp - (1. / (2. * K));
        tmp_tp_d0 = tmp_tp_d0*(-1); //fix d0 sign

        static double pi = 4.0*atan(1.0);
        float delphi = tmp_tp_phi-atan2(-K*tmp_tp_x0p,K*tmp_tp_y0p);
        if (delphi<-pi) delphi+=2.0*pi;
        if (delphi>pi) delphi-=2.0*pi;
        tmp_tp_z0 = tmp_tp_vz+tmp_tp_t*delphi/(2.0*K);

        if (fabs(tmp_tp_z0) > TP_maxZ0) continue;

        float dxy = sqrt(tmp_tp_vx*tmp_tp_vx + tmp_tp_vy*tmp_tp_vy);
        if (dxy > 1.0) continue;

        // only consider TPs associated with >= 1 cluster, or >= X stubs, or have stubs in >= X layers (configurable options)
        if (MCTruthTTClusterHandle.isValid()){
            if (MCTruthTTClusterHandle->findTTClusterRefs(tp_ptr).size() < 1) continue;
        }

        std::vector< edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > > theStubRefs = MCTruthTTStubHandle->findTTStubRefs(tp_ptr);
        int nStubTP = (int) theStubRefs.size();
        tmp_tp_nstub = nStubTP;
        tmp_tp_chi2dof = 0;
        tmp_tp_bendchi2 = 0;
        tmp_tp_DeltaZ = fabs(tmp_tp_z0);
        tmp_tp_seed = 0;

        // how many layers/disks have stubs?
        int hasStubInLayer[11] = {0};
        for (unsigned int is=0; is<theStubRefs.size(); is++) {
            DetId detid( theStubRefs.at(is)->getDetId() );
            int layer = -1;
            if ( detid.subdetId()==StripSubdetector::TOB ) {
                layer = static_cast<int>(tTopo->layer(detid)) - 1; //fill in array as entries 0-5
            }
            else if ( detid.subdetId()==StripSubdetector::TID ) {
                layer = static_cast<int>(tTopo->layer(detid)) + 5; //fill in array as entries 6-10
            }

            //treat genuine stubs separately (==2 is genuine, ==1 is not)
            if (MCTruthTTStubHandle->findTrackingParticlePtr(theStubRefs.at(is)).isNull() && hasStubInLayer[layer]<2) hasStubInLayer[layer] = 1;
            else hasStubInLayer[layer] = 2;
        }

        int nStubLayerTP = 0;
        for (int isum=0; isum<11; isum++) {
            if ( hasStubInLayer[isum] >= 1) nStubLayerTP   += 1;
        }

        if (TP_minNStubs>0 && nStubTP<TP_minNStubs) continue;
        if (TP_minNLayers>0 && nStubLayerTP<TP_minNLayers) continue;
        
        if(BDTMethod){
        BDT_response2 = reader2->EvaluateMVA("BDT method2");
        BDT_response3 = reader3->EvaluateMVA("BDT method3");
        if (BDT_response2<0&&BDT_response3<0) continue;
        }
        else{
        if (tmp_tp_signal!=1) continue;
//        float TP_maxDeltaZ=1;
//        if      ( fabs(tmp_tp_eta)>=0   &&  fabs(tmp_tp_eta)<0.7)  TP_maxDeltaZ = 0.4;
//        else if ( fabs(tmp_tp_eta)>=0.7 &&  fabs(tmp_tp_eta)<1.0)  TP_maxDeltaZ = 0.6;
//        else if ( fabs(tmp_tp_eta)>=1.0 &&  fabs(tmp_tp_eta)<1.2)  TP_maxDeltaZ = 0.76;
//        else if ( fabs(tmp_tp_eta)>=1.2 &&  fabs(tmp_tp_eta)<1.6)  TP_maxDeltaZ = 1.0;
//        else if ( fabs(tmp_tp_eta)>=1.6 &&  fabs(tmp_tp_eta)<2.0)  TP_maxDeltaZ = 1.7;
//        else if ( fabs(tmp_tp_eta)>=2.0 &&  fabs(tmp_tp_eta)<=2.4) TP_maxDeltaZ = 2.2;
//
//        if ( tmp_tp_DeltaZ > TP_maxDeltaZ) continue;
            
        }

        tp_sumPx += tmp_tp_pt*cos(tmp_tp_phi);
        tp_sumPx += tmp_tp_pt*sin(tmp_tp_phi);
    } //end loop tracking particles
  
  
    
  
  // calculate trkMET
  float et = sqrt( sumPx*sumPx + sumPy*sumPy );
  double etmiss_PU = sqrt( sumPx_PU*sumPx_PU + sumPy_PU*sumPy_PU );

  math::XYZTLorentzVector missingEt( -sumPx, -sumPy, 0, et);

  int ibx = 0;
  METCollection->push_back( L1TkEtMissParticle( missingEt,
    L1TkEtMissParticle::kMET,
    etTot,
    etmiss_PU,
    etTot_PU,
    ibx ) );

  iEvent.put( std::move(METCollection), "trkMET");
    
  // calculate TrueMET
  float tp_et = sqrt( tp_sumPx*tp_sumPx + tp_sumPy*tp_sumPy );
  double tp_etmiss_PU = sqrt( tp_sumPx_PU*tp_sumPx_PU + tp_sumPy_PU*tp_sumPy_PU );
    
  math::XYZTLorentzVector tp_missingEt( -tp_sumPx, -tp_sumPy, 0, tp_et);
    
  tp_METCollection->push_back( L1TkEtMissParticle( tp_missingEt,
    L1TkEtMissParticle::kMET,
    tp_etTot,
    tp_etmiss_PU,
    tp_etTot_PU,
    ibx ) );
    
  iEvent.put( std::move(tp_METCollection), "trueMET");
    
} // end producer

void L1TrackerEtMissProducer::beginJob() {
}

void L1TrackerEtMissProducer::endJob() {
}

DEFINE_FWK_MODULE(L1TrackerEtMissProducer);

