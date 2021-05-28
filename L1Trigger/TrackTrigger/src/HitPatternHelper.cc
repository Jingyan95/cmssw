//
//  Created by J.Li on 1/23/21.
//

#include "L1Trigger/TrackTrigger/interface/HitPatternHelper.h"
#include <algorithm>
#include <cmath>

namespace HPH {

    SensorModule::SensorModule(bool isbarrel,
                               bool isPS,
                               int numColumns,
                               int layerid,
                               double r,
                               double z,
                               double pitchCol,
                               double tilt)
    :isbarrel_(isbarrel),
    isPS_(isPS),
    numColumns_(numColumns),
    layerid_(layerid),
    r_(r),
    z_(z),
    pitchCol_(pitchCol),
    tilt_(tilt)
    {
      sin_ = std::sin(tilt_);
      cos_ = std::cos(tilt_);
    }

    Setup::Setup(const edm::ParameterSet& iConfig,
                 const TrackerGeometry& trackerGeometry,
                 const TrackerTopology& trackerTopology)
    :trackerGeometry_(&trackerGeometry),
    trackerTopology_(&trackerTopology)
    {
     
        iConfig_ = iConfig;
        
        for (const auto& gd : trackerGeometry_->dets()) {
          DetId detid = (*gd).geographicalId();
          if (detid.subdetId() != StripSubdetector::TOB && detid.subdetId() != StripSubdetector::TID)
            continue;  // only run on OT
          if (!trackerTopology_->isLower(detid))
            continue;                              // loop on the stacks: choose the lower arbitrarily
    
          // Get the DetSets of the Clusters
          const GeomDetUnit* det0 = trackerGeometry_->idToDetUnit(detid);
          const auto* theGeomDet = dynamic_cast<const PixelGeomDetUnit*>(det0);
          const PixelTopology* topol = dynamic_cast<const PixelTopology*>(&(theGeomDet->specificTopology()));
          const GlobalPoint pos0 = det0->position();
          const GlobalPoint pos1 = trackerGeometry_->idToDetUnit(trackerTopology_->partnerDetId(detid))->position();
          double r = pos0.perp();
          double z = pos0.z();
    
          bool flipped = pos0.mag() > pos1.mag();
          bool isbarrel = detid.subdetId() == StripSubdetector::TOB;
          bool isPS = trackerGeometry_->getDetectorType(detid) == TrackerGeometry::ModuleType::Ph2PSP;
          double tilt = flipped ? atan2(pos1.z() - pos0.z(), pos0.perp() - pos1.perp())
            : atan2(pos0.z() - pos1.z(), pos1.perp() - pos0.perp());
    
          int layerid = isbarrel ? trackerTopology_->layer(detid) : trackerTopology_->layer(detid)+10;
          int numColumns = topol->ncolumns();
          double pitchCol = topol->pitch().second;
    
          SensorModules_.emplace_back(isbarrel, isPS, numColumns, layerid, r, z, pitchCol, tilt);
    
        }
    
        sort(SensorModules_.begin(), SensorModules_.end(), smallerR);
        sort(SensorModules_.begin(), SensorModules_.end(), smallerZ);
        SensorModules_.erase(unique(SensorModules_.begin(), SensorModules_.end(), equalRZ), SensorModules_.end());
        
    }

    HitPatternHelper::HitPatternHelper(const Setup* setup,
                                       int hitpattern,
                                       double cot,
                                       double z0)
    :hitpattern_(hitpattern),
    cot_(cot),
    z0_(z0),
    numExpLayer_(0),
    numMissingLayer_(0),
    numMissingPS_(0),
    numMissing2S_(0),
    numPS_(0),
    num2S_(0),
    numMissingInterior1_(0),
    numMissingInterior2_(0),
    Setup_(setup),
    layers_(),
    binary_(11,0),
    HPHdebug_(Setup_->HPHdebug()),
    chosenRofZ_(Setup_->chosenRofZ())
    {
        
     float kfzRef = z0_ + chosenRofZ_ * cot_;
     int kf_eta_reg = 0;
     for (int iEtaSec = 1; iEtaSec < ((int)etaRegions_.size() - 1); iEtaSec++) {  // Doesn't apply eta < 2.4 cut.
         float etaMax = etaRegions_[iEtaSec];
         float zRefMax = chosenRofZ_ / tan(2. * atan(exp(-etaMax)));
         if (kfzRef > zRefMax){
          kf_eta_reg = iEtaSec;
         }
     }
     etaSector_ = kf_eta_reg;
     if (kf_eta_reg < ((int)etaRegions_.size()-1) / 2) {
         kf_eta_reg = ((int)etaRegions_.size()-1) / 2 - 1 - kf_eta_reg;
     }
     else {
         kf_eta_reg = kf_eta_reg - (int) (etaRegions_.size()-1) / 2;
     }
       
     for (SensorModule sm : Setup_->SensorModules()) {
         double d = (z0_ - sm.z() + sm.r() * cot_) / (sm.cos() - sm.sin() * cot_);
         if (abs(d) < sm.numColumns() * sm.pitchCol() / 2.){
             layers_.push_back(sm);
         }
     }
        
     sort( layers_.begin(), layers_.end(), smallerID );
     layers_.erase( unique( layers_.begin(), layers_.end(), equalID ), layers_.end() );
        
     numExpLayer_= layers_.size();
        
    int nbits = floor(log2(hitpattern_)) + 1;
    int lay_i = 0;
    bool seq = false;
    for (int i = 0; i < nbits; i++) {
        lay_i = ((1 << i) & hitpattern_) >> i;  //0 or 1 in ith bit (right to left)

        if (lay_i && !seq)
            seq = true;  //sequence starts when first 1 found
        if (!lay_i && seq){
            numMissingInterior1_++;
            bool realhit=false;
            for (int j : hitmap_[kf_eta_reg][i]){
                if (j<1) continue;
                int k = findLayer(j);
                if (k>0) realhit=true;
            }
            if (realhit) numMissingInterior2_++;
        }
    }
        
     for (int i=0;i<7;i++){
         if (HPHdebug_){
         std::cout<<"KF Layer = "<<i<<std::endl;
         }
         for (int j : hitmap_[kf_eta_reg][i]){
             if (HPHdebug_){
             std::cout<<" KF expects = "<<j<<std::endl;
             }
             if (j<1) continue;
             int k = findLayer(j);
             if (k<0) continue;
             if (HPHdebug_){
             std::cout<<" Confirmed by layermap"<<std::endl;
             }
             if (((1<<i)&hitpattern_)>>i){
                 if (HPHdebug_){
                 std::cout<<" Layer Found in HitPattern"<<std::endl;
                 }
                 binary_[ReducedId(j)] = 1;
                 if (layers_[k].isPS()){
                     numPS_++;
                 }
                 else{
                     num2S_++;
                 }
             }
             else{
                 if (HPHdebug_){
                 std::cout<<" Layer Missing in HitPattern"<<std::endl;
                 }
                 if (layers_[k].isPS()){
                     numMissingPS_++;
                 }
                 else{
                     numMissing2S_++;
                 }
             }
         }
     }
        if (HPHdebug_){
            if (HPHdebug_){
            std::cout<<"=================================================="<<std::endl;
            std::cout<<"hitpattern = "<<hitpattern_<<", numPS = "<<numPS_<<", num2S = "<<num2S_<<", missingPS = "<<numMissingPS_<<", missing2S = "<<numMissing2S_<<std::endl;
            std::cout<<"=================================================="<<std::endl;
            }
        std::cout<<"////////////////////////////////////////////////////////////////////////////////////////////////////////////////////"<<std::endl;
        }
    }

    int HitPatternHelper::ReducedId(int layerid){
        if (layerid<=6){
            layerid=layerid-1;
            return layerid;
        }
        else{
            layerid=layerid-5;
            return layerid;
        }
    };

    int HitPatternHelper::findLayer(int layerid){
        for (int i=0; i<(int)layers_.size();i++) {
            if (layerid == (int)layers_[i].layerid()){
                return i;
            }
         }
        return -1;
    }

}
