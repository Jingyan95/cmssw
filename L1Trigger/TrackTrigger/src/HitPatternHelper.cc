//
//  Created by J.Li on 1/23/21.
//

#include "L1Trigger/TrackTrigger/interface/HitPatternHelper.h"
#include <algorithm>
#include <cmath>

HitPatternHelper::HitPatternHelper(std::vector<trackerDTC::SensorModule> const& sensorModules, int hitpattern, double cot, double z0)
:hitpattern_(hitpattern),
cot_(cot),
z0_(z0),
numExpLayer_(0),
numMissingLayer_(0),
numMissingPS_(0),
numMissing2S_(0),
numPS_(0),
num2S_(0),
sensorModules_(sensorModules),
layers_(),
binary_(11,0)
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
   
 for (trackerDTC::SensorModule sm : sensorModules_) {
     double d = (z0_ - sm.z() + sm.r() * cot_) / (sm.cos() - sm.sin() * cot_);
     if (abs(d) < sm.numColumns() * sm.pitchCol() / 2.){
         layers_.push_back(sm);
     }
 }
    
 sort( layers_.begin(), layers_.end(), smallerID );
 layers_.erase( unique( layers_.begin(), layers_.end(), equalID ), layers_.end() );
    
 numExpLayer_= layers_.size();
    
// if (kf_eta_reg==5) debug=true;
    
 for (int i=0;i<7;i++){
     if (debug){
     std::cout<<"KF Layer = "<<i<<std::endl;
     }
     for (int j : hitmap_[kf_eta_reg][i]){
         if (debug){
         std::cout<<" KF expects = "<<j<<std::endl;
         }
         if (j<1) continue;
         int k = findLayer(j);
         if (k<0) continue;
         if (debug){
         std::cout<<" Confirmed by layermap"<<std::endl;
         }
         if (((1<<i)&hitpattern_)>>i){
             if (debug){
             std::cout<<" Layer Found in HitPattern"<<std::endl;
             }
             binary_[ReducedId(j)] = 1;
             if (layers_[k].psModule()){
                 numPS_++;
             }
             else{
                 num2S_++;
             }
         }
         else{
             if (debug){
             std::cout<<" Layer Missing in HitPattern"<<std::endl;
             }
             if (layers_[k].psModule()){
                 numMissingPS_++;
             }
             else{
                 numMissing2S_++;
             }
         }
     }
 }
    if (debug){
        if (debug){
        std::cout<<"=================================================="<<std::endl;
        std::cout<<"hitpattern = "<<hitpattern_<<", numPS = "<<numPS_<<", num2S = "<<num2S_<<", missingPS = "<<numMissingPS_<<", missing2S = "<<numMissing2S_<<std::endl;
        std::cout<<"=================================================="<<std::endl;
        }
    std::cout<<"////////////////////////////////////////////////////////////////////////////////////////////////////////////////////"<<std::endl;
    }
}

int HitPatternHelper::ReducedId(int layerId){
    if (layerId<=6){
        layerId=layerId-1;
        return layerId;
    }
    else{
        layerId=layerId-5;
        return layerId;
    }
};

int HitPatternHelper::findLayer(int layerId){
    for (int i=0; i<(int)layers_.size();i++) {
        if (layerId == (int)layers_[i].layerId()){
            return i;
        }
     }
    return -1;
}
