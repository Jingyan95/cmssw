//
//  Created by J.Li on 1/23/21.
//

#ifndef L1Trigger_TrackTrigger_interface_HitPatternHelper_h
#define L1Trigger_TrackTrigger_interface_HitPatternHelper_h

#include <stdio.h>
#include <iostream>
#include <vector>

#include "L1Trigger/TrackerDTC/interface/SensorModule.h"

class HitPatternHelper{
    
   public:
    
    HitPatternHelper(std::vector<trackerDTC::SensorModule> const& sensorModules, int hitpattern, double cot, double z0);
    ~HitPatternHelper(){}
    
    int getetaSector(){return etaSector_;}
    int getnumExpLayer(){return numExpLayer_;} //how many layers KF expects?
    int getnumMissingPS(){return numMissingPS_;}// how many PS layers are missing?
    int getnumMissing2S(){return numMissing2S_;}// how many 2S layers are missing?
    int getnumPS(){return numPS_;}
    int getnum2S(){return num2S_;}
    std::vector<int> getbinary(){return binary_;}
    static auto smallerID(trackerDTC::SensorModule lhs, trackerDTC::SensorModule rhs){ return lhs.layerId() < rhs.layerId(); }
    static auto equalID(trackerDTC::SensorModule lhs, trackerDTC::SensorModule rhs){ return lhs.layerId() == rhs.layerId(); }
    
    int ReducedId(int layerId);
    int findLayer(int layerId);
    
   private:
   
    int etaSector_;
    bool debug = false;
    int hitpattern_;
    double cot_;
    double z0_;
    int numExpLayer_;
    int numMissingLayer_;
    int numMissingPS_;
    int numMissing2S_;
    int numPS_;
    int num2S_;
    std::vector<trackerDTC::SensorModule> sensorModules_;
    std::vector<trackerDTC::SensorModule> layers_;
    std::vector<int> binary_;
    std::vector<float> etaRegions_ = {-2.4, -2.08, -1.68, -1.26, -0.90, -0.62, -0.41, -0.20, 0.0, 0.20, 0.41, 0.62, 0.90, 1.26, 1.68, 2.08, 2.4};
    std::vector<int> hitmap_[8][7] = {
        {{1}, {2},  {3},    {4},  {5},     {6},     {0}},
        {{1}, {2},  {3},    {4},  {5},     {6},     {0}},
        {{1}, {2},  {3},    {4},  {5},     {6},     {0}},
        {{1}, {2},  {3},    {4},  {5},     {6},     {0}},
        {{1}, {2},  {3},    {4},  {5, 11}, {6, 12}, {13}},
        {{1}, {2},  {3, 4}, {11}, {12},    {13},    {14, 15}},
        {{1}, {2},  {11},   {12}, {13},    {14},    {15}},
        {{1}, {11}, {12},   {13}, {14},    {15},    {0}},
    };
    float chosenRofZ_ = 50.0;
    
};

#endif
