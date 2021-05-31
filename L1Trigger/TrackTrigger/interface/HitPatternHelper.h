//
//  Created by J.Li on 1/23/21.
//

#ifndef L1Trigger_TrackTrigger_interface_HitPatternHelper_h
#define L1Trigger_TrackTrigger_interface_HitPatternHelper_h

#include "FWCore/Framework/interface/data_default_record_trait.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "L1Trigger/TrackTrigger/interface/HitPatternHelperRcd.h"

#include <stdio.h>
#include <iostream>
#include <vector>

using namespace std;
using namespace edm;

namespace HPH {

    class SensorModule{
        
       public:
        
        SensorModule(){}
        SensorModule(bool isbarrel,
                     bool isPS,
                     int numColumns,
                     int layerid,
                     double r,
                     double z,
                     double pitchCol,
                     double tilt);
        ~SensorModule() {}
        
        bool isbarrel() const { return isbarrel_; }
        bool isPS() const { return isPS_; }
        int numColumns() const { return numColumns_; }
        int layerid() const { return layerid_; }
        double r() const { return r_; }
        double z() const { return z_; }
        double pitchCol() const { return pitchCol_; }
        double tilt() const { return tilt_; }
        double sin() const { return sin_; }
        double cos() const { return cos_; }
        
       private:
        
        bool isbarrel_;
        bool isPS_;
        int numColumns_;
        int layerid_;
        double r_;
        double z_;
        double pitchCol_;
        double tilt_;
        double sin_;
        double cos_;

    };

    class Setup{
        
       public:
       
        Setup() {}
        Setup(const edm::ParameterSet& iConfig,
              const TrackerGeometry& trackerGeometry,
              const TrackerTopology& trackerTopology);
        ~Setup() {}
        
        static auto smallerR(SensorModule lhs, SensorModule rhs){ return lhs.r() < rhs.r(); }
        static auto smallerZ(SensorModule lhs, SensorModule rhs){ return lhs.z() < rhs.z(); }
        static auto equalRZ(SensorModule lhs, SensorModule rhs){ return abs(lhs.r() - rhs.r()) < delta_ && abs(lhs.z() - rhs.z()) < delta_; }
        std::vector<SensorModule> SensorModules() const {return SensorModules_;}
        
        bool HPHdebug() const { return iConfig_.getParameter<bool>("HPHdebug"); }
        bool useNewKF() const { return iConfig_.getParameter<bool>("useNewKF"); }
        double chosenRofZ() const { return iConfig_.getParameter<double>("chosenRofZ"); }
       
       private:
        
        edm::ParameterSet iConfig_;
        const TrackerGeometry* trackerGeometry_;
        const TrackerTopology* trackerTopology_;
        static constexpr double delta_ = 1.e-3;
        std::vector<SensorModule> SensorModules_;
        
     
    };

    class HitPatternHelper{
        
       public:
        
        HitPatternHelper(){}
        HitPatternHelper(const Setup* setup,
                         int hitpattern,
                         double cot,
                         double z0);
        ~HitPatternHelper(){}
        
        int getetaSector(){return etaSector_;}
        int getnumExpLayer(){return numExpLayer_;} //how many layers KF expects?
        int getnumMissingPS(){return numMissingPS_;}// how many PS layers are missing?
        int getnumMissing2S(){return numMissing2S_;}// how many 2S layers are missing?
        int getnumPS(){return numPS_;}
        int getnum2S(){return num2S_;}
        int getnumMissingInterior1(){return numMissingInterior1_;}
        int getnumMissingInterior2(){return numMissingInterior2_;}
        std::vector<int> getbinary(){return binary_;}
        static auto smallerID(SensorModule lhs, SensorModule rhs){ return lhs.layerid() < rhs.layerid(); }
        static auto equalID(SensorModule lhs, SensorModule rhs){ return lhs.layerid() == rhs.layerid(); }
        
        int ReducedId(int layerId);
        int findLayer(int layerId);
        
       private:
       
        int etaSector_;
        int hitpattern_;
        double cot_;
        double z0_;
        int numExpLayer_;
        int numMissingLayer_;
        int numMissingPS_;
        int numMissing2S_;
        int numPS_;
        int num2S_;
        int numMissingInterior1_;
        int numMissingInterior2_;
        const Setup* Setup_;
        std::vector<SensorModule> layers_;
        std::vector<int> binary_;
        bool HPHdebug_;
        bool useNewKF_;
        float chosenRofZ_;
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
        
    };

}

EVENTSETUP_DATA_DEFAULT_RECORD(HPH::Setup, HPH::SetupRcd);

#endif
