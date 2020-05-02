#ifndef L1Trigger_TrackFindingTracklet_interface_CleanTrackMemory_h
#define L1Trigger_TrackFindingTracklet_interface_CleanTrackMemory_h

#include "L1Trigger/TrackFindingTracklet/interface/MemoryBase.h"

#include <vector>

namespace Trklet {

  class Settings;
  class Tracklet;
  class L1SimTrack;

  class CleanTrackMemory : public MemoryBase {
  public:
    CleanTrackMemory(std::string name, const Settings* const settings, unsigned int iSector, double phimin, double phimax);

    void addTrack(Tracklet* tracklet) { tracks_.push_back(tracklet); }
    
    unsigned int nTracks() const { return tracks_.size(); }
    
    void clean() { tracks_.clear(); }
    
    bool foundTrack(std::ofstream& outres, L1SimTrack simtrk);
    
    void writeCT(bool first);

    
  private:
    double phimin_;
    double phimax_;
    std::vector<Tracklet*> tracks_;
  };
  
};
#endif
