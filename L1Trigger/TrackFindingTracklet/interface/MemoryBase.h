//Base class for memory modules
#ifndef L1Trigger_TrackFindingTracklet_interface_MemoryBase_h
#define L1Trigger_TrackFindingTracklet_interface_MemoryBase_h

#include "Settings.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <fstream>
#include <sstream>
#include <cassert>
#include <bitset>

using namespace Trklet;
using namespace std;

class MemoryBase{

public:

  MemoryBase(string name, const Settings* const settings, unsigned int iSector):
  settings_(settings){
    name_=name;
    iSector_=iSector;
    bx_=0;
    event_=0;
  }

  virtual ~MemoryBase(){}

  string getName() const {return name_;}
  string getLastPartOfName() const {return name_.substr(name_.find_last_of('_')+1);}

  virtual void clean()=0;

  //method sets the layer and disk based on the name. pos is the position in the
  //memory name where the layer or disk is specified
  void initLayerDisk(unsigned int pos, int& layer, int& disk){

    string subname=name_.substr(pos,2);
    layer=0;
    disk=0;

    if (subname=="L1") layer=1;
    if (subname=="L2") layer=2;
    if (subname=="L3") layer=3;
    if (subname=="L4") layer=4;
    if (subname=="L5") layer=5;
    if (subname=="L6") layer=6;
    if (subname=="D1") disk=1;
    if (subname=="D2") disk=2;
    if (subname=="D3") disk=3;
    if (subname=="D4") disk=4;
    if (subname=="D5") disk=5;
    if (layer==0&&disk==0) {
      edm::LogPrint("Tracklet") << "Memoryname = "<<name_<<" subname = "<<subname
				<<" layer "<<layer<<" disk "<<disk;
    }
    assert((layer!=0)||(disk!=0));
  }

  unsigned int initLayerDisk(unsigned int pos){

    int layer,disk;
    initLayerDisk(pos,layer,disk);

    if (disk>0) return 5+disk;
    return layer-1;

  }


  
  // Based on memory name check if this memory is used for special seeding:
  // overlap is layer-disk seeding
  // extra is the L2L3 seeding
  // extended is the seeding for displaced tracks
  void initSpecialSeeding(unsigned int pos, bool& overlap, bool& extra, bool& extended) {
    
    overlap=false;
    extra=false;
    extended=false;

    string subname=name_.substr(pos,1);

    static const std::set<std::string> overlapset = {"X","Y","W","Q","R","S","T","Z","x","y","w","q","r","s","t","z"};
    overlap=overlapset.find(subname)!=overlapset.end();
    
    static const std::set<std::string> extraset = {"I","J","K","L"};
    extra=extraset.find(subname)!=extraset.end();

    static const std::set<std::string> extendedset = {"a","b","c","d","e","f","g","h","x","y","z","w","q","r","s","t"};
    extended=extendedset.find(subname)!=extendedset.end();
    
  }

  //Used for a hack below due to MAC OS case sensitiviy problem for files
  void findAndReplaceAll(std::string & data, std::string toSearch, std::string replaceStr){

    // Get the first occurrence
    size_t pos = data.find(toSearch);
    
    // Repeat till end is reached
    while( pos != std::string::npos){
      // Replace this occurrence of Sub String
      data.replace(pos, toSearch.size(), replaceStr);
      // Get the next occurrence from the current position
      pos =data.find(toSearch, pos + replaceStr.size());
    }
  }
  
  void openFile(bool first, std::string filebase){
    
    std::string fname=filebase;
    fname+=getName();

    findAndReplaceAll(fname,"PHIa","PHIaa");
    findAndReplaceAll(fname,"PHIb","PHIbb");
    findAndReplaceAll(fname,"PHIc","PHIcc");
    findAndReplaceAll(fname,"PHId","PHIdd");

    findAndReplaceAll(fname,"PHIx","PHIxx");
    findAndReplaceAll(fname,"PHIy","PHIyy");
    findAndReplaceAll(fname,"PHIz","PHIzz");
    findAndReplaceAll(fname,"PHIw","PHIww");

    
    fname+="_";
    ostringstream oss;
    oss << iSector_+1;
    if (iSector_+1<10) fname+="0";
    fname+=oss.str();
    fname+=".dat";

    if (first) {
      bx_ = 0;
      event_ = 1;
      out_.open(fname.c_str());
    } else {
      out_.open(fname.c_str(),std::ofstream::app);
    }
      
    out_ << "BX = "<<(bitset<3>)bx_ << " Event : " << event_ << endl;

    bx_++;
    event_++;
    if (bx_>7) bx_=0;
  }

  
protected:

  string name_;
  unsigned int iSector_;
  
  ofstream out_;
  int bx_;
  int event_;

  const Settings* const settings_;

};

#endif
