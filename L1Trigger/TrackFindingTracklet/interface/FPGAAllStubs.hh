// This class holds all the stubs in a DTC region for a give layer
#ifndef FPGAALLSTUBS_H
#define FPGAALLSTUBS_H

#include "L1TStub.hh"
#include "FPGAStub.hh"
#include "FPGAMemoryBase.hh"

#include <ctype.h>

using namespace std;

class FPGAAllStubs:public FPGAMemoryBase{

public:

  FPGAAllStubs(string name, unsigned int iSector, 
	       double phimin, double phimax):
    FPGAMemoryBase(name,iSector){
    phimin_=phimin;
    phimax_=phimax;
    string subname=name.substr(3,2);
    if (subname[0]=='_') subname=name.substr(3,2);
    
    layer_ = 0; 
    disk_  = 0;
 
    if (subname=="L1") layer_=1;
    if (subname=="L2") layer_=2;
    if (subname=="L3") layer_=3;
    if (subname=="L4") layer_=4;
    if (subname=="L5") layer_=5;
    if (subname=="L6") layer_=6;
    if (subname=="D1") disk_=1;
    if (subname=="D2") disk_=2;
    if (subname=="D3") disk_=3;
    if (subname=="D4") disk_=4;
    if (subname=="D5") disk_=5;
    if (subname=="F1") disk_=1;
    if (subname=="F2") disk_=2;
    if (subname=="F3") disk_=3;
    if (subname=="F4") disk_=4;
    if (subname=="F5") disk_=5;
    if (subname=="B1") disk_=-1;
    if (subname=="B2") disk_=-2;
    if (subname=="B3") disk_=-3;
    if (subname=="B4") disk_=-4;
    if (subname=="B5") disk_=-5;
    if (layer_==0&&disk_==0) {
      cout << name<<" subname = "<<subname<<" "<<layer_<<" "<<disk_<<endl;
    }   
    assert((layer_!=0)||(disk_!=0));

    assert(name.substr(5,3)=="PHI");
  }

  void addStub(std::pair<FPGAStub*,L1TStub*> stub) {
    stubs_.push_back(stub);
  }

  unsigned int nStubs() const {return stubs_.size();}

  FPGAStub* getFPGAStub(unsigned int i) const {return stubs_[i].first;}
  L1TStub* getL1TStub(unsigned int i) const {return stubs_[i].second;}
  std::pair<FPGAStub*,L1TStub*> getStub(unsigned int i) const {return stubs_[i];}

  void clean() {
    stubs_.clear();
  }

  void writeStubs(bool first) {

    std::string fname="MemPrints/Stubs/AllStubs_";
    fname+=getName();
    fname+="_";
    ostringstream oss;
    oss << iSector_+1;
    if (iSector_+1<10) fname+="0";
    fname+=oss.str();
    fname+=".dat";
    if (first) {
      bx_=0;
      event_=1;
      out_.open(fname.c_str());
    }
    else
      out_.open(fname.c_str(),std::ofstream::app);

    out_ << "BX = "<<(bitset<3>)bx_ << " Event : " << event_ << endl;

    for (unsigned int j=0;j<stubs_.size();j++){
      string stub= (layer_>0)? stubs_[j].first->str()
      : stubs_[j].first->strdisk();

      if (j<16) out_ <<"0";
      out_ << hex << j << dec ;
      out_ <<" "<<stub << endl;
    }
    out_.close();

    bx_++;
    event_++;
    if (bx_>7) bx_=0;


  }

  int layer() const { return layer_;}
  int disk() const { return disk_;}

private:

  double phimin_;
  double phimax_;
  std::vector<std::pair<FPGAStub*,L1TStub*> > stubs_;

  int layer_;
  int disk_;

};

#endif
