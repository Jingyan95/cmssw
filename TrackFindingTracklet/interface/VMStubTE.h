#ifndef VMSTUBTE_H
#define VMSTUBTE_H

#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <vector>

#include "FPGAWord.h"
#include "Util.h"
#include "Constants.h"
#include "Stub.h"
#include "L1TStub.h"

using namespace std;

class VMStubTE{

public:

  VMStubTE() {}
  

  VMStubTE(std::pair<Stub*, L1TStub*> stub, unsigned int finephi, FPGAWord finerz, FPGAWord bend, unsigned int vmbits, FPGAWord allstubindex) {
    stub_=stub;
    finephi_=finephi;
    finerz_=finerz;
    bend_=bend;
    vmbits_=vmbits;
    allStubIndex_=allstubindex;
  }

 
  ~VMStubTE() {

  }

  unsigned int finephi() const {
    return finephi_;
  }
  
  FPGAWord finerz() const {
    return finerz_;
  }

  FPGAWord bend() const {
    return bend_;
  }
  
  unsigned int vmbits() const {
    return vmbits_;
  }

  std::pair<Stub*, L1TStub*> stub() const {
    return stub_;
  }

  bool isPSmodule() const {
    return stub_.first->isPSmodule();
  }
  
private:

  unsigned int finephi_;
  FPGAWord finerz_;
  FPGAWord bend_;
  unsigned int vmbits_;
  FPGAWord allStubIndex_;
  std::pair<Stub*, L1TStub*> stub_;

};



#endif



