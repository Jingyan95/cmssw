#ifndef L1Trigger_TrackFindingTracklet_interface_TETableBase_h
#define L1Trigger_TrackFindingTracklet_interface_TETableBase_h

#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <vector>

using namespace std;

class TETableBase {
public:
  TETableBase() {}

  virtual ~TETableBase() {}

  virtual int lookup(int, int) {
    assert(0);  //Should never get here - this should be a pure virtual fcn
  };

  void writeVMTable(std::string name, bool positive = true) {
    ofstream out;
    out.open(name.c_str());
    out << "{" << endl;
    for (unsigned int i = 0; i < table_.size(); i++) {
      if (i != 0)
        out << "," << endl;

      assert(nbits_ > 0);

      int itable = table_[i];
      if (positive) {
        if (table_[i] < 0)
          itable = (1 << nbits_) - 1;
      }

      //FPGAWord tmp;
      //tmp.set(itable, nbits_,positive,__LINE__,__FILE__);
      //out << tmp.str() << endl;

      out << itable;
    }
    out << endl << "};" << endl;
    out.close();
  }

protected:
  vector<int> table_;
  int nbits_;
};

#endif
