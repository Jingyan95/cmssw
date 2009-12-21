// -*- C++ -*
/* 
Functions to return total energy contained in NxN (3x3/5x5/7x7)
Hcal towers aroud a given DetId. 

Inputs : 
1. HcalTopology, 
2. DetId around whic NxN is to be formed, 
3. HcalRecHitCollection,
4. Number of towers to be navigated along eta and phi along 
   one direction (navigation is done alone +-deta and +-dphi).
5. option to include HO

Authors:  Seema Sharma, Sunanda Banerjee
Created: August 2009
*/

#ifndef CalibrationIsolatedParticleseHCALMatrix_h
#define CalibrationIsolatedParticleseHCALMatrix_h

// system include files
#include <memory>
#include <map>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/DetId/interface/DetId.h"

#include "Geometry/CaloTopology/interface/HcalTopology.h"

namespace spr{

  template< typename T>
  double eHCALmatrix(const HcalTopology* topology, const DetId& det, edm::Handle<T>& hits, int ieta, int iphi, bool includeHO=false, bool algoNew=true, bool debug=false);

  template< typename T>
  std::vector< std::pair< DetId,double> > eHCALmatrix(const HcalTopology* topology, const DetId& det, edm::Handle<T>& hits, int ieta, int iphi, bool includeHO=false, bool debug=false);

  template< typename T>
  double eHCALmatrix(const HcalTopology* topology, const DetId& det, edm::Handle<T>& hits, int ietaE, int ietaW, int iphiN, int iphiS, bool includeHO=false, bool debug=false);

  template< typename T>
  std::pair<double,int> eHCALmatrixTotal(const HcalTopology* topology, const DetId& det, edm::Handle<T>& hits, int ieta, int iphi, bool includeHO=false, bool debug=false);

  template <typename T>
  std::vector<typename T::const_iterator> hitHCALmatrix(const HcalTopology* topology, const DetId& det, edm::Handle<T>& hits, int ieta, int iphi, bool includeHO=false, bool debug=false);

  template <typename T>
  std::vector<typename T::const_iterator> hitHCALmatrixNew(const HcalTopology* topology, const DetId& det, edm::Handle<T>& hits, int ieta, int iphi, bool includeHO=false, bool debug=false);

  template <typename T>
  std::vector<typename T::const_iterator> hitHCALmatrixTotal(const HcalTopology* topology, const DetId& det, edm::Handle<T>& hits, int ietaE, int ietaW, int iphiN, int iphiS, bool includeHO=false, bool debug=false);

  template <typename T>
  std::vector<typename T::const_iterator> hitHCALneighbours(std::vector<DetId>& vNeighboursDetId, std::vector<DetId>& dets, const HcalTopology* topology, edm::Handle<T>& hits, bool includeHO=false, bool debug=false);

  template <typename T>
  std::vector<typename T::const_iterator> hitHCALneighbours(std::vector<DetId>& vdets, const HcalTopology* topology, edm::Handle<T>& hits, bool includeHO=false, bool debug=false);

  template <typename T>
  std::vector<std::pair<DetId,double> > energyDetIdHCALneighbours(std::vector<DetId>& vdets, const HcalTopology* topology, edm::Handle<T>& hits, bool includeHO=false, bool debug=false);
}

#include "Calibration/IsolatedParticles/interface/eHCALMatrix.icc"
#endif
