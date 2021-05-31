import FWCore.ParameterSet.Config as cms

HitPatternHelper_params = cms.PSet (

  HPHdebug   = cms.bool(True),   
  useNewKF   = cms.bool(False),
  chosenRofZ = cms.double(50.0)

)
