import FWCore.ParameterSet.Config as cms

HitPatternHelper_params = cms.PSet (

  HPHdebug   = cms.bool(True),                     # checks if input sample production is configured as current process
  chosenRofZ = cms.double(50.0)                      # enable emulation of truncation, lost stubs are filled in BranchLost

)
