import FWCore.ParameterSet.Config as cms

L1TrackerEtMiss = cms.EDProducer('L1TrackerEtMissProducer',
    # L1TrackInputTag = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),
     L1TrackInputTag = cms.InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks"),
    # L1VertexInputTag = cms.InputTag("VertexProducer", "l1vertextdr"),
     L1VertexInputTag = cms.InputTag("L1TkPrimaryVertex", "l1vertextdr"),
     L1Tk_maxZ0 = cms.double ( 15. ) ,    # in cm
     L1Tk_maxEta = cms.double ( 2.4 ) ,
     L1Tk_maxChi2dof = cms.double( 50. ),
     L1Tk_maxBendchi2 = cms.double( 1.75 ),
     L1Tk_minPt = cms.double( 2. ),       # in GeV
     L1Tk_maxDeltaZ = cms.double( 3. ),      # in cm
     L1Tk_minNStubs = cms.int32( 4 ),     # min number of stubs for the tracks to enter in TrkMET calculation
     L1Tk_minNStubsPS = cms.int32( 2 ),   # min number of stubs in the PS Modules
     L1Tk_maxPt = cms.double( 200. ),	 # in GeV. When maxPt > 0, tracks with PT above maxPt are considered as
                                     # mismeasured and are treated according to HighPtTracks below.
                                     # When maxPt < 0, no special treatment is done for high PT tracks.
     L1Tk_HighPtTracks = cms.int32( 1 ),  # when = 0 : truncation. Tracks with PT above maxPt are ignored
                                     # when = 1 : saturation. Tracks with PT above maxPt are set to PT=maxPt.
                                     # When maxPt < 0, no special treatment is done for high PT tracks.
)

