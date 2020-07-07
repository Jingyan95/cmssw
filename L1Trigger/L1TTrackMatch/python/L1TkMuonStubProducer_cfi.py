import FWCore.ParameterSet.Config as cms

L1TkMuonStub = cms.EDProducer("L1TkMuonStubProducer",
    ###############################################
    ############################################### common stuff    
    L1EMTFTrackCollectionInputTag = cms.InputTag("simEmtfDigis"),
    L1EMTFHitCollectionInputTag = cms.InputTag("simEmtfDigis"),
    L1TrackInputTag = cms.InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks"),
    ###############################################
    ############################################### TP algo
    emtfMatchAlgoVersion = cms.string( 'DynamicWindows' ), # version of matching with Trackes (string ID) :  DynamicWindows
    ###############################################
    ############################################### DynamicWindows algo
    mu_stub_station = cms.int32(12),

    emtfcorr_boundaries     = cms.FileInPath('L1Trigger/L1TMuon/data/emtf_luts/matching_windows_boundaries.root'),
    emtfcorr_S1_theta_windows  = cms.FileInPath('L1Trigger/L1TMuon/data/emtf_luts/matching_windows_station1_theta_q99.root'),
    emtfcorr_S1_phi_windows    = cms.FileInPath('L1Trigger/L1TMuon/data/emtf_luts/matching_windows_station1_phi_q99.root'),
    emtfcorr_theta_windows  = cms.FileInPath('L1Trigger/L1TMuon/data/emtf_luts/matching_windows_theta_q99.root'),
    emtfcorr_phi_windows    = cms.FileInPath('L1Trigger/L1TMuon/data/emtf_luts/matching_windows_phi_q99.root'),
    ## block to control the evolution of the matching window vs pt
    ## if do_relax_factors = False ==> global scale of upper and lower boundaries by "final_window_factor"
    ## if do_relax_factors = True  ==> progressive linear scale, the factor is
    ##      - initial_window_factor for pt <= pt_start_relax
    ##      - final_window_factor for pt >= pt_end_relax
    ##      - and a linear interpolation in the middle
    ## facror = 0 --> no changes to the window size
    initial_window_factor   = cms.double(0.0),
    final_window_factor     = cms.double(0.0),
    pt_start_relax          = cms.double(2.0),
    pt_end_relax            = cms.double(6.0),
    do_relax_factors        = cms.bool(False),
    ##
    n_trk_par      = cms.int32(4), # 4 or 5
    min_trk_p      = cms.double(3.5),
    max_trk_aeta   = cms.double(2.5),
    max_trk_chi2   = cms.double(100.0),
    min_trk_nstubs = cms.int32(4),
    require_BX0    = cms.bool(False),

)

# Station = 1, windows 90
L1TkMuonStubS1 = L1TkMuonStub.clone(

    mu_stub_station = cms.int32(1),

)

L1TkMuonStubS1_relaxed = L1TkMuonStub.clone(

    mu_stub_station = cms.int32(1),
    final_window_factor     = cms.double(0.5),

)

# Station = 2
L1TkMuonStubS2 = L1TkMuonStub.clone(

    mu_stub_station = cms.int32(2),

)

# Station = 12, 
L1TkMuonStubS12 = L1TkMuonStub.clone()

L1TkMuonStubS12_relaxed = L1TkMuonStub.clone(

    mu_stub_station = cms.int32(12),
    final_window_factor     = cms.double(0.5),
    do_relax_factors        = cms.bool(True),

)

L1TkMuonStubS12q90q99 = L1TkMuonStub.clone(

    mu_stub_station = cms.int32(12),
    emtfcorr_S1_theta_windows  = cms.FileInPath('L1Trigger/L1TMuon/data/emtf_luts/matching_windows_station1_theta_q90.root'),
    emtfcorr_S1_phi_windows    = cms.FileInPath('L1Trigger/L1TMuon/data/emtf_luts/matching_windows_station1_phi_q90.root'),
    emtfcorr_theta_windows  = cms.FileInPath('L1Trigger/L1TMuon/data/emtf_luts/matching_windows_theta_q99.root'),
    emtfcorr_phi_windows    = cms.FileInPath('L1Trigger/L1TMuon/data/emtf_luts/matching_windows_phi_q99.root'),

)

# Station = 1234, windows 90
L1TkMuonStubS1234 = L1TkMuonStub.clone(

    mu_stub_station = cms.int32(1234),

)
