VLQSemiLepPresel
================

Vector-Like-Quark analysis code. Preselection for various analyses. Common code specific to VLQ analyses.


The latest pre-selected events reside on dust:
/nfs/dust/cms/user/tholenhe/VLQSemiLepPreSel/PHYS14-ntuple2-v2


Selections
==========

The selection requirements are defined in ``include/VLQSLPS_selectionItems.h``. 
Git tags are used for versioning:

v1 
--

- leading ak4Jet pT > 200 GeV 
- primary lepton pT > 50 GeV
- ST > 500 GeV


v2
--

- leading ak4Jet pT > 100 GeV 
- primary lepton pT > 50 GeV
- ST > 400 GeV
- number of loose csvv2 btags (ak4Jets) >= 1
- 2D-Cut (dR(l, j) > 0.2 OR dpt(l, j) > 10.)


v3
--

- IDs: 
  - ElectronID_Spring15_50ns_medium_noIso, PtEtaCut(20.0, 2.4)
  - MuonIDTight, PtEtaCut(20.0, 2.4)
  - JetPFID(JetPFID::WP_LOOSE), PtEtaCut(30.0,3.6)

- leading ak4Jet pT > 100 GeV 
- primary lepton pT > 25 GeV
- ST > 400 GeV
- number of loose csvv2 btags (ak4Jets) >= 1
- 2D-Cut (dR(l, j) > 0.25 OR dpt(l, j) > 40.)
- any of these triggers fired (OR combination):
  - "HLT_Ele105_CaloIdVT_GsfTrkIdT_v*"
  - "HLT_Ele32_eta2p1_WP75_Gsf_v*"
  - "HLT_Mu45_eta2p1_v*"
  - "HLT_Mu50_v*"
  - "HLT_IsoMu24_eta2p1_v*"
  - "HLT_IsoMu27_v*"
  - "HLT_PFHT800_v*"


v4, v5
------

- IDs: 
  - ElectronID_Spring15_50ns_medium_noIso, PtEtaCut(20.0, 2.4)
  - MuonIDTight, PtEtaCut(20.0, 2.4)
  - JetPFID(JetPFID::WP_LOOSE), PtEtaCut(30.0,3.6)

- leading ak4Jet pT > 100 GeV 
- primary lepton pT > 50 GeV (Muon), > 115 GeV (Electron)
- ST > 400 GeV
- 2D-Cut (dR(l, j) > 0.4 OR dpt(l, j) > 40.)
- any of these triggers fired (OR combination):
  - "HLT_Ele105_CaloIdVT_GsfTrkIdT_v*",
  - "HLT_Ele32_eta2p1_WPLoose_Gsf_v*",
  - "HLT_Mu45_eta2p1_v*",
  - "HLT_IsoMu24_eta2p1_v*",
  - "HLT_PFHT800_v*",


v6 (for 25ns runs)
------------------

- IDs: 
  - ElectronID_Spring15_25ns_medium_noIso, PtEtaCut(20.0, 2.4)
  - MuonIDMedium, PtEtaCut(20.0, 2.4)
  - JetPFID(JetPFID::WP_LOOSE), PtEtaCut(30.0,3.6)

- leading ak4Jet pT > 100 GeV 
- primary lepton pT > 50 GeV (Muon), > 115 GeV (Electron)
- ST > 400 GeV
- 2D-Cut (dR(l, j) > 0.4 OR dpt(l, j) > 40.)
- number of loose csvv2 btags (ak4Jets) >= 1
- any of these triggers fired (OR combination):
  - "HLT_Ele105_CaloIdVT_GsfTrkIdT_v*",
  - "HLT_Mu45_eta2p1_v*",
  - "HLT_PFHT800_v*",

