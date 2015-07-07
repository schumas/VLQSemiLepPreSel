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



