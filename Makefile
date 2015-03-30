LIBRARY := SUHH2VLQSemiLepPreSel
# enable par creation; this is necessary for all packages containing AnalysisModules
# to be loaded from by AnalysisModuleRunner.
PAR := 1
USERLDFLAGS := -lSUHH2core -lSUHH2common -lGenVector
include ../Makefile.common
