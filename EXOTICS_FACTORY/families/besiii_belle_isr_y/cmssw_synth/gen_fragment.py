import FWCore.ParameterSet.Config as cms

from Configuration.Generator.Pythia8CommonSettings_cfi import pythia8CommonSettingsBlock
from Configuration.Generator.Pythia8CP5Settings_cfi import pythia8CP5SettingsBlock

# SYNTHETIC template: Y(4260)/Y(4360)/Y(4660) proxy for ISR cross sections.

generator = cms.EDFilter(
    "Pythia8GeneratorFilter",
    maxEventsToPrint=cms.untracked.int32(1),
    pythiaPylistVerbosity=cms.untracked.int32(1),
    filterEfficiency=cms.untracked.double(1.0),
    pythiaHepMCVerbosity=cms.untracked.bool(False),
    comEnergy=cms.double(13000.0),
    PythiaParameters=cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CP5SettingsBlock,
        processParameters=cms.vstring(
            "SoftQCD:nonDiffractive = on",
            "9900426:new = Y1 Y1 3 0 0 4.260 0.060 0.010 0.010",
            "9900436:new = Y2 Y2 3 0 0 4.360 0.060 0.010 0.010",
            "9900466:new = Y3 Y3 3 0 0 4.660 0.060 0.010 0.010",
            "9900426:oneChannel = 1 0.6 0 443 211 -211",
            "9900426:oneChannel = 1 0.4 0 100443 211 -211",
            "9900436:oneChannel = 1 0.6 0 443 211 -211",
            "9900436:oneChannel = 1 0.4 0 100443 211 -211",
            "9900466:oneChannel = 1 0.6 0 443 211 -211",
            "9900466:oneChannel = 1 0.4 0 100443 211 -211",
        ),
        parameterSets=cms.vstring("pythia8CommonSettings", "pythia8CP5Settings", "processParameters"),
    ),
)
