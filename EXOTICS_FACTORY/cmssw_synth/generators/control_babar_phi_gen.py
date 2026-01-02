import FWCore.ParameterSet.Config as cms

from Configuration.Generator.Pythia8CommonSettings_cfi import pythia8CommonSettingsBlock
from Configuration.Generator.Pythia8CP5Settings_cfi import pythia8CP5SettingsBlock

# SYNTHETIC template: phi(1680)/phi(2170) proxy for e+e- cross sections.

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
            "9902168:new = Phi1 Phi1 3 0 0 1.680 0.050 0.010 0.010",
            "9902217:new = Phi2 Phi2 3 0 0 2.170 0.070 0.010 0.010",
            "9902168:oneChannel = 1 0.5 0 333 9010221",
            "9902168:oneChannel = 1 0.5 0 333 211 -211",
            "9902217:oneChannel = 1 0.5 0 333 9010221",
            "9902217:oneChannel = 1 0.5 0 333 211 -211",
        ),
        parameterSets=cms.vstring("pythia8CommonSettings", "pythia8CP5Settings", "processParameters"),
    ),
)
