import FWCore.ParameterSet.Config as cms

from Configuration.Generator.Pythia8CommonSettings_cfi import pythia8CommonSettingsBlock
from Configuration.Generator.Pythia8CP5Settings_cfi import pythia8CP5SettingsBlock

# SYNTHETIC template: omega(1420)/omega(1650) proxy for e+e- cross sections.

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
            "9901420:new = Omega1 Omega1 3 0 0 1.420 0.040 0.010 0.010",
            "9901650:new = Omega2 Omega2 3 0 0 1.650 0.050 0.010 0.010",
            "9901420:oneChannel = 1 0.5 0 223 9010221",
            "9901420:oneChannel = 1 0.5 0 223 211 -211",
            "9901650:oneChannel = 1 0.5 0 223 9010221",
            "9901650:oneChannel = 1 0.5 0 223 211 -211",
        ),
        parameterSets=cms.vstring("pythia8CommonSettings", "pythia8CP5Settings", "processParameters"),
    ),
)
