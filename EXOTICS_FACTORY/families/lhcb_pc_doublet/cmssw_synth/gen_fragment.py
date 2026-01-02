import FWCore.ParameterSet.Config as cms

from Configuration.Generator.Pythia8CommonSettings_cfi import pythia8CommonSettingsBlock
from Configuration.Generator.Pythia8CP5Settings_cfi import pythia8CP5SettingsBlock

# SYNTHETIC template: Pc(4440)/Pc(4457) proxy (not a physics model).

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
            "9900443:new = Pc1 Pc1 3 0 0 4.440 0.020 0.010 0.010",
            "9900444:new = Pc2 Pc2 3 0 0 4.457 0.020 0.010 0.010",
            "9900443:oneChannel = 1 1.0 0 443 2212",
            "9900444:oneChannel = 1 1.0 0 443 2212",
        ),
        parameterSets=cms.vstring("pythia8CommonSettings", "pythia8CP5Settings", "processParameters"),
    ),
)
