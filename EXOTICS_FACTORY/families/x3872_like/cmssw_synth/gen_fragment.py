import FWCore.ParameterSet.Config as cms

from Configuration.Generator.Pythia8CommonSettings_cfi import pythia8CommonSettingsBlock
from Configuration.Generator.Pythia8CP5Settings_cfi import pythia8CP5SettingsBlock

# SYNTHETIC template: X(3872)-like proxy (not a physics model).

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
            "9903872:new = X1 X1 3 0 0 3.872 0.005 0.001 0.001",
            "9903873:new = X2 X2 3 0 0 3.910 0.010 0.001 0.001",
            "9903872:oneChannel = 1 0.6 0 443 211 -211",
            "9903872:oneChannel = 1 0.4 0 443 211 -211 111",
            "9903873:oneChannel = 1 0.6 0 443 211 -211",
            "9903873:oneChannel = 1 0.4 0 443 211 -211 111",
        ),
        parameterSets=cms.vstring("pythia8CommonSettings", "pythia8CP5Settings", "processParameters"),
    ),
)
