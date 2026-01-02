import FWCore.ParameterSet.Config as cms

from Configuration.Generator.Pythia8CommonSettings_cfi import pythia8CommonSettingsBlock
from Configuration.Generator.Pythia8CP5Settings_cfi import pythia8CP5SettingsBlock

# SYNTHETIC template: Pcs(4459) proxy (not a physics model).

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
            "9904459:new = Pcs1 Pcs1 3 0 0 4.459 0.020 0.010 0.010",
            "9904458:new = Pcs2 Pcs2 3 0 0 4.443 0.020 0.010 0.010",
            "9904459:oneChannel = 1 1.0 0 443 3122",
            "9904458:oneChannel = 1 1.0 0 443 3122",
        ),
        parameterSets=cms.vstring("pythia8CommonSettings", "pythia8CP5Settings", "processParameters"),
    ),
)
