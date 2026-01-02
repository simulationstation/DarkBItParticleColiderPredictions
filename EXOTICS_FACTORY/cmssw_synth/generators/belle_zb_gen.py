import FWCore.ParameterSet.Config as cms

from Configuration.Generator.Pythia8CommonSettings_cfi import pythia8CommonSettingsBlock
from Configuration.Generator.Pythia8CP5Settings_cfi import pythia8CP5SettingsBlock

# SYNTHETIC template: Zb(10610)/Zb(10650) proxy (not a physics model).

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
            # Define synthetic Zb states.
            "9900553:new = Zb1 Zb1 3 0 0 10.610 0.018 0.010 0.010",
            "9900554:new = Zb2 Zb2 3 0 0 10.650 0.018 0.010 0.010",
            # Force decays to Upsilon(nS) pi; hb PDGs included as placeholders.
            "9900553:oneChannel = 1 0.5 0 553 211",
            "9900553:oneChannel = 1 0.3 0 100553 211",
            "9900553:oneChannel = 1 0.2 0 200553 211",
            "9900553:oneChannel = 1 0.1 0 10513 211",
            "9900553:oneChannel = 1 0.1 0 20513 211",
            "9900554:oneChannel = 1 0.5 0 553 211",
            "9900554:oneChannel = 1 0.3 0 100553 211",
            "9900554:oneChannel = 1 0.2 0 200553 211",
            "9900554:oneChannel = 1 0.1 0 10513 211",
            "9900554:oneChannel = 1 0.1 0 20513 211",
        ),
        parameterSets=cms.vstring("pythia8CommonSettings", "pythia8CP5Settings", "processParameters"),
    ),
)
