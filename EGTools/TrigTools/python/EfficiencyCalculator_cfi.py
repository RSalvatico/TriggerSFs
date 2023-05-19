import FWCore.ParameterSet.Config as cms

EfficiencyCalculator = cms.EDAnalyzer('EfficiencyCalculator',
                                      stageL1Trigger = cms.uint32(2),
                                      runningOnData      = cms.bool(False),
                                      runningEra         = cms.int32(0), # One of the possible python types for the C++ type "int" 
)
