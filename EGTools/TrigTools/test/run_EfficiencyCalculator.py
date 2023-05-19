import FWCore.ParameterSet.Config as cms
import glob
import os

process = cms.Process("USER")

#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimL1EmulatorRepack_Full_cff')
process.load('HLTrigger.Configuration.HLT_GRun_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v37') 

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# dirName = "/eos/cms/store/group/phys_egamma/ec/rsalvatico/HCALIsoTuning/126X_afterECALFix/"
# #dirName = "/afs/cern.ch/user/r/rselvati/work/private/Effi/HCALIsoTuning/CMSSW_13_0_0/src/EGTools/TrigTools/test/Ele32_Nominal/"
# fileList = filter(os.path.isfile, glob.glob(dirName + "*.root"))
# fList = []
# for f in fileList:
#     fs = str(f).replace("/eos/","file:/eos/")
#     #fs = str(f).replace("/afs/","file:/afs/")
#     fList.append(fs)
# print(fList)
#from rootfile_list import theList
process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring("root://cms-xrd-global.cern.ch//store/data/Run2018D/EGamma/MINIAOD/UL2018_MiniAODv2-v2/120000/003D1380-15BF-534F-9DF4-3E7D7743D4F4.root")
                            )

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing()
options.register('runningOnData',
                 True, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Data/MC config flag")

options.register('runningEra',
                 0, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Year configuration")

#Output file
process.TFileService = cms.Service("TFileService",
   fileName = cms.string("EfficiencyTest.root")
)

process.load("EGTools.TrigTools.EfficiencyCalculator_cfi")
process.EfficiencyCalculator.runningOnData = options.runningOnData
process.EfficiencyCalculator.runningEra    = options.runningEra


#process.EfficiencyCalculator = cms.EDAnalyzer('EfficiencyCalculator',
#                                              stageL1Trigger = cms.uint32(2)
#)

#from EgammaUser.EgammaPostRecoTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
#setupEgammaPostRecoSeq(process,era='2018-UL')  

process.p = cms.Path(process.EfficiencyCalculator)
process.schedule = cms.Schedule(process.p)
