import FWCore.ParameterSet.Config as cms

_correctionFile2016Legacy    = "EgammaAnalysis/ElectronTools/data/ScalesSmearings/Legacy2016_07Aug2017_FineEtaR9_v3_ele_unc"
_correctionFile2017Nov17     = "EgammaAnalysis/ElectronTools/data/ScalesSmearings/Run2017_17Nov2017_v1_ele_unc"
_correctionFile2016ULpreVFP  = "EgammaAnalysis/ElectronTools/data/ScalesSmearings/Run2016_UltraLegacy_preVFP_RunFineEtaR9Gain"
_correctionFile2016ULpostVFP = "EgammaAnalysis/ElectronTools/data/ScalesSmearings/Run2016_UltraLegacy_postVFP_RunFineEtaR9Gain"
_correctionFile2017UL        = "EgammaAnalysis/ElectronTools/data/ScalesSmearings/Run2017_24Feb2020_runEtaR9Gain_v2"
_correctionFile2018UL        = "EgammaAnalysis/ElectronTools/data/ScalesSmearings/Run2018_29Sep2020_RunFineEtaR9Gain"

calibratedEgammaSettings = cms.PSet(minEtToCalibrate = cms.double(5.0),
                                    semiDeterministic = cms.bool(True),
                                    correctionFile = cms.string(_correctionFile2017Nov17),
                                    recHitCollectionEB = cms.InputTag('reducedEcalRecHitsEB'),
                                    recHitCollectionEE = cms.InputTag('reducedEcalRecHitsEE'),
                                    produceCalibratedObjs = cms.bool(True)
                                   )

from Configuration.Eras.Modifier_run2_miniAOD_devel_cff import run2_miniAOD_devel

from Configuration.Eras.Modifier_run2_egamma_2016_cff import run2_egamma_2016
from Configuration.Eras.Modifier_tracker_apv_vfp30_2016_cff import tracker_apv_vfp30_2016
(run2_miniAOD_devel & run2_egamma_2016 & tracker_apv_vfp30_2016).toModify(calibratedEgammaSettings,correctionFile = _correctionFile2016ULpreVFP)
(run2_miniAOD_devel & run2_egamma_2016 & ~tracker_apv_vfp30_2016).toModify(calibratedEgammaSettings,correctionFile = _correctionFile2016ULpostVFP)

from Configuration.ProcessModifiers.run2_miniAOD_UL_cff import run2_miniAOD_UL

from Configuration.Eras.Modifier_run2_egamma_2017_cff import run2_egamma_2017
(run2_miniAOD_UL & run2_egamma_2017).toModify(calibratedEgammaSettings,correctionFile = _correctionFile2017UL)

from Configuration.Eras.Modifier_run2_egamma_2018_cff import run2_egamma_2018
(run2_miniAOD_UL & run2_egamma_2018).toModify(calibratedEgammaSettings,correctionFile = _correctionFile2018UL)

from Configuration.Eras.Modifier_run2_miniAOD_80XLegacy_cff import run2_miniAOD_80XLegacy
run2_miniAOD_80XLegacy.toModify(calibratedEgammaSettings,correctionFile = _correctionFile2016Legacy)

calibratedEgammaPatSettings = calibratedEgammaSettings.clone(
    recHitCollectionEB = cms.InputTag('reducedEgamma','reducedEBRecHits'),
    recHitCollectionEE = cms.InputTag('reducedEgamma','reducedEERecHits')
    )

ecalTrkCombinationRegression = cms.PSet(
    ecalTrkRegressionConfig = cms.PSet(
        rangeMinLowEt = cms.double(-1.),
        rangeMaxLowEt = cms.double(3.0),
        rangeMinHighEt = cms.double(-1.),
        rangeMaxHighEt = cms.double(3.0),
        lowEtHighEtBoundary = cms.double(50.),
        forceHighEnergyTrainingIfSaturated = cms.bool(False),
        ebLowEtForestName = cms.string('electron_eb_ECALTRK_lowpt'),
        ebHighEtForestName = cms.string('electron_eb_ECALTRK'),
        eeLowEtForestName = cms.string('electron_ee_ECALTRK_lowpt'),
        eeHighEtForestName = cms.string('electron_ee_ECALTRK')
        ),
    ecalTrkRegressionUncertConfig = cms.PSet(
        rangeMinLowEt = cms.double(0.0002),
        rangeMaxLowEt = cms.double(0.5),
        rangeMinHighEt = cms.double(0.0002),
        rangeMaxHighEt = cms.double(0.5),
        lowEtHighEtBoundary = cms.double(50.),  
        forceHighEnergyTrainingIfSaturated = cms.bool(False),
        ebLowEtForestName = cms.string('electron_eb_ECALTRK_lowpt_var'),
        ebHighEtForestName = cms.string('electron_eb_ECALTRK_var'),
        eeLowEtForestName = cms.string('electron_ee_ECALTRK_lowpt_var'),
        eeHighEtForestName = cms.string('electron_ee_ECALTRK_var')
        ),
    maxEcalEnergyForComb=cms.double(200.),
    minEOverPForComb=cms.double(0.025),
    maxEPDiffInSigmaForComb=cms.double(15.),
    maxRelTrkMomErrForComb=cms.double(10.),
    
)

calibratedElectrons = cms.EDProducer("CalibratedElectronProducer",
                                     calibratedEgammaSettings,                                   
                                     epCombConfig = ecalTrkCombinationRegression,
                                     src = cms.InputTag('gedGsfElectrons'),
                                     )

calibratedPatElectrons = cms.EDProducer("CalibratedPatElectronProducer",
                                        calibratedEgammaPatSettings,
                                        epCombConfig = ecalTrkCombinationRegression,
                                        src = cms.InputTag('slimmedElectrons'), 
                                        )

calibratedPhotons = cms.EDProducer("CalibratedPhotonProducer",
                                   calibratedEgammaSettings,
                                   src = cms.InputTag('gedPhotons'),    
                                  )
calibratedPatPhotons = cms.EDProducer("CalibratedPatPhotonProducer",
                                      calibratedEgammaPatSettings,
                                      src = cms.InputTag('slimmedPhotons'),
                                      )

def prefixName(prefix,name):
    return prefix+name[0].upper()+name[1:]



