import FWCore.ParameterSet.Config as cms

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask, addToProcessAndTask

def miniAOD_customizeCommon(process):
    process.patMuons.isoDeposits = cms.PSet()
    process.patElectrons.isoDeposits = cms.PSet()
    process.patTaus.isoDeposits = cms.PSet()
    process.patPhotons.isoDeposits = cms.PSet()
    #
    process.patMuons.embedTrack         = True  # used for IDs
    process.patMuons.embedCombinedMuon  = True  # used for IDs
    process.patMuons.embedMuonBestTrack = True  # used for IDs
    process.patMuons.embedStandAloneMuon = True # maybe?
    process.patMuons.embedPickyMuon = False   # no, use best track
    process.patMuons.embedTpfmsMuon = False   # no, use best track
    process.patMuons.embedDytMuon   = False   # no, use best track
    process.patMuons.addPuppiIsolation = cms.bool(True)
    process.patMuons.puppiIsolationChargedHadrons = cms.InputTag("muonPUPPIIsolation","h+-DR040-ThresholdVeto000-ConeVeto000")
    process.patMuons.puppiIsolationNeutralHadrons = cms.InputTag("muonPUPPIIsolation","h0-DR040-ThresholdVeto000-ConeVeto001")
    process.patMuons.puppiIsolationPhotons        = cms.InputTag("muonPUPPIIsolation","gamma-DR040-ThresholdVeto000-ConeVeto001")
    process.patMuons.puppiNoLeptonsIsolationChargedHadrons = cms.InputTag("muonPUPPINoLeptonsIsolation","h+-DR040-ThresholdVeto000-ConeVeto000")
    process.patMuons.puppiNoLeptonsIsolationNeutralHadrons = cms.InputTag("muonPUPPINoLeptonsIsolation","h0-DR040-ThresholdVeto000-ConeVeto001")
    process.patMuons.puppiNoLeptonsIsolationPhotons        = cms.InputTag("muonPUPPINoLeptonsIsolation","gamma-DR040-ThresholdVeto000-ConeVeto001")

    process.patMuons.computeMiniIso = True
    process.patMuons.computeMuonMVA = True
    process.patMuons.computeSoftMuonMVA = True

    process.patMuons.addTriggerMatching = True
    from Configuration.ProcessModifiers.run2_miniAOD_UL_cff import run2_miniAOD_UL
    from Configuration.Eras.Modifier_run2_muon_2016_cff import run2_muon_2016
    from Configuration.Eras.Modifier_run2_muon_2017_cff import run2_muon_2017
    from Configuration.Eras.Modifier_run2_muon_2018_cff import run2_muon_2018
    run2_muon_2016.toModify( process.patMuons, effectiveAreaVec = [0.0735,0.0619,0.0465,0.0433,0.0577])
    run2_muon_2017.toModify( process.patMuons, effectiveAreaVec = [0.0566, 0.0562, 0.0363, 0.0119, 0.0064])
    run2_muon_2018.toModify( process.patMuons, effectiveAreaVec = [0.0566, 0.0562, 0.0363, 0.0119, 0.0064])
    run2_muon_2016.toModify( process.patMuons, mvaTrainingFile = "RecoMuon/MuonIdentification/data/mu_2016_BDTG.weights.xml")
    run2_miniAOD_UL.toModify( process.patMuons, getdBFromTrack = True)
    process.patMuons.computePuppiCombinedIso = True
    #
    # disable embedding of electron and photon associated objects already stored by the ReducedEGProducer
    process.patElectrons.embedGsfElectronCore = False  ## process.patElectrons.embed in AOD externally stored gsf electron core
    process.patElectrons.embedSuperCluster    = False  ## process.patElectrons.embed in AOD externally stored supercluster
    process.patElectrons.embedPflowSuperCluster         = False  ## process.patElectrons.embed in AOD externally stored supercluster
    process.patElectrons.embedSeedCluster               = False  ## process.patElectrons.embed in AOD externally stored the electron's seedcluster
    process.patElectrons.embedBasicClusters             = False  ## process.patElectrons.embed in AOD externally stored the electron's basic clusters
    process.patElectrons.embedPreshowerClusters         = False  ## process.patElectrons.embed in AOD externally stored the electron's preshower clusters
    process.patElectrons.embedPflowBasicClusters        = False  ## process.patElectrons.embed in AOD externally stored the electron's pflow basic clusters
    process.patElectrons.embedPflowPreshowerClusters    = False  ## process.patElectrons.embed in AOD externally stored the electron's pflow preshower clusters
    process.patElectrons.embedRecHits         = False  ## process.patElectrons.embed in AOD externally stored the RecHits - can be called from the PATElectronProducer
    process.patElectrons.electronSource = cms.InputTag("reducedEgamma","reducedGedGsfElectrons")
    process.patElectrons.usePfCandidateMultiMap = True
    process.patElectrons.pfCandidateMultiMap    = cms.InputTag("reducedEgamma","reducedGsfElectronPfCandMap")
    process.patElectrons.electronIDSources = cms.PSet()
    run2_miniAOD_UL.toModify( process.patElectrons, getdBFromTrack = True)
    from Configuration.Eras.Modifier_run2_miniAOD_80XLegacy_cff import run2_miniAOD_80XLegacy
    run2_miniAOD_80XLegacy.toModify(process.patElectrons,
                                    addPFClusterIso = cms.bool(True),
                                    ecalPFClusterIsoMap = cms.InputTag("reducedEgamma", "eleEcalPFClusIso"),
                                    hcalPFClusterIsoMap = cms.InputTag("reducedEgamma", "eleHcalPFClusIso"))
    from Configuration.Eras.Modifier_run2_miniAOD_94XFall17_cff import run2_miniAOD_94XFall17
    run2_miniAOD_94XFall17.toModify(process.patElectrons,
                                    addPFClusterIso = cms.bool(True),
                                    ecalPFClusterIsoMap = cms.InputTag("reducedEgamma", "eleEcalPFClusIso"),
                                    hcalPFClusterIsoMap = cms.InputTag("reducedEgamma", "eleHcalPFClusIso"))

    #add puppi isolation in miniAOD
    process.patElectrons.addPuppiIsolation = cms.bool(True)
    process.patElectrons.puppiIsolationChargedHadrons = cms.InputTag("egmElectronPUPPIIsolation","h+-DR030-BarVeto000-EndVeto001")
    process.patElectrons.puppiIsolationNeutralHadrons = cms.InputTag("egmElectronPUPPIIsolation","h0-DR030-BarVeto000-EndVeto000")
    process.patElectrons.puppiIsolationPhotons        = cms.InputTag("egmElectronPUPPIIsolation","gamma-DR030-BarVeto000-EndVeto008")
    process.patElectrons.puppiNoLeptonsIsolationChargedHadrons = cms.InputTag("egmElectronPUPPINoLeptonsIsolation","h+-DR030-BarVeto000-EndVeto001")
    process.patElectrons.puppiNoLeptonsIsolationNeutralHadrons = cms.InputTag("egmElectronPUPPINoLeptonsIsolation","h0-DR030-BarVeto000-EndVeto000")
    process.patElectrons.puppiNoLeptonsIsolationPhotons        = cms.InputTag("egmElectronPUPPINoLeptonsIsolation","gamma-DR030-BarVeto000-EndVeto008")

    process.patElectrons.computeMiniIso = cms.bool(True)

    process.elPFIsoDepositChargedPAT.src = cms.InputTag("reducedEgamma","reducedGedGsfElectrons")
    process.elPFIsoDepositChargedAllPAT.src = cms.InputTag("reducedEgamma","reducedGedGsfElectrons")
    process.elPFIsoDepositNeutralPAT.src = cms.InputTag("reducedEgamma","reducedGedGsfElectrons")
    process.elPFIsoDepositGammaPAT.src = cms.InputTag("reducedEgamma","reducedGedGsfElectrons")
    process.elPFIsoDepositPUPAT.src = cms.InputTag("reducedEgamma","reducedGedGsfElectrons")
    #
    process.patPhotons.embedSuperCluster = False ## whether to process.patPhotons.embed in AOD externally stored supercluster
    process.patPhotons.embedSeedCluster               = False  ## process.patPhotons.embed in AOD externally stored the photon's seedcluster
    process.patPhotons.embedBasicClusters             = False  ## process.patPhotons.embed in AOD externally stored the photon's basic clusters
    process.patPhotons.embedPreshowerClusters         = False  ## process.patPhotons.embed in AOD externally stored the photon's preshower clusters
    process.patPhotons.embedRecHits         = False  ## process.patPhotons.embed in AOD externally stored the RecHits - can be called from the PATPhotonProducer

    #add puppi isolation in miniAOD
    process.patPhotons.addPuppiIsolation = cms.bool(True)
    process.patPhotons.puppiIsolationChargedHadrons = cms.InputTag("egmPhotonPUPPIIsolation","h+-DR030-")
    process.patPhotons.puppiIsolationNeutralHadrons = cms.InputTag("egmPhotonPUPPIIsolation","h0-DR030-")
    process.patPhotons.puppiIsolationPhotons        = cms.InputTag("egmPhotonPUPPIIsolation","gamma-DR030-")

    from Configuration.Eras.Modifier_run2_miniAOD_80XLegacy_cff import run2_miniAOD_80XLegacy
    run2_miniAOD_80XLegacy.toModify(process.patPhotons,
                                    addPFClusterIso = cms.bool(True),
                                    ecalPFClusterIsoMap = cms.InputTag("reducedEgamma", "phoEcalPFClusIso"),
                                    hcalPFClusterIsoMap = cms.InputTag("reducedEgamma", "phoHcalPFClusIso"))
    from Configuration.Eras.Modifier_run2_miniAOD_94XFall17_cff import run2_miniAOD_94XFall17
    run2_miniAOD_94XFall17.toModify(process.patPhotons,
                                    addPFClusterIso = cms.bool(True),
                                    ecalPFClusterIsoMap = cms.InputTag("reducedEgamma", "phoEcalPFClusIso"),
                                    hcalPFClusterIsoMap = cms.InputTag("reducedEgamma", "phoHcalPFClusIso"))
    #the 80X legacy customsations are done in ootPhotonProducer for OOT photons
    run2_miniAOD_94XFall17.toModify(process.patOOTPhotons,
                                    addPFClusterIso = cms.bool(True),
                                    ecalPFClusterIsoMap = cms.InputTag("reducedEgamma", "ootPhoEcalPFClusIso"),
                                    hcalPFClusterIsoMap = cms.InputTag("reducedEgamma", "ootPhoHcalPFClusIso"))


    process.patPhotons.photonSource = cms.InputTag("reducedEgamma","reducedGedPhotons")
    process.patPhotons.electronSource = cms.InputTag("reducedEgamma","reducedGedGsfElectrons")

    process.phPFIsoDepositChargedPAT.src = cms.InputTag("reducedEgamma","reducedGedPhotons")
    process.phPFIsoDepositChargedAllPAT.src = cms.InputTag("reducedEgamma","reducedGedPhotons")
    process.phPFIsoDepositNeutralPAT.src = cms.InputTag("reducedEgamma","reducedGedPhotons")
    process.phPFIsoDepositGammaPAT.src = cms.InputTag("reducedEgamma","reducedGedPhotons")
    process.phPFIsoDepositPUPAT.src = cms.InputTag("reducedEgamma","reducedGedPhotons")
    #
    process.patOOTPhotons.photonSource = cms.InputTag("reducedEgamma","reducedOOTPhotons")
    process.patOOTPhotons.electronSource = cms.InputTag("reducedEgamma","reducedGedGsfElectrons")
    #
    process.selectedPatJets.cut = cms.string("pt > 10")
    process.selectedPatMuons.cut = cms.string("pt > 5 || isPFMuon || (pt > 3 && (isGlobalMuon || isStandAloneMuon || numberOfMatches > 0 || muonID('RPCMuLoose')))")

    from Configuration.Eras.Modifier_phase2_muon_cff import phase2_muon
    phase2_muon.toModify(process.selectedPatMuons, cut = "pt > 5 || isPFMuon || (pt > 3 && (isGlobalMuon || isStandAloneMuon || numberOfMatches > 0 || muonID('RPCMuLoose') || muonID('ME0MuonArbitrated') || muonID('GEMMuonArbitrated')) )")

    process.selectedPatElectrons.cut = cms.string("")
    process.selectedPatTaus.cut = cms.string("pt > 18. && tauID('decayModeFindingNewDMs')> 0.5")
    process.selectedPatPhotons.cut = cms.string("")

    from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection

    from PhysicsTools.PatAlgos.slimming.applySubstructure_cff import applySubstructure
    applySubstructure( process )
    
    #
    from PhysicsTools.PatAlgos.tools.trigTools import switchOnTriggerStandAlone
    switchOnTriggerStandAlone( process, outputModule = '' )
    process.patTrigger.packTriggerPathNames = cms.bool(True)
    #
    # apply type I + other PFMEt corrections to pat::MET object
    # and estimate systematic uncertainties on MET

    from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncForMiniAODProduction
    runMetCorAndUncForMiniAODProduction(process, metType="PF",
                                        jetCollUnskimmed="patJets")

    #caloMET computation
    from PhysicsTools.PatAlgos.tools.metTools import addMETCollection
    addMETCollection(process,
                     labelName = "patCaloMet",
                     metSource = "caloMetM"
                     )

    #noHF pfMET =========

    task = getPatAlgosToolsTask(process)

    process.noHFCands = cms.EDFilter("GenericPFCandidateSelector",
                                     src=cms.InputTag("particleFlow"),
                                     cut=cms.string("abs(pdgId)!=1 && abs(pdgId)!=2 && abs(eta)<3.0")
                                     )
    task.add(process.noHFCands)

    runMetCorAndUncForMiniAODProduction(process,
                                        pfCandColl=cms.InputTag("noHFCands"),
                                        recoMetFromPFCs=True, #needed for HF removal
                                        jetSelection="pt>15 && abs(eta)<3.",
                                        postfix="NoHF"
                                        )

    process.load('PhysicsTools.PatAlgos.slimming.slimmedMETs_cfi')
    task.add(process.slimmedMETs)

    addToProcessAndTask('slimmedMETsNoHF', process.slimmedMETs.clone(), process, task)
    process.slimmedMETsNoHF.src = cms.InputTag("patMETsNoHF")
    process.slimmedMETsNoHF.rawVariation =  cms.InputTag("patPFMetNoHF")
    process.slimmedMETsNoHF.t1Uncertainties = cms.InputTag("patPFMetT1%sNoHF")
    process.slimmedMETsNoHF.t01Variation = cms.InputTag("patPFMetT0pcT1NoHF")
    process.slimmedMETsNoHF.t1SmearedVarsAndUncs = cms.InputTag("patPFMetT1Smear%sNoHF")
    process.slimmedMETsNoHF.tXYUncForRaw = cms.InputTag("patPFMetTxyNoHF")
    process.slimmedMETsNoHF.tXYUncForT1 = cms.InputTag("patPFMetT1TxyNoHF")
    process.slimmedMETsNoHF.tXYUncForT01 = cms.InputTag("patPFMetT0pcT1TxyNoHF")
    process.slimmedMETsNoHF.tXYUncForT1Smear = cms.InputTag("patPFMetT1SmearTxyNoHF")
    process.slimmedMETsNoHF.tXYUncForT01Smear = cms.InputTag("patPFMetT0pcT1SmearTxyNoHF")
    del process.slimmedMETsNoHF.caloMET
    # ================== NoHF pfMET

    #  ==================  CHSMET
    process.CHSCands = cms.EDFilter("CandPtrSelector",
                                    src=cms.InputTag("packedPFCandidates"),
                                    cut=cms.string("fromPV(0) > 0")
                                    )
    task.add(process.CHSCands)

    process.pfMetCHS = cms.EDProducer("PFMETProducer",
                                      src = cms.InputTag("CHSCands"),
                                      alias = cms.string('pfMet'),
                                      globalThreshold = cms.double(0.0),
                                      calculateSignificance = cms.bool(False),
                                      )
    task.add(process.pfMetCHS)

    addMETCollection(process,
                     labelName = "patCHSMet",
                     metSource = "pfMetCHS"
                     )

    process.patCHSMet.computeMETSignificance = cms.bool(False)

    #  ==================  CHSMET

    #  ==================  TrkMET
    process.TrkCands = cms.EDFilter("CandPtrSelector",
                                    src=cms.InputTag("packedPFCandidates"),
                                    cut=cms.string("charge()!=0 && pvAssociationQuality()>=4 && vertexRef().key()==0")
                                    )
    task.add(process.TrkCands)

    process.pfMetTrk = cms.EDProducer("PFMETProducer",
                                      src = cms.InputTag("TrkCands"),
                                      alias = cms.string('pfMet'),
                                      globalThreshold = cms.double(0.0),
                                      calculateSignificance = cms.bool(False),
                                      )

    task.add(process.pfMetTrk)

    addMETCollection(process,
                     labelName = "patTrkMet",
                     metSource = "pfMetTrk"
                     )

    process.patTrkMet.computeMETSignificance = cms.bool(False)

    #  ==================  TrkMET


    ## PU JetID
    process.load("RecoJets.JetProducers.PileupJetID_cfi")
    task.add(process.pileUpJetIDTask)

    process.patJets.userData.userFloats.src = [ cms.InputTag("pileupJetId:fullDiscriminant"), ]
    process.patJets.userData.userInts.src = [ cms.InputTag("pileupJetId:fullId"), ]

    ## Quark Gluon Likelihood
    process.load('RecoJets.JetProducers.QGTagger_cfi')
    task.add(process.QGTaggerTask)

    process.patJets.userData.userFloats.src += [ 'QGTagger:qgLikelihood', ]

    #HF jet shower shape
    process.load('RecoJets.JetProducers.hfJetShowerShape_cfi')
    task.add(process.hfJetShowerShape)

    run2_miniAOD_UL.toModify(process.patJets.userData.userFloats,
                             src = process.patJets.userData.userFloats.src + ['hfJetShowerShape:sigmaEtaEta', 'hfJetShowerShape:sigmaPhiPhi'])
    run2_miniAOD_UL.toModify(process.patJets.userData.userInts,
                             src = process.patJets.userData.userInts.src + ['hfJetShowerShape:centralEtaStripSize', 'hfJetShowerShape:adjacentEtaStripsSize'])

    ## DeepCSV meta discriminators (simple arithmethic on output probabilities)
    process.load('RecoBTag.Combined.deepFlavour_cff')
    task.add(process.pfDeepCSVDiscriminatorsJetTags)
    process.patJets.discriminatorSources.extend([
            cms.InputTag('pfDeepCSVDiscriminatorsJetTags:BvsAll' ),
            cms.InputTag('pfDeepCSVDiscriminatorsJetTags:CvsB'   ),
            cms.InputTag('pfDeepCSVDiscriminatorsJetTags:CvsL'   ),
            ])

    ## CaloJets
    process.caloJetMap = cms.EDProducer("RecoJetDeltaRValueMapProducer",
         src = process.patJets.jetSource,
         matched = cms.InputTag("ak4CaloJets"),
         distMax = cms.double(0.4),
         values = cms.vstring('pt','emEnergyFraction'),
	 valueLabels = cms.vstring('pt','emEnergyFraction'),
	 lazyParser = cms.bool(True) )
    task.add(process.caloJetMap)
    process.patJets.userData.userFloats.src += [ 'caloJetMap:pt', 'caloJetMap:emEnergyFraction' ]

    #Muon object modifications
    from PhysicsTools.PatAlgos.slimming.muonIsolationsPUPPI_cfi import makeInputForPUPPIIsolationMuon
    makeInputForPUPPIIsolationMuon(process)

    #EGM object modifications
    from PhysicsTools.PatAlgos.slimming.egmIsolationsPUPPI_cfi import makeInputForPUPPIIsolationEgm
    makeInputForPUPPIIsolationEgm(process)
    from RecoEgamma.EgammaTools.egammaObjectModificationsInMiniAOD_cff import egamma_modifications
    process.slimmedElectrons.modifierConfig.modifications = egamma_modifications
    process.slimmedPhotons.modifierConfig.modifications   = egamma_modifications

    #VID Electron IDs
    process.patElectrons.addElectronID = cms.bool(True)
    electron_ids = ['RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff',
                    'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV71_cff',
                    'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V1_cff',
                    'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V2_cff',
                    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V1_cff',
                    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V1_cff',
                    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V2_cff',
                    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V2_cff',
                    'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff',
                    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff',
                    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_HZZ_V1_cff',
                    ]
    switchOnVIDElectronIdProducer(process,DataFormat.MiniAOD, task)
    process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag("reducedEgamma","reducedGedGsfElectrons")
    process.electronMVAValueMapProducer.src = cms.InputTag('reducedEgamma','reducedGedGsfElectrons')

    # To use older DataFormats, the electronMVAValueMapProducer MUST take a updated electron collection
    # such that the conversion variables are filled correctly.
    process.load("RecoEgamma.EgammaTools.gedGsfElectronsTo106X_cff")
    run2_miniAOD_80XLegacy.toModify(task, func=lambda t: t.add(process.gedGsfElectronsFrom80XTo106XTask))
    run2_miniAOD_80XLegacy.toModify(process.electronMVAValueMapProducer,
                                     keysForValueMaps = cms.InputTag('reducedEgamma','reducedGedGsfElectrons'),
                                     src = cms.InputTag("gedGsfElectronsFrom80XTo106X"))

    run2_miniAOD_94XFall17.toModify(task, func=lambda t: t.add(process.gedGsfElectronsFrom94XTo106XTask))
    run2_miniAOD_94XFall17.toModify(process.electronMVAValueMapProducer,
                                     keysForValueMaps = cms.InputTag('reducedEgamma','reducedGedGsfElectrons'),
                                     src = cms.InputTag("gedGsfElectronsFrom94XTo106X"))

    for idmod in electron_ids:
        setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection,None,False,task)

    #VID Photon IDs
    process.patPhotons.addPhotonID = cms.bool(True)
    photon_ids = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_V1_TrueVtx_cff',
                  'RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_V2_cff',
                  'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V1p1_cff',
                  'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V2_cff',
                  'RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring16_V2p2_cff',
                  'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring16_nonTrig_V1_cff']
    switchOnVIDPhotonIdProducer(process,DataFormat.AOD, task)
    process.egmPhotonIDs.physicsObjectSrc = cms.InputTag("reducedEgamma","reducedGedPhotons")
    process.photonMVAValueMapProducer.src = cms.InputTag('reducedEgamma','reducedGedPhotons')
    for idmod in photon_ids:
        setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection,None,False,task)

    #add the cut base IDs bitmaps of which cuts passed
    from RecoEgamma.EgammaTools.egammaObjectModifications_tools import makeVIDBitsModifier
    egamma_modifications.append(makeVIDBitsModifier(process,"egmGsfElectronIDs","egmPhotonIDs"))

    #-- Adding boosted taus
    from RecoTauTag.Configuration.boostedHPSPFTaus_cfi import addBoostedTaus
    addBoostedTaus(process)
    process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
    process.load("RecoTauTag.Configuration.HPSPFTaus_cff")
    #-- Adding customization for 94X 2017 legacy reMiniAOD
    from Configuration.Eras.Modifier_run2_miniAOD_94XFall17_cff import run2_miniAOD_94XFall17
    _makePatTausTaskWithRetrainedMVATauID = process.makePatTausTask.copy()
    _makePatTausTaskWithRetrainedMVATauID.add(process.hpsPFTauDiscriminationByIsolationMVArun2v1DBoldDMwLTTask,
                                              process.hpsPFTauDiscriminationByIsolationMVArun2v1DBnewDMwLTTask,
                                              process.hpsPFTauIsolationSums03Task,
                                              process.hpsPFTauDiscriminationByIsolationMVArun2v1DBdR03oldDMwLTTask)
    run2_miniAOD_94XFall17.toReplaceWith(
        process.makePatTausTask, _makePatTausTaskWithRetrainedMVATauID
        )
    #-- Adding customization for UL reMiniAOD:
    # running retrained anti-e MVA6 discriminants as setup in the HPS PFTau
    # master configuration for this era.
    _makePatTausTaskWithMVA6ElectronRejection = process.makePatTausTask.copy()
    _makePatTausTaskWithMVA6ElectronRejection.add(
        process.hpsPFTauDiscriminationByMVA6ElectronRejectionTask
    )
    run2_miniAOD_UL.toReplaceWith(
        process.makePatTausTask, _makePatTausTaskWithMVA6ElectronRejection
    )
    #-- Adding DeepTauID
    # deepTau v2
    _updatedTauName = 'slimmedTausDeepIDsv2'
    _noUpdatedTauName = 'slimmedTausNoDeepIDs'
    import RecoTauTag.RecoTau.tools.runTauIdMVA as tauIdConfig
    tauIdEmbedder = tauIdConfig.TauIDEmbedder(
        process, cms, debug = False,
        originalTauName = _noUpdatedTauName,
        updatedTauName = _updatedTauName,
        toKeep = ['deepTau2017v2']
    )
    tauIdEmbedder.runTauID()
    addToProcessAndTask(_noUpdatedTauName, process.slimmedTaus.clone(),process,task)
    delattr(process, 'slimmedTaus')
    process.slimmedTaus = getattr(process, _updatedTauName).clone()
    process.deepTauIDTask = cms.Task(process.deepTau2017v2, process.slimmedTaus)
    task.add(process.deepTauIDTask)

    # deepTau v2p1
    _updatedTauNameNew = 'slimmedTausDeepIDsv2p1'
    tauIdEmbedderNew = tauIdConfig.TauIDEmbedder(
        process, cms, debug = False,
        originalTauName = _noUpdatedTauName,
        updatedTauName = _updatedTauNameNew,
        toKeep = ['deepTau2017v2p1']
    )
    tauIdEmbedderNew.runTauID()
    deepTauIDTaskNew_ = cms.Task(process.deepTau2017v2p1,process.slimmedTaus)

    from Configuration.Eras.Modifier_run2_tau_ul_2016_cff import run2_tau_ul_2016
    from Configuration.Eras.Modifier_run2_tau_ul_2018_cff import run2_tau_ul_2018
    for era in [run2_miniAOD_UL,run2_tau_ul_2016,run2_tau_ul_2018]:
        era.toReplaceWith(process.slimmedTaus,
                          getattr(process, _updatedTauNameNew))
        era.toReplaceWith(process.deepTauIDTask,
                          deepTauIDTaskNew_)

    #-- Adding tauID against dead ECal towers to taus
    # note: default (AOD) behoviour of the tauID modified (leading track
    #       extrapolation to ECAL enabled) in the HPS PFTau master
    #       configuration for this era.
    _makePatTausTaskWithDeadECalVeto = process.makePatTausTask.copy()
    _makePatTausTaskWithDeadECalVeto.add(
        process.hpsPFTauDiscriminationByDeadECALElectronRejection
    )
    run2_miniAOD_UL.toReplaceWith(
        process.makePatTausTask, _makePatTausTaskWithDeadECalVeto
    )
    _withDeadEcalTauIDPs = cms.PSet(
        process.patTaus.tauIDSources,
        againstElectronDeadECAL = cms.InputTag("hpsPFTauDiscriminationByDeadECALElectronRejection")
    )
    run2_miniAOD_UL.toModify(process.patTaus,
                             tauIDSources = _withDeadEcalTauIDPs)
    #... and to boosted taus
    _withDeadEcalTauIDBoostedPs = cms.PSet(
        process.patTausBoosted.tauIDSources,
        againstElectronDeadECAL = cms.InputTag("hpsPFTauDiscriminationByDeadECALElectronRejectionBoosted")
    )
    run2_miniAOD_UL.toModify(process.patTausBoosted,
                             tauIDSources = _withDeadEcalTauIDBoostedPs)

    #-- Adding customization for 80X 2016 legacy reMiniAOD and 2018 heavy ions
    from Configuration.Eras.Modifier_run2_miniAOD_80XLegacy_cff import run2_miniAOD_80XLegacy
    from Configuration.Eras.Modifier_pp_on_AA_2018_cff import pp_on_AA_2018
    _makePatTausTaskWithTauReReco = process.makePatTausTask.copy()
    _makePatTausTaskWithTauReReco.add(process.PFTauTask)
    (run2_miniAOD_80XLegacy | pp_on_AA_2018).toReplaceWith(
        process.makePatTausTask, _makePatTausTaskWithTauReReco
        )

    # Adding puppi jets
    if not hasattr(process, 'ak4PFJetsPuppi'): #MM: avoid confilct with substructure call
        process.load('RecoJets.JetProducers.ak4PFJetsPuppi_cfi')
        task.add(process.ak4PFJets)
        task.add(process.ak4PFJetsPuppi)
    process.ak4PFJetsPuppi.doAreaFastjet = True # even for standard ak4PFJets this is overwritten in RecoJets/Configuration/python/RecoPFJets_cff
    from RecoJets.JetAssociationProducers.j2tParametersVX_cfi import j2tParametersVX
    process.ak4PFJetsPuppiTracksAssociatorAtVertex = cms.EDProducer("JetTracksAssociatorAtVertex",
        j2tParametersVX,
        jets = cms.InputTag("ak4PFJetsPuppi")
    )
    task.add(process.ak4PFJetsPuppiTracksAssociatorAtVertex)
    process.patJetPuppiCharge = cms.EDProducer("JetChargeProducer",
        src = cms.InputTag("ak4PFJetsPuppiTracksAssociatorAtVertex"),
        var = cms.string('Pt'),
        exp = cms.double(1.0)
    )
    task.add(process.patJetPuppiCharge)

    noDeepFlavourDiscriminators = [x.value() for x in process.patJets.discriminatorSources if not "DeepFlavour" in x.value()]
    addJetCollection(process, postfix   = "", labelName = 'Puppi', jetSource = cms.InputTag('ak4PFJetsPuppi'),
                    jetCorrections = ('AK4PFPuppi', ['L2Relative', 'L3Absolute'], ''),
                    pfCandidates = cms.InputTag("particleFlow"),
                    algo= 'AK', rParam = 0.4, btagDiscriminators = noDeepFlavourDiscriminators
                    )

    process.patJetGenJetMatchPuppi.matched = 'slimmedGenJets'

    process.patJetsPuppi.jetChargeSource = cms.InputTag("patJetPuppiCharge")

    process.selectedPatJetsPuppi.cut = cms.string("pt > 15")

    from PhysicsTools.PatAlgos.slimming.applyDeepBtagging_cff import applyDeepBtagging
    applyDeepBtagging( process )

    addToProcessAndTask('slimmedJetsPuppiNoMultiplicities', process.slimmedJetsNoDeepFlavour.clone(), process, task)
    process.slimmedJetsPuppiNoMultiplicities.src = cms.InputTag("selectedPatJetsPuppi")
    process.slimmedJetsPuppiNoMultiplicities.packedPFCandidates = cms.InputTag("packedPFCandidates")

    from PhysicsTools.PatAlgos.patPuppiJetSpecificProducer_cfi import patPuppiJetSpecificProducer
    from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
    process.patPuppiJetSpecificProducer = patPuppiJetSpecificProducer.clone(
      src=cms.InputTag("slimmedJetsPuppiNoMultiplicities"),
    )
    task.add(process.patPuppiJetSpecificProducer)
    updateJetCollection(
       process,
       labelName = 'PuppiJetSpecific',
       jetSource = cms.InputTag('slimmedJetsPuppiNoMultiplicities'),
    )
    process.updatedPatJetsPuppiJetSpecific.userData.userFloats.src = ['patPuppiJetSpecificProducer:puppiMultiplicity', 'patPuppiJetSpecificProducer:neutralPuppiMultiplicity', 'patPuppiJetSpecificProducer:neutralHadronPuppiMultiplicity', 'patPuppiJetSpecificProducer:photonPuppiMultiplicity', 'patPuppiJetSpecificProducer:HFHadronPuppiMultiplicity', 'patPuppiJetSpecificProducer:HFEMPuppiMultiplicity' ]
    process.slimmedJetsPuppi = process.selectedUpdatedPatJetsPuppiJetSpecific.clone()
    delattr(process, 'selectedUpdatedPatJetsPuppiJetSpecific')

    task.add(process.slimmedJetsPuppi)

    ## puppi met
    from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppies
    makePuppies( process );

    runMetCorAndUncForMiniAODProduction(process, metType="Puppi",
                                        pfCandColl=cms.InputTag("puppiForMET"),
                                        jetCollUnskimmed="slimmedJetsPuppi",
                                        recoMetFromPFCs=True,
                                        jetFlavor="AK4PFPuppi",
                                        postfix="Puppi"
                                        )

    process.load('PhysicsTools.PatAlgos.slimming.slimmedMETs_cfi')
    task.add(process.slimmedMETs)
    addToProcessAndTask('slimmedMETsPuppi', process.slimmedMETs.clone(), process, task)
    process.slimmedMETsPuppi.src = cms.InputTag("patMETsPuppi")
    process.slimmedMETsPuppi.rawVariation =  cms.InputTag("patPFMetPuppi")
    process.slimmedMETsPuppi.t1Uncertainties = cms.InputTag("patPFMetT1%sPuppi")
    process.slimmedMETsPuppi.t01Variation = cms.InputTag("patPFMetT0pcT1Puppi")
    process.slimmedMETsPuppi.t1SmearedVarsAndUncs = cms.InputTag("patPFMetT1Smear%sPuppi")
    process.slimmedMETsPuppi.tXYUncForRaw = cms.InputTag("patPFMetTxyPuppi")
    process.slimmedMETsPuppi.tXYUncForT1 = cms.InputTag("patPFMetT1TxyPuppi")
    process.slimmedMETsPuppi.tXYUncForT01 = cms.InputTag("patPFMetT0pcT1TxyPuppi")
    process.slimmedMETsPuppi.tXYUncForT1Smear = cms.InputTag("patPFMetT1SmearTxyPuppi")
    process.slimmedMETsPuppi.tXYUncForT01Smear = cms.InputTag("patPFMetT0pcT1SmearTxyPuppi")
    del process.slimmedMETsPuppi.caloMET

    process.load('RecoMET.METPUSubtraction.deepMETProducer_cfi')

    process.deepMETsResolutionTune = process.deepMETProducer.clone()
    process.deepMETsResponseTune = process.deepMETProducer.clone()
    process.deepMETsResponseTune.graph_path = 'RecoMET/METPUSubtraction/data/deepmet/deepmet_resp_v1_2018.pb'

    from Configuration.Eras.Modifier_run2_jme_2016_cff import run2_jme_2016
    run2_jme_2016.toModify(
        process.deepMETsResponseTune,
        graph_path="RecoMET/METPUSubtraction/data/deepmet/deepmet_resp_v1_2016.pb"
    )

    run2_miniAOD_UL.toModify(task, func=lambda t: t.add(process.deepMETsResolutionTune, process.deepMETsResponseTune))
    run2_miniAOD_UL.toModify(process.slimmedMETs, addDeepMETs = True)

    # add DetIdAssociatorRecords to EventSetup (for isolatedTracks)
    process.load("TrackingTools.TrackAssociator.DetIdAssociatorESProducer_cff")

    # EGamma objects from HGCal are not yet in GED
    # so add companion collections for Phase-II MiniAOD production
    from Configuration.Eras.Modifier_phase2_hgcal_cff import phase2_hgcal
    process.load("RecoEgamma.EgammaTools.slimmedEgammaFromMultiCl_cff")
    phase2_hgcal.toModify(task, func=lambda t: t.add(process.slimmedEgammaFromMultiClTask))

    # L1 pre-firing weights for 2016 and 2017
    from Configuration.Eras.Modifier_run2_L1prefiring_cff import run2_L1prefiring
    from Configuration.Eras.Modifier_stage2L1Trigger_cff import stage2L1Trigger
    from Configuration.Eras.Modifier_stage2L1Trigger_2017_cff import stage2L1Trigger_2017
    process.load("PhysicsTools.PatUtils.L1PrefiringWeightProducer_cff")
    stage2L1Trigger.toModify(process.prefiringweight, DataEraECAL = "2017BtoF", JetMaxMuonFraction = -1, DoMuons = cms.bool(False) )
    stage2L1Trigger_2017.toModify(process.prefiringweight, DataEraECAL = "2017BtoF", JetMaxMuonFraction = -1, DoMuons = cms.bool(False))
    run2_L1prefiring.toModify(task, func=lambda t: t.add(process.prefiringweight))

def miniAOD_customizeMC(process):
    task = getPatAlgosToolsTask(process)
    #GenJetFlavourInfos
    process.load("PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi")
    task.add(process.selectedHadronsAndPartons)
    task.add(process.selectedHadronsAndPartonsForGenJetsFlavourInfos)

    process.load("PhysicsTools.JetMCAlgos.AK4GenJetFlavourInfos_cfi")
    task.add(process.ak4GenJetFlavourInfos)

    process.load('PhysicsTools.PatAlgos.slimming.slimmedGenJetsFlavourInfos_cfi')
    task.add(process.slimmedGenJetsFlavourInfos)

    #slimmed pileup information
    process.load('PhysicsTools.PatAlgos.slimming.slimmedAddPileupInfo_cfi')
    task.add(process.slimmedAddPileupInfo)

    process.muonMatch.matched = "prunedGenParticles"
    process.electronMatch.matched = "prunedGenParticles"
    process.electronMatch.src = cms.InputTag("reducedEgamma","reducedGedGsfElectrons")
    process.photonMatch.matched = "prunedGenParticles"
    process.photonMatch.src = cms.InputTag("reducedEgamma","reducedGedPhotons")
    process.ootPhotonMatch.matched = "prunedGenParticles"
    process.ootPhotonMatch.src = cms.InputTag("reducedEgamma","reducedOOTPhotons")
    process.tauMatch.matched = "prunedGenParticles"
    process.tauGenJets.GenParticles = "prunedGenParticles"
    #Boosted taus
    process.tauMatchBoosted.matched = "prunedGenParticles"
    process.tauGenJetsBoosted.GenParticles = "prunedGenParticles"
    process.patJetPartons.particles = "genParticles"
    process.patJetPartonMatch.matched = "prunedGenParticles"
    process.patJetPartonMatch.mcStatus = [ 3, 23 ]
    process.patJetGenJetMatch.matched = "slimmedGenJets"
    process.patJetGenJetMatchAK8.matched =  "slimmedGenJetsAK8"
    process.patMuons.embedGenMatch = False
    process.patElectrons.embedGenMatch = False
    process.patPhotons.embedGenMatch = False
    process.patOOTPhotons.embedGenMatch = False
    process.patTaus.embedGenMatch = False
    process.patTausBoosted.embedGenMatch = False
    process.patJets.embedGenPartonMatch = False
    #also jet flavour must be switched
    process.patJetFlavourAssociation.rParam = 0.4

def miniAOD_customizeOutput(out):
    from PhysicsTools.PatAlgos.slimming.MicroEventContent_cff import MiniAODOverrideBranchesSplitLevel

    MiniAODOverrideBranchesSplitLevelWithBs = MiniAODOverrideBranchesSplitLevel.copy()
    MiniAODOverrideBranchesSplitLevelWithBs.extend([cms.untracked.PSet(branch = cms.untracked.string("recoVertexs_offlineSlimmedPrimaryVerticesWithBS__*"),splitLevel=cms.untracked.int32(99))])

    from Configuration.Eras.Modifier_bParking_cff import bParking
    from Configuration.Eras.Modifier_run2_miniAOD_devel_cff import run2_miniAOD_devel

    out.overrideBranchesSplitLevel = MiniAODOverrideBranchesSplitLevel
    (bParking | run2_miniAOD_devel).toModify(out, overrideBranchesSplitLevel = MiniAODOverrideBranchesSplitLevelWithBs)

    out.splitLevel = cms.untracked.int32(0)
    out.dropMetaData = cms.untracked.string('ALL')
    out.fastCloning= cms.untracked.bool(False)
    out.overrideInputFileSplitLevels = cms.untracked.bool(True)
    out.compressionAlgorithm = cms.untracked.string('LZMA')

def miniAOD_customizeData(process):
    from PhysicsTools.PatAlgos.tools.coreTools import runOnData
    runOnData( process, outputModules = [] )
    process.load("RecoCTPPS.Configuration.recoCTPPS_cff")
    process.load('L1Trigger.L1TGlobal.simGtExtFakeProd_cfi')
    task = getPatAlgosToolsTask(process)
    from Configuration.Eras.Modifier_ctpps_2016_cff import ctpps_2016
    from Configuration.ProcessModifiers.run2_miniAOD_UL_cff import run2_miniAOD_UL
    (ctpps_2016 & ~run2_miniAOD_UL).toModify(task, func=lambda t: t.add(process.ctppsLocalTrackLiteProducer, process.ctppsProtons))
    (ctpps_2016 & run2_miniAOD_UL).toModify(task, func=lambda t: t.add(process.recoCTPPSTask))
    run2_miniAOD_UL.toModify(task, func=lambda t: t.add(process.simGtExtUnprefireable))

def miniAOD_customizeAllData(process):
    miniAOD_customizeCommon(process)
    miniAOD_customizeData(process)
    return process

def miniAOD_customizeAllMC(process):
    miniAOD_customizeCommon(process)
    miniAOD_customizeMC(process)
    return process

def miniAOD_customizeAllMCFastSim(process):
    miniAOD_customizeCommon(process)
    miniAOD_customizeMC(process)
    from PhysicsTools.PatAlgos.slimming.metFilterPaths_cff import miniAOD_customizeMETFiltersFastSim
    process = miniAOD_customizeMETFiltersFastSim(process)
    from PhysicsTools.PatAlgos.slimming.isolatedTracks_cfi import miniAOD_customizeIsolatedTracksFastSim
    process = miniAOD_customizeIsolatedTracksFastSim(process)
    process.patMuons.addTriggerMatching = False

    return process
