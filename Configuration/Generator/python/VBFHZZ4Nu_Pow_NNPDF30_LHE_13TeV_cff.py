import FWCore.ParameterSet.Config as cms

# link to card:
# https://github.com/cms-sw/genproductions/tree/master/bin/Powheg/production/V2/13TeV/Higgs/VBF_H_NNPDF30_13TeV/VBF_H_NNPDF30_13TeV_M-125.input

externalLHEProducer = cms.EDProducer("ExternalLHEProducer",
    args = cms.vstring('/cvmfs/cms.cern.ch/phys_generator/gridpacks/slc6_amd64_gcc630/13TeV/Powheg/V2/RelValidation/VBFH/VBF_H_slc6_amd64_gcc630_CMSSW_9_3_9_patch1_VBF_new.tgz'),
    nEvents = cms.untracked.uint32(5000),
    numberOfParameters = cms.uint32(1),
    outputFile = cms.string('cmsgrid_final.lhe'),
    scriptName = cms.FileInPath('GeneratorInterface/LHEInterface/data/run_generic_tarball_cvmfs.sh'),
    generateConcurrently = cms.untracked.bool(True)
)
