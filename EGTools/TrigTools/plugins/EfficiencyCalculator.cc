#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
//#include "DataFormats/HLTReco/interface/EgammaObject.h"
#include "L1Trigger/L1TGlobal/plugins/L1TGlobalProducer.h"
//#include "HLTrigger/Egamma/plugins/HLTEgammaL1TMatchFilterRegional.cc"
#include "CondFormats/DataRecord/interface/L1TGlobalParametersRcd.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include <vector>
#include <string>
#include <iostream>
#include <TH1F.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include "Math/VectorUtil.h"
#include <stdlib.h>

#define TWOPI 6.283185308
/****************************************************************************
Plugin meant to calculate trigger efficiencies.
 

Author Riccardo Salvatico (KU), 2023
*****************************************************************************/

class EfficiencyCalculator : public edm::one::EDAnalyzer<edm::one::WatchRuns> {
 
private:
  
  edm::EDGetTokenT<std::vector<pat::Electron> > eleToken_;
  edm::EDGetTokenT<edm::TriggerResults > triggerResultsToken_;
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone> > triggerObjectsToken_;
  edm::EDGetTokenT<std::vector<reco::GenParticle> > genToken_;
  //edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupSummaryToken_;
  edm::Service<TFileService> fs;
  float barrel_end_ = 1.4442;
  float endcap_end_ = 2.5;
  bool  runningOnData_;
  int   runningEra_;

  TH2D* num_ele;
  TH2D* den_ele;

  double binsX[9] = {-2.5,-2.0,-1.566,-1.4442, 0.0, 1.4442, 1.566, 2.0, 2.5};
  double binsY[14] = {10,20,30,32,34,36,38,40,45,50,75,100,200,500};

  //TEfficiency *efficiency;

public:
  explicit EfficiencyCalculator(const edm::ParameterSet& iConfig);
  ~EfficiencyCalculator(){}
  
 private:
  virtual void beginRun(const edm::Run& run,const edm::EventSetup& iSetup);
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override{}
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endJob(){}
};


EfficiencyCalculator::EfficiencyCalculator(const edm::ParameterSet& iConfig):
  runningOnData_(iConfig.getParameter<bool>("runningOnData")),
  runningEra_(iConfig.getParameter<int>("runningEra"))
{
  //It is important to pick the "HLT" label of TriggerResults, because the RECO and PAT collections do not contain the HLT path and filter names 
  triggerResultsToken_ = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults","","HLT")); 
  triggerObjectsToken_ = consumes<std::vector<pat::TriggerObjectStandAlone> >(edm::InputTag("slimmedPatTrigger"));
  eleToken_            = consumes<std::vector<pat::Electron> >(edm::InputTag("slimmedElectrons"));
  if(!runningOnData_){
    genToken_          = consumes<std::vector<reco::GenParticle> >(edm::InputTag("prunedGenParticles"));
  }
  //pileupSummaryToken_ = consumes<std::vector<PileupSummaryInfo> >(edm::InputTag("addPileupInfo"));

  //pT(y) VS eta(x)
  num_ele = fs->make<TH2D>("num_ele",";#eta_{SC};#it{p}_{T} (GeV)",8,binsX,13,binsY);
  den_ele = fs->make<TH2D>("den_ele",";#eta_{SC};#it{p}_{T} (GeV)",8,binsX,13,binsY);

  //efficiency = fs->make<TEfficiency>("efficiency",";#eta_{SC};#it{p}_{T} (GeV)",8,binsX,13,binsY);
  //efficiency = fs->make<TEfficiency>(*num_ele,*den_ele);

}

// std::vector<const pat::TriggerObjectStandAlone*> matchTrigObjs(const float eta,const float phi,const float pT,const std::vector<pat::TriggerObjectStandAlone>& trigObjs,const float maxDeltaR=0.1, const float maxDpT=1.) {

//   std::vector<const pat::TriggerObjectStandAlone*> matchedObjs;
//   const float maxDR2 = maxDeltaR*maxDeltaR;
//   for(const auto trigObj : trigObjs){
//     const float deltaR2 = reco::deltaR2(eta, phi, trigObj.eta(), trigObj.phi());
//     if(deltaR2 < maxDR2 && (fabs(pT - trigObj.pt())/pT) < maxDpT) matchedObjs.push_back(&trigObj);
//   }
//   return matchedObjs;
// }


std::vector<const pat::TriggerObjectStandAlone*> matchTrigObjs(const float eta,const float phi,const float pT,const std::vector<pat::TriggerObjectStandAlone>& trigObjs,const float maxDeltaR=0.1, const float maxDpT=1.) {

  std::vector<const pat::TriggerObjectStandAlone*> matchedObjs;
  const float maxDR2 = maxDeltaR*maxDeltaR;
  for(auto& trigObj : trigObjs){
    const float deltaR2 = reco::deltaR2(eta, phi, trigObj.eta(), trigObj.phi());
    if(deltaR2 < maxDR2 && (fabs(pT - trigObj.pt())/pT) < maxDpT) matchedObjs.push_back(&trigObj);
  }
  return matchedObjs;
}


std::vector<const reco::GenParticle*> getGenparts(const std::vector<reco::GenParticle>& genparts,const int pid=11, bool antipart=true, const int status=1){

  std::vector<const reco::GenParticle*> selected;
  if(genparts.empty()) return selected;

  for(auto& part : genparts){
    const int pdg_id = part.pdgId();
    if(pdg_id == pid || (antipart && abs(pdg_id) == abs(pid))){
      if(part.isHardProcess() && status == 1){
	selected.push_back(&part);
      }
    }
  }
  return selected;
}


std::vector<pat::Electron>* matchToGen(const std::vector<pat::Electron>& recoElectrons, const std::vector<reco::GenParticle>& genparts,const int pid=11,bool antipart=true,const float max_dr=0.1,const int status=1){

  const float max_dr2 = max_dr*max_dr;
  bool isMatched = false;
  auto selectedParts = getGenparts(genparts,pid,antipart,status);
  std::vector<pat::Electron>* matchedRecoElectrons = nullptr;

  for(auto recoEle : recoElectrons){
    
    for(auto genPart : selectedParts){
      const float dr2 = reco::deltaR2(recoEle.eta(), recoEle.phi(), genPart->eta(), genPart->phi());
      if(dr2 < max_dr2){
	isMatched = true;
	break;
      }
    }
    if(isMatched) matchedRecoElectrons->push_back(recoEle);
  }
  return matchedRecoElectrons;
}


bool isGoodEvent(const std::vector<pat::Electron>& recoElectrons) {

  if(recoElectrons.size() < 2) return false;

  for(auto recoEle : recoElectrons){
    if(recoEle.electronID("cutBasedElectronID-Fall17-94X-V2-tight")) return true;
  }

  return false;
}


float calculateInvMass(const pat::Electron tagElectron, const pat::Electron probeCandidate) {
  
  TLorentzVector tag;
  TLorentzVector probe;

  tag.SetPxPyPzE(tagElectron.px(),tagElectron.py(),tagElectron.pz(),tagElectron.energy());
  probe.SetPxPyPzE(probeCandidate.px(),probeCandidate.py(),probeCandidate.pz(),probeCandidate.energy());

  float invMass = (tag + probe).M();

  return invMass;
  
}


//we need to initalise the menu each run (menu can and will change on run boundaries)
void EfficiencyCalculator::beginRun(const edm::Run& run,const edm::EventSetup& setup)
{
  //bool changed=false;
  //hltPSProv_.init(run,setup,hltProcess_,changed);
}

void EfficiencyCalculator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle<std::vector<pat::Electron> > ele;
  iEvent.getByToken(eleToken_,ele);

  edm::Handle<edm::TriggerResults > triggerResults;
  iEvent.getByToken(triggerResultsToken_,triggerResults);

  edm::Handle<std::vector<pat::TriggerObjectStandAlone> > triggerObjects;
  iEvent.getByToken(triggerObjectsToken_,triggerObjects);

  edm::Handle<std::vector<reco::GenParticle> > gen;
  if(!runningOnData_) iEvent.getByToken(genToken_,gen);
  

  // edm::Handle<std::vector<PileupSummaryInfo> >  PUInfo;
  // iEvent.getByToken(pileupSummaryToken_, PUInfo);
  
  // float nPU = -1;
 
  // std::vector<PileupSummaryInfo>::const_iterator PUi; 
 
  // for(PUi = PUInfo->begin(); PUi != PUInfo->end(); ++PUi) {
  //   const int BX = PUi->getBunchCrossing();
  //   if(BX == 0){
  //     nPU = PUi->getTrueNumInteractions(); //Corresponding to the in-time PU
  //   }
  // }

  auto electrons = ele.product();

  auto trigRes = triggerResults.product();

  if(!runningOnData_) electrons = matchToGen(*electrons,*gen.product());
  if(electrons->size() == 0 || electrons == nullptr) return;
  if(!isGoodEvent(*electrons)) return;

  std::vector<pat::Electron> listOfTags;
  std::vector<pat::TriggerObjectStandAlone> unpackedTrigObjs;


  //const edm::TriggerNames &names = iEvent.triggerNames(*trigRes);

  // for (unsigned int i = 0, n = trigRes->size(); i < n; ++i){
  //   if(!trigRes->accept(i)) continue;
  //   std::string tmp_triggername = names.triggerName(i);
  //   std::cout << "name: " << tmp_triggername << std::endl;
  // } 

  for(auto& trigObj : *triggerObjects){
      unpackedTrigObjs.push_back(trigObj);
      unpackedTrigObjs.back().unpackFilterLabels(iEvent,*trigRes);
  }

  //Fill in a list of Tags
  for(auto& el : *electrons){

    if(fabs(el.superCluster()->eta()) > barrel_end_ || el.pt() < 32.) continue;

    bool isTagTriggerMatched = false;
    //std::vector<const pat::TriggerObjectStandAlone*> matchedTrigObjs = matchTrigObjs(el.superCluster()->eta(),el.superCluster()->phi(),el.pt(),unpackedTrigObjs);

    auto matchedTrigObjs = matchTrigObjs(el.superCluster()->eta(),el.superCluster()->phi(),el.pt(),unpackedTrigObjs);
    //std::cout << "matchedTrigObjs size: " << matchedTrigObjs.size() << std::endl;

    for(const auto trigObj : matchedTrigObjs) {
      if(trigObj->hasFilterLabel("hltEle32WPTightGsfTrackIsoFilter")) isTagTriggerMatched = true;
    }

    if(isTagTriggerMatched) listOfTags.push_back(el);
  }

  //std::cout << "listOfTags size: " << listOfTags.size() << std::endl;

  //Selection on Probes
  for(auto& el : *electrons){

    bool isProbeTriggerMatched = false;    

    for(auto tag : listOfTags){
      
      float invMass = calculateInvMass(tag, el);

      double eventWeight = 1.;
      
      if(tag.pt() != el.pt() && invMass > 70. && invMass < 110. && el.electronID("cutBasedElectronID-Fall17-94X-V2-tight")){
	//std::cout << "SC: " << el.superCluster()->eta() << "  pT: " << el.pt() << std::endl;
	den_ele->Fill(el.superCluster()->eta(),el.pt(),eventWeight);
	//efficiency->Fill(false,el.superCluster()->eta(),el.pt());

	//std::vector<const pat::TriggerObjectStandAlone*> matchedTrigObjs = matchTrigObjs(el.superCluster()->eta(),el.superCluster()->phi(),el.pt(),unpackedTrigObjs);
	auto matchedTrigObjs = matchTrigObjs(el.superCluster()->eta(),el.superCluster()->phi(),el.pt(),unpackedTrigObjs);

	for(const auto trigObj : matchedTrigObjs) {
	  if(trigObj->hasFilterLabel("hltEle32WPTightGsfTrackIsoFilter")) isProbeTriggerMatched = true;
	}

	if(isProbeTriggerMatched){
	  num_ele->Fill(el.superCluster()->eta(),el.pt(),eventWeight);
	  //efficiency->Fill(true,el.superCluster()->eta(),el.pt());
	}
      }
    }
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(EfficiencyCalculator);
