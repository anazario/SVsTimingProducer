// -*- C++ -*-
//
// Package:    TimingWithSVs/SVsTimingProducer
// Class:      SVsTimingProducer
//
/**\class SVsTimingProducer SVsTimingProducer.cc TimingWithSVs/SVsTimingProducer/plugins/SVsTimingProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Andres Abreu
//         Created:  Sat, 21 Oct 2023 20:06:20 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"


//
// class declaration
//

class LostTrackProducer : public edm::stream::EDProducer<> {
public:
  explicit LostTrackProducer(const edm::ParameterSet&);
  ~LostTrackProducer() = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:

  enum class TrkStatus { NOTUSED = 0, PFCAND, PFCANDNOTRKPROPS, PFELECTRON, PFPOSITRON, VTX };

  void produce(edm::Event&, const edm::EventSetup&) override;
  bool passTrkCuts(const reco::Track& tr) const;

  // ----------member data ---------------------------
  const edm::EDGetTokenT<reco::PFCandidateCollection> pfCandsToken_;
  const edm::EDGetTokenT<reco::TrackCollection> generalTrackToken_;
  const edm::EDGetTokenT<reco::VertexCollection> svsToken_;
  const double minPt_;
  const double minHits_;
  const double minPixelHits_;
};

LostTrackProducer::LostTrackProducer(const edm::ParameterSet& iConfig) 
  : pfCandsToken_(consumes<reco::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("inputCandidates"))),
    generalTrackToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("inputTracks"))),
    svsToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("secondaryVertices"))),
    minPt_(iConfig.getParameter<double>("minPt")),
    minHits_(iConfig.getParameter<uint32_t>("minHits")),
    minPixelHits_(iConfig.getParameter<uint32_t>("minPixelHits")) {

  produces<std::vector<reco::Track>>("lostTracks");
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void LostTrackProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<reco::PFCandidateCollection> pfCandsHandle;
  iEvent.getByToken(pfCandsToken_, pfCandsHandle);

  edm::Handle<reco::TrackCollection> generalTracksHandle;
  iEvent.getByToken(generalTrackToken_, generalTracksHandle);

  edm::Handle<reco::VertexCollection> svsHandle;
  iEvent.getByToken(svsToken_, svsHandle);

  //std::cout << "Total general tracks: " << generalTracksHandle->size() << std::endl;

  auto lostTracks = std::make_unique<reco::TrackCollection>();

  std::vector<TrkStatus> trackStatus(generalTracksHandle->size(), TrkStatus::NOTUSED);

  for (unsigned int ic = 0, nc = pfCandsHandle->size(); ic < nc; ++ic) {
    //edm::Ref<reco::PFCandidateCollection> r(pfCandsHandle, ic);
    const reco::PFCandidate& cand = (*pfCandsHandle)[ic];
    if (cand.charge() && cand.trackRef().isNonnull() && cand.trackRef().id() == generalTracksHandle.id()) {
      if (cand.pdgId() == 11)
	trackStatus[cand.trackRef().key()] = TrkStatus::PFELECTRON;
      else if (cand.pdgId() == -11)
	trackStatus[cand.trackRef().key()] = TrkStatus::PFPOSITRON;
      //else if ((*pf2pc)[r]->numberOfHits() > 0)
      //trkStatus[cand.trackRef().key()] = TrkStatus::PFCAND;
      else
	trackStatus[cand.trackRef().key()] = TrkStatus::PFCANDNOTRKPROPS;
    }
  }

  for (const auto& secVert : *svsHandle) {
    for (auto trkIt = secVert.tracks_begin(); trkIt != secVert.tracks_end(); trkIt++) {
      if (trackStatus[trkIt->key()] == TrkStatus::NOTUSED)
	trackStatus[trkIt->key()] = TrkStatus::VTX;
    }
  }

  for (unsigned int trkIndx = 0; trkIndx < generalTracksHandle->size(); trkIndx++) {
    reco::TrackRef trk(generalTracksHandle, trkIndx);
    if (trackStatus[trkIndx] == TrkStatus::VTX || (trackStatus[trkIndx] == TrkStatus::NOTUSED && passTrkCuts(*trk))) {
      lostTracks->emplace_back(*trk);
    }
  }

  //std::cout << "total lost tracks: " << lostTracks->size() << std::endl;

  iEvent.put(std::move(lostTracks), "lostTracks");  
}

bool LostTrackProducer::passTrkCuts(const reco::Track& tr) const {
  const bool passTrkHits = tr.pt() > minPt_ && tr.numberOfValidHits() >= minHits_ &&
    tr.hitPattern().numberOfValidPixelHits() >= minPixelHits_;
  //const bool passTrkQual = passesQuality(tr, qualsToAutoAccept_);
  
  //return passTrkHits || passTrkQual || passThroughCut_(tr);
  return passTrkHits;  
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void LostTrackProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(LostTrackProducer);
