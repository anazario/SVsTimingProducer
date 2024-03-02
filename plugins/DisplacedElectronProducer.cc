// -*- C++ -*-
//
// Package:    TimingWithSVs/SVsTimingProducer
// Class:      DisplacedElectronProducer
//
/**\class DisplacedElectronProducer DisplacedElectronProducer.cc TimingWithSVs/SVsTimingProducer/plugins/DisplacedElectronProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Andres Abreu
//         Created:  Wed, 10 Jan 2024 12:11:43 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/Point3D.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/JetReco/interface/TrackExtrapolation.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
//#include "RecoJets/JetAssociationProducers/interface/TrackExtrapolator.h"

#include "TrackingTools/Records/interface/DetIdAssociatorRecord.h"
#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"
#include "TrackingTools/TrackAssociator/interface/TrackAssociatorParameters.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/GeomPropagators/interface/HelixExtrapolatorToLine2Order.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"

//#include "CondFormats/HcalObjects/interface/HcalChannelQuality.h"
//#include "CondFormats/HcalObjects/interface/HcalChannelStatus.h"
//#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
//#include "CondFormats/DataRecord/interface/HcalChannelQualityRcd.h"
//#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>

//Local includes
#include "TimingWithSVs/SVsTimingProducer/interface/Hungarian.h"
#include "TimingWithSVs/SVsTimingProducer/interface/DeltaRMatch.h"
#include "TimingWithSVs/SVsTimingProducer/interface/MatchTracksToSC.h"
#include "GenElectronClassifier.h"
//
// class declaration
//

bool IsIndexMatched(const std::vector<int>& matchedIndexes, int index);

typedef ROOT::Math::PtEtaPhiMVector LorentzVec;

class DisplacedElectronProducer : public edm::stream::EDProducer<> {
public:
  explicit DisplacedElectronProducer(const edm::ParameterSet&);
  ~DisplacedElectronProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  typedef edm::Handle<reco::TrackCollection> TracksHandle;
  typedef edm::Handle<reco::GsfTrackCollection> GsfTracksHandle;
  typedef edm::Handle<reco::SuperClusterCollection> SCsHandle;
  
  void produce(edm::Event&, const edm::EventSetup&) override;

  void GetSignalGenElectrons(edm::Event& iEvent);

  template <typename T>
  int FindIndex(const std::vector<T>& vec, const T& value) const;
  void AssignTrack(reco::Electron &electron, const reco::TrackRef &trackRef) const;
  void AssignTrack(reco::Electron &electron, const reco::GsfTrackRef &trackRef) const;
  template <typename T>
  reco::ElectronCollection MakeElectrons(edm::Handle<std::vector<T>> &tracks, 
					 SCsHandle &superClusters, 
					 const std::vector<int> &trackIndex, 
					 const std::vector<int> &scIndex) const;
  template <typename T>
  reco::ElectronCollection MatchTracksToSuperClusters(edm::Handle<std::vector<T>> &tracks, 
						      SCsHandle &superClusters, 
						      edm::Event& iEvent,
						      const edm::EventSetup &iSetup,
						      float minDeltaR = 1.);

  // ----------member data ---------------------------
  //edm::EDGetTokenT<reco::PFCandidateCollection> pfCollectionToken_;
  edm::EDGetTokenT<reco::TrackCollection> generalTrackToken_;
  edm::EDGetTokenT<reco::GsfTrackCollection> gsfElectronTrackToken_;
  //edm::EDGetTokenT<reco::GsfTrackCollection> lowPtGsfElectronTrackToken_;
  //edm::EDGetTokenT<reco::GsfElectronCollection> electronsToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  edm::EDGetTokenT<reco::SuperClusterCollection> superClusterToken_;
  edm::EDGetTokenT<reco::SuperClusterCollection> ootSuperClusterToken_;

  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;
  edm::EDGetTokenT<reco::TrackCollection> displacedTrackToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magneticFieldToken_;
  //edm::EDGetTokenT<std::vector<reco::TrackExtrapolation> > trackExtrapolationToken_;

  TrackDetectorAssociator trackAssociator_;
  TrackAssociatorParameters trackAssocParameters_;

  reco::GenParticleCollection signalGenElectrons_;

};

DisplacedElectronProducer::DisplacedElectronProducer(const edm::ParameterSet& iConfig) :
  //pfCollectionToken_(consumes<reco::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCollectionSrc")) ),
  generalTrackToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("generalTracksSrc")) ),
  gsfElectronTrackToken_(consumes<reco::GsfTrackCollection>(iConfig.getParameter<edm::InputTag>("gsfElectronTracksSrc")) ),
  //lowPtGsfElectronTrackToken_(consumes<reco::GsfTrackCollection>(iConfig.getParameter<edm::InputTag>("lowPtGsfElectronTracksSrc")) ),
  //electronsToken_(consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("electronSrc")) ),
  genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticleSrc")) ),
  superClusterToken_(consumes<reco::SuperClusterCollection>(iConfig.getParameter<edm::InputTag>("superClusters")) ),
  ootSuperClusterToken_(consumes<reco::SuperClusterCollection>(iConfig.getParameter<edm::InputTag>("ootSuperClusters")) ),

  caloGeometryToken_(esConsumes()),
  displacedTrackToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("displacedTracksSrc")) ),
  magneticFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>())
  //trackExtrapolationToken_(consumes<std::vector<reco::TrackExtrapolation> >(iConfig.getParameter<edm::InputTag>("trackExtrapolatorSrc")) ) 
{
  // TrackAssociator parameters
  edm::ParameterSet parameters = iConfig.getParameter<edm::ParameterSet>("TrackAssociatorParameters");
  edm::ConsumesCollector iC = consumesCollector();
  trackAssocParameters_.loadParameters(parameters, iC);
  trackAssociator_.useDefaultPropagator();

  produces<reco::ElectronCollection>("displacedElectrons").setBranchAlias("displacedElectrons");
}

DisplacedElectronProducer::~DisplacedElectronProducer() {}

//
// member functions
//
// ------------ method called to produce the data  ------------
void DisplacedElectronProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace std;

  typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double> > Point;
  typedef ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double> > Vector;

  auto magfield = iSetup.getTransientHandle(magneticFieldToken_);
  CaloGeometry ecalGeometry = iSetup.getData(caloGeometryToken_);

  //edm::Handle<reco::PFCandidateCollection> pfCandidates;
  //iEvent.getByToken(pfCollectionToken_, pfCandidates);

  edm::Handle<reco::TrackCollection> generalTracks;
  iEvent.getByToken(generalTrackToken_, generalTracks);

  //edm::Handle<reco::GsfElectronCollection> electrons;
  //iEvent.getByToken(electronsToken_, electrons);
  
  edm::Handle<reco::TrackCollection> displacedTracks;
  iEvent.getByToken(displacedTrackToken_, displacedTracks);

  edm::Handle<reco::GsfTrackCollection> gsfTracks;
  iEvent.getByToken(gsfElectronTrackToken_, gsfTracks);

  //edm::Handle<vector<reco::TrackExtrapolation> > extrapolatedTracks;
  //iEvent.getByToken(trackExtrapolationToken_, extrapolatedTracks);

  edm::Handle<reco::SuperClusterCollection> superClusterCollection;
  iEvent.getByToken(superClusterToken_, superClusterCollection);

  edm::Handle<reco::SuperClusterCollection> ootSuperClusterCollection;
  iEvent.getByToken(ootSuperClusterToken_, ootSuperClusterCollection);

  // Smart pointers to containers of output collections
  std::unique_ptr<reco::ElectronCollection> displacedElectrons = std::make_unique<reco::ElectronCollection>();

  reco::TrackCollection allTracks(*generalTracks), matchedTracks;

  GetSignalGenElectrons(iEvent); 

  // Add tracks from other collections to generalTracks that don't overlap
  allTracks.insert(allTracks.end(), displacedTracks->begin(), displacedTracks->end());
  allTracks.insert(allTracks.end(), gsfTracks->begin(), gsfTracks->end());

  if (signalGenElectrons_.size() == 0) {
    iEvent.put(std::move(displacedElectrons), "displacedElectrons");
    return;
  }
  /* Requires TrackExtra!!!
  const size_t ngsf = gsfTracks->size();
  for (unsigned int igsf = 0; igsf < ngsf; ++igsf) {
    reco::GsfTrackRef gref(gsfTracks, igsf);
    reco::TrackRef trk = gref->seedRef().castTo<reco::ElectronSeedRef>()->ctfTrack();

    if (trk.id() != generalTracks.id()) 
      std::cout << "tracks could not be matched!" << std::endl; 
  }
  */
  //double cost = assigner.Solve(mat, assignments);

  //reco::SuperClusterCollection allSCs(*superClusterCollection);
  //allSCs.insert(allSCs.end(), ootSuperClusterCollection->begin(), ootSuperClusterCollection->end());

  DeltaRMatch<reco::Track, reco::GenParticle> match(allTracks, signalGenElectrons_, 1);
  
  //const std::vector<int> matchedTrackIndexes(match.GetMatchedIndexesA());
  //const std::vector<int> matchedGenIndexes(match.GetMatchedIndexesB());
  //const std::vector<double> matchedDeltaRs(match.GetMatchedDeltaRs());
  /*
  for(int i = 0; i < match.GetNMatches(); i++) {
    PrintGenElectronInfo(signalGenElectrons_[matchedGenIndexes[i]]);
    //std::cout << "DeltaR: " << matchedDeltaRs[i] << ", pt = " << allTracks[matchedTrackIndexes[i]].pt() 
    //<< ", Charge: " << allTracks[matchedTrackIndexes[i]].charge() << std::endl;
  }
  std::cout << std::endl;
  */
  
  std::vector<reco::ElectronCollection> temp_electrons = {MatchTracksToSuperClusters(generalTracks, superClusterCollection, iEvent, iSetup),
							  MatchTracksToSuperClusters(generalTracks, ootSuperClusterCollection, iEvent, iSetup),
							  MatchTracksToSuperClusters(displacedTracks, superClusterCollection, iEvent, iSetup),
							  MatchTracksToSuperClusters(displacedTracks, ootSuperClusterCollection, iEvent, iSetup),	
							  MatchTracksToSuperClusters(gsfTracks, superClusterCollection, iEvent, iSetup),
							  MatchTracksToSuperClusters(gsfTracks, ootSuperClusterCollection, iEvent, iSetup)};
  
  reco::TrackCollection generalTrackVec;
  for(auto const &track : *generalTracks)
    if(track.pt() > 2)
      generalTrackVec.emplace_back(track);

  reco::SuperClusterCollection superClusters;
  for(auto const &sc : *superClusterCollection) 
    superClusters.emplace_back(sc);
  

  MatchTracksToSC<reco::Track> assigner(iEvent, iSetup, magfield, ecalGeometry, trackAssocParameters_, generalTrackVec, superClusters);
  /*
  int electronSetIndex(0);
  for(auto const &electrons : temp_electrons) {
    std::cout << "Electron set " << electronSetIndex << " with " << electrons.size() << " electrons: " << std::endl;
    displacedElectrons->insert(displacedElectrons->end(), electrons.begin(), electrons.end());
    electronSetIndex++;
  }
  std::cout << "Saved " << displacedElectrons->size() << " electrons" << std::endl;
  */
  for(auto const &electron : *displacedElectrons) {
    electron.gsfTrack();
  }
  //HungarianAlgorithm assigner;
  
  
  iEvent.put(std::move(displacedElectrons), "displacedElectrons"); 
}// Producer end

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DisplacedElectronProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  /*
  desc.add<edm::InputTag>("generalTracksSrc", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("gsfElectronTracksSrc", edm::InputTag("electronGsfTracks"));
  desc.add<edm::InputTag>("displacedTracksSrc", edm::InputTag("displacedTracks"));
  desc.add<edm::InputTag>("superClusters", edm::InputTag("particleFlowEGamma"));
  desc.add<edm::InputTag>("ootSuperClusters", edm::InputTag("particleFlowSuperClusterOOTECAL"));
  desc.add<edm::InputTag>("genParticleSrc", edm::InputTag("genParticles"));
  desc.add<edm::InputTag>("displacedElectrons", edm::InputTag("displacedElectrons"));
  desc.add("tkAssocParamBlock");
  */
  descriptions.addDefault(desc);
  
}

void DisplacedElectronProducer::GetSignalGenElectrons(edm::Event& iEvent) {
  
  signalGenElectrons_.clear();

  edm::Handle<reco::GenParticleCollection> genCandidates;
  iEvent.getByToken(genParticlesToken_, genCandidates);

  for(const auto &genCandidate : *genCandidates) {

    if(genCandidate.status() != 1)
      continue;

    if(abs(genCandidate.pdgId()) == 11 && isSignalGenElectron(genCandidate)) {
      signalGenElectrons_.push_back(genCandidate);
    }
     
  }
}

template <typename T>
int DisplacedElectronProducer::FindIndex(const std::vector<T>& vec, const T& value) const {
  auto it = std::find(vec.begin(), vec.end(), value);
  if (it != vec.end()) {
    return std::distance(vec.begin(), it);
  } else {
    // Return -1 or some other value to indicate that the element was not found
    return -1;
  }
}

void DisplacedElectronProducer::AssignTrack(reco::Electron &electron, const reco::TrackRef &trackRef) const {
  electron.setTrack(trackRef);
}

void DisplacedElectronProducer::AssignTrack(reco::Electron &electron, const reco::GsfTrackRef &trackRef) const {
  electron.setGsfTrack(trackRef);
}

template <typename T>
reco::ElectronCollection DisplacedElectronProducer::MakeElectrons(edm::Handle<std::vector<T>> &tracks,
								  SCsHandle &superClusters,
								  const std::vector<int> &trackIndex,
								  const std::vector<int> &scIndex) const {

  reco::ElectronCollection electrons;

  for(size_t i = 0; i < trackIndex.size(); i++) {
    const T track(tracks->at(i));
    const reco::SuperCluster superCluster(superClusters->at(i));

    const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > electron4Vec(track.px(), track.py(), track.pz(), superCluster.rawEnergy());
    reco::Electron matchedElectron(track.charge(), electron4Vec, math::XYZPoint(track.vx(), track.vy(), track.vz()));
    matchedElectron.setSuperCluster(reco::SuperClusterRef(superClusters, scIndex[i]));
    AssignTrack(matchedElectron, edm::Ref<std::vector<T>>(tracks, trackIndex[i]));

    electrons.push_back(matchedElectron);
  }

  return electrons;
}

template <typename T>
reco::ElectronCollection DisplacedElectronProducer::MatchTracksToSuperClusters(edm::Handle<std::vector<T>> &tracks, 
									       SCsHandle &superClusters, 
									       edm::Event& iEvent,
									       const edm::EventSetup &iSetup,
									       float minDeltaR) {

  reco::ElectronCollection electrons;

  auto magfield = iSetup.getTransientHandle(magneticFieldToken_);
  CaloGeometry ecalGeometry = iSetup.getData(caloGeometryToken_);
  
  int trackIndex(0), scMatchedIndex(-1), matchedTrackIndex(-1);
  std::vector<int> matchedSCIndexes, matchedTrackIndexes;
  std::vector<float> matchedDeltaR;
  
  // Loop over track collection
  for(auto const& track : *tracks) {

    bool isMatched(false);
    float deltaR(-1.);
    int minMatchedDeltaR(minDeltaR);

    reco::SuperCluster matchedSC;

    FreeTrajectoryState initialState = trajectoryStateTransform::initialFreeState(track, magfield.product());
    TrackDetMatchInfo trackDetInfo = trackAssociator_.associate(iEvent, iSetup, trackAssocParameters_, &initialState);

    // Loop over ECAL det IDs were tracks are propagated
    for(const auto &id : trackDetInfo.crossedEcalIds) {
      GlobalPoint trackHitAtEcal(ecalGeometry.getGeometry(id)->getPosition());
      const float ecalHit_eta(trackHitAtEcal.eta());
      const float ecalHit_phi(trackHitAtEcal.phi());

      int scIndex(0);
      for(const auto &sc : *superClusters) {

        const float superCluster_eta(sc.eta());
        const float superCluster_phi(sc.phi());

        deltaR = sqrt(reco::deltaR2(ecalHit_eta, ecalHit_phi, superCluster_eta, superCluster_phi) );

	// Make sure supercluster isn't already matched to a track                                                                                                               
        if(IsIndexMatched(matchedSCIndexes,scIndex)) 
	  continue;

        if(deltaR < minMatchedDeltaR) {
	  minMatchedDeltaR = deltaR;
	  isMatched = true;
          matchedSC = sc;

	  // Save index of matched track and supercluster
          scMatchedIndex = scIndex;
	  matchedTrackIndex = trackIndex;
        }
        scIndex++;
      }// End SuperCluster Loop
    }// End DetID Loop 
    if(isMatched) {

      // Get match list index of previous track that best-matched this same supercluster 
      const int matchIndex(FindIndex(matchedSCIndexes, scMatchedIndex) ) ;
      if(matchIndex > 0) {
	// if the new matched track has a smaller deltaR than the old match replace it
	if(matchedDeltaR[matchIndex] < minMatchedDeltaR) {
	  matchedTrackIndexes[matchIndex] = matchedTrackIndex;
	  continue;
	}
      }
      
      // Save match information (track and supercluster index and match DeltaR)
      matchedSCIndexes.push_back(scMatchedIndex);
      matchedTrackIndexes.push_back(matchedTrackIndex);
      matchedDeltaR.push_back(minMatchedDeltaR);

      // Construct displaced electron object
      electrons = MakeElectrons(tracks, superClusters, matchedTrackIndexes, matchedSCIndexes);

      //std::cout << "track " << trackIndex;
      //std::cout << " matched to supercluster " << scMatchedIndex << " with track: deltaR = " << deltaR << ", raw energy: " << matchedSC.rawEnergy() << std::endl;
    }

    trackIndex++;
  }//End track loop
  //std::cout << std::endl;

  return electrons;
}

bool IsIndexMatched(const std::vector<int>& matchedIndexes, int index) {
  // Check if the index is present in the matchedIndexes vector
  return std::find(matchedIndexes.begin(), matchedIndexes.end(), index) != matchedIndexes.end();
}

//define this as a plug-in
DEFINE_FWK_MODULE(DisplacedElectronProducer);
