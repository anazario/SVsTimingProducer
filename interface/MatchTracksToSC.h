#ifndef TimingWithSVs_SVsTimingProducer_MatchTracksToSC_h
#define TimingWithSVs_SVsTimingProducer_MatchTracksToSC_h

// -*- C++ -*-
//
// Package:    TimingWithSVs/SVsTimingProducer
// Class:      MatchTracksToSC
//
/**\class MatchTracksToSC MatchTracksToSC.h TimingWithSVs/SVsTimingProducer/interface/MatchTracksToSC.h

   Description: Class for matching any set of tracks to a set of reco::SuperCluster. 

   Implementation:                             

   This class uses the Hungarian algorithm for solving the assignment problem using 
   a cost matrix of the deltaR values between each possible pairing. The goal of the 
   algorithm is to obtain the assignments for which the total sum of deltaR values 
   is the smallest.   

   To run the algorithm all that is needed is to call the constructor with all the 
   properly defined inputs (some of which require cmssw).  

*/
//
// Original Author:  Andres Abreu
//         Created:  Fri, 26 Jan 2024 12:11:43 GMT
//
//

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "TrackingTools/Records/interface/DetIdAssociatorRecord.h"
#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"
#include "TrackingTools/TrackAssociator/interface/TrackAssociatorParameters.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "TimingWithSVs/SVsTimingProducer/interface/Hungarian.h"
#include "TimingWithSVs/SVsTimingProducer/interface/DeltaRMatch.h"


//===============================================================================================//
//                                  class: MatchedTrackSCPair                                    //
//-----------------------------------------------------------------------------------------------// 
// Matched track-to-supercluster pair class used internally by the MatchTracksToSC class. Aside  //
// from the matched track and supercluster it stores the deltaR of the match and the location of //
// the extrapolated track at the face of the ECAL. Can take any type of track as input.          //
//-----------------------------------------------------------------------------------------------//
//===============================================================================================//
template <class T> class MatchedTrackSCPair {

 public:
  MatchedTrackSCPair() = default;
  virtual ~MatchedTrackSCPair() = default;

 MatchedTrackSCPair(const T &track, const reco::SuperCluster &superCluster, const double deltaR, const GlobalPoint &trackLocationAtECAL)
   : track_(track), superCluster_(superCluster), deltaR_(deltaR), trackLocationAtECAL_(trackLocationAtECAL) {}


  // Getter functions
  double GetDeltaR() const { return deltaR_; }
  T GetTrack() const { return track_; }
  reco::SuperCluster GetSuperCluster() const { return superCluster_; }
  GlobalPoint GetTrackLocationAtECAL() const { return trackLocationAtECAL_; }

 private:

  T track_;
  reco::SuperCluster superCluster_;
  double deltaR_;
  GlobalPoint trackLocationAtECAL_;

};

template <class T> class MatchedTrackSCPairs : public std::vector<MatchedTrackSCPair<T>> {
 public:

 MatchedTrackSCPairs() 
   : cost_(-1) {}

 MatchedTrackSCPairs(const double &cost) 
   : cost_(cost) {}

  virtual ~MatchedTrackSCPairs() = default;

  double GetCost() const {return cost_;}

 private:

  double cost_;  
  
};

//===============================================================================================//
//                                     class: MatchTrackToSC                                     //
//-----------------------------------------------------------------------------------------------//
// Class for solving the assignment problem between a set of tracks (any type) and a set of reco //
// superclusters (reco::SuperCluster). The class has only one constructor which takes a total of //
// 7 arguments: the edm::Event and edm::EventSetup, accessible from any CMSSW analyzer/producer, //
// an edm::ESTransientHandle to the magnetic field in that event, a CaloGeometry object, the     //
// track association parameters needed for the propagation tool and the collection of tracks and //
// superclusters to be matched. This implementation makes uses the Hungarian (aka Kuhn-Munkres)  //
// algorithm to solve the assignment problem, using deltaR between the supercluster location and //
// the propagated track location at the ECAL. The following methods are available after calling  //
// the constructor:                                                                              //
//                                                                                               // 
//-----------------------------------------------------------------------------------------------//
//                                   method: double GetCost()                                    //
// Description:                                                                                  // 
// Returns the total cost associated with the solved assignment problem. This value is simply    //
// the net sum of all the match deltaR's in the event.                                           //
//                                                                                               //
//-----------------------------------------------------------------------------------------------//
//                         method: std::vector<T> GetMatchedTracks()                             //
// Description:                                                                                  //
// Returns the complete set of matched tracks in order of ascending deltaR. T indicates the type //
// of track used as input.                                                                       //
//                                                                                               //
//-----------------------------------------------------------------------------------------------//
//                    method: reco::SuperClusterCollection GetMatchedSCs()                       //
// Description:                                                                                  //
// Returns the set of matched superclusters in an event in order of ascending values of deltaR.  //
//                                                                                               //
//-----------------------------------------------------------------------------------------------//  
//             method: std::vector<MatchedTrackSCPair<T>> GetMatchedTrackSCPairs()               //
// Description:                                                                                  //
// Returns the set of track-to-supercluster pairs. See the MatchedTrackSCPair class description  //
// for more details.                                                                             //
//                                                                                               //
//-----------------------------------------------------------------------------------------------//
//===============================================================================================//
template <class T> class MatchTracksToSC {

 public:

  MatchTracksToSC(const edm::Event &iEvent, 
		  const edm::EventSetup &iSetup, 
		  const edm::ESTransientHandle<MagneticField> &magneticField,
		  const CaloGeometry &caloGeometry,
		  const TrackAssociatorParameters &parameters,
		  const std::vector<T> &tracks, 
		  const reco::SuperClusterCollection &superClusters) {

    trackAssociator_.useDefaultPropagator();

    const TrackRecHitLocations propagatedTrackLocationsAtECAL = GetPropagatedTrackAtECAL(iEvent, iSetup, magneticField, caloGeometry, parameters, tracks);

    Matrix<double> costMatrix = TrackToSuperClusterCostMatrix(propagatedTrackLocationsAtECAL, superClusters);
    HungarianAlgorithm assigner;

    assigner.Solve(costMatrix, matchedIndexes_);

    matchedPairs_ = ConstructMatchedPairs(costMatrix);

    cost_ = 0;
    for(auto const pair : matchedPairs_)
      cost_ += pair.GetMatchCriteria(); 

    FillMatchContainers(tracks, superClusters, propagatedTrackLocationsAtECAL);    

    if(superClusters.size() == 0)
      PrintMatchSummary();

  }

  // Getter Methods
  double GetCost() const { return cost_; }
  std::vector<T> GetMatchedTracks() const { return matchedTracks_; }
  reco::SuperClusterCollection GetMatchedSCs() const { return matchedSCs_; }
  std::vector<GlobalPoint> GetMatchedTrackLocationAtECAL() const { return matchedDetIdLocations_; }
  MatchedTrackSCPairs<T> GetMatchedTrackSCPairs() const { return matchedTrackSCPairs_; }

 private:

  template <typename S> 
  using Matrix = std::vector<std::vector<S> >;

  typedef std::vector<std::vector<GlobalPoint> > TrackRecHitLocations;

  double cost_;
  Matrix<int> minMatchedDetIndex_;
  std::vector<int> matchedIndexes_;
  std::vector<MatchedPair> matchedPairs_;
  std::vector<T> matchedTracks_;
  std::vector<reco::SuperCluster> matchedSCs_;
  std::vector<GlobalPoint> matchedDetIdLocations_;
  MatchedTrackSCPairs<T> matchedTrackSCPairs_;

  TrackDetectorAssociator trackAssociator_;

  std::vector<MatchedPair> GetMatchedPairs() const { return matchedPairs_; }

  //===============================================================================================//
  //-----------------------------------------------------------------------------------------------//
  // Propagate all input tracks to the ECAL by using their initial free state. The propagation is  //
  // handled by the TrackDetectorAssociator class from which one can get the full list of detector //
  // components that the propagated track was estimated to have crossed. This function returns the //
  // set of GlobalPoint objects associated to the propagated track location at the ECAL.           //
  //-----------------------------------------------------------------------------------------------//
  //===============================================================================================//
  TrackRecHitLocations GetPropagatedTrackAtECAL(const edm::Event& iEvent, 
						const edm::EventSetup& iSetup, 
						const edm::ESTransientHandle<MagneticField> &magneticField,
						const CaloGeometry &caloGeometry,
						const TrackAssociatorParameters &parameters,
						const std::vector<T> &tracks) {

    TrackRecHitLocations locationCollection;

    for(const auto &track : tracks) {

      FreeTrajectoryState initialState = trajectoryStateTransform::initialFreeState(track, magneticField.product());

      TrackDetMatchInfo trackDetInfo = trackAssociator_.associate(iEvent, iSetup, parameters, &initialState);

      std::vector<GlobalPoint> trackLocations;

      for(const auto &id : trackDetInfo.crossedEcalIds) {
        GlobalPoint trackHitAtEcal(caloGeometry.getGeometry(id)->getPosition());
        trackLocations.emplace_back(trackHitAtEcal);
      }
      locationCollection.emplace_back(trackLocations);
    }
    
    return locationCollection;
  }

  //===============================================================================================//
  //-----------------------------------------------------------------------------------------------//
  // Find the minimum deltaR value between the crossed ECAL components by the propagated track and //
  // a given super cluster. The index of the matched detID can be passed by reference.             //                                
  //-----------------------------------------------------------------------------------------------//
  //===============================================================================================//
  double GetMinDeltaRTrackSC(const std::vector<GlobalPoint> &crossedLocationECAL, const reco::SuperCluster &superCluster, int &index = -1) const {

    // Index should be undefined when no match is possible
    index = -1;

    // Initialize deltaR with a high value when unmatched 
    double deltaR(999.);
 
    const double etaSC(superCluster.eta());
    const double phiSC(superCluster.phi());
    std::vector<double> deltaRs;

    for(const auto &location : crossedLocationECAL) {
      const double etaTrack(location.eta());
      const double phiTrack(location.phi());
      
      deltaRs.push_back( sqrt(reco::deltaR2(etaSC, phiSC, etaTrack, phiTrack)) );
    }

    if(deltaRs.size() > 0) {
      auto minIt = std::min_element(deltaRs.begin(), deltaRs.end());
      deltaR = *minIt;
      index = std::distance(deltaRs.begin(), minIt);
    }

    return deltaR;
  }

  //===============================================================================================// 
  //-----------------------------------------------------------------------------------------------//
  // Private method for determining the cost matrix. The method takes two inputs: the collection   //
  // of detector component locations that the propagated track is predicted to have passed through //
  // and the collection of super clusters. Only the smallest deltaR value of the available ECAL    //
  // detIDs is used in the construction of the cost matrix.                                        //
  //-----------------------------------------------------------------------------------------------//
  //===============================================================================================// 
  Matrix<double> TrackToSuperClusterCostMatrix(const TrackRecHitLocations &trackLocationsAtECAL,
					       const reco::SuperClusterCollection &superClusters) {
    
    const int trackDimension(trackLocationsAtECAL.size());
    const int scDimension(superClusters.size());
    Matrix<double> costMatrix(trackDimension, std::vector<double>(scDimension));

    for(int ti = 0; ti < trackDimension; ti++) {
      const std::vector<GlobalPoint> trackLocations = trackLocationsAtECAL[ti];

      std::vector<int> trackToSCsIndex;
      for(int si = 0; si < scDimension; si++) {
	int detIndex;
	const reco::SuperCluster superCluster(superClusters[si]);
	costMatrix[ti][si] = GetMinDeltaRTrackSC(trackLocations, superCluster, detIndex);
	trackToSCsIndex.push_back(detIndex);
      }
      minMatchedDetIndex_.emplace_back(trackToSCsIndex);
    }
    return costMatrix;
  }

  //===============================================================================================// 
  //-----------------------------------------------------------------------------------------------// 
  // Save the indexes of matched tracks to super clusters, along with the matched DeltaR, as       //
  // MatchedPair objects. The function returns a vector of all matched pairs in the current event. //
  // The vector is sorted in order of best matched tracks to super clusters (smallest deltaR).     //
  //-----------------------------------------------------------------------------------------------// 
  //===============================================================================================// 
  std::vector<MatchedPair> ConstructMatchedPairs(const Matrix<double> &costMatrix) const {

    std::vector<MatchedPair> matchedPairs;
    for(unsigned int mi = 0; mi < costMatrix.size(); mi++) {
      if(matchedIndexes_[mi] == -1)
        continue;

      const double deltaR = costMatrix[mi][matchedIndexes_[mi]];
      matchedPairs.emplace_back(MatchedPair(mi, matchedIndexes_[mi], deltaR));
    }

    // Sort matched pairs in ascending order (smallest deltaR values first).
    std::sort(matchedPairs.begin(), matchedPairs.end(), [](const MatchedPair& a, const MatchedPair& b) {
	return a.GetMatchCriteria() < b.GetMatchCriteria(); 
      });

    return matchedPairs;
  }

  //===============================================================================================// 
  //-----------------------------------------------------------------------------------------------// 
  // Fill the vectors of matched tracks, super clusters and the propagated track location at the   //
  // ECAL. The order of the matches is the same for all containers (meaning that for any arbitrary //
  // index i all of the objects of every container at i are related to each other).                //
  //-----------------------------------------------------------------------------------------------//
  //===============================================================================================//  
  void FillMatchContainers(const std::vector<T> &tracks, const reco::SuperClusterCollection &superClusters, const TrackRecHitLocations &locations) {

    matchedTracks_.clear();
    matchedSCs_.clear();
    matchedDetIdLocations_.clear();
    matchedTrackSCPairs_ = MatchedTrackSCPairs<T>(cost_);

    for(auto const &matches : matchedPairs_) {
      int matchedDetIndex = minMatchedDetIndex_[matches.GetIndexA()][matches.GetIndexB()];
      
      // If the algorithm had no choice but to match a track that could not be propagated to
      // the ECAL, ignore it  
      if(matchedDetIndex < 0)
	continue;

      matchedTracks_.emplace_back(tracks[matches.GetIndexA()]);
      matchedSCs_.emplace_back(superClusters[matches.GetIndexB()]);
      matchedDetIdLocations_.emplace_back(locations[matches.GetIndexA()][matchedDetIndex]);

      matchedTrackSCPairs_.emplace_back(MatchedTrackSCPair<T>(tracks[matches.GetIndexA()], 
							      superClusters[matches.GetIndexB()],
							      matches.GetMatchCriteria(),
							      locations[matches.GetIndexA()][matchedDetIndex]));
    }
  }

  void PrintMatchSummary() {

    std::cout<< "There were " << matchedTrackSCPairs_.size() << " matches in this event:" << std::endl;
    
    for(auto const &match : matchedTrackSCPairs_) {
       
      std::cout << "  DeltaR = " << match.GetDeltaR() << std::endl;
      std::cout << "  Track Pt = " << match.GetTrack().pt() << std::endl;
      std::cout << "  SuperCluster Energy = " << match.GetSuperCluster().rawEnergy() << std::endl;
      
      GlobalPoint locationAtECAL = match.GetTrackLocationAtECAL();
      std::cout << "  Location at ECAL: (" << locationAtECAL.x() << ", " << locationAtECAL.y() << ", " 
		<< locationAtECAL.z() << ") " << std::endl;
      std::cout << std::endl;
    }
  }
  
};

#endif
