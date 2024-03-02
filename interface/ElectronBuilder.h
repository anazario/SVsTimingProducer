#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

#include "TimingWithSVs/SVsTimingProducer/interface/Hungarian.h"
#include "TimingWithSVs/SVsTimingProducer/interface/DeltaRMatch.h"

template <class T> class ElectronBuilder {

 public:

  ElectronBuilder(const std::vector<T> &tracks, const SuperClusterCollection &superClusters) {
    
  }

  

 private:

  double cost_;
  std::vector<int> matchedIndexes_;
  std::vector<MatchedPair> matchedPairs_;

  std::vector<std::vector<double>> CalculateCostMatrix(const std::vector<A> &objectsA, const std::vector<B> &objectsB) const {

    const int dimensionA = objectsA.size();
    const int dimensionB = objectsB.size();
    std::vector<std::vector<double>> costMatrix(dimensionA, std::vector<double>(dimensionB) );

    for(int ai = 0; ai < dimensionA; ai++) {

      const double etaA = objectsA[ai].eta();
      const double phiA = objectsA[ai].phi();

      for(int bi = 0; bi < dimensionB; bi++) {

        const double etaB = objectsB[bi].eta();
        const double phiB = objectsB[bi].phi();

        costMatrix[ai][bi] = sqrt(reco::deltaR2(etaA, phiA, etaB, phiB));

      }// End B loop
    }// End A loop

    return costMatrix;
  }
  
}
