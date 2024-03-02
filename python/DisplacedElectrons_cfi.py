import FWCore.ParameterSet.Config as cms

from TrackingTools.TrackAssociator.default_cfi import TrackAssociatorParameterBlock

tkAssocParamBlock = TrackAssociatorParameterBlock.clone()
tkAssocParamBlock.TrackAssociatorParameters.useMuon = cms.bool(False)
tkAssocParamBlock.TrackAssociatorParameters.useCalo = cms.bool(False)
tkAssocParamBlock.TrackAssociatorParameters.useHO = cms.bool(False)
tkAssocParamBlock.TrackAssociatorParameters.usePreshower = cms.bool(False)
tkAssocParamBlock.TrackAssociatorParameters.EERecHitCollectionLabel = cms.InputTag("reducedEcalRecHitsEE")
tkAssocParamBlock.TrackAssociatorParameters.EBRecHitCollectionLabel = cms.InputTag("reducedEcalRecHitsEB")
tkAssocParamBlock.TrackAssociatorParameters.HBHERecHitCollectionLabel = cms.InputTag("reducedHcalRecHits","hbhereco")
tkAssocParamBlock.TrackAssociatorParameters.HORecHitCollectionLabel = cms.InputTag("reducedHcalRecHits","horeco")

displacedElectrons = cms.EDProducer("DisplacedElectronProducer",
                                    tkAssocParamBlock,
                                    generalTracksSrc = cms.InputTag("generalTracks"),
                                    gsfElectronTracksSrc = cms.InputTag("electronGsfTracks"),
                                    #lowPtGsfElectronTracksSrc = cms.InputTag("lowPtGsfEleGsfTracks"),
                                    #pfCollectionSrc = cms.InputTag("particleFlow"),
                                    #electronSrc = cms.InputTag("gedGsfElectrons"),
                                    displacedTracksSrc = cms.InputTag("displacedTracks"),
                                    #trackExtrapolatorSrc = cms.InputTag("trackExtrapolator"),
                                    genParticleSrc = cms.InputTag("genParticles"),
                                    superClusters = cms.InputTag("particleFlowEGamma"),
                                    ootSuperClusters = cms.InputTag("particleFlowSuperClusterOOTECAL", "particleFlowSuperClusterOOTECALBarrel"),
                                    displacedElectrons = cms.InputTag("displacedElectrons")
)

