import FWCore.ParameterSet.Config as cms

lostTracks = cms.EDProducer("LostTrackProducer",
                            inputCandidates = cms.InputTag("particleFlow"),
                            packedPFCandidates  = cms.InputTag("packedPFCandidates"),
                            inputTracks = cms.InputTag("generalTracks"),
                            secondaryVertices = cms.InputTag("inclusiveSecondaryVertices"),
                            minPt = cms.double(0.95),
                            minHits = cms.uint32(8),
                            minPixelHits = cms.uint32(1),
                            qualsToAutoAccept = cms.vstring("highPurity"),
)

