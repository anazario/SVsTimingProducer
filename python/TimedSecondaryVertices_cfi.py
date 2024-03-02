import FWCore.ParameterSet.Config as cms

timedSecondaryVertices = cms.EDProducer('SVsTimingProducer',
                                        electronSrc = cms.InputTag("gedGsfElectrons"),#"lowPtGsfElectrons"),
                                        timedVertices = cms.InputTag("timedVertices"),
                                        vertexLabel = cms.string("timedVertices"),
)

