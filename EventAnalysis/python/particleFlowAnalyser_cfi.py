import FWCore.ParameterSet.Config as cms

particleFlowAnalyser = cms.EDAnalyzer(
    "ParticleFlowAnalyser",
    pfCandidateSrc = cms.InputTag('packedPFCandidates'),
    ptMin = cms.double(0.01),
    absEtaMax = cms.double(5.),
    )
