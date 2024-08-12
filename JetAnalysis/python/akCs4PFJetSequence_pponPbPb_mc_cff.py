import FWCore.ParameterSet.Config as cms

from HeavyIonsAnalysis.JetAnalysis.inclusiveJetAnalyzer_cff import *

akCs4PFJetAnalyzer = inclusiveJetAnalyzer.clone(
    jetTag = cms.InputTag("slimmedJets"),
    genjetTag = 'slimmedGenJets',
    rParam = 0.4,
    fillGenJets = True,
    isMC = True,
    genParticles = cms.untracked.InputTag("prunedGenParticles"),
    eventInfoTag = cms.InputTag("generator"),
    jetName = cms.untracked.string("akCs4PF"),
    genPtMin = cms.untracked.double(5),
    hltTrgResults = cms.untracked.string('TriggerResults::'+'HISIGNAL'),
    )



ak4PFJetAnalyzer = inclusiveJetAnalyzer.clone(
    jetTag = cms.InputTag("slimmedJets"),
    genjetTag = 'slimmedGenJets',
    rParam = 0.4,
    fillGenJets = True,
    isMC = True,
    genParticles = cms.untracked.InputTag("prunedGenParticles"),
    eventInfoTag = cms.InputTag("generator"),
    jetName = cms.untracked.string("ak4PF"),
    genPtMin = cms.untracked.double(5),
    hltTrgResults = cms.untracked.string('TriggerResults::'+'HISIGNAL'),
    )

ak2PFJetAnalyzer = inclusiveJetAnalyzer.clone(
    jetTag = cms.InputTag("slimmedJets"),
    genjetTag = 'slimmedGenJets',
    rParam = 0.2,
    fillGenJets = True,
    isMC = True,
    #doSubEvent = True,
    #useHepMC = cms.untracked.bool(False),
    genParticles = cms.untracked.InputTag("prunedGenParticles"),
    eventInfoTag = cms.InputTag("generator"),
    jetName = cms.untracked.string("ak2PF"),
    genPtMin = cms.untracked.double(5),
    hltTrgResults = cms.untracked.string('TriggerResults::'+'HISIGNAL'),
    #genDroppedBranches = cms.InputTag("ak2GenJets","droppedBranches")
    )
 

ak6PFJetAnalyzer = inclusiveJetAnalyzer.clone(
    jetTag = cms.InputTag("slimmedJets"),
    genjetTag = 'slimmedGenJets',
    rParam = 0.6,
    fillGenJets = True,
    isMC = True,
    #doSubEvent = True,
    #useHepMC = cms.untracked.bool(False),
    genParticles = cms.untracked.InputTag("prunedGenParticles"),
    eventInfoTag = cms.InputTag("generator"),
    jetName = cms.untracked.string("ak6PF"),
    genPtMin = cms.untracked.double(5),
    hltTrgResults = cms.untracked.string('TriggerResults::'+'HISIGNAL'),
    #genDroppedBranches = cms.InputTag("ak6GenJets","droppedBranches")
    )

ak8PFJetAnalyzer = inclusiveJetAnalyzer.clone(
    jetTag = cms.InputTag("slimmedJets"),
    genjetTag = 'slimmedGenJets',
    rParam = 0.8,
    fillGenJets = True,
    isMC = True,
    #doSubEvent = True,
    #useHepMC = cms.untracked.bool(False),
    genParticles = cms.untracked.InputTag("prunedGenParticles"),
    eventInfoTag = cms.InputTag("generator"),
    jetName = cms.untracked.string("ak8PF"),
    genPtMin = cms.untracked.double(5),
    hltTrgResults = cms.untracked.string('TriggerResults::'+'HISIGNAL'),
    #genDroppedBranches = cms.InputTag("ak8GenJets","droppedBranches")
    ) 
