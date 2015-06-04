import FWCore.ParameterSet.Config as cms

# JETS  CA8 ----------------------------

from RecoJets.JetProducers.ak5PFJets_cfi import ak5PFJets
ca8PFJetsCHS = ak5PFJets.clone(
    src = 'pfNoPileUp',
    jetPtMin = cms.double(25.0),
    doAreaFastjet = cms.bool(True),
    rParam = cms.double(0.8),
    jetAlgorithm = cms.string("CambridgeAachen"),
    )

jetSource = 'ca8PFJetsCHS'

# corrections 
from PhysicsTools.PatAlgos.recoLayer0.jetCorrFactors_cfi import *
patJetCorrFactorsCA8CHS = patJetCorrFactors.clone()
patJetCorrFactorsCA8CHS.src = jetSource
# will need to add L2L3 corrections in the cfg
patJetCorrFactorsCA8CHS.levels = ['L1FastJet', 'L2Relative', 'L3Absolute']
patJetCorrFactorsCA8CHS.payload = 'AK7PFchs'
patJetCorrFactorsCA8CHS.useRho = True

# parton and gen jet matching

from PhysicsTools.PatAlgos.mcMatchLayer0.jetMatch_cfi import *
patJetPartonMatchCA8CHS = patJetPartonMatch.clone()
patJetPartonMatchCA8CHS.src = jetSource
patJetGenJetMatchCA8CHS = patJetGenJetMatch.clone()
patJetGenJetMatchCA8CHS.src = jetSource
patJetGenJetMatchCA8CHS.matched = 'ca8GenJetsNoNu'

from PhysicsTools.PatAlgos.mcMatchLayer0.jetFlavourId_cff import *
patJetPartonAssociationCA8CHS = patJetPartonAssociation.clone()
patJetPartonAssociationCA8CHS.jets = jetSource

# pat jets

from RecoJets.JetAssociationProducers.ak5JTA_cff import *
ca8CHSJetTracksAssociatorAtVertex=ak5JetTracksAssociatorAtVertex.clone()
ca8CHSJetTracksAssociatorAtVertex.jets=jetSource
from RecoBTag.Configuration.RecoBTag_cff import * # btagging sequence
# basic b-tagging
impactParameterTagInfosCA8CHS=impactParameterTagInfos.clone()
impactParameterTagInfosCA8CHS.jetTracks='ca8CHSJetTracksAssociatorAtVertex'
secondaryVertexTagInfosCA8CHS=secondaryVertexTagInfos.clone()
secondaryVertexTagInfosCA8CHS.trackIPTagInfos='impactParameterTagInfosCA8CHS'
## secondary vertex based taggers
# combined secondary vertex tagger
combinedSecondaryVertexBJetTagsCA8CHS=combinedSecondaryVertexBJetTags.clone()
combinedSecondaryVertexBJetTagsCA8CHS.tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosCA8CHS"),
    cms.InputTag("secondaryVertexTagInfosCA8CHS"))
# simple SV taggers
simpleSecondaryVertexHighEffBJetTagsCA8CHS = simpleSecondaryVertexHighEffBJetTags.clone(
  tagInfos = cms.VInputTag(cms.InputTag("secondaryVertexTagInfosCA8CHS"))
)
simpleSecondaryVertexHighPurBJetTagsCA8CHS = simpleSecondaryVertexHighPurBJetTags.clone(
  tagInfos = cms.VInputTag(cms.InputTag("secondaryVertexTagInfosCA8CHS"))
)
## impact parameter based taggers
# jet probability taggers
jetProbabilityBJetTagsCA8CHS = jetProbabilityBJetTags.clone(
  tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosCA8CHS"))
)
jetBProbabilityBJetTagsCA8CHS = jetBProbabilityBJetTags.clone(
  tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosCA8CHS"))
)
# track counting taggers
trackCountingHighEffBJetTagsCA8CHS = trackCountingHighEffBJetTags.clone(
  tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosCA8CHS"))
)
trackCountingHighPurBJetTagsCA8CHS = trackCountingHighPurBJetTags.clone(
  tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosCA8CHS"))
)
# b-tagging sequence
btaggingCA8CHS=cms.Sequence(ca8CHSJetTracksAssociatorAtVertex+
                            impactParameterTagInfosCA8CHS+
                            secondaryVertexTagInfosCA8CHS+
                            combinedSecondaryVertexBJetTagsCA8CHS+
                            simpleSecondaryVertexHighEffBJetTagsCA8CHS+
                            simpleSecondaryVertexHighPurBJetTagsCA8CHS+
                            jetProbabilityBJetTagsCA8CHS+
                            jetBProbabilityBJetTagsCA8CHS+
                            trackCountingHighEffBJetTagsCA8CHS+
                            trackCountingHighPurBJetTagsCA8CHS)

from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cfi import *
patJetsCA8CHS = patJets.clone()
patJetsCA8CHS.jetSource = jetSource
patJetsCA8CHS.addJetCharge = False
patJetsCA8CHS.embedCaloTowers = False
patJetsCA8CHS.embedPFCandidates = False
patJetsCA8CHS.addAssociatedTracks = False
patJetsCA8CHS.addBTagInfo = True
patJetsCA8CHS.addDiscriminators = True
patJetsCA8CHS.addJetID = False
patJetsCA8CHS.tagInfoSources = cms.VInputTag(cms.InputTag("secondaryVertexTagInfosCA8CHS"))
patJetsCA8CHS.trackAssociationSource = cms.InputTag("ca8CHSJetTracksAssociatorAtVertex")
patJetsCA8CHS.discriminatorSources = cms.VInputTag(cms.InputTag("combinedSecondaryVertexBJetTagsCA8CHS"),
                                                   cms.InputTag("simpleSecondaryVertexHighEffBJetTagsCA8CHS"),
                                                   cms.InputTag("simpleSecondaryVertexHighPurBJetTagsCA8CHS"),
                                                   cms.InputTag("jetProbabilityBJetTagsCA8CHS"),
                                                   cms.InputTag("jetBProbabilityBJetTagsCA8CHS"),
                                                   cms.InputTag("trackCountingHighEffBJetTagsCA8CHS"),
                                                   cms.InputTag("trackCountingHighPurBJetTagsCA8CHS"))
patJetsCA8CHS.getJetMCFlavour = False
patJetsCA8CHS.jetCorrFactorsSource = cms.VInputTag(cms.InputTag('patJetCorrFactorsCA8CHS'))
patJetsCA8CHS.genPartonMatch = cms.InputTag('patJetPartonMatchCA8CHS')
patJetsCA8CHS.genJetMatch = cms.InputTag('patJetGenJetMatchCA8CHS')

from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import *
selectedPatJetsCA8CHS = selectedPatJets.clone()
selectedPatJetsCA8CHS.src = 'patJetsCA8CHS'
selectedPatJetsCA8CHS.cut = 'pt()>30'

from RecoJets.Configuration.RecoGenJets_cff import ak7GenJetsNoNu
ca8GenJetsNoNu = ak7GenJetsNoNu.clone()
ca8GenJetsNoNu.rParam = 0.8
ca8GenJetsNoNu.jetAlgorithm = "CambridgeAachen"
ca8GenJetsNoNu.jetPtMin = 25
ca8GenJetsNoNu.doAreaFastjet = True


from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cfi import *
patGenJetsCA8CHS = patJets.clone()
patGenJetsCA8CHS.jetSource = 'ca8GenJetsNoNu'
patGenJetsCA8CHS.addGenJetMatch = False
patGenJetsCA8CHS.addGenPartonMatch = False
patGenJetsCA8CHS.addJetCharge = False
patGenJetsCA8CHS.embedCaloTowers = False
patGenJetsCA8CHS.embedPFCandidates = False
patGenJetsCA8CHS.addAssociatedTracks = False
patGenJetsCA8CHS.addBTagInfo = False
patGenJetsCA8CHS.addDiscriminators = False
patGenJetsCA8CHS.addJetID = False
patGenJetsCA8CHS.getJetMCFlavour = False
patGenJetsCA8CHS.addJetCorrFactors = False

from RecoJets.Configuration.GenJetParticles_cff import genParticlesForJetsNoNu
jetMCSequenceCA8CHS = cms.Sequence(
    patJetPartonMatchCA8CHS +
    genParticlesForJetsNoNu +
    ca8GenJetsNoNu +
    patGenJetsCA8CHS +
    patJetGenJetMatchCA8CHS
    )

PATCMGJetSequenceCA8CHS = cms.Sequence(
    ca8PFJetsCHS +
    jetMCSequenceCA8CHS +
    patJetCorrFactorsCA8CHS +
    btaggingCA8CHS + # needs generalTracks
    patJetsCA8CHS +
    selectedPatJetsCA8CHS
    )

# JETS PRUNED CA8 ----------------------------

from RecoJets.JetProducers.ak5PFJetsPruned_cfi import ak5PFJetsPruned
ca8PFJetsCHSpruned = ak5PFJetsPruned.clone(
    src = 'pfNoPileUp',
    jetPtMin = cms.double(25.0),
    doAreaFastjet = cms.bool(True),
    rParam = cms.double(0.8),
    jetAlgorithm = cms.string("CambridgeAachen"),
    )

jetSource = 'ca8PFJetsCHSpruned'

# corrections 
from PhysicsTools.PatAlgos.recoLayer0.jetCorrFactors_cfi import *
patJetCorrFactorsCA8CHSpruned = patJetCorrFactors.clone()
patJetCorrFactorsCA8CHSpruned.src = jetSource
# will need to add L2L3 corrections in the cfg
patJetCorrFactorsCA8CHSpruned.levels = ['L1FastJet', 'L2Relative', 'L3Absolute']
patJetCorrFactorsCA8CHSpruned.payload = 'AK7PFchs'
patJetCorrFactorsCA8CHSpruned.useRho = True

# parton and gen jet matching

from PhysicsTools.PatAlgos.mcMatchLayer0.jetMatch_cfi import *
patJetPartonMatchCA8CHSpruned = patJetPartonMatch.clone()
patJetPartonMatchCA8CHSpruned.src = jetSource
patJetGenJetMatchCA8CHSpruned = patJetGenJetMatch.clone()
patJetGenJetMatchCA8CHSpruned.src = jetSource
patJetGenJetMatchCA8CHSpruned.matched = 'ca8GenJetsNoNu'

from PhysicsTools.PatAlgos.mcMatchLayer0.jetFlavourId_cff import *
patJetPartonAssociationCA8CHSpruned = patJetPartonAssociation.clone()
patJetPartonAssociationCA8CHSpruned.jets = jetSource

# pat jets

from RecoJets.JetAssociationProducers.ak5JTA_cff import *
ca8CHSprunedJetTracksAssociatorAtVertex=ak5JetTracksAssociatorAtVertex.clone()
ca8CHSprunedJetTracksAssociatorAtVertex.jets=jetSource
from RecoBTag.Configuration.RecoBTag_cff import * # btagging sequence
# basic b-tagging
impactParameterTagInfosCA8CHSpruned=impactParameterTagInfos.clone()
impactParameterTagInfosCA8CHSpruned.jetTracks='ca8CHSprunedJetTracksAssociatorAtVertex'
secondaryVertexTagInfosCA8CHSpruned=secondaryVertexTagInfos.clone()
secondaryVertexTagInfosCA8CHSpruned.trackIPTagInfos='impactParameterTagInfosCA8CHSpruned'
## secondary vertex based taggers
# combined secondary vertex tagger
combinedSecondaryVertexBJetTagsCA8CHSpruned=combinedSecondaryVertexBJetTags.clone()
combinedSecondaryVertexBJetTagsCA8CHSpruned.tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosCA8CHSpruned"),
    cms.InputTag("secondaryVertexTagInfosCA8CHSpruned"))
# simple SV taggers
simpleSecondaryVertexHighEffBJetTagsCA8CHSpruned = simpleSecondaryVertexHighEffBJetTags.clone(
  tagInfos = cms.VInputTag(cms.InputTag("secondaryVertexTagInfosCA8CHSpruned"))
)
simpleSecondaryVertexHighPurBJetTagsCA8CHSpruned = simpleSecondaryVertexHighPurBJetTags.clone(
  tagInfos = cms.VInputTag(cms.InputTag("secondaryVertexTagInfosCA8CHSpruned"))
)
## impact parameter based taggers
# jet probability taggers
jetProbabilityBJetTagsCA8CHSpruned = jetProbabilityBJetTags.clone(
  tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosCA8CHSpruned"))
)
jetBProbabilityBJetTagsCA8CHSpruned = jetBProbabilityBJetTags.clone(
  tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosCA8CHSpruned"))
)
# track counting taggers
trackCountingHighEffBJetTagsCA8CHSpruned = trackCountingHighEffBJetTags.clone(
  tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosCA8CHSpruned"))
)
trackCountingHighPurBJetTagsCA8CHSpruned = trackCountingHighPurBJetTags.clone(
  tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosCA8CHSpruned"))
)
# b-tagging sequence
btaggingCA8CHSpruned=cms.Sequence(ca8CHSprunedJetTracksAssociatorAtVertex+
                            impactParameterTagInfosCA8CHSpruned+
                            secondaryVertexTagInfosCA8CHSpruned+
                            combinedSecondaryVertexBJetTagsCA8CHSpruned+
                            simpleSecondaryVertexHighEffBJetTagsCA8CHSpruned+
                            simpleSecondaryVertexHighPurBJetTagsCA8CHSpruned+
                            jetProbabilityBJetTagsCA8CHSpruned+
                            jetBProbabilityBJetTagsCA8CHSpruned+
                            trackCountingHighEffBJetTagsCA8CHSpruned+
                            trackCountingHighPurBJetTagsCA8CHSpruned)

patJetsCA8CHSpruned = patJets.clone()
patJetsCA8CHSpruned.jetSource = jetSource
patJetsCA8CHSpruned.addJetCharge = False
patJetsCA8CHSpruned.embedCaloTowers = False
patJetsCA8CHSpruned.embedPFCandidates = False
patJetsCA8CHSpruned.addAssociatedTracks = False
patJetsCA8CHSpruned.addBTagInfo = True
patJetsCA8CHSpruned.addDiscriminators = True
patJetsCA8CHSpruned.addJetID = False
patJetsCA8CHSpruned.tagInfoSources = cms.VInputTag(cms.InputTag("secondaryVertexTagInfosCA8CHSpruned"))
patJetsCA8CHSpruned.trackAssociationSource = cms.InputTag("ca8CHSprunedJetTracksAssociatorAtVertex")
patJetsCA8CHSpruned.discriminatorSources = cms.VInputTag(cms.InputTag("combinedSecondaryVertexBJetTagsCA8CHSpruned"),
                                                   cms.InputTag("simpleSecondaryVertexHighEffBJetTagsCA8CHSpruned"),
                                                   cms.InputTag("simpleSecondaryVertexHighPurBJetTagsCA8CHSpruned"),
                                                   cms.InputTag("jetProbabilityBJetTagsCA8CHSpruned"),
                                                   cms.InputTag("jetBProbabilityBJetTagsCA8CHSpruned"),
                                                   cms.InputTag("trackCountingHighEffBJetTagsCA8CHSpruned"),
                                                   cms.InputTag("trackCountingHighPurBJetTagsCA8CHSpruned"))
patJetsCA8CHSpruned.getJetMCFlavour = False
patJetsCA8CHSpruned.jetCorrFactorsSource = cms.VInputTag(cms.InputTag('patJetCorrFactorsCA8CHSpruned'))
patJetsCA8CHSpruned.genPartonMatch = cms.InputTag('patJetPartonMatchCA8CHSpruned')
patJetsCA8CHSpruned.genJetMatch = cms.InputTag('patJetGenJetMatchCA8CHSpruned')

selectedPatJetsCA8CHSprunedPre = selectedPatJets.clone()
selectedPatJetsCA8CHSprunedPre.src = 'patJetsCA8CHSpruned'
selectedPatJetsCA8CHSprunedPre.cut = 'pt()>30'

ca8PrunedGenJetsNoNu = ak7GenJetsNoNu.clone()
ca8PrunedGenJetsNoNu.rParam = 0.8
ca8PrunedGenJetsNoNu.jetAlgorithm = "CambridgeAachen"
ca8PrunedGenJetsNoNu.doAreaFastjet = True
ca8PrunedGenJetsNoNu.usePruning = cms.bool(True)
ca8PrunedGenJetsNoNu.useExplicitGhosts = cms.bool(True)
ca8PrunedGenJetsNoNu.writeCompound = cms.bool(True)
ca8PrunedGenJetsNoNu.jetCollInstanceName = cms.string("SubJets")
ca8PrunedGenJetsNoNu.nFilt = cms.int32(2)
ca8PrunedGenJetsNoNu.zcut = cms.double(0.1)
ca8PrunedGenJetsNoNu.rcut_factor = cms.double(0.5)
ca8PrunedGenJetsNoNu.jetPtMin = 25

patGenJetsCA8CHSpruned = patJets.clone()
patGenJetsCA8CHSpruned.jetSource = 'ca8PrunedGenJetsNoNu'
patGenJetsCA8CHSpruned.addGenJetMatch = False
patGenJetsCA8CHSpruned.addGenPartonMatch = False
patGenJetsCA8CHSpruned.addJetCharge = False
patGenJetsCA8CHSpruned.embedCaloTowers = False
patGenJetsCA8CHSpruned.embedPFCandidates = False
patGenJetsCA8CHSpruned.addAssociatedTracks = False
patGenJetsCA8CHSpruned.addBTagInfo = False
patGenJetsCA8CHSpruned.addDiscriminators = False
patGenJetsCA8CHSpruned.addJetID = False
patGenJetsCA8CHSpruned.getJetMCFlavour = False
patGenJetsCA8CHSpruned.addJetCorrFactors = False

#### Adding subjet b-tagging

from RecoJets.JetAssociationProducers.ak5JTA_cff import *
ca8CHSprunedSubjetsJetTracksAssociatorAtVertex=ak5JetTracksAssociatorAtVertex.clone()
ca8CHSprunedSubjetsJetTracksAssociatorAtVertex.jets=cms.InputTag('ca8PFJetsCHSpruned','SubJets')
from RecoBTag.Configuration.RecoBTag_cff import * # btagging sequence
impactParameterTagInfosCA8CHSprunedSubjets=impactParameterTagInfos.clone()
impactParameterTagInfosCA8CHSprunedSubjets.jetTracks='ca8CHSprunedSubjetsJetTracksAssociatorAtVertex'
secondaryVertexTagInfosCA8CHSprunedSubjets=secondaryVertexTagInfos.clone()
secondaryVertexTagInfosCA8CHSprunedSubjets.trackIPTagInfos='impactParameterTagInfosCA8CHSprunedSubjets'
combinedSecondaryVertexBJetTagsCA8CHSprunedSubjets=combinedSecondaryVertexBJetTags.clone()
combinedSecondaryVertexBJetTagsCA8CHSprunedSubjets.tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosCA8CHSprunedSubjets"),
    cms.InputTag("secondaryVertexTagInfosCA8CHSprunedSubjets"))

# simple SV taggers
simpleSecondaryVertexHighEffBJetTagsCA8CHSprunedSubjets = simpleSecondaryVertexHighEffBJetTags.clone(
  tagInfos = cms.VInputTag(cms.InputTag("secondaryVertexTagInfosCA8CHSprunedSubjets"))
)
simpleSecondaryVertexHighPurBJetTagsCA8CHSprunedSubjets = simpleSecondaryVertexHighPurBJetTags.clone(
  tagInfos = cms.VInputTag(cms.InputTag("secondaryVertexTagInfosCA8CHSprunedSubjets"))
)
## impact parameter based taggers
# jet probability taggers
jetProbabilityBJetTagsCA8CHSprunedSubjets = jetProbabilityBJetTags.clone(
  tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosCA8CHSprunedSubjets"))
)
jetBProbabilityBJetTagsCA8CHSprunedSubjets = jetBProbabilityBJetTags.clone(
  tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosCA8CHSprunedSubjets"))
)
# track counting taggers
trackCountingHighEffBJetTagsCA8CHSprunedSubjets = trackCountingHighEffBJetTags.clone(
  tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosCA8CHSprunedSubjets"))
)
trackCountingHighPurBJetTagsCA8CHSprunedSubjets = trackCountingHighPurBJetTags.clone(
  tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosCA8CHSprunedSubjets"))
)
# b-tagging sequence
btaggingCA8CHSprunedSubjets=cms.Sequence(ca8CHSprunedSubjetsJetTracksAssociatorAtVertex+
                            impactParameterTagInfosCA8CHSprunedSubjets+
                            secondaryVertexTagInfosCA8CHSprunedSubjets+
                            combinedSecondaryVertexBJetTagsCA8CHSprunedSubjets+
                            simpleSecondaryVertexHighEffBJetTagsCA8CHSprunedSubjets+
                            simpleSecondaryVertexHighPurBJetTagsCA8CHSprunedSubjets+
                            jetProbabilityBJetTagsCA8CHSprunedSubjets+
                            jetBProbabilityBJetTagsCA8CHSprunedSubjets+
                            trackCountingHighEffBJetTagsCA8CHSprunedSubjets+
                            trackCountingHighPurBJetTagsCA8CHSprunedSubjets)

patJetsCA8CHSprunedSubjets = patJets.clone()
patJetsCA8CHSprunedSubjets.jetSource = cms.InputTag('ca8PFJetsCHSpruned','SubJets')
patJetsCA8CHSprunedSubjets.addGenJetMatch = False
patJetsCA8CHSprunedSubjets.addGenPartonMatch = False
patJetsCA8CHSprunedSubjets.addJetCharge = False
patJetsCA8CHSprunedSubjets.embedCaloTowers = False
patJetsCA8CHSprunedSubjets.embedPFCandidates = False
patJetsCA8CHSprunedSubjets.addAssociatedTracks = False
patJetsCA8CHSprunedSubjets.addBTagInfo = True
patJetsCA8CHSprunedSubjets.addDiscriminators = True
patJetsCA8CHSprunedSubjets.addJetID = False
patJetsCA8CHSprunedSubjets.tagInfoSources = cms.VInputTag(cms.InputTag("secondaryVertexTagInfosCA8CHSprunedSubjets"))
patJetsCA8CHSprunedSubjets.trackAssociationSource = cms.InputTag("ca8CHSprunedSubjetsJetTracksAssociatorAtVertex")
patJetsCA8CHSprunedSubjets.discriminatorSources = cms.VInputTag(cms.InputTag("combinedSecondaryVertexBJetTagsCA8CHSprunedSubjets"))
patJetsCA8CHSprunedSubjets.getJetMCFlavour = False
patJetsCA8CHSprunedSubjets.addJetCorrFactors = False

selectedPatJetsCA8CHSpruned = cms.EDProducer("BoostedJetMerger",
                                                      jetSrc=cms.InputTag("selectedPatJetsCA8CHSprunedPre"),
                                                      subjetSrc=cms.InputTag("patJetsCA8CHSprunedSubjets")
    )

jetMCSequenceCA8CHSpruned = cms.Sequence(
    patJetPartonMatchCA8CHSpruned +
    genParticlesForJetsNoNu +
    ca8PrunedGenJetsNoNu +
    patGenJetsCA8CHSpruned +
    patJetGenJetMatchCA8CHSpruned
    )

PATCMGJetSequenceCA8CHSpruned = cms.Sequence(
    ca8PFJetsCHSpruned +
    jetMCSequenceCA8CHSpruned +
    patJetCorrFactorsCA8CHSpruned +
    btaggingCA8CHSpruned +
    patJetsCA8CHSpruned +
    selectedPatJetsCA8CHSprunedPre +
    btaggingCA8CHSprunedSubjets +
    patJetsCA8CHSprunedSubjets +
    selectedPatJetsCA8CHSpruned
    )


# JETS TRIMMED CA8 ----------------------------

from RecoJets.JetProducers.ak5PFJetsTrimmed_cfi import ak5PFJetsTrimmed
ca8PFJetsCHStrimmed = ak5PFJetsTrimmed.clone(
    src = 'pfNoPileUp',
    jetPtMin = cms.double(25.0),
    doAreaFastjet = cms.bool(True),
    rParam = cms.double(0.8),
    jetAlgorithm = cms.string("CambridgeAachen"),
    rFilt = cms.double(0.1),
    trimPtFracMin = cms.double(0.03)
    )

jetSource = 'ca8PFJetsCHStrimmed'

# corrections 
from PhysicsTools.PatAlgos.recoLayer0.jetCorrFactors_cfi import *
patJetCorrFactorsCA8CHStrimmed = patJetCorrFactors.clone()
patJetCorrFactorsCA8CHStrimmed.src = jetSource
# will need to add L2L3 corrections in the cfg
patJetCorrFactorsCA8CHStrimmed.levels = ['L1FastJet', 'L2Relative', 'L3Absolute']
patJetCorrFactorsCA8CHStrimmed.payload = 'AK7PFchs'
patJetCorrFactorsCA8CHStrimmed.useRho = True

# parton and gen jet matching

from PhysicsTools.PatAlgos.mcMatchLayer0.jetMatch_cfi import *
patJetPartonMatchCA8CHStrimmed = patJetPartonMatch.clone()
patJetPartonMatchCA8CHStrimmed.src = jetSource
patJetGenJetMatchCA8CHStrimmed = patJetGenJetMatch.clone()
patJetGenJetMatchCA8CHStrimmed.src = jetSource
patJetGenJetMatchCA8CHStrimmed.matched = 'ca8GenJetsNoNu'

from PhysicsTools.PatAlgos.mcMatchLayer0.jetFlavourId_cff import *
patJetPartonAssociationCA8CHStrimmed = patJetPartonAssociation.clone()
patJetPartonAssociationCA8CHStrimmed.jets = jetSource

# pat jets

from RecoJets.JetAssociationProducers.ak5JTA_cff import *
ca8CHStrimmedJetTracksAssociatorAtVertex=ak5JetTracksAssociatorAtVertex.clone()
ca8CHStrimmedJetTracksAssociatorAtVertex.jets=jetSource
from RecoBTag.Configuration.RecoBTag_cff import * # btagging sequence
impactParameterTagInfosCA8CHStrimmed=impactParameterTagInfos.clone()
impactParameterTagInfosCA8CHStrimmed.jetTracks='ca8CHStrimmedJetTracksAssociatorAtVertex'
secondaryVertexTagInfosCA8CHStrimmed=secondaryVertexTagInfos.clone()
secondaryVertexTagInfosCA8CHStrimmed.trackIPTagInfos='impactParameterTagInfosCA8CHStrimmed'
combinedSecondaryVertexBJetTagsCA8CHStrimmed=combinedSecondaryVertexBJetTags.clone()
combinedSecondaryVertexBJetTagsCA8CHStrimmed.tagInfos = cms.VInputTag(cms.InputTag("impactParameterTagInfosCA8CHStrimmed"),
    cms.InputTag("secondaryVertexTagInfosCA8CHStrimmed"))
btaggingCA8CHStrimmed=cms.Sequence(ca8CHStrimmedJetTracksAssociatorAtVertex+impactParameterTagInfosCA8CHStrimmed+secondaryVertexTagInfosCA8CHStrimmed+combinedSecondaryVertexBJetTagsCA8CHStrimmed)

patJetsCA8CHStrimmed = patJets.clone()
patJetsCA8CHStrimmed.jetSource = jetSource
patJetsCA8CHStrimmed.addJetCharge = False
patJetsCA8CHStrimmed.embedCaloTowers = False
patJetsCA8CHStrimmed.embedPFCandidates = False
patJetsCA8CHStrimmed.addAssociatedTracks = False
patJetsCA8CHStrimmed.addBTagInfo = True
patJetsCA8CHStrimmed.addDiscriminators = True
patJetsCA8CHStrimmed.addJetID = False
patJetsCA8CHStrimmed.tagInfoSources = cms.VInputTag(cms.InputTag("secondaryVertexTagInfosCA8CHStrimmed"))
patJetsCA8CHStrimmed.trackAssociationSource = cms.InputTag("ca8CHStrimmedJetTracksAssociatorAtVertex")
patJetsCA8CHStrimmed.discriminatorSources = cms.VInputTag(cms.InputTag("combinedSecondaryVertexBJetTagsCA8CHStrimmed"))
patJetsCA8CHStrimmed.getJetMCFlavour = False
patJetsCA8CHStrimmed.jetCorrFactorsSource = cms.VInputTag(cms.InputTag('patJetCorrFactorsCA8CHStrimmed'))
patJetsCA8CHStrimmed.genPartonMatch = cms.InputTag('patJetPartonMatchCA8CHStrimmed')
patJetsCA8CHStrimmed.genJetMatch = cms.InputTag('patJetGenJetMatchCA8CHStrimmed')

selectedPatJetsCA8CHStrimmed = selectedPatJets.clone()
selectedPatJetsCA8CHStrimmed.src = 'patJetsCA8CHStrimmed'
selectedPatJetsCA8CHStrimmed.cut = 'pt()>30'

ca8TrimmedGenJetsNoNu = ak7GenJetsNoNu.clone()
ca8TrimmedGenJetsNoNu.rParam = 0.8
ca8TrimmedGenJetsNoNu.jetAlgorithm = "CambridgeAachen"
ca8TrimmedGenJetsNoNu.doAreaFastjet = True
ca8TrimmedGenJetsNoNu.useTrimming = cms.bool(True)
ca8TrimmedGenJetsNoNu.useExplicitGhosts = cms.bool(True)
ca8TrimmedGenJetsNoNu.writeCompound = cms.bool(True)
ca8TrimmedGenJetsNoNu.rFilt = cms.double(0.1)
ca8TrimmedGenJetsNoNu.trimPtFracMin = cms.double(0.03)
ca8TrimmedGenJetsNoNu.jetPtMin = 25

patGenJetsCA8CHStrimmed = patJets.clone()
patGenJetsCA8CHStrimmed.jetSource = 'ca8TrimmedGenJetsNoNu'
patGenJetsCA8CHStrimmed.addGenJetMatch = False
patGenJetsCA8CHStrimmed.addGenPartonMatch = False
patGenJetsCA8CHStrimmed.addJetCharge = False
patGenJetsCA8CHStrimmed.embedCaloTowers = False
patGenJetsCA8CHStrimmed.embedPFCandidates = False
patGenJetsCA8CHStrimmed.addAssociatedTracks = False
patGenJetsCA8CHStrimmed.addBTagInfo = False
patGenJetsCA8CHStrimmed.addDiscriminators = False
patGenJetsCA8CHStrimmed.addJetID = False
patGenJetsCA8CHStrimmed.getJetMCFlavour = False
patGenJetsCA8CHStrimmed.addJetCorrFactors = False

jetMCSequenceCA8CHStrimmed = cms.Sequence(
    patJetPartonMatchCA8CHStrimmed +
    genParticlesForJetsNoNu +
    ca8TrimmedGenJetsNoNu +
    patGenJetsCA8CHStrimmed +
    patJetGenJetMatchCA8CHStrimmed
    )

PATCMGJetSequenceCA8CHStrimmed = cms.Sequence(
    ca8PFJetsCHStrimmed +
    jetMCSequenceCA8CHStrimmed +
    patJetCorrFactorsCA8CHStrimmed +
    btaggingCA8CHStrimmed +
    patJetsCA8CHStrimmed +
    selectedPatJetsCA8CHStrimmed
    )




#### Adding Nsubjetiness

selectedPatJetsCA8CHSwithNsub = cms.EDProducer("NjettinessAdder",
                              src=cms.InputTag("selectedPatJetsCA8CHS"),
                              cone=cms.double(0.8)
                              )

#### Adding QJets

selectedPatJetsCA8CHSwithQjets = cms.EDProducer("QjetsAdder",
                               src=cms.InputTag("selectedPatJetsCA8CHSwithNsub"),
                               zcut=cms.double(0.1),
                               dcutfctr=cms.double(0.5),
                               expmin=cms.double(0.0),
                               expmax=cms.double(0.0),
                               rigidity=cms.double(0.1),
                               ntrial = cms.int32(50),
                               cutoff=cms.double(200.0),
                               jetRad= cms.double(2.0),
                               jetAlgo=cms.string("CA"),
                               preclustering = cms.int32(50),
                              )

ca8Jets = cms.Sequence( PATCMGJetSequenceCA8CHS + PATCMGJetSequenceCA8CHSpruned + selectedPatJetsCA8CHSwithNsub + selectedPatJetsCA8CHSwithQjets )
