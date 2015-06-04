import FWCore.ParameterSet.Config as cms

from CMGTools.Common.Tools.cmsswRelease import isNewerThan

sep_line = '-'*67
print sep_line
print 'CMG PAT-tuplizer, contact Colin before any modification'
print sep_line

process = cms.Process("PAT")


print 'querying database for source files'


########################################################
## User options
########################################################

## Basic options
createPATtuple = False
createCMGtuple = False
skimEvents = True
runOnMC    = False
runOld5XGT = False
runOnFastSim = False
runQJets = True

## Z-->ee/mumu filter
EE = ("patElectronsWithTrigger patElectronsWithTrigger")
MUMU = ("patMuonsWithTrigger patMuonsWithTrigger")
DILEPTON_KINCUT = ("70.0 < mass < 110.0 && pt > 80.0")

## W-->enu/munu filter
ENU = ("patMETs patElectronsWithTrigger")
MUNU = ("patMETs patMuonsWithTrigger")
METLEPTON_KINCUT = ("pt > 80.0")

## JetMet filter
JET_KINCUT = ("pt > 100")
MET_KINCUT = ("et > 50")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

## Input files
from CMGTools.Production.datasetToSource import *

datasetInfo = (
    'CMS',
    '/DYJetsToLL_PtZ-100_TuneZ2star_8TeV_ext-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM',
    '.*root')
process.source = datasetToSource(
    *datasetInfo
    )



process.source.fileNames = process.source.fileNames[:20]
# If you want you can overwrite the previous input filenames like this:
#####process.source.fileNames = ['file:root://eoscms//eos/cms'/store/cmst3/group/cmgtools/CMG/DY2JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/V5_B/PFAOD_0.root']
#process.source.fileNames = ['file:root://eoscms//eos/cms'/store/data/Run2012C/DoublePhotonHighPt/AOD/PromptReco-v2/000/200/600/0AA4E2BB-2BE5-E111-B3A8-00237DDC5C24.root']
##process.source.fileNames = ['file:root://eoscms//eos/cms'/store/cmst3/user/bonato/patTuple/2012/EXOVVtest/BulkG_ZZllqq_M1000_c0p2_STEP3_AODSIM_24_1_fvt.root']
process.source.fileNames = [
'/store/data/Run2012D/JetHT/AOD/22Jan2013-v1/10002/E841A6AE-2798-E211-870D-E0CB4E1A118F.root',
'/store/data/Run2012D/JetHT/AOD/22Jan2013-v1/10002/E8416FFA-8197-E211-9ABE-90E6BA19A214.root',
'/store/data/Run2012D/JetHT/AOD/22Jan2013-v1/10002/E816B170-E597-E211-8C3E-90E6BA0D09CE.root',
'/store/data/Run2012D/JetHT/AOD/22Jan2013-v1/10002/E6BCE747-9897-E211-BBD3-20CF305616D0.root',
'/store/data/Run2012D/JetHT/AOD/22Jan2013-v1/10002/E8F31E8D-0A98-E211-981C-20CF30561711.root',
'/store/data/Run2012D/JetHT/AOD/22Jan2013-v1/10002/E89A7E5D-9D97-E211-A76F-90E6BA19A205.root',
'/store/data/Run2012D/JetHT/AOD/22Jan2013-v1/10002/E888E66D-4B97-E211-B24B-E0CB4E29C4DE.root',
'/store/data/Run2012D/JetHT/AOD/22Jan2013-v1/10002/E853050A-3D98-E211-9E26-00259073E398.root',
'/store/data/Run2012D/JetHT/AOD/22Jan2013-v1/10002/EC72543F-0E98-E211-86AE-485B39800C1D.root',
'/store/data/Run2012D/JetHT/AOD/22Jan2013-v1/10002/EAF16CB5-1B98-E211-8551-485B39800BB9.root',
'/store/data/Run2012D/JetHT/AOD/22Jan2013-v1/10002/EAB0443F-D897-E211-8D15-E0CB4E1A117E.root',
'/store/data/Run2012D/JetHT/AOD/22Jan2013-v1/10002/E8F6EC47-9897-E211-A47D-20CF3058709A.root',
'/store/data/Run2012D/JetHT/AOD/22Jan2013-v1/10002/ECCBD06E-4B97-E211-A0FB-0025907277BE.root',
'/store/data/Run2012D/JetHT/AOD/22Jan2013-v1/10002/ECC36573-4B97-E211-AFCE-001E4F3F28DC.root',
'/store/data/Run2012D/JetHT/AOD/22Jan2013-v1/10002/ECC0B73B-0698-E211-8E7C-E0CB4E29C4D5.root',
'/store/data/Run2012D/JetHT/AOD/22Jan2013-v1/10002/ECA3DB0E-D297-E211-8FE7-20CF3056171C.root',
'/store/data/Run2012D/JetHT/AOD/22Jan2013-v1/10002/EEEB563D-1A98-E211-B610-001EC9D4A1FD.root',
'/store/data/Run2012D/JetHT/AOD/22Jan2013-v1/10002/EE46C244-9D97-E211-9713-E0CB4E29C4FD.root',
'/store/data/Run2012D/JetHT/AOD/22Jan2013-v1/10002/EE4594A2-5097-E211-9C63-0025907277BE.root',
'/store/data/Run2012D/JetHT/AOD/22Jan2013-v1/10002/EE2AA631-0198-E211-AF5A-485B3989725B.root',
'/store/data/Run2012D/JetHT/AOD/22Jan2013-v1/10002/F0A85055-0D98-E211-B422-E0CB4E4408DE.root',
'/store/data/Run2012D/JetHT/AOD/22Jan2013-v1/10002/F0854BE9-1E98-E211-B66E-E0CB4E19F95B.root',
'/store/data/Run2012D/JetHT/AOD/22Jan2013-v1/10002/F069364C-2998-E211-B939-E0CB4E19F9A1.root',
'/store/data/Run2012D/JetHT/AOD/22Jan2013-v1/10002/F03BEF05-3098-E211-9FF6-002618FDA1BD.root',
'/store/data/Run2012D/JetHT/AOD/22Jan2013-v1/10002/F2451EA6-4697-E211-A5EF-E0CB4E4408DF.root',
'/store/data/Run2012D/JetHT/AOD/22Jan2013-v1/10002/F0EF9210-C497-E211-8E26-00259073E3D6.root',
'/store/data/Run2012D/JetHT/AOD/22Jan2013-v1/10002/F0D1FB55-9D97-E211-9C5C-E0CB4E1A117B.root',
'/store/data/Run2012D/JetHT/AOD/22Jan2013-v1/10002/F0AAF115-D297-E211-A22D-E0CB4E1A116A.root',
'/store/data/Run2012D/JetHT/AOD/22Jan2013-v1/10002/F2983694-B597-E211-8686-BCAEC54B304A.root',
'/store/data/Run2012D/JetHT/AOD/22Jan2013-v1/10002/F288EDE1-9397-E211-B204-002590747DEC.root',
'/store/data/Run2012D/JetHT/AOD/22Jan2013-v1/10002/F288751A-8897-E211-AA73-001EC9D81A3E.root',
]

## Maximal Number of Events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

print sep_line
print process.source.fileNames
print sep_line 


########################################################
## Setup PAT Sequence
########################################################

print 'loading the main CMG sequence'

process.load('CMGTools.Common.PAT.PATCMG_cff')
from CMGTools.Common.eventContent.patEventContentCMG_cff import patEventContentCMG

#### Adding AK7 jets

process.load("ExoDiBosonResonances.PATtupleProduction.PAT_ak7jets_cff")
process.PATCMGSequence += process.PATCMGJetSequenceAK7CHS

#### Adding AK7 pruned jets

process.load("CMGTools.Common.PAT.jetSubstructure_cff")
process.PATCMGSequence.remove(process.PATCMGJetSequenceCHSpruned) # don't produce the AK5 pruned collections
process.jetMCSequenceAK7CHSpruned.remove(process.ak7GenJetsNoNu) # don't cluster the ak7GenJetsNoNu twice
process.selectedPatJetsAK7CHSpruned.cut = 'pt()>20'
process.PATCMGSequence += process.PATCMGJetSequenceAK7CHSpruned
patEventContentCMG+=['keep *_ak7PFJetsCHSpruned_SubJets_*']
patEventContentCMG+=['keep *_ak7GenJetsNoNu_*_*']


#### Adding CA8 jets and CA8 pruned jets

process.load("ExoDiBosonResonances.PATtupleProduction.PAT_ca8jets_cff")
process.PATCMGSequence += process.PATCMGJetSequenceCA8CHS
process.PATCMGSequence += process.PATCMGJetSequenceCA8CHSpruned
patEventContentCMG+=['keep *_ca8PFJetsCHSpruned_SubJets_*']
patEventContentCMG+=['keep *_ca8GenJetsNoNu_*_*']


#### Adding Nsubjetiness

process.selectedPatJetsAK7CHSwithNsub = cms.EDProducer("NjettinessAdder",
                              src=cms.InputTag("selectedPatJetsAK7CHS"),
                              cone=cms.double(0.7)
                              )
process.PATCMGSequence += process.selectedPatJetsAK7CHSwithNsub
patEventContentCMG+=['drop patJets_selectedPatJetsAK7CHS_*_*']
process.PATCMGSequence += process.selectedPatJetsCA8CHSwithNsub
patEventContentCMG+=['drop patJets_selectedPatJetsCA8CHS_*_*']

#### Adding QJets

process.selectedPatJetsAK7CHSwithQjets = cms.EDProducer("QjetsAdder",
			   src=cms.InputTag("selectedPatJetsAK7CHSwithNsub"),
			   zcut=cms.double(0.1),
			   dcutfctr=cms.double(0.5),
			   expmin=cms.double(0.0),
			   expmax=cms.double(0.0),
			   rigidity=cms.double(0.1),
			   ntrial = cms.int32(50),
			   cutoff=cms.double(100.0),
			   jetRad= cms.double(0.7),
			   jetAlgo=cms.string("AK"),
			   preclustering = cms.int32(30),
			  )
if not runQJets:
    process.selectedPatJetsAK7CHSwithQjets.cutoff=100000.0
process.PATCMGSequence += process.selectedPatJetsAK7CHSwithQjets
patEventContentCMG+=['drop patJets_selectedPatJetsAK7CHSwithNsub_*_*']
process.PATCMGSequence += process.selectedPatJetsCA8CHSwithQjets
patEventContentCMG+=['drop patJets_selectedPatJetsCA8CHSwithNsub_*_*']

######ADD PU JET ID

from  CMGTools.External.pujetidsequence_cff import puJetId, puJetMva
process.puJetIdAK7CHS = puJetId.clone(
    jets ='selectedPatJetsAK7CHSwithQjets',
    jec = 'AK7chs'
    )
process.PATCMGSequence += process.puJetIdAK7CHS
patEventContentCMG+=['keep *_puJetIdAK7CHS_*_*']
process.puJetIdCA8CHS = puJetId.clone(
    jets ='selectedPatJetsCA8CHSwithQjets',
    jec = 'AK7chs'
    )

from CMGTools.External.pujetidproducer_cfi import  stdalgos_4x, stdalgos_5x, stdalgos, cutbased, chsalgos_4x, chsalgos_5x, chsalgos

process.puJetMvaAK5CHS= puJetMva.clone(
    jetids = cms.InputTag("puJetIdCHS"),
####    jets ='selectedPatJetsCHS',
    jets ='patJetsWithVarCHS',
    algos =  chsalgos
    )

process.puJetMvaAK7CHS= puJetMva.clone(
    jetids = cms.InputTag("puJetIdAK7CHS"),
    jets ='selectedPatJetsAK7CHSwithQjets',
    algos =  chsalgos
    )

process.puJetMvaCA8CHS= puJetMva.clone(
    jetids = cms.InputTag("puJetIdCA8CHS"),
    jets ='selectedPatJetsCA8CHSwithQjets',
    algos =  chsalgos
    )

process.puJetIdAK5Sequence = cms.Sequence(                      process.puJetMvaAK5CHS)
process.puJetIdAK7Sequence = cms.Sequence(process.puJetIdAK7CHS+process.puJetMvaAK7CHS)
process.puJetIdCA8Sequence = cms.Sequence(process.puJetIdCA8CHS+process.puJetMvaCA8CHS)
#### these are moved down below this same cfg, after having built AK5 CHS jets
#process.PATCMGSequence += process.puJetIdAK5Sequence
#process.PATCMGSequence += process.puJetIdAK7Sequence
#process.PATCMGSequence += process.puJetIdCA8Sequence
patEventContentCMG+=['keep *_puJetIdAK7CHS_*_*']
patEventContentCMG+=['keep *_puJetIdCA8CHS_*_*']
patEventContentCMG+=['keep *_puJetMvaAK5CHS_*_*']
patEventContentCMG+=['keep *_puJetMvaAK7CHS_*_*']
patEventContentCMG+=['keep *_puJetMvaCA8CHS_*_*']

if runOnMC is False:
    # removing MC stuff
    print 'removing MC stuff, as we are running on Data'

    process.patElectrons.addGenMatch = False
    process.makePatElectrons.remove( process.electronMatch )
    
    process.patMuons.addGenMatch = False
    process.makePatMuons.remove( process.muonMatch )
    
    process.PATCMGSequence.remove( process.PATCMGGenSequence )
    process.PATCMGJetSequence.remove( process.jetMCSequence )
    process.PATCMGJetSequence.remove( process.patJetFlavourId )
    process.patJets.addGenJetMatch = False
    process.patJets.addGenPartonMatch = False

   #### if isNewerThan('CMSSW_5_2_0'):
    process.PATCMGJetSequenceCHSpruned.remove( process.jetMCSequenceCHSpruned )
    process.patJetsCHSpruned.addGenJetMatch = False
    process.patJetsCHSpruned.addGenPartonMatch = False
    process.PATCMGJetSequenceAK7CHS.remove( process.jetMCSequenceAK7CHS )
    process.patJetsAK7CHS.addGenJetMatch = False
    process.patJetsAK7CHS.addGenPartonMatch = False
    process.PATCMGJetSequenceAK7CHSpruned.remove( process.jetMCSequenceAK7CHSpruned )
    process.patJetsAK7CHSpruned.addGenJetMatch = False
    process.patJetsAK7CHSpruned.addGenPartonMatch = False
    process.PATCMGJetSequenceCA8CHS.remove( process.jetMCSequenceCA8CHS )
    process.patJetsCA8CHS.addGenJetMatch = False
    process.patJetsCA8CHS.addGenPartonMatch = False
    process.PATCMGJetSequenceCA8CHSpruned.remove( process.jetMCSequenceCA8CHSpruned )
    process.patJetsCA8CHSpruned.addGenJetMatch = False
    process.patJetsCA8CHSpruned.addGenPartonMatch = False
    
    process.PATCMGTauSequence.remove( process.tauGenJets )
    process.PATCMGTauSequence.remove( process.tauGenJetsSelectorAllHadrons )
    process.PATCMGTauSequence.remove( process.tauGenJetMatch )
    process.PATCMGTauSequence.remove( process.tauMatch )
    process.patTaus.addGenJetMatch = False
    process.patTaus.addGenMatch = False

    process.patMETs.addGenMET = False 
    process.patMETsRaw.addGenMET = False 

    # setting up JSON file
    # json = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/DCSOnly/json_DCSONLY.txt'
    # print 'using json file: ', json
    # from CMGTools.Common.Tools.applyJSON_cff import *
    # applyJSON(process, json )

    # adding L2L3Residual corrections
    process.patJetCorrFactors.levels.append('L2L3Residual')

     ###   if isNewerThan('CMSSW_5_2_0'):
    process.patJetCorrFactorsCHSpruned.levels.append('L2L3Residual')
    process.patJetCorrFactorsAK7CHS.levels.append('L2L3Residual')
    process.patJetCorrFactorsAK7CHSpruned.levels.append('L2L3Residual')
    process.patJetCorrFactorsCA8CHS.levels.append('L2L3Residual')
    process.patJetCorrFactorsCA8CHSpruned.levels.append('L2L3Residual')
    

#### Adding HEEP and modified isolation
### Boosted electrons isolation
### Remake the HEEP ID with no isolation cuts
print "Adding HEEP and modified isolation..."

process.load("RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi")
from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import *

process.heepIdNoIso = cms.EDProducer("HEEPIdValueMapProducer",
                                     eleLabel = cms.InputTag("gsfElectrons"),
                                     barrelCuts = cms.PSet(heepBarrelCuts),
                                     endcapCuts = cms.PSet(heepEndcapCuts),
                                     eleIsolEffectiveAreas = cms.PSet(heepEffectiveAreas),
                                     eleRhoCorrLabel = cms.InputTag("kt6PFJets", "rho"),
                                     verticesLabel = cms.InputTag("offlinePrimaryVertices"),
                                     applyRhoCorrToEleIsol = cms.bool(True),
                                     writeIdAsInt = cms.bool(True)
                                     )

process.heepIdNoIso.barrelCuts.cuts=cms.string("et:detEta:ecalDriven:dEtaIn:dPhiIn:hadem:e2x5Over5x5:nrMissHits:dxy")
process.heepIdNoIso.endcapCuts.cuts=cms.string("et:detEta:ecalDriven:dEtaIn:dPhiIn:hadem:sigmaIEtaIEta:nrMissHits:dxy")

process.heepIdNoIsoEles = cms.EDProducer("tsw::HEEPGsfProducer",
                                         cutValueMap = cms.InputTag("heepIdNoIso"),
                                         inputGsfEles = cms.InputTag("gsfElectrons")
                                         )

from TSWilliams.BstdZeeTools.bstdzeemodisolproducer_cff import *

process.modElectronIso = cms.EDProducer("BstdZeeModIsolProducer",
                                        bstdZeeModIsolParams,
                                        vetoGsfEles = cms.InputTag("heepIdNoIsoEles")
                                        )

# Adding these new value maps to the patElectrons
process.patElectrons.userIsolation.user = cms.VPSet(
    cms.PSet(src = cms.InputTag("modElectronIso","track")),
    cms.PSet(src = cms.InputTag("modElectronIso","ecal")),
    cms.PSet(src = cms.InputTag("modElectronIso","hcalDepth1")),
    )

# And now redo the PAT electrons with Mathias' module
process.heepPatElectrons = cms.EDProducer("HEEPAttStatusToPAT",
                                          eleLabel = cms.InputTag("patElectrons"),
                                          barrelCuts = cms.PSet(heepBarrelCuts),
                                          endcapCuts = cms.PSet(heepEndcapCuts),
                                          applyRhoCorrToEleIsol = cms.bool(True),
                                          eleIsolEffectiveAreas = cms.PSet (heepEffectiveAreas),
                                          eleRhoCorrLabel = cms.InputTag("kt6PFJets","rho"),
                                          verticesLabel = cms.InputTag("offlinePrimaryVerticesWithBS"),
                                          )

# With no isolation cuts. In this way, the newly created PAT electrons have a complete HEEP ID, minus isolation cuts,
# while they ALSO have three userIsolations corresponding to the modified isolations.
process.heepPatElectrons.barrelCuts.cuts=cms.string("et:detEta:ecalDriven:dEtaIn:dPhiIn:hadem:e2x5Over5x5:nrMissHits:dxy")
process.heepPatElectrons.endcapCuts.cuts=cms.string("et:detEta:ecalDriven:dEtaIn:dPhiIn:hadem:sigmaIEtaIEta:nrMissHits:dxy")

process.PATCMGSequence.replace( process.patElectrons,
                                process.heepIdNoIso +
                                process.heepIdNoIsoEles +
                                process.modElectronIso +
                                process.patElectrons +
                                process.heepPatElectrons )

# Finally, change the source of the selectedPatElectrons.
# The changes trickle down from here.
process.selectedPatElectrons.src = cms.InputTag("heepPatElectrons")

#### Adding new TuneP muons and new Track Errors
# Load the tune P muons
print "Adding new TuneP muons and Track Errors..."

# This adds process.tunePmuons and process.muonTrackError
process.load("ExoDiBosonResonances.EDBRMuon.newTuneP_cff")

# Change the source of the patMuons
process.patMuons.muonSource = "tunePmuons"
# Add the user float
process.patMuons.userData.userFloats.src = ['muonTrackError']

# Since we don't use these isoDeposits, we might as well take them out.
process.patMuons.isoDeposits = cms.PSet()
process.patMuons.isolationValues = cms.PSet()
process.patMuons.embedCaloMETMuonCorrs = False # Don't use
process.patMuons.embedTcMETMuonCorrs = False # Don't use

# Put the new modules in the sequence
if runOnMC is False:
    process.PATCMGSequence.replace( process.patMuons,
                                    process.tunePmuons +
                                    process.muonTrackError +
                                    process.patMuons )

if runOnMC is True:
        process.PATCMGSequence.replace( process.muonMatch,
                                        process.tunePmuons +
                                        process.muonTrackError +
                                        process.muonMatch )
        process.muonMatch.src = "tunePmuons"

print 'cloning the jet sequence to build PU chs jets'

from PhysicsTools.PatAlgos.tools.helpers import cloneProcessingSnippet
process.PATCMGJetCHSSequence = cloneProcessingSnippet(process, process.PATCMGJetSequence, 'CHS')
process.PATCMGJetCHSSequence.insert( 0, process.ak5PFJetsCHS )
from CMGTools.Common.Tools.visitorUtils import replaceSrc
replaceSrc( process.PATCMGJetCHSSequence, 'ak5PFJets', 'ak5PFJetsCHS')
replaceSrc( process.PATCMGJetCHSSequence, 'particleFlow', 'pfNoPileUp')
jecPayload = 'AK5PFchs'
process.patJetsWithVarCHS.payload = jecPayload
process.patJetCorrFactorsCHS.payload = jecPayload
process.puJetIdCHS.jec = jecPayload
process.cmgPUJetMvaCHS.jec = jecPayload
process.selectedPatJetsCHS.cut = 'pt()>10'

# Change the soft muons? Change the MET?
###
### WW ANALYSIS - PAY ATTENTION TO THIS
###

#process.softMuonTagInfos.leptons    = "tunePmuons"
#process.softMuonTagInfosCHS.leptons = "tunePmuons"

###
### WW ANALYSIS - PAY ATTENTION TO THIS
###
#

########################################################
## Path definition
########################################################

process.dump = cms.EDAnalyzer('EventContentAnalyzer')

process.load('CMGTools.Common.PAT.addFilterPaths_cff')
process.p = cms.Path(
    process.prePathCounter + 
    process.PATCMGSequence +
    process.PATCMGJetCHSSequence+
    process.puJetIdAK5Sequence+process.puJetIdAK7Sequence+process.puJetIdCA8Sequence
    )

process.p += process.postPathCounter

# For testing, you can remove some of the objects:
# NOTE: there are a few dependencies between these sequences
# process.PATCMGSequence.remove(process.PATCMGPileUpSubtractionSequence)
# process.PATCMGSequence.remove(process.PATCMGRhoSequence)
# process.PATCMGSequence.remove(process.PATCMGMuonSequence)
# process.PATCMGSequence.remove(process.PATCMGElectronSequence)
# process.PATCMGSequence.remove(process.PATCMGGenSequence)
# process.PATCMGSequence.remove(process.PATCMGJetSequence)
# process.PATCMGSequence.remove(process.PATCMGTauSequence)
# process.PATCMGSequence.remove(process.PATCMGMetSequence)
# process.p.remove(process.PATCMGJetCHSSequence)
# process.p.remove(process.PATCMGTriggerSequence)
# process.p.remove(process.PATCMGPhotonSequence)
# process.p.remove(process.PATCMGVertexSequence)
# process.p.remove(process.PATCMGPhotonSequence)
# process.p.remove(process.MetSignificanceSequence)
# process.p.remove(process.PATCMGMetRegressionSequence)
# process.p.remove(process.PATCMGJetSequenceCHSpruned)

if runOnFastSim :
    process.vertexWeightSequence.remove(process.vertexWeight3DMay10ReReco)
    process.vertexWeightSequence.remove(process.vertexWeight3DMay10ReReco)
    process.vertexWeightSequence.remove(process.vertexWeight3DPromptRecov4)
    process.vertexWeightSequence.remove(process.vertexWeight3D05AugReReco)
    process.vertexWeightSequence.remove(process.vertexWeight3DPromptRecov6)
    process.vertexWeightSequence.remove(process.vertexWeight3D2invfb)
    process.vertexWeightSequence.remove(process.vertexWeight3D2011B)
    process.vertexWeightSequence.remove(process.vertexWeight3D2011AB)
    process.vertexWeightSequence.remove(process.vertexWeight3DFall11May10ReReco)
    process.vertexWeightSequence.remove(process.vertexWeight3DFall11PromptRecov4)
    process.vertexWeightSequence.remove(process.vertexWeight3DFall1105AugReReco)
    process.vertexWeightSequence.remove(process.vertexWeight3DFall11PromptRecov6)
    process.vertexWeightSequence.remove(process.vertexWeight3DFall112invfb)
    process.vertexWeightSequence.remove(process.vertexWeight3DFall112011B)
    process.vertexWeightSequence.remove(process.vertexWeight3DFall112011AB)


########################################################
## Skim events
########################################################


## ZToEE
process.ZToEEcand = cms.EDProducer("CandViewShallowCloneCombiner",
                                   decay = cms.string(EE),
                                   checkCharge = cms.bool(False),
                                   cut = cms.string(DILEPTON_KINCUT)
                                   )

process.ZToEEfilter = cms.EDFilter("CandViewCountFilter",
                                   src = cms.InputTag("ZToEEcand"),
                                   minNumber = cms.uint32(1) )

process.ZToEEskimSequence = cms.Sequence( process.ZToEEcand * process.ZToEEfilter )
process.ZToEEskimPath = cms.Path( process.ZToEEskimSequence )


## ZToMUMU
process.ZToMUMUcand = cms.EDProducer("CandViewShallowCloneCombiner", 
                                     decay = cms.string(MUMU),
                                     checkCharge = cms.bool(False),
                                     cut = cms.string(DILEPTON_KINCUT)  
                                     )

process.ZToMUMUfilter = cms.EDFilter("CandViewCountFilter",
                                     src = cms.InputTag("ZToMUMUcand"),
                                     minNumber = cms.uint32(1) )

process.ZToMUMUskimSequence = cms.Sequence( process.ZToMUMUcand * process.ZToMUMUfilter )
process.ZToMUMUskimPath = cms.Path( process.ZToMUMUskimSequence )


## WToENU
#MT="sqrt(2*daughter(0).pt*daughter(1).pt*(1 - cos(daughter(0).phi - daughter(1).phi)))"
process.WToENUcand = cms.EDProducer("CandViewShallowCloneCombiner",
                                    decay = cms.string(ENU),
                                    checkCharge = cms.bool(False),
                                    cut = cms.string(METLEPTON_KINCUT)  
                                    )

process.WToENUfilter = cms.EDFilter("CandViewCountFilter",
                                    src = cms.InputTag("WToENUcand"),
                                    minNumber = cms.uint32(1) )

process.WToENUskimSequence = cms.Sequence( process.WToENUcand * process.WToENUfilter )
process.WToENUskimPath = cms.Path( process.WToENUskimSequence )


## WToMUNU
#MT="sqrt(2*daughter(0).pt*daughter(1).pt*(1 - cos(daughter(0).phi - daughter(1).phi)))"
process.WToMUNUcand = cms.EDProducer("CandViewShallowCloneCombiner",
                                     decay = cms.string(MUNU),
                                     checkCharge = cms.bool(False),
                                     cut = cms.string(METLEPTON_KINCUT)  
                                     )

process.WToMUNUfilter = cms.EDFilter("CandViewCountFilter",
                                     src = cms.InputTag("WToMUNUcand"),
                                     minNumber = cms.uint32(1) )

process.WToMUNUskimSequence = cms.Sequence( process.WToMUNUcand * process.WToMUNUfilter )
process.WToMUNUskimPath = cms.Path( process.WToMUNUskimSequence )

## JetMET
### Add cuts in jet and MET
process.selectedPatJetsForCut = process.selectedPatJetsCA8CHS.clone(
                                    cut = cms.string(JET_KINCUT),
                                    filter = cms.bool(True)
                                    )

process.selectedPatMETs = cms.EDFilter( "PATMETSelector",
                        		       src = cms.InputTag("patMETs"),
                                       cut = cms.string(MET_KINCUT),
                                       filter = cms.bool(True)
                                       )

process.JetMETskimSequence = cms.Sequence( process.selectedPatJetsForCut * process.selectedPatMETs )
process.JetMETskimPath = cms.Path( process.JetMETskimSequence )

########################################################
## PAT output definition
########################################################

## Define event selection
EventSelection = cms.vstring('ZToEEskimPath','ZToMUMUskimPath','WToENUskimPath','WToMUNUskimPath','JetMETskimPath')
if not skimEvents:
    EventSelection.append('p')

## Output Module Configuration (expects a path 'p')
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string("out.root"),
                               # save only events passing the full path
                               SelectEvents   = cms.untracked.PSet( SelectEvents = EventSelection ),
                               # save PAT Layer 1 output; you need a '*' to
                               # unpack the list of commands 'patEventContent'
                               outputCommands = patEventContentCMG
                               )
# needed to override the CMG format, which drops the pat taus
process.out.outputCommands.append('keep patTaus_selectedPatTaus_*_*')

#### drop collections not used by EXO-VV analysis
process.out.outputCommands.append('drop *_cmg*_*_*')
process.out.outputCommands.append('drop *_particleFlow*_*_*')
process.out.outputCommands.append('drop *_pfNoPileUp_*_*')
process.out.outputCommands.append('drop *_pfSelectedPhotons_*_*')
process.out.outputCommands.append('drop *_phPFIsoDeposit*_*_*')

#FIXME now keeping the whole event content...
# process.out.outputCommands.append('keep *_*_*_*')

process.outpath = cms.EndPath()

if createPATtuple:
    process.outpath += process.out


########################################################
## CMG output definition
########################################################

from CMGTools.Common.eventContent.patEventContentCMG_cff import everything
process.outcmg = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('cmgTuple.root'),
    SelectEvents   = cms.untracked.PSet( SelectEvents = EventSelection ),
    outputCommands = everything,
    dropMetaData = cms.untracked.string('PRIOR')
    )

if createCMGtuple:
    process.outpath += process.outcmg


########################################################
## Conditions 
########################################################

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")

### Set the global tag from
### https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions#Winter13_2012_A_B_C_D_datasets_r
EXOVV_GT= 'START53_V7G::All'
if runOnMC:
    if 'START53' in datasetInfo[1]:
        EXOVV_GT = 'START53_V23::All' 
    elif 'START52' in datasetInfo[1]:
        EXOVV_GT = 'START52_V9F::All'
else :
    if 'Run2012D' in datasetInfo[1]:
        EXOVV_GT = 'FT_53_V21_AN4::All'
    else:
        EXOVV_GT = 'FT_53_V21_AN4::All'
        
process.GlobalTag.globaltag = EXOVV_GT  ###cms.string("START53_V21::All")

### do it like CMGTools (but tags in getGlobalTag are not the latest for Jan22nd)
#from CMGTools.Common.Tools.getGlobalTag import getGlobalTag
### OLD way: process.GlobalTag.globaltag = getGlobalTag( runOnMC, runOld5XGT )
#process.GlobalTag.globaltag = getGlobalTagByDataset( runOnMC, datasetInfo[1])

print 'Global tag       : ', process.GlobalTag.globaltag
###

########################################################
## Below, stuff that you probably don't want to modify
########################################################

## Geometry and Detector Conditions (needed for a few patTuple production steps)

## Get the schedule
from CMGTools.Common.PAT.patCMGSchedule_cff import getSchedule
process.schedule = getSchedule(process, runOnMC, runOnFastSim)
## Add filters to the schedule
process.schedule.append( process.ZToEEskimPath )
process.schedule.append( process.ZToMUMUskimPath )
process.schedule.append( process.WToENUskimPath )
process.schedule.append( process.WToMUNUskimPath )
process.schedule.append( process.JetMETskimPath )


## Also add the TOBTEC Fakes Filter
process.load("KStenson.TrackingFilters.tobtecfakesfilter_cfi")
print process.tobtecfakesfilter
process.tobtecfakesfilterPath = cms.Path(~process.tobtecfakesfilter)
process.schedule.append( process.tobtecfakesfilterPath )

## Close the schedule
process.schedule.append( process.outpath )

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.MessageLogger.suppressWarning = cms.untracked.vstring('ecalLaserCorrFilter')
## Options and Output Report
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


## Print the schedule
print process.schedule

print sep_line

print 'Fastjet instances (dominating our processing time...):'
from CMGTools.Common.Tools.visitorUtils import SeqVisitor
v = SeqVisitor('FastjetJetProducer')
process.p.visit(v)

print sep_line

print 'starting CMSSW'

if not runOnMC and isNewerThan('CMSSW_5_2_0'):
    process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3Residual")
