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
createPATtuple = True
createCMGtuple = False
skimEvents = True
runOnMC    = True
runOld5XGT = False
runOnFastSim = False

## Z-->ee/mumu filter
EE = ("patElectronsWithTrigger patElectronsWithTrigger")
MUMU = ("patMuonsWithTrigger patMuonsWithTrigger")
DILEPTON_KINCUT = ("70.0 < mass < 110.0 && pt > 5.0")

## W-->enu/munu filter
ENU = ("patMETs patElectronsWithTrigger")
MUNU = ("patMETs patMuonsWithTrigger")
METLEPTON_KINCUT = ("pt > 10.0")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

## Options and Output Report
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

## Input files
process.load("ExoDiBosonResonances.PATtupleProduction.JHU_Bulk600_ZZ_c1_cff")
#process.source.fileNames = ["file:/home/trtomei/shared/QstarGIToQZee_M2000_AOD53X_1kEvts_20121113.root"]
outputFileName = "newPatTuple_ZZ_600_HEEP.root"

# If you want you can overwrite the previous input filenames like this:
#process.source.fileNames = ['file:root://eoscms//eos/cms/store/cmst3/group/cmgtools/CMG/DY2JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/V5_B/PFAOD_0.root']

## Maximal Number of Events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

print sep_line
print process.source.fileNames
print sep_line 


########################################################
## Setup PAT Sequence
########################################################

print 'loading the main CMG sequence'

process.load('CMGTools.Common.PAT.PATCMG_cff')

#### Adding AK7 pruned jets

from CMGTools.Common.PAT.jetSubstructure_cff import *
process.PATCMGSequence += PATCMGJetSequenceAK7CHSpruned

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

    if isNewerThan('CMSSW_5_2_0'):
        process.PATCMGJetSequenceCHSpruned.remove( process.jetMCSequenceCHSpruned )
        process.patJetsCHSpruned.addGenJetMatch = False
        process.patJetsCHSpruned.addGenPartonMatch = False
        process.PATCMGJetSequenceAK7CHSpruned.remove( process.jetMCSequenceAK7CHSpruned )
        process.patJetsAK7CHSpruned.addGenJetMatch = False
        process.patJetsAK7CHSpruned.addGenPartonMatch = False

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
    if isNewerThan('CMSSW_5_2_0'):
        process.patJetCorrFactorsCHSpruned.levels.append('L2L3Residual')
        process.patJetCorrFactorsAK7CHSpruned.levels.append('L2L3Residual')


print 'cloning the jet sequence to build PU chs jets'

from PhysicsTools.PatAlgos.tools.helpers import cloneProcessingSnippet
process.PATCMGJetCHSSequence = cloneProcessingSnippet(process, process.PATCMGJetSequence, 'CHS')
process.PATCMGJetCHSSequence.insert( 0, process.ak5PFJetsCHS )
from CMGTools.Common.Tools.visitorUtils import replaceSrc
replaceSrc( process.PATCMGJetCHSSequence, 'ak5PFJets', 'ak5PFJetsCHS')
replaceSrc( process.PATCMGJetCHSSequence, 'particleFlow', 'pfNoPileUp')
process.patJetCorrFactorsCHS.payload = 'AK5PFchs'
process.selectedPatJetsCHS.cut = 'pt()>10'

### Boosted electrons isolation
### Remake the HEEP ID with no isolation cuts
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


########################################################
## Path definition
########################################################

process.dump = cms.EDAnalyzer('EventContentAnalyzer')

process.load('CMGTools.Common.PAT.addFilterPaths_cff')
process.p = cms.Path(
    process.prePathCounter + 
    process.PATCMGSequence +
    process.PATCMGJetCHSSequence
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


########################################################
## PAT output definition
########################################################

## Define event selection
EventSelection = cms.vstring('ZToEEskimPath','ZToMUMUskimPath','WToENUskimPath','WToMUNUskimPath')
if not skimEvents:
    EventSelection.append('p')

## Output Module Configuration (expects a path 'p')
from CMGTools.Common.eventContent.patEventContentCMG_cff import patEventContentCMG
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string(outputFileName),
                               # save only events passing the full path
                               SelectEvents   = cms.untracked.PSet( SelectEvents = EventSelection ),
                               # save PAT Layer 1 output; you need a '*' to
                               # unpack the list of commands 'patEventContent'
                               outputCommands = patEventContentCMG
                               )
# needed to override the CMG format, which drops the pat taus
process.out.outputCommands.append('keep patTaus_selectedPatTaus_*_*')

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

from CMGTools.Common.Tools.getGlobalTag import getGlobalTag

process.GlobalTag.globaltag = getGlobalTag( runOnMC, runOld5XGT )
print 'Global tag       : ', process.GlobalTag.globaltag


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
## Close the schedule
process.schedule.append( process.outpath )
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
