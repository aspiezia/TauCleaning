import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.tools.helpers import *


process = cms.Process("CMG")
###########
# Options #
###########
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents))
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 5000


#process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#from Configuration.PyReleaseValidation.autoCond import autoCond
#process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )
process.GlobalTag.globaltag = cms.string("START53_V7A::All")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff")

readFiles = cms.untracked.vstring()
process.source = cms.Source("PoolSource",
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
                            fileNames = readFiles
                            )

readFiles.extend(['file:/afs/cern.ch/user/t/tomei/public/patTuple_JHU_Bulk600_ZZ_c1_small.root'])


process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string("test.root"),
                               outputCommands = cms.untracked.vstring('drop *',
                                                                      'keep *_*_*_*'
                                                                      )
                               )

process.outpath = cms.EndPath(process.out)


################
# Ele Sequence #
################
              

process.test = cms.EDProducer("NjettinessAdder",
                              src=cms.InputTag("selectedPatJets"),
                              cone=cms.double(0.5)
                              )

process.test2 = cms.EDProducer("QjetsAdder",
                               src=cms.InputTag("test"),
                               zcut=cms.double(0.1),
                               dcutfctr=cms.double(0.5),
                               expmin=cms.double(0.0),
                               expmax=cms.double(0.0),
                               rigidity=cms.double(0.1),
                               ntrial = cms.int32(50),
                               cutoff=cms.double(10.0),
                               jetRad= cms.double(0.5),
                               jetAlgo=cms.string("AK"),
                               preclustering = cms.int32(50),
                              )

process.testp = cms.Path(process.test+process.test2 )
