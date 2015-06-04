import FWCore.ParameterSet.Config as cms
import sys

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.TFileService = cms.Service("TFileService", fileName = cms.string("histo.root") )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
"/store/group/phys_exotica/leptonsPlusJets/ExoDiBosonResonances/tomei/store/user/tomei/BulkG_ZZ_lljj_M2000_G120-JHU/BulkG_ZZ_lljj_M2000_G120-JHU/c8f8ed334db8a7d6f56c62266b1dfa5b/Bulk_AODSIM_6_1_GAz.root",
"/store/group/phys_exotica/leptonsPlusJets/ExoDiBosonResonances/tomei/store/user/tomei/BulkG_ZZ_lljj_M2000_G120-JHU/BulkG_ZZ_lljj_M2000_G120-JHU/c8f8ed334db8a7d6f56c62266b1dfa5b/Bulk_AODSIM_7_1_NaL.root",
"/store/group/phys_exotica/leptonsPlusJets/ExoDiBosonResonances/tomei/store/user/tomei/BulkG_ZZ_lljj_M2000_G120-JHU/BulkG_ZZ_lljj_M2000_G120-JHU/c8f8ed334db8a7d6f56c62266b1dfa5b/Bulk_AODSIM_8_1_vxY.root",
"/store/group/phys_exotica/leptonsPlusJets/ExoDiBosonResonances/tomei/store/user/tomei/BulkG_ZZ_lljj_M2000_G120-JHU/BulkG_ZZ_lljj_M2000_G120-JHU/c8f8ed334db8a7d6f56c62266b1dfa5b/Bulk_AODSIM_9_1_xA2.root",
"/store/group/phys_exotica/leptonsPlusJets/ExoDiBosonResonances/tomei/store/user/tomei/BulkG_ZZ_lljj_M2000_G120-JHU/BulkG_ZZ_lljj_M2000_G120-JHU/c8f8ed334db8a7d6f56c62266b1dfa5b/Bulk_AODSIM_5_1_Y18.root",
"/store/group/phys_exotica/leptonsPlusJets/ExoDiBosonResonances/tomei/store/user/tomei/BulkG_ZZ_lljj_M2000_G120-JHU/BulkG_ZZ_lljj_M2000_G120-JHU/c8f8ed334db8a7d6f56c62266b1dfa5b/Bulk_AODSIM_4_1_ql5.root",
"/store/group/phys_exotica/leptonsPlusJets/ExoDiBosonResonances/tomei/store/user/tomei/BulkG_ZZ_lljj_M2000_G120-JHU/BulkG_ZZ_lljj_M2000_G120-JHU/c8f8ed334db8a7d6f56c62266b1dfa5b/Bulk_AODSIM_3_1_RPy.root",
"/store/group/phys_exotica/leptonsPlusJets/ExoDiBosonResonances/tomei/store/user/tomei/BulkG_ZZ_lljj_M2000_G120-JHU/BulkG_ZZ_lljj_M2000_G120-JHU/c8f8ed334db8a7d6f56c62266b1dfa5b/Bulk_AODSIM_2_1_SYw.root",
"/store/group/phys_exotica/leptonsPlusJets/ExoDiBosonResonances/tomei/store/user/tomei/BulkG_ZZ_lljj_M2000_G120-JHU/BulkG_ZZ_lljj_M2000_G120-JHU/c8f8ed334db8a7d6f56c62266b1dfa5b/Bulk_AODSIM_21_1_lTl.root",
"/store/group/phys_exotica/leptonsPlusJets/ExoDiBosonResonances/tomei/store/user/tomei/BulkG_ZZ_lljj_M2000_G120-JHU/BulkG_ZZ_lljj_M2000_G120-JHU/c8f8ed334db8a7d6f56c62266b1dfa5b/Bulk_AODSIM_20_1_xRz.root",
"/store/group/phys_exotica/leptonsPlusJets/ExoDiBosonResonances/tomei/store/user/tomei/BulkG_ZZ_lljj_M2000_G120-JHU/BulkG_ZZ_lljj_M2000_G120-JHU/c8f8ed334db8a7d6f56c62266b1dfa5b/Bulk_AODSIM_1_1_cfP.root",
"/store/group/phys_exotica/leptonsPlusJets/ExoDiBosonResonances/tomei/store/user/tomei/BulkG_ZZ_lljj_M2000_G120-JHU/BulkG_ZZ_lljj_M2000_G120-JHU/c8f8ed334db8a7d6f56c62266b1dfa5b/Bulk_AODSIM_19_1_3NN.root",
"/store/group/phys_exotica/leptonsPlusJets/ExoDiBosonResonances/tomei/store/user/tomei/BulkG_ZZ_lljj_M2000_G120-JHU/BulkG_ZZ_lljj_M2000_G120-JHU/c8f8ed334db8a7d6f56c62266b1dfa5b/Bulk_AODSIM_18_1_1wJ.root",
"/store/group/phys_exotica/leptonsPlusJets/ExoDiBosonResonances/tomei/store/user/tomei/BulkG_ZZ_lljj_M2000_G120-JHU/BulkG_ZZ_lljj_M2000_G120-JHU/c8f8ed334db8a7d6f56c62266b1dfa5b/Bulk_AODSIM_17_1_OF1.root",
"/store/group/phys_exotica/leptonsPlusJets/ExoDiBosonResonances/tomei/store/user/tomei/BulkG_ZZ_lljj_M2000_G120-JHU/BulkG_ZZ_lljj_M2000_G120-JHU/c8f8ed334db8a7d6f56c62266b1dfa5b/Bulk_AODSIM_16_1_2VM.root",
"/store/group/phys_exotica/leptonsPlusJets/ExoDiBosonResonances/tomei/store/user/tomei/BulkG_ZZ_lljj_M2000_G120-JHU/BulkG_ZZ_lljj_M2000_G120-JHU/c8f8ed334db8a7d6f56c62266b1dfa5b/Bulk_AODSIM_15_1_Kji.root",
"/store/group/phys_exotica/leptonsPlusJets/ExoDiBosonResonances/tomei/store/user/tomei/BulkG_ZZ_lljj_M2000_G120-JHU/BulkG_ZZ_lljj_M2000_G120-JHU/c8f8ed334db8a7d6f56c62266b1dfa5b/Bulk_AODSIM_14_1_kjk.root",
"/store/group/phys_exotica/leptonsPlusJets/ExoDiBosonResonances/tomei/store/user/tomei/BulkG_ZZ_lljj_M2000_G120-JHU/BulkG_ZZ_lljj_M2000_G120-JHU/c8f8ed334db8a7d6f56c62266b1dfa5b/Bulk_AODSIM_13_1_ICX.root",
"/store/group/phys_exotica/leptonsPlusJets/ExoDiBosonResonances/tomei/store/user/tomei/BulkG_ZZ_lljj_M2000_G120-JHU/BulkG_ZZ_lljj_M2000_G120-JHU/c8f8ed334db8a7d6f56c62266b1dfa5b/Bulk_AODSIM_12_1_HGm.root",
"/store/group/phys_exotica/leptonsPlusJets/ExoDiBosonResonances/tomei/store/user/tomei/BulkG_ZZ_lljj_M2000_G120-JHU/BulkG_ZZ_lljj_M2000_G120-JHU/c8f8ed334db8a7d6f56c62266b1dfa5b/Bulk_AODSIM_11_1_eNC.root",
"/store/group/phys_exotica/leptonsPlusJets/ExoDiBosonResonances/tomei/store/user/tomei/BulkG_ZZ_lljj_M2000_G120-JHU/BulkG_ZZ_lljj_M2000_G120-JHU/c8f8ed334db8a7d6f56c62266b1dfa5b/Bulk_AODSIM_10_1_r0R.root"
    )
                            )
process.maxEvents= cms.untracked.PSet(
    input = cms.untracked.int32(-1),
    )
################################
### NEEDED FOR THE PAT MUONS ###
################################
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")

### Set the global tag from
### https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions#Winter13_2012_A_B_C_D_datasets_r
process.GlobalTag.globaltag = 'START53_V21::All'

process.TransientTrackBuilderESProducer = cms.ESProducer("TransientTrackBuilderESProducer",
                                                         ComponentName = cms.string('TransientTrackBuilder')
                                                         )

### Load the tune P muons
process.load("ExoDiBosonResonances.EDBRMuon.newTuneP_cff")

### Load the PAT muons
process.load("PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi")

### Change their source
process.patMuons.muonSource = "tunePmuons"
process.patMuons.addGenMatch = False
process.patMuons.embedGenMatch = False
process.patMuons.embedCaloMETMuonCorrs = False
process.patMuons.embedTcMETMuonCorrs = False
process.patMuons.userData.userFloats.src = ['muonTrackError']

### Validation
process.validation    = cms.EDAnalyzer("HPTMNewTunePAnalyzer")
process.PATvalidation = cms.EDAnalyzer("HPTMNewTunePPATAnalyzer") # Careful, this spits a lot of printfs.

process.p = cms.Path(process.tunePmuons + process.muonTrackError + process.patMuons + process.validation)
                                              
process.out = cms.OutputModule("PoolOutputModule",
     outputCommands = cms.untracked.vstring('drop *',
                                            'keep *_muonTrackError_*_*',
                                            'keep *_tunePmuons_*_*',
                                            'keep *_muons_*_*',
                                            'keep *_patMuons_*_*'),
                               fileName = cms.untracked.string("BulkG_ZZ_lljj_M2000_G120-JHU.root")
                               )

#process.e = cms.EndPath(process.out)
