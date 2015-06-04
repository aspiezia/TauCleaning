import FWCore.ParameterSet.Config as cms

BTagWeightsKinFitMu = cms.EDProducer("BTagWeightProducerWmunuDiJet",
                                    src = cms.InputTag("HLTWeightsKinFitMu") ,
									BTagJets=cms.InputTag("jetAK5"),
									EffmapFilename = cms.FileInPath("ExoDiBosonResonances/EDBRCommon/data/BtaggingEffs/WW_TuneZ2star_8TeV_pythia6_tauola_AK5PF_CSVM_bTaggingEfficiencyMap.root"),
									scale_b = cms.double(0),#0 noscale, 1 scale up, -1 scale down
									scale_light = cms.double(0),#0 noscale, 1 scale up, -1 scale down
                                    )

    
    
BTagWeightsMu = BTagWeightsKinFitMu.clone()
BTagWeightsMu.src = "HLTWeightsMu"
    
    
BTagWeightsMergedMu =cms.EDProducer("BTagWeightProducerWmunuSingleJet",
                                   src = cms.InputTag( "HLTWeightsMergedMu") ,
								   BTagJets=cms.InputTag("jetAK5"),
								   EffmapFilename = cms.FileInPath("ExoDiBosonResonances/EDBRCommon/data/BtaggingEffs/WW_TuneZ2star_8TeV_pythia6_tauola_AK5PF_CSVM_bTaggingEfficiencyMap.root"),
								   scale_b = cms.double(0),#0 noscale, 1 scale up, -1 scale down
								   scale_light = cms.double(0),#0 noscale, 1 scale up, -1 scale down
                                   )


