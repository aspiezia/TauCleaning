import FWCore.ParameterSet.Config as cms


### values used for 2012 analysis, to be used with HLTWeightProducer2012
dummyHLTMuEffWeights =  cms.PSet( lumi = cms.double(1.),#dummy
                               lep1 = cms.VPSet( cms.PSet(bin = cms.vdouble(0.0,3.0,0.0,7000. ,1.0))),
                               lep2 = cms.VPSet( cms.PSet(bin = cms.vdouble(0.0,3.0,0.0,7000. ,1.0))),
                               hlt  = cms.VPSet( cms.PSet(bin = cms.vdouble(0.0,999.0,0.0,999. ,1.0))),
                               )

### put the data/MC SF for the lepton ID and the data/MC SF for the trigger
### (this is a change for V2 of the analysis)
globMuIdSF_Run2012= cms.VPSet(
    cms.PSet(bin = cms.vdouble(0.0, 0.9, 20.000, 40.000, 0.994607)),
    cms.PSet(bin = cms.vdouble(0.9, 1.2, 20.000, 40.000, 0.994208)),
    cms.PSet(bin = cms.vdouble(1.2, 2.1, 20.000, 40.000, 0.9962)),
    cms.PSet(bin = cms.vdouble(2.1, 2.4, 20.000, 40.000, 0.986778)),
    cms.PSet(bin = cms.vdouble(0.0, 0.9, 40.000, 60.000, 0.992214)),
    cms.PSet(bin = cms.vdouble(0.9, 1.2, 40.000, 60.000, 0.993411)),
    cms.PSet(bin = cms.vdouble(1.2, 2.1, 40.000, 60.000, 0.994307)),
    cms.PSet(bin = cms.vdouble(2.1, 2.4, 40.000, 60.000, 0.980866)),
    cms.PSet(bin = cms.vdouble(0.0, 0.9, 60.000, 80.000, 0.999195)),
    cms.PSet(bin = cms.vdouble(0.9, 1.2, 60.000, 80.000, 0.998496)),
    cms.PSet(bin = cms.vdouble(1.2, 2.1, 60.000, 80.000, 0.999985)),
    cms.PSet(bin = cms.vdouble(2.1, 2.4, 60.000, 80.000, 0.978304)),
    cms.PSet(bin = cms.vdouble(0.0, 0.9, 80.000, 100.000, 1.00498)),
    cms.PSet(bin = cms.vdouble(0.9, 1.2, 80.000, 100.000, 1.00308)),
    cms.PSet(bin = cms.vdouble(1.2, 2.1, 80.000, 100.000, 1.00895)),
    cms.PSet(bin = cms.vdouble(2.1, 2.4, 80.000, 100.000, 0.968944)),
    cms.PSet(bin = cms.vdouble(0.0, 0.9, 100.000, 5000.000, 1.00418)),
    cms.PSet(bin = cms.vdouble(0.9, 1.2, 100.000, 5000.000, 1.00288)),
    cms.PSet(bin = cms.vdouble(1.2, 2.1, 100.000, 5000.000, 1.01453)),
    cms.PSet(bin = cms.vdouble(2.1, 2.4, 100.000, 5000.000, 1.02846)),
)

tkMuIdSF_Run2012= cms.VPSet(
    cms.PSet(bin = cms.vdouble(0.0, 0.9, 20.000, 40.000, 0.996303)),
    cms.PSet(bin = cms.vdouble(0.9, 1.2, 20.000, 40.000, 0.995903)),
    cms.PSet(bin = cms.vdouble(1.2, 2.1, 20.000, 40.000, 0.996)),
    cms.PSet(bin = cms.vdouble(2.1, 2.4, 20.000, 40.000, 0.988059)),
    cms.PSet(bin = cms.vdouble(0.0, 0.9, 40.000, 60.000, 0.994707)),
    cms.PSet(bin = cms.vdouble(0.9, 1.2, 40.000, 60.000, 0.995604)),
    cms.PSet(bin = cms.vdouble(1.2, 2.1, 40.000, 60.000, 0.994606)),
    cms.PSet(bin = cms.vdouble(2.1, 2.4, 40.000, 60.000, 0.983822)),
    cms.PSet(bin = cms.vdouble(0.0, 0.9, 60.000, 80.000, 1.00029)),
    cms.PSet(bin = cms.vdouble(0.9, 1.2, 60.000, 80.000, 0.999892)),
    cms.PSet(bin = cms.vdouble(1.2, 2.1, 60.000, 80.000, 0.999985)),
    cms.PSet(bin = cms.vdouble(2.1, 2.4, 60.000, 80.000, 0.978107)),
    cms.PSet(bin = cms.vdouble(0.0, 0.9, 80.000, 100.000, 1.00698)),
    cms.PSet(bin = cms.vdouble(0.9, 1.2, 80.000, 100.000, 1.00308)),
    cms.PSet(bin = cms.vdouble(1.2, 2.1, 80.000, 100.000, 1.00417)),
    cms.PSet(bin = cms.vdouble(2.1, 2.4, 80.000, 100.000, 0.969437)),
    cms.PSet(bin = cms.vdouble(0.0, 0.9, 100.000, 5000.000, 1.00488)),
    cms.PSet(bin = cms.vdouble(0.9, 1.2, 100.000, 5000.000, 1.00209)),
    cms.PSet(bin = cms.vdouble(1.2, 2.1, 100.000, 5000.000, 1.01403)),
    cms.PSet(bin = cms.vdouble(2.1, 2.4, 100.000, 5000.000, 0.993478)),
)

doubleMuTrigEff_Run2012= cms.VPSet(
    cms.PSet(bin = cms.vdouble(0.0,0.9, 0.0,0.9 ,0.983689)),
    cms.PSet(bin = cms.vdouble(0.0,0.9, 0.9,1.2 ,0.988776)),
    cms.PSet(bin = cms.vdouble(0.0,0.9, 1.2,2.1 ,1.01425)),
    cms.PSet(bin = cms.vdouble(0.0,0.9, 2.1,2.4 ,0.925808)),

    cms.PSet(bin = cms.vdouble(0.9, 1.2, 0.0,0.9 ,0.971349)),
    cms.PSet(bin = cms.vdouble(0.9, 1.2, 0.9,1.2 ,1.03600)),
    cms.PSet(bin = cms.vdouble(0.9, 1.2, 1.2,2.1 ,1.00243)),
    cms.PSet(bin = cms.vdouble(0.9, 1.2, 2.1,2.4 ,0.981699)),

    cms.PSet(bin = cms.vdouble(1.2, 2.1, 0.0,0.9 ,1.00282)),
    cms.PSet(bin = cms.vdouble(1.2, 2.1, 0.9,1.2 ,1.01962)),
    cms.PSet(bin = cms.vdouble(1.2, 2.1, 1.2,2.1 ,1.00827)),
    cms.PSet(bin = cms.vdouble(1.2, 2.1, 2.1,2.4 ,0.983534)),

    cms.PSet(bin = cms.vdouble(2.1, 2.4, 0.0,0.9 ,1.03590)),
    cms.PSet(bin = cms.vdouble(2.1, 2.4, 0.9,1.2 ,1.08973)),
    cms.PSet(bin = cms.vdouble(2.1, 2.4, 1.2,2.1 ,0.971467)),
    cms.PSet(bin = cms.vdouble(2.1, 2.4, 2.1,2.4 ,0.981039))
    )

#### to be expanded according to the modifications of the HLT triggers that we use
#### (a set of efficieincies for each of them, function of etan and pt, calculated with T&P)

HLTWeightsKinFitMu = cms.EDProducer("HLTWeightProducer2012Mu",
                                    src = cms.InputTag("cmgDiMuonDiJetKinFitEDBR") ,
                                    isMuChannel = cms.bool(True),
                                    ranges = cms.VPSet(
                                          cms.PSet( lumi = cms.double(19747.0),#Run2012
                                                    lep1 = globMuIdSF_Run2012,
                                                    lep2 = tkMuIdSF_Run2012,
                                                    hlt = doubleMuTrigEff_Run2012 
                                                    )
                                          #  dummyHLTMuEffWeights
                                                                
                                  
                                          )# end ranges
                                    )       
HLTWeightsMu = HLTWeightsKinFitMu.clone()
HLTWeightsMu.src = "cmgDiMuonDiJetEDBR"


HLTWeightsMergedMu =cms.EDProducer("HLTWeightProducer2012MuVJet",
                                   src = cms.InputTag( "cmgDiMuonVJetEDBR") ,
                                   isMuChannel = cms.bool(True),
                                   ranges = cms.VPSet(
                                            cms.PSet(lumi = cms.double(19747.0),
                                            lep1 = globMuIdSF_Run2012,
                                            lep2 = tkMuIdSF_Run2012,
                                            hlt = doubleMuTrigEff_Run2012
                                                     )
                                            ) # end ranges
                                   ### dummyHLTMuEffWeights
                                   )



###
#######################################################
###
### values used for 2011 analysis EXO-11-102, to be used with HLTWeightProducerMu
doubleMuLeg1Eff_Run2011A= cms.VPSet(
    cms.PSet(bin = cms.vdouble(0.0,0.8,10.,20. ,0.975)),
    cms.PSet(bin = cms.vdouble(0.0,0.8,20.,40. ,0.975)),
    cms.PSet(bin = cms.vdouble(0.0,0.8,40.,150.,0.973)),
    cms.PSet(bin = cms.vdouble(0.8,2.1,10.,20. ,0.955)),
    cms.PSet(bin = cms.vdouble(0.8,2.1,20.,40. ,0.950)),
    cms.PSet(bin = cms.vdouble(0.8,2.1,40.,150.,0.946)),
    cms.PSet(bin = cms.vdouble(2.1,2.4,10.,20. ,0.958)),
    cms.PSet(bin = cms.vdouble(2.1,2.4,20.,40. ,0.922)),
    cms.PSet(bin = cms.vdouble(2.1,2.4,40.,150.,0.908))
    )                                                                                

doubleMuLeg2Eff_Run2011A= cms.VPSet(
    cms.PSet(bin = cms.vdouble(0.0,0.8,10.,20. ,0.975)),
    cms.PSet(bin = cms.vdouble(0.0,0.8,20.,40. ,0.975)),
    cms.PSet(bin = cms.vdouble(0.0,0.8,40.,150.,0.973)),
    cms.PSet(bin = cms.vdouble(0.8,2.1,10.,20. ,0.955)),
    cms.PSet(bin = cms.vdouble(0.8,2.1,20.,40. ,0.950)),
    cms.PSet(bin = cms.vdouble(0.8,2.1,40.,150.,0.946)),
    cms.PSet(bin = cms.vdouble(2.1,2.4,10.,20. ,0.958)),
    cms.PSet(bin = cms.vdouble(2.1,2.4,20.,40. ,0.922)),
    cms.PSet(bin = cms.vdouble(2.1,2.4,40.,150.,0.908))
    )                                                                                
singleMuEff_Run2011A= cms.VPSet(
    cms.PSet(bin = cms.vdouble(0.0,0.8,10.,20. ,0.975)),
    cms.PSet(bin = cms.vdouble(0.0,0.8,20.,40. ,0.975)),
    cms.PSet(bin = cms.vdouble(0.0,0.8,40.,150.,0.973)),
    cms.PSet(bin = cms.vdouble(0.8,2.1,10.,20. ,0.955)),
    cms.PSet(bin = cms.vdouble(0.8,2.1,20.,40. ,0.950)),
    cms.PSet(bin = cms.vdouble(0.8,2.1,40.,150.,0.946)),
    cms.PSet(bin = cms.vdouble(2.1,2.4,10.,20. ,0.958)),
    cms.PSet(bin = cms.vdouble(2.1,2.4,20.,40. ,0.922)),
    cms.PSet(bin = cms.vdouble(2.1,2.4,40.,150.,0.908))
    )    
myOLDDummyHLTMuEffWeights =  cms.PSet( lumi = cms.double(1.),#dummy
                               double1 = cms.VPSet( cms.PSet(bin = cms.vdouble(0.0,3.0,0.0,7000. ,1.0))),
                               double2 = cms.VPSet( cms.PSet(bin = cms.vdouble(0.0,3.0,0.0,7000. ,1.0))),
                               single  = cms.VPSet( cms.PSet(bin = cms.vdouble(0.0,3.0,0.0,7000. ,1.0))),
                               )
HLTWeightsKinFitMuOLD = cms.EDProducer("HLTWeightProducerMu",
                                  src = cms.InputTag("cmgDiMuonDiJetKinFitEDBR") ,
                                  ranges = cms.VPSet(
    #  cms.PSet( lumi = cms.double(217.),#Run2011A1
    #                                                                 double1 = doubleMuLeg1Eff_Run2011A,
    #                                                                 double2 = doubleMuLeg2Eff_Run2011A,
    #                                                                 single = singleMuEff_Run2011A
    #                                                                 ),
    #                                                       cms.PSet( lumi = cms.double(920.),#Run2011A2
    #                                                                 double1 = doubleMuLeg1Eff_Run2011A,
    #                                                                 double2 = doubleMuLeg2Eff_Run2011A,
    #                                                                 single = singleMuEff_Run2011A
    #                                                                 ),
                                                      myOLDDummyHLTMuEffWeights
                                                                
                                  
                                                      )# end ranges
                                                      )          


 
