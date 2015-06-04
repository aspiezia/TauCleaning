export SCRAM_ARCH=slc5_amd64_gcc462
wget -O installCMGTools.py "http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/CMG/CMGTools/Common/bin/installCMGTools.py?view=co"
python installCMGTools.py -r CMSSW_5_3_9 --recipe=rel_V5_15_0
cd CMGTools/CMSSW_5_3_9/src
cvs up -r V00-02-06 CalibCalorimetry/EcalLaserCorrection/src/EcalLaserDbService.cc
cmsenv
chmod +x CMGTools/Production/scripts/*.py
chmod +x CMGTools/Common/scripts/*.py
cd EGamma/EGammaAnalysisTools/data
cat download.url | xargs wget
cd -
cvs up -r 1.14 DataFormats/CaloRecHit/src/classes_def.xml 
rm -r CMGTools/RootTools 
rm -r CMGTools/ZJetsTutorial 
rm -r CMGTools/H2TauTau 
rm -r CMGTools/Utilities 
rm -r TauAnalysis/SVFitStandAlone 
rm -r TauAnalysis/CandidateTools
cvs co -r V00-09-03 -d SHarper/HEEPAnalyzer UserCode/SHarper/HEEPAnalyzer
cvs co -r V00-02-03  -d TSWilliams/BstdZeeTools  UserCode/TSWilliams/BstdZee/BstdZeeTools
cvs co -A -d KStenson/TrackingFilters UserCode/KStenson/TrackingFilters

git clone git://github.com/cms-edbr/ExoDiBosonResonances -b edbr_tautau
mv ExoDiBosonResonances/PATtupleProduction/data/Buildfile_KStenson KStenson/TrackingFilters/BuildFile.xml
rm -rf KStenson/TrackingFilters/src
mv KStenson/TrackingFilters/plugins KStenson/TrackingFilters/src 
cvs co -r V00-02-03s TauAnalysis/CandidateTools

cd ../../../

cmsrel CMSSW_5_3_13
cd CMSSW_5_3_13/src/
cmsenv
git cms-addpkg PhysicsTools/IsolationAlgos
git cherry-pick 1bccce7edaeea05a60ef1b0db7e44652321d027c
git cms-merge-topic -u veelken:from-CMSSW_5_3_13
git cms-merge-topic -u veelken:5_3_X_boostedTaus
git cms-merge-topic -u cms-tau-pog:CMSSW_5_3_X_boostedTaus

cp -r ../../CMGTools/CMSSW_5_3_9/src/AnalysisDataFormats .
cp -r ../../CMGTools/CMSSW_5_3_9/src/CMGTools .
cp -r ../../CMGTools/CMSSW_5_3_9/src/CalibCalorimetry .
cp -r ../../CMGTools/CMSSW_5_3_9/src/CondFormats .
cp -r ../../CMGTools/CMSSW_5_3_9/src/DPGAnalysis .
cp -r ../../CMGTools/CMSSW_5_3_9/src/DataFormats/CaloRecHit DataFormats/
cp -r ../../CMGTools/CMSSW_5_3_9/src/DataFormats/StdDictionaries DataFormats/
cp -r ../../CMGTools/CMSSW_5_3_9/src/EGamma .
cp -r ../../CMGTools/CMSSW_5_3_9/src/EgammaAnalysis .
cp -r ../../CMGTools/CMSSW_5_3_9/src/ExoDiBosonResonances .
cp -r ../../CMGTools/CMSSW_5_3_9/src/FWCore/GuiBrowsers FWCore/
cp -r ../../CMGTools/CMSSW_5_3_9/src/KStenson .
cp -r ../../CMGTools/CMSSW_5_3_9/src/Muon .
cp -r ../../CMGTools/CMSSW_5_3_9/src/PhysicsTools/PatUtils PhysicsTools/
cp -r ../../CMGTools/CMSSW_5_3_9/src/RecoEgamma .
cp -r ../../CMGTools/CMSSW_5_3_9/src/RecoLuminosity .
cp -r ../../CMGTools/CMSSW_5_3_9/src/RecoMET .
cp -r ../../CMGTools/CMSSW_5_3_9/src/RecoParticleFlow/PFProducer RecoParticleFlow/
cp -r ../../CMGTools/CMSSW_5_3_9/src/RecoParticleFlow/PostProcessing RecoParticleFlow/
cp -r ../../CMGTools/CMSSW_5_3_9/src/SHarper .
cp -r ../../CMGTools/CMSSW_5_3_9/src/TSWilliams .
cp -r ../../CMGTools/CMSSW_5_3_9/src/TauAnalysis .

cp /afs/cern.ch/user/v/veelken/public/boostedTaus/2013Nov04/BoostedTauSeedsProducer.cc RecoJets/JetProducers/plugins/
cp /afs/cern.ch/user/v/veelken/public/boostedTaus/2013Dec02/VirtualJetProducer.cc RecoJets/JetProducers/plugins/
cp /afs/cern.ch/user/v/veelken/public/boostedTaus/2013Dec02/FastjetJetProducer.cc RecoJets/JetProducers/plugins/
cp /afs/cern.ch/work/c/caber/public/ForNewPatTuplizer08Dec2013/MuonIso/boostedMuonIsolation_cff.py PhysicsTools/IsolationAlgos/python/
cp /afs/cern.ch/work/c/caber/public/ForNewPatTuplizer08Dec2013/ElectronIso/boostedElectronIsolation_cff.py PhysicsTools/IsolationAlgos/python/
cp /afs/cern.ch/work/c/caber/public/ForNewPatTuplizer08Dec2013/PATTaus/PATTaus_cff.py CMGTools/Common/python/PAT/
cp /afs/cern.ch/work/a/aspiezia/public/PattuplesTau/PATJets_cff.py CMGTools/Common/python/PAT/
cp /afs/cern.ch/work/a/aspiezia/public/PattuplesTau/boostedHPSPFTaus_cff.py RecoTauTag/Configuration/python/
cp EGamma/EGammaAnalysisTools/data/Electrons_BDTG_NonTrigV0_Cat* ExoDiBosonResonances/PATtupleProduction/prod/submitJobsWithCrab/

emacs CMGTools/Common/src/TauFactory.cc &
# change every PFCandidateRef -> PFCandidatePtr

emacs AnalysisDataFormats/CMGTools/interface/Tau.h &
#edit (because we have more than 45 DMs): #define NCMGTAUIDS 45 -> #define NCMGTAUIDS 200

scram b -j 8


