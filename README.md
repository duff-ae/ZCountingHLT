# create a project area
# Any pre-release of (CMSSW_15_0_X) >CMSSW_15_0_6 can be used
cmsrel CMSSW_15_0_6 

cd CMSSW_15_0_6/src
cmsenv

mkdir Demo
cd Demo

git clone https://github.com/duff-ae/ZCountingHLT HLTScouting

cd HLTScouting
scram b

# init proxy
voms-proxy-init --voms cms

# process test file
cmsRun python/zcounting_cfg.py inFile=root://cms-xrd-global.cern.ch//store/data/Run2025C/ScoutingPFRun3/HLTSCOUT/v1/000/392/542/00000/5a05e12c-b12d-43b2-9869-859cb8ade2ec.root outFile=test_physics.root

