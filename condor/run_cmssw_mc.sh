#!/bin/bash

cd /afs/cern.ch/user/a/alshevel/private/CMSSW_15_0_6/src/Demo/HLTScouting
cmsenv

input_file="inFile=$1"
output_file="outFile=$2"
echo "Starting cmsRun with input file $input_file"

cmsRun python/mc_cfg.py $input_file $output_file
