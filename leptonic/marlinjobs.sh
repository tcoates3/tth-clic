#!/bin/bash

source /cvmfs/clicdp.cern.ch/iLCSoft/builds/2016-09-27/x86_64-slc6-gcc48-opt/init_ilcsoft.sh

cd /eos/user/y/yixuan/clic/mymarlin/outputs/

for i in $( ls ); do
    inputFile=$i
    echo $i
    inputFilePath=/eos/user/y/yixuan/clic/mymarlin/outputs/$inputFile
    echo $inputFilePath
    cd /eos/user/y/yixuan/clic/mymarlin/output/
    Marlin --global.LCIOInputFiles=$inputFilePath /eos/user/y/yixuan/clic/mymarlin/myLeptonsJetsExtra.xml
    mv testVertexed1.slcio ${inputFile}_extra.slcio
done