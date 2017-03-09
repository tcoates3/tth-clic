#!/bin/bash

source /cvmfs/clicdp.cern.ch/iLCSoft/builds/2016-09-27/x86_64-slc6-gcc48-opt/init_ilcsoft.sh

export MARLIN_DLL=/eos/user/y/yixuan/clic/mymarlin/lib/libmymarlin.so

cd /eos/user/y/yixuan/clic/mymarlin/outputs/
n=0
for i in $( ls ); do
    echo $n
    inputFile=$i
    echo $i
    inputFilePath=/eos/user/y/yixuan/clic/mymarlin/outputs/$inputFile
    echo $inputFilePath
    cd /eos/user/y/yixuan/clic/mymarlin/roots/
    Marlin --global.LCIOInputFiles=$inputFilePath /eos/user/y/yixuan/clic/mymarlin/MyHiggsHadronic.xml
    mv output.root $n.root
    n=$((n+1))
done