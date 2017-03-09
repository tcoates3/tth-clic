#!/bin/bash

#source /eos/user/y/yixuan/clic/mymarlin/init_ilcsoft.sh
#export MARLIN_DLL=/eos/user/y/yixuan/clic/mymarlin/lib/libmymarlin.so
cd /eos/user/y/yixuan/clic/mymarlin/outputs/
n=0
for i in $( ls ); do
    inputFile=$i
    echo $i
    inputFilePath=/eos/user/y/yixuan/clic/mymarlin/outputs/$inputFile
    echo $inputFilePath
    cd /eos/user/y/yixuan/clic/mymarlin/flavourtagged/
    Marlin --global.LCIOInputFiles=$inputFilePath /eos/user/y/yixuan/clic/mymarlin/myFlavourTagging_sl.xml
    #flavortagged=/eos/user/y/yixuan/clic/mymarlin/flavourtagged/flavortagged.slcio
    #Marlin --global.LCIOInputFiles=$flavortagged /eos/user/y/yixuan/clic/mymarlin/MyHiggsHadronic.xml
    mv flavortagged.slcio $n.slcio
    n=$((n+1))
done