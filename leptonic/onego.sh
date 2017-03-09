#!/bin/bash


cd /eos/experiment/clicdp/grid/ilc/prod/clic/1.4tev/tt/SID/REC/00002417/001/

n=4

for i in $( ls -I tt_rec_2417_1000.slcio -I tt_rec_2417_1001.slcio -I tt_rec_2417_1002.slcio -I tt_rec_2417_1003.slcio -I tt_rec_2417_1004.slcio -I tt_rec_2417_1005.slcio -I tt_rec_2417_1006.slcio -I tt_rec_2417_1007.slcio -I tt_rec_2417_1008.slcio -I tt_rec_2417_1009.slcio -I tt_rec_2417_1010.slcio -I tt_rec_2417_1011.slcio -I tt_rec_2417_1012.slcio -I tt_rec_2417_1013.slcio -I tt_rec_2417_1014.slcio -I tt_rec_2417_1015.slcio -I tt_rec_2417_1016.slcio -I tt_rec_2417_1017.slcio ); do
    inputFile=$i
    echo $i
    inputFilePath=/eos/experiment/clicdp/grid/ilc/prod/clic/1.4tev/tt/SID/REC/00002417/001/$inputFile
    echo $inputFilePath
    cd /eos/user/y/yixuan/clic/mymarlin
    java -Xmx1024m -Xms128m -jar /afs/cern.ch/eng/clic/software/lcsim/lcsim-2_5/target/lcsim-2.5-bin.jar -DinputFile=${inputFilePath} /eos/user/y/yixuan/clic/mymarlin/myPreMarlinDrivers.xml
    Marlin /eos/user/y/yixuan/clic/mymarlin/myLeptonsJetsVertexing_sl.xml
    
    echo $n
    cd /eos/user/y/yixuan/clic/mymarlin/TreeMaker/
    ./tth.exe
    cd /eos/user/y/yixuan/clic/mymarlin/TreeMaker/root/bkg/tt/
    mv tt.root $n.root
    n=$((n+1))
done
