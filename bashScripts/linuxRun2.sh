#!/bin/bash
#Notes about script:
#this script can be run from anywhere, it will save the output into the same folder
#it needs to have a smarticlesMoves.csv file in the same directory as this file.
#to create the csv, run matlab file in smarticles folder matlabScripts/GenerateSmarticleCsv.m
#and direct the dialog window to the directory containing this file.
#finally you must fill in the else part of the if statement below
if hostname=='phys42232.physics.gatech.edu'; then
smartRunFile="/home/wsavoie/Documents/ChronoPkgs/SmarticlesBuild/SmarticlesSystem_Dynamic"
dataDir="/home/wsavoie/Documents/ChronoPkgs/Smarticles/data";
runDir="/home/wsavoie/Documents/SmarticleResults"
else
smartRunFile="/path/to/SmarticleSystem_Dynamic/here"
dataDir="/your/data/directory/here"
runDir="directory you want to run from here"
fi

if diff ../Smarticles/data/ $dataDir;then #if no differences
echo "found directory";
else
echo "copying directory";
rm -r data
cp -R $dataDir ../Smarticles/data/
fi


# #read from params file
# lw=$(grep -E "lw" params | awk '{print $2}');
# dt=$(grep -E "dt" params | awk '{print $2}');
# nl=$(grep -E "numlayers" params | awk '{print $2}');
# re=$(grep -E "read" params | awk '{print $2}');
# pa=$(grep -E "pa" params | awk '{print $2}');
# echo "run vars!: $lw $dt $nl $re $pa"

cd $runDir
foldName=Active20%
lwArr=(0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8);
dtArr=(0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025);
nlArr=(1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1);
reArr=(0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1); 
paArr=(.2 .2 .2 .2 .2 1 1 1 1 1 1 1 1 1 1);
ang1Arr=(0 30 60 90 120 90 90 90 90 90)
ang2Arr=(0 30 60 90 120 90 90 90 90 90)
#boxangArr=(-20 -25 -30 -35 -40)
boxangArr=(-45 -45 -45 -45 -45)
numPerLay=(35 35 35 35 35 5 5 5 5 5 5 5 5 5 5 5 5 5)
saveFrame=1
#0.02 0.002
changeToStressPerc=(0.3 0.001 0.005 0.0025 0.00125)
mkdir $foldName
for laz in `seq 0 0`; do
	for repeats in `seq 0 0`; do
		for angs in `seq 0 4`; do
		  a=./$foldName/${boxangArr[$angs]}_lw${lwArr[$angs]}_a1=${ang1Arr[$angs]}_a2=${ang2Arr[$angs]}_$(date '+%Y%m%d_%H%M%S')
		  mkdir $a
		  cp smarticleMoves.csv $a
		  cp ./checkpoints/${boxangArr[$angs]}.csv ./$a/smarticles.csv
		  cd $a
		  $smartRunFile ${lwArr[$angs]} 0.00025 ${nlArr[$angs]} ${reArr[$angs]} ${paArr[$angs]} ${ang1Arr[$angs]} ${ang2Arr[$angs]} ${boxangArr[$angs]} ${numPerLay[$angs]} ${changeToStressPerc[$angs]} $saveFrame;
		  cd ../..;
		  # cp -r ./PostProcess $a
		done;
	done;
done;



echo "\\n blah \\n";
echo ${lwArr[*]};
# $smartRunFile $lw $dt $nl $re $pa
