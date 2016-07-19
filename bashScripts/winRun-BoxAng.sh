#!/bin/bash
smartRunFile="D:/ChronoCode/chronoPkgs/SmarticlesBuild/Release/SmarticlesSystem_Dynamic.exe"
dataDir="D:/ChronoCode/chronoPkgs/Smarticles/data"
runDir="A:/SmarticleRun"
if diff ../Smarticles/data $dataDir;then #if no differences
echo "found directory";
else
echo "copying directory";
rm -r data
cp -R $dataDir data
fi


#read from params file
# lw=$(grep -E "lw" params | awk '{print $2}');
# dt=$(grep -E "dt" params | awk '{print $2}');
# nl=$(grep -E "numlayers" params | awk '{print $2}');
# re=$(grep -E "read" params | awk '{print $2}');
# pa=$(grep -E "pa" params | awk '{print $2}');
echo "run vars!: $lw $dt $nl $re $pa"

# lwArr=(0.1 0.2 0.3 0.5 0.7 0.8 0.9 1.0);
# dtArr=(0.0005 0.0005 0.0005 0.0005 0.0005 0.0005 0.0005 0.0005);
# nlArr=(100 90 80 70 60 50 40 40);
# reArr=(0 0 0 0 0);
# paArr=(0 0 0 0 0);
# for i in `seq 0 7`; do
  # echo $i
  # a=./volFrac2/${lwArr[$i]}-$(date '+%Y%m%d-%H%M%S')
  # mkdir $a
  # # cd $a
  # $smartRunFile ${lwArr[$i]} ${dtArr[$i]} ${nlArr[$i]} ${reArr[$i]} ${paArr[$i]};
  # cp -r ./PostProcess $a
  # # cd ../..;
# done;
# echo "\\n blah \\n";
# echo ${lwArr[*]};
# # $smartRunFile $lw $dt $nl $re $pa
cd $runDir
foldName=TestActive
lwArr=(0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8);
dtArr=(0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025);
nlArr=(1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1);
reArr=(0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1); 
paArr=(1 .98 .75 .75 .75 .75 .75 .75 .75 .75 1 1 1 1 1 1 1 1 1 1);
ang1Arr=(120 120 30 45 60 75 90 105 120 90 90 90 90 90)
ang2Arr=(-120 -120 30 45 60 75 90 105 120 90 90 90)
#boxangArr=(-20 -25 -30 -35 -40)
boxangArr=(-40 -40 -45 -45 -45 -45 -45 -45 -45 -45)
numPerLay=(100 100 35 35 35 35 35 35 35 35 5 5 5 5 5 5)
saveFrame=0
#0.02 0.002
changeToStressPerc=(0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3)
mkdir $foldName
for laz in `seq 0 0`; do
	for repeats in `seq 0 0`; do
		for angs in `seq 0 0`; do
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