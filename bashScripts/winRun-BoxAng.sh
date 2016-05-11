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
foldName=BoxAngChangeTorPct30v4
lwArr=(0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8);
dtArr=(0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025);
nlArr=(1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1);
reArr=(1 1 1 1 1 1 1 1 1 1 1 1);
paArr=(1 1 1 1 1 1 1 1 1 1 1 1);
ang1Arr=(90 90 90 90 90 90 90 90 90 90)
ang2Arr=(90 90 90 90 90 90 90 90 90 90)
#boxangArr=(-20 -25 -30 -35 -40)
boxangArr=(-20 -30 -40)
numPerLay=(5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5)
#0.02 0.002
lazy=(0.01 0.001 0.005 0.0025 0.00125)
mkdir $foldName
for laz in `seq 0 4`; do
	for repeats in `seq 0 3`; do
		for angs in `seq 0 2`; do
		  a=./$foldName/${boxangArr[$angs]}_${lazy[$laz]}_$(date '+%Y%m%d_%H%M%S')
		  mkdir $a
		  cp smarticleMoves.csv $a
		  cp ./checkpoints/${boxangArr[$angs]}.csv ./$a/smarticles.csv
		  cd $a
		  $smartRunFile 0.8 0.00025 6 1 1 0 0 ${boxangArr[$angs]} 5 ${lazy[$laz]};
		  cd ../..;
		  # cp -r ./PostProcess $a
		done;
	done;
done;