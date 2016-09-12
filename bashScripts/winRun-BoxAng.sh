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

lwArr=(0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8);
dtArr=(.001 .001 .001 .001 .001 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025);
nlArr=(11 11 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1);
reArr=(0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1); 
paArr=(1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1);
ang1Arr=(90 90 90 90 90 90 90 90 90 90 90 90 90 90)
ang2Arr=(90 90 90 90 90 90 90 90 90 90 90 90)
#boxangArr=(-20 -25 -30 -35 -40)
boxangArr=(-30 -30 -30 -30 -30 -30 -30 -30 -30 -30)
numPerLay=(5 5 20 20 20 20 20 20 20 20 5 5 5 5 5 5)
saveFrame=0
#0.02 0.002
changeToStressPerc=(0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3)


########regular run
# foldName='hopper_no_OT_r=pi-4'
# mkdir $foldName
# for laz in `seq 0 0`; do
	# for repeats in `seq 0 0`; do
		# for angs in `seq 0 4`; do
		  # a=./$foldName/${boxangArr[$angs]}_lw${lwArr[$angs]}_a1=${ang1Arr[$angs]}_a2=${ang2Arr[$angs]}_$(date '+%Y%m%d_%H%M%S')
		  # mkdir $a
		  # cp smarticleMoves.csv $a
		  # #cp ./checkpoints/${boxangArr[$angs]}.csv ./$a/smarticles.csv
		  # cd $a
		  # $smartRunFile ${lwArr[$angs]} 0.00025 ${nlArr[$angs]} ${reArr[$angs]} ${paArr[$angs]} ${ang1Arr[$angs]} ${ang2Arr[$angs]} ${boxangArr[$angs]} ${numPerLay[$angs]} ${changeToStressPerc[$angs]} $saveFrame;
		  # cd ../..;
		  # # cp -r ./PostProcess $a
		# done;
	# done;
# done;

##if running many runs with different smarticle.csv files
# foldName='LineGait_r=1.0/hopper_no_OT_r='
 foldName='Opening=1.25w_No_OT/hopper_no_OT_r='
# csvNames=(pi2 pi3 pi4 pi5 pi6 pi7 pi8);
# csvNames=(0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6);
csvNames=(1.0);
# csvNames=(0.2a45 0.4a45 0.6a45 0.8a45 1.0a45 1.2a45 1.4a45 1.6a45);
# csvNames=( 1.0a0 1.0a15 1.0a30 1.0a45 1.0a60 1.0a75 1.0a90);
for i in "${csvNames[@]}"; do	tempFoldname=$foldName$i
	mkdir $tempFoldname
	for angs in `seq 0 1`; do
		a=./$tempFoldname/${boxangArr[$angs]}_lw${lwArr[$angs]}_a1=${ang1Arr[$angs]}_a2=${ang2Arr[$angs]}_$(date '+%Y%m%d_%H%M%S')
		mkdir $a
		cp ./smarticleMovesCSV/$i".csv" ./$a/smarticleMoves.csv
		echo $i".csv"
		# cp ./checkpoints/${boxangArr[$angs]}.csv ./$a/smarticles.csv
		cd $a
		$smartRunFile ${lwArr[$angs]} 0.001 ${nlArr[$angs]} ${reArr[$angs]} ${paArr[$angs]} ${ang1Arr[$angs]} ${ang2Arr[$angs]} ${boxangArr[$angs]} ${numPerLay[$angs]} ${changeToStressPerc[$angs]} $saveFrame;
	cd $runDir
	done;
done;