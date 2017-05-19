#!/bin/bash

if [hostname=='phys42232.physics.gatech.edu']; then
smartRunFile="/home/ws/SmartSim/SmarticlesBuild/SmarticlesSystem_Dynamic"
dataDir="/home/ws/SmartSim/Smarticles/data"
runDir="/home/ws/SmartSim/Results"
sys="/home/ws/SmartSim/Results"
csvFold="/home/ws/SmartSim/Results"
else
smartRunFile="D:/ChronoCode/chronoPkgs/SmarticlesBuild/Release/SmarticlesSystem_Dynamic.exe"
dataDir="D:/ChronoCode/chronoPkgs/Smarticles/data"
runDir="A:/SmarticleRun"
sys="/cygdrive/a/SmarticleRun"
csvFold="/cygdrive/d/SimResults/Chrono/SmarticleU/tests"
fi


if diff ../Smarticles/data $dataDir;then #if no differences
echo "found directory";
else
echo "copying directory";
rm -r data
cp -R $dataDir data
fi



lwArr=(0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8 0.8);
dtArr=(.001 .001 .001 .001 .001 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025);
nlArr=(1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1);
reArr=(0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1); 
paArr=(1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1);
ang1Arr=(90 90 90 90 90 90 90 90 90 90 90 90 90 90)
ang2Arr=(90 90 90 90 90 90 90 90 90 90 90 90)
#boxangArr=(-20 -25 -30 -35 -40)
boxangArr=(-30 -30 -30 -30 -30 -30 -30 -30 -30 -30)
numPerLay=(2 5 20 20 20 20 20 20 20 20 5 5 5 5 5 5)
saveFrame=1
#0.02 0.002
changeToStressPerc=(0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)
#+x +y -x -y
dirs=(0 1 2 3)
inactiveParticle=0;
#######regular run
foldName=Amoeba_Ring_t_60_Square_pi_m_.055_d_$inactiveParticle
# foldName='Amoeba_Ring_t_60_Square_s_pi_d_1'
mkdir $sys/$foldName
for fric in `seq 4 4`; do
	for deadDir in `seq 1 3`; do
		# for repeats in `seq 0 20`; do
			for repeats in `seq 0 0`; do
			for robs in `seq 1 1`; do
			  a=$sys/$foldName/f_${changeToStressPerc[$fric]}_rob_${numPerLay[$robs]}_v_${repeats}_dir_${deadDir}
			  mkdir $a
			  cp $csvFold/smarticleMoves.csv $a/smarticleMoves.csv
			  cd $a
			  # # $smartRunFile ${lwArr[$angs]} 0.00025 ${nlArr[$angs]} ${reArr[$angs]} ${paArr[$angs]} ${ang1Arr[$angs]} ${ang2Arr[$angs]} ${boxangArr[$angs]} ${numPerLay[$angs]} ${changeToStressPerc[$angs]} $saveFrame;
			  $smartRunFile 0.8 0.001 1 0 1 0 0 0 ${numPerLay[$robs]} ${changeToStressPerc[$fric]} $saveFrame $deadDir 801 1 $inactiveParticle
			  cd ..;
			  # cp -r ./PostProcess $a
			  done;
		done;
	done;
done;

##if running many runs with different smarticle.csv files
# foldName='LineGait_r=1.0/hopper_no_OT_r='

# /cygdrive/d/SimResults/Chrono/SmarticleU/tests ./bash_ringSmarticle.sh 

