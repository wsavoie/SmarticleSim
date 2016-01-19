#!/bin/bash
smartRunFile="D:/ChronoCode/chronoPkgs/SmarticlesBuild/Release/SmarticlesSystem_Dynamic.exe"
dataDir="D:/ChronoCode/chronoPkgs/Smarticles/data";

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

foldName=StressArmLwSweep
lwArr=(0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1);
dtArr=(0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025 0.00025);
nlArr=(40 40 40 40 40 40 40 40 40);
reArr=(0 0 0 0 0 0 0 0 0 0 0 0);
paArr=(1 1 1 1 1 1 1 1 1 0 0 0);
ang1Arr=(90 90 90 90 90 90 90 90 90 90)
ang2Arr=(90 90 90 90 90 90 90 90 90 90)
mkdir $foldName
for i in `seq 0 8`; do
  echo $i
  a=./$foldName/${lwArr[$i]}-${nlArr[$i]}-${ang1Arr[$i]}-$(date '+%Y%m%d-%H%M%S')
  mkdir $a
  cp smarticleMoves.csv $a
  cd $a
  $smartRunFile ${lwArr[$i]} ${dtArr[$i]} ${nlArr[$i]} ${reArr[$i]} ${paArr[$i]} ${ang1Arr[$i]} ${ang2Arr[$i]};
  cd ../..;
  # cp -r ./PostProcess $a
done;