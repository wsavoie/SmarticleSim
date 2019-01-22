#!/bin/bash
smartRunFile="D:/ChronoCode/chronoPkgs/SmarticlesBuild/Release/SmarticlesSystem_Dynamic.exe"
dataDir="D:/ChronoCode/chronoPkgs/Smarticles/data";
##mf="B:\SmartSimResults\12_5"
mf=/cygdrive/b/SmartS*/12*
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
# for i in `seq 0 7`; dosg
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
#cd ../../../../../b/SmartS*/12*
# foldName=RemoveBucketAt1sNoGravUShape
# foldName=RemoveBucketAt1s
# foldName=RemoveBucketStraightShape


# lwArr=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9);
# dtVal=0.0002
# nl=(60 58 56 54 52 50 48 45 44)
# npl=(15 14 13 12 11 10 9 8 8)
# re=-1
# reArr=(0 0 0 0 0 0 0 0 0 0 0 0);

# lwArr=(0.4 0.5 0.6 0.8 0.9);
# dtVal=0.0002
# nl=(35 30 22 20 20)
# npl=(14 13 11 10 10)
foldName=SCRemoveBucket
procedure="SCRemoveBucket.txt"

lwArr=(0.4 0.5 0.6 0.7 0.8);
nl=(38 30 28 27 26)
npl=(13 13 11 10 10)

# lwArr=(0.4 0.6 0.7 0.8);
# nl=(35 22 21 20)
# npl=(14 11 11 10)

# lwArr=(0.6);
# nl=(22)
# npl=(11)

dtVal=0.0002
vibAmp=(5)
re=-1
pa=1
ang1=90
ang2=90
vers=$1
windx=0
windy=0
ss=$2
mkdir $foldName
for j in `seq 0 4`; do
	for i in `seq 0 0`; do
	  echo $i
	  mkdir ./$foldName/lw_${lwArr[$j]}_nl_${nl[$j]}_npl_${npl[$j]}_vib_${vibAmp[$i]}
	  a=$foldName/lw_${lwArr[$j]}_nl_${nl[$j]}_npl_${npl[$j]}_vib_${vibAmp[$i]}/v${vers}
	  mkdir $a
	  cp smarticleMoves.csv $a
	  cp ${procedure} $a
	  cd $a
	  $smartRunFile ${lwArr[$j]} ${dtVal} ${nl[$j]} ${re} ${pa} ${ang1} ${ang2} 0 ${npl[$j]} 0.2 ${ss} 0 ${windx} ${windy} 1 ${vibAmp[$i]} ${procedure};
	  cd $mf
	  # cp -r ./PostProcess $a
	done;
done;

powershell.exe -File "sendemail.ps1"