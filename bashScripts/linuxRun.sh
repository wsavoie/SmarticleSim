#!/bin/bash
smartRunFile="/home/wsavoie/Documents/ChronoPkgs/SmarticlesBuild/SmarticlesSystem_Dynamic"
dataDir="/home/wsavoie/Documents/ChronoPkgs/Smarticles/data";

if diff ../Smarticles/data $dataDir;then #if no differences
echo "found directory";
else
echo "copying directory";
rm -r data
cp -R $dataDir data
fi


#read from params file
lw=$(grep -E "lw" params | awk '{print $2}');
dt=$(grep -E "dt" params | awk '{print $2}');
nl=$(grep -E "numlayers" params | awk '{print $2}');
re=$(grep -E "read" params | awk '{print $2}');
pa=$(grep -E "pa" params | awk '{print $2}');
echo "run vars!: $lw $dt $nl $re $pa"

lwArr=(0.1 0.3 0.5 0.7);
dtArr=(0.00025 0.00025 0.00025 0.00025);
nlArr=(40 40 40 40);
reArr=(0 0 0 0);
paArr=(0.0 0.0 0.0 0.0);

# for i in `seq 0 4`; do
for i in `seq 2 3`; do
  echo $i
  $smartRunFile ${lwArr[$i]} ${dtArr[$i]} ${nlArr[$i]} ${reArr[$i]} ${paArr[$i]};
  mkdir ./volFracFilling/${lwArr[$i]}/
  cp -R PostProcess ./volFracFilling/${lwArr[$i]}/
done
#
# for i in `seq 0 4`; do
#   echo $i
#   $smartRunFile ${lwArr[$i]} ${dtArr[$i]} ${nlArr[$i]} ${reArr[$i]} ${paArr[$i]};
#   mkdir ./MixedSmartWithOT10/r2/${paArr[$i]}/
#   cp -R PostProcess ./MixedSmartWithOT10/r2/${paArr[$i]}
# done

echo "\\n blah \\n";
echo ${lwArr[*]};
# $smartRunFile $lw $dt $nl $re $pa
