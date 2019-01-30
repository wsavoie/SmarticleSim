#!/bin/bash
smartRunFile="/home/wsavoie/Documents/ChronoPkgs/SmarticlesBuild/SmarticlesSystem_Dynamic"
dataDir="/home/wsavoie/Documents/ChronoPkgs/Smarticles/data";

if diff ../Smarticles/data/ $dataDir;then #if no differences
echo "found directory";
else
echo "copying directory";
rm -r data
cp -R $dataDir ../Smarticles/data/
fi


#read from params file
lw=$(grep -E "lw" params | awk '{print $2}');
dt=$(grep -E "dt" params | awk '{print $2}');
nl=$(grep -E "numlayers" params | awk '{print $2}');
re=$(grep -E "read" params | awk '{print $2}');
pa=$(grep -E "pa" params | awk '{print $2}');
echo "run vars!: $lw $dt $nl $re $pa"

lwArr=(0.1 0.3 0.5 0.7 0.9);
dtArr=(0.0005 0.0005 0.0005 0.0005 0.0005);
nlArr=(60 55 45 40 30);
reArr=(0 0 0 0 0);
paArr=(0 0 0 0 0);

# for i in `seq 0 4`; do
for i in `seq 4 -1 4`; do
  echo $i
  $smartRunFile ${lwArr[$i]} ${dtArr[$i]} ${nlArr[$i]} ${reArr[$i]} ${paArr[$i]};
  mkdir ./test/${lwArr[$i]}-$(date '+%Y%m%d_%H%M%S')/
  cp -R PostProcess ./diffWidthVolFracFilling/${lwArr[$i]}-$(date '+%Y%m%d_%H%M%S')/
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
