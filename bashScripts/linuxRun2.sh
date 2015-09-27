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

# lwArr=(0.01 1.1 1.3 1.4 0.8);
# dtArr=(0.0005 0.0005 0.0005 0.0005 0.0005);
# nlArr=(80 30 30 25 40);
# lwArr=(0.1 0.3 0.5 0.7 0.9);
# dtArr=(0.0005 0.0005 0.0005 0.0005 0.0005);
# nlArr=(75 60 45 40 35);
lwArr=(0.1 0.2 0.3 0.5 0.7 0.8 0.9 1.0);
dtArr=(0.0005 0.0005 0.0005 0.0005 0.0005 0.0005 0.0005 0.0005);
nlArr=(100 90 80 70 60 50 40 40);
reArr=(0 0 0 0 0);
paArr=(0 0 0 0 0);

# for i in `seq 0 4`; do
for i in `seq 0 7`; do
  echo $i
  a=./Cyl4.4VolFracFilling/${lwArr[$i]}-$(date '+%Y%m%d-%H%M%S')
  mkdir $a
  cp smarticleMoves.csv $a
  cd $a
  $smartRunFile ${lwArr[$i]} ${dtArr[$i]} ${nlArr[$i]} ${reArr[$i]} ${paArr[$i]};
  cd ../..;
done;

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
