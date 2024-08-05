#!/bin/bash

#T='1 2 3 4 5 6 7 8 9 10'
#for j in $T; do
for j in {1..20}; do
 echo "t = " ${j}
 rm run2.in
 cat temp1.in >> run2.in
 echo "/random/setSeeds ${j} ${j}" >>run2.in
 cat temp3.in >> run2.in
# cp RunGun ./RunGun1
 ./RunGun run2.in
 sleep 10
 cd ../AnalysisTools
 root -b mytest9.C
# sleep 10
# root -b mytest5.C
 sleep 10
 cd ../build
# rm RunGun1
done

