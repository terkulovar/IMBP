#!/bin/bash

#T='1 2 3 4 5 6 7 8 9 10'
#for j in $T; do
for j in {1..20}; do
 echo "t = " ${j}
 rm run1.in
 cat temp1.in >> run1.in
 echo "/random/setSeeds ${j} ${j}" >>run1.in
 cat temp2.in >> run1.in
# cp RunGun ./RunGun1
 ./RunGun run1.in
 sleep 10
 cd ../AnalysisTools
 root -b mytest8.C
# sleep 10
# root -b mytest5.C
 sleep 10
 cd ../build
# rm RunGun1
done

