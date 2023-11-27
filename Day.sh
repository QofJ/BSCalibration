#!/bin/bash
source ./setup.bsh
month=12
for i in {1..249}
do
   day=`expr 2022000 + $i`
   echo $day
   sed -e 's/day/'$day'/g' day.sh > day$day.sh
   chmod 777 day$day.sh
   hep_sub day$day.sh -g lhaaso
done