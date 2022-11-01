#!/bin/bash

rm avg_runtimes.csv

while read matrix
do
  name=`echo $matrix | cut -d',' -f2`
  avg_runtime=`grep $name runtimes | awk '{sum += $3} END {print sum /=10}'`

  echo $name,$avg_runtime >> avg_runtimes.csv
done < matrix.ID.name

