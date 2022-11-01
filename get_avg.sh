#!/bin/bash

# $1: ORI/TR

rm avg.csv

if [ "$1" == "ORI" ]
then
  echo "matrix,build DAG,build matrix" >> avg.csv
else
echo "matrix,build DAG,transform DAG,build matrix" >> avg.csv
fi

while read matrix
do
  name=`echo $matrix | cut -d',' -f2`
  build_dag=`cat chrono.avg."$name" | head -1 | cut -d: -f2` >> avg.csv
  if [ "$1" == "TR" ]
  then
    transform_dag=`cat chrono.avg."$name" | head -2 | tail -1 | cut -d: -f2` >> avg.csv
  fi

  build_matrix=`cat chrono.avg."$name" | tail -1 | cut -d: -f2` >> avg.csv
  if [ "$1" == "ORI" ]
  then
    echo "${name},${build_dag},${build_matrix}" >> avg.csv
  else
    echo "${name},${build_dag},${transform_dag},${build_matrix}" >> avg.csv
  fi
#done < ../matrix.ID.name
done < matrix.ID.name

