#!/bin/bash

# $1: chrono.matrixNameA
# $2: transformed/not (TR/ORI)

`grep "build DAG" chrono."$1" >>  chrono.build_DAG`
awk -F: '{sum += $2} END{avg=sum/10 ; print "avg. build DAG: " avg}' chrono.build_DAG >> chrono.avg."$1"

if [ "$2" == "TR" ]
then
  `grep "transform DAG" chrono."$1" >> chrono.transform_DAG`
  awk -F: '{sum += $2} END{avg=sum/10 ; print "avg. transform DAG: " avg}' chrono.transform_DAG >> chrono.avg."$1"

  `grep "matrix" chrono."$1" >> chrono.build_transformed_matrix`
  awk -F: '{sum += $2} END{avg=sum/10 ; print "avg. build transformed matrix: " avg}' chrono.build_transformed_matrix >> chrono.avg."$1"
else
  `grep "matrix" chrono."$1" >> chrono.build_matrix`
  awk -F: '{sum += $2} END{avg=sum/10 ; print "avg. build matrix: " avg}' chrono.build_matrix >> chrono.avg."$1"
fi

[ -e file ] && rm chrono.build_DAG  chrono.build_matrix
[ -e file ] && rm chrono.transform_DAG chrono.build_transformed_matrix
