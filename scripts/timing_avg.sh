#!/bin/bash

# $1: chrono.matrixNameA
# $2: transformed/not (TR/ORI)

rm chrono.build_DAG chrono.transform_DAG chrono.dump chrono.dump.matrix chrono.rewrite chrono.bldLvls chrono.calcFLOPS chrono.analyze chrono.rewrite_module chrono.policy chrono.without_dumps &> /dev/null

grep "build DAG" chrono."$1" >>  chrono.build_DAG
awk -F: '{sum += $2} END{avg=sum/7 ; print "avg. build DAG:," avg}' chrono.build_DAG >> chrono.avg."$1".csv

grep "buildLevels" chrono."$1" >> chrono.bldLvls
awk -F: '{sum += $2} END{avg=sum/7 ; print "avg. build levels -ARL; AIR:," avg}' chrono.bldLvls >> chrono.avg."$1".csv

grep "calculateFLOPS" chrono."$1" >> chrono.calcFLOPS
awk -F: '{sum += $2} END{avg=sum/7 ; print "avg. calculateFLOPS:," avg}' chrono.calcFLOPS >> chrono.avg."$1".csv

grep "matrix & code" chrono."$1" >> chrono.dump
awk -F, '{sum += $2} END{avg=sum/7 ; print "avg. dumping matrix & code:," avg}' chrono.dump >> chrono.avg."$1".csv

grep "dump data" chrono."$1" >> chrono.dump.matrix
awk -F:, '{sum += $2} END{avg=sum/7 ; print "avg. dump matrix:," avg}' chrono.dump.matrix >> chrono.avg."$1".csv

grep "rewrite module" chrono."$1" >> chrono.rewrite_module
awk -F: '{sum += $2} END{avg=sum/7 ; print "avg. rewrite module:," avg}' chrono.rewrite_module >> chrono.avg."$1".csv

grep "without dumps" chrono."$1" >> chrono.without_dumps
awk -F, '{sum += $2} END{avg=sum/7 ; print "avg. rewrite without dumps:," avg}' chrono.without_dumps >> chrono.avg."$1".csv

if [ "$2" == "TR" ]
then
  grep "transform DAG" chrono."$1" >> chrono.transform_DAG
  awk -F: '{sum += $2} END{avg=sum/7 ; print "avg. transform DAG:," avg}' chrono.transform_DAG >> chrono.avg."$1".csv

  grep "total rewrite" chrono."$1" >> chrono.rewrite
  awk -F, '{sum += $2} END{avg=sum/7 ; print "avg. rewrite:," avg}' chrono.rewrite >> chrono.avg."$1".csv

  grep "analyzeForCriteria" chrono."$1" >> chrono.analyze
  awk -F: '{sum += $2} END{avg=sum/7 ; print "avg. analyzeForCriteria:," avg}' chrono.analyze >> chrono.avg."$1".csv

  grep "policy" chrono."$1" >> chrono.policy
  awk -F: '{sum += $2} END{avg=sum/7 ; print "avg. policy(strategy):," avg}' chrono.policy >> chrono.avg."$1".csv
fi

