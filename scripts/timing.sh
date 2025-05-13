#!/bin/bash

# $1: matrix name
# $2: transformed/not: TR/ORI

rm chrono chrono.dump."$1"  &> /dev/null

`grep "chrono" report."$1" >> chrono`
`grep "dump" chrono >> chrono.dump."$1"`

if [[ "$2" == "TR" ]]; then
  awk -F':' '{sum += $2} END{print "chrono dump transformed matrix & code:," sum}' chrono.dump."$1" >> chrono."$1"
else
  awk -F':' '{sum += $2} END{print "chrono dump matrix & code:," sum}' chrono.dump."$1" >> chrono."$1"
fi

rewrite_exec=`grep "calculate rewriteExecutor" chrono | cut -d' ' -f4`
dump_runXcalcX=`grep -e "dump runX.c" -e "dump calculateX.c" chrono | awk -F: '{sum += $2} END {print sum}'`
rewrite_exec_nodump=`echo "$rewrite_exec - $dump_runXcalcX" | bc`
echo "rewriteExecutor without dump: $rewrite_exec_nodump" >> chrono."$1"

dump_data_main=`grep -e "dump main" -e "dump data" chrono | awk -F: '{sum += $2} END {print sum}'`
calculateB=`grep "calculate B" chrono | cut -d: -f2`
rewrite_total=`echo "$calculateB + $rewrite_exec_nodump" | bc`
if [[ "$2" == "TR" ]]; then
  echo "total rewrite without dumps(build transformed matrix):, $rewrite_total" >> chrono."$1"
else
  echo "total rewrite without dumps(build matrix):, $rewrite_total" >> chrono."$1"
fi

echo "dump data:, `grep "dump data" chrono | cut -d: -f2 `" >> chrono."$1"


grep "build DAG" chrono >> chrono."$1"
grep "rewrite module" chrono >> chrono."$1"
grep "buildLevels" chrono >> chrono."$1"
grep "calculateFLOPS" chrono >> chrono."$1"

if [[ "$2" == "TR" ]]; then
  grep "transform DAG" chrono >> chrono."$1"
  grep "analyzeForCriteria" chrono >> chrono."$1"
  grep "policy" chrono >> chrono."$1"
fi

