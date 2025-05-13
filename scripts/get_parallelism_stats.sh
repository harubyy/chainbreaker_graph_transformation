#!/bin/bash

if [ $# -ne 1 ]; then
  echo -e "enter the path\ne.g. results/merge_off/rewrite/16_20"
  exit
fi

rm thread_stats.csv &> /dev/null
echo ",max th.,single th. levels,multi th. levels,multi th. funs,single th. levels ptg.,single th. funs ptg" >> thread_stats.csv

while read matrix
do
  ID=`echo $matrix | cut -d',' -f1`
  name=`echo $matrix | cut -d',' -f2`

  echo $name

  files=`ls "$1"/${name}/run*.c`
  max_threads=`grep "parallel" $files | cut -d'(' -f2 | cut -d')' -f1 | sort | tail -1`

#  echo $max_threads

  # levels
  num_single_threaded=`grep single $files | wc -l`
  num_multi_threaded=`grep sections $files | wc -l`

  # functions

  tot=`wc -l --total=only "$1"/${name}/calculators*.h`
  cnt=`ls "$1"/${name}/calculators* | wc -l`
  fun_all=$((tot-cnt))
  num_multi_threaded_fun=$(($fun_all - $num_single_threaded))

#  echo $num_single_threaded
#  echo $num_multi_threaded
#  echo $num_multi_threaded_fun

  level_ptg=`echo "scale=2; $num_single_threaded / ($num_single_threaded + $num_multi_threaded)" | bc `

  fun_ptg=`echo "scale=2; $num_single_threaded / $fun_all" | bc `

  echo $name","$max_threads","$num_single_threaded","$num_multi_threaded","$num_multi_threaded_fun","$level_ptg","$fun_ptg >> thread_stats.csv
#done < matrix.ID.name
done < matrix.ID.name_
