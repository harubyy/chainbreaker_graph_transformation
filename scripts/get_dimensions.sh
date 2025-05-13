#!/bin/bash

if [ $# -ne 1 ]; then
  echo "usage: ./get_dimensions.sh results/merge_on/rewrite/3cri/16_20/"
  exit
fi

while read matrix
do
  ID=`echo $matrix | cut -d',' -f1`
  name=`echo $matrix | cut -d',' -f2`

  path="$1"/report."$name"
  line_begin=`grep -n dimen "$path" | cut -d: -f1`
  line_end=$(($line_begin+3))

  echo "$name": >> dim.csv
  sed -ne "${line_begin},${line_end}p" "$path" >> dim.csv
  echo "" >> dim.csv

done < matrix.ID.name
