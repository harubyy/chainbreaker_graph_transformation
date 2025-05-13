#!/bin/bash

if [ $# -ne 1 ] ; then
  echo "usage: ./collect_errCnt.sh results/merge_on/rewrite/16_20"
  exit
fi

merge=`echo "$1" | cut -d/ -f2`
mode=`echo "$1" | cut -d/ -f3`
par=`echo "$1" | cut  -d/ -f4`
file_name=errCnt_"$merge"_"$mode"_"$par".csv
rm $file_name

echo $file_name

cd "$1"/
echo "BEFORE," >> $file_name
grep --color errCnt report.* | cut -d. -f2 | cut -d: -f1,3 --output-delimiter=, >> $file_name
echo "AFTER," >> $file_name

while read matrix
do
  name=`echo $matrix | cut -d',' -f2`

  errcnt=`grep --color errCnt "$name"/out | cut -d':' -f2`
  echo "$name","$errcnt" >> $file_name
done < ~/repos/chainbreaker/matrix.ID.name

cd -
