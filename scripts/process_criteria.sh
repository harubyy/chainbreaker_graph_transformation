#!/bin/bash

if [ $# -ne 1 ]; then
  echo "enter the path"
  echo "e.g.  results/merge_off/rewrite/16_20/"
  exit
fi

grep "num. of levels:" -m 3 $1/report.* | sed '1~3d' | cut -d. -f2-4 >> num_levels.csv

grep ARL $1/report.*  >> ARL.csv
grep AIR $1/report.* | awk -F"." '{print $2}'  >> AIR.csv
grep ALC $1/report.* | awk -F"." '{print $2}'  >> ALC.csv
