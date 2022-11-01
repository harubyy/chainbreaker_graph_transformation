#!/bin/bash

# $1: matrix name
# $2: transformed/not: TR/ORI

rm chrono

`grep "chrono" report."$1" >> chrono`

if [ "$2" == "ORI" ]
then
  awk -F':' 'NR<2 {print $0} NR >= 2 {sum += $2} END{print "chrono build matrix: " sum}' chrono >> chrono."$1"
else
  awk -F':' 'NR<3 {print $0} NR >= 3 {sum += $2} END{print "chrono build transformed matrix: " sum}' >> chrono."$1"
fi

