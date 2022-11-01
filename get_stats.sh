#!/bin/bash

if [ "$1" == "--help" ] || [ "$#" -ne 1 ]
then
  echo -e "parameter #1:"
  echo -e "report name"
  exit
fi

matrix=`echo $1 | cut -d. -f2`

start=`grep -n BEFORE $1 | head -1 | cut -d':' -f1`
end=`grep -n BEFORE $1 | tail -1 | cut -d':' -f1`
sed -ne "${start},${end}p" $1 >> "$matrix".before.csv 

start=`grep -n AFTER $1 | head -1 | cut -d':' -f1`
end=`grep -n AFTER $1 | tail -1 | cut -d':' -f1`
sed -ne "${start},${end}p" $1 >> "$matrix".after.csv 
