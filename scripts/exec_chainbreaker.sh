#!/bin/bash

# $1: path to folders of the matrices

if [ $# -ne 1 ]; then
  echo "usage: ./exec_chainbreaker.sh par_8_20_merge_on_rewrite/new3cri"
  exit
fi

num_cores=`nproc --all`
num_cores=$(($num_cores/2))

rm "$1"/comp.report.csv
echo ", num. files, total size in MB" >> "$1"/comp.report.csv

while read matrix
do
  ID=`echo $matrix | cut -d',' -f1`
  name=`echo $matrix | cut -d',' -f2`

  cd "$1"/"$name"/
  echo "building chainbreaker for $name with ID: $ID"

  num_files=`ls *.c *.h | wc -l`
  echo -n "${name},${num_files},"  >> "${HOME}"/repos/chainbreaker/"$1"/comp.report.csv
  #soetimes ls fails, then manually used:
  #find -type f -name '*.c'  | wc -l
  #find -type f -name '*.h'  | wc -l
  #summed the two using echo $(())

  if [ "$num_files" -gt 60000 ]; then
    echo "too many files to build. aborting..."
  else 
#    sed -i "s|printf(\"x|//printf(\"/x|" main.c
    sed -i 's|/tmp/|'$HOME'/repos/chainbreaker/binary/new_3cri/|' main.c

    total_bytes=`du -ch *.c *.h | tail -1 | cut -f1`
    #sometimes du fails, then manually used:
    #du -x -d 1
    echo "$num_files with a total of $total_bytes"
    echo "$total_bytes" >> "${HOME}"/repos/chainbreaker/"$1"/comp.report.csv

    make -j $num_cores &> /dev/null
  fi
  rm *.o  &> /dev/null
  cd - &> /dev/null
done < matrix.ID.name

