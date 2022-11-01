#!/bin/bash

# $1: matrix name
# $2: transformed/not: TR/ORI
# $3: strategy: original/costMap/threeCriteria/threeCriteria_rDist1

rowPtr=`grep rowPtr ../reports_"$3"/report."$1" | cut -d' ' -f2`
colIdx=`grep colIdx ../reports_"$3"/report."$1" | cut -d' ' -f2`
values=`grep values ../reports_"$3"/report."$1" | cut -d' ' -f2`

if [ "$2" == ORI ]
then 
  echo "original matrix:"
  echo "rowPtr colIdx values:"
  echo "$rowPtr $colIdx $values $1"
  echo "  echo \"running for $1\"" >> run 
  echo "  ./sptrsv $(( $rowPtr - 2 )) $colIdx $1 false >> runtimes" >> run
else
  echo "transformed matrix:"
  echo "rowPtr colIdx values:"
  echo "$rowPtr $colIdx $values $1"
  echo "echo \"running for $1\"" >> run 
  echo "  ./sptrsv $(( $rowPtr - 2 )) $colIdx $1 true >> runtimes " >> run
fi

