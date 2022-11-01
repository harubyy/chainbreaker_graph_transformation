#!/bin/bash

# $1: transformed/not: ORI/TR
# $2: strategy: original/costMap/threeCriteria/threeCriteria_rDist1

rm run

echo "#!/bin/bash" >> run
echo "rm runtimes" >> run

while read matrix
do
  name=`echo $matrix | cut -d',' -f2`

  echo "for i in {1..15}" >> run
  echo "do" >> run
  echo "  echo \"iteration: $i\""

  ./extract_sizes.sh "$name" "$1" "$2"

  echo "done" >> run
done < matrix.ID.name


chmod ug+x run
