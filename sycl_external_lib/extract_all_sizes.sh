#!/bin/bash

# $1: transformed/not: ORI/TR
# $2: strategy: original/costMap/threeCriteria/threeCriteria_rDist1

rm run

echo "#!/bin/bash" >> run
echo "rm runtimes" >> run
echo "for i in {1..10}" >> run
echo "do" >> run
echo "  echo \"iteration: $i\""

while read matrix
do
  name=`echo $matrix | cut -d',' -f2`

  ./extract_sizes.sh "$name" "$1" "$2"
done < matrix.ID.name

echo "done" >> run

chmod ug+x run
