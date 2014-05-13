#!/bin/bash
#./bin/atoms -v -O 1 -i 10 -n 1k |&  grep PERF | head -2 | tail -1 | cut -c 10-200 | tr -d ' '*

#./bin/atoms -v -O 1 -i 10 -n 1k |&  grep PERF | head -2 | tail -1 | cut -c 20-30 | tr -d ' '
#| cut -f 2

#./bin/atoms -v -O 1 -i 10 -n 1k |&  sed "s/^$ .+\s*([0-9]+)\s*([0-9]+)\s*([0-9].?[0-9]+)"

#./bin/atoms -v -O 1 -i 10 -n 1k |&  grep "/s ^\[PERF]\s*([0-9]\+)\s*([0-9]\+)\s*([0-9]\.\?[0-9]\+)"

# ok mais ne marche pas 
#^\[PERF]\s*([0-9]+)\s*([0-9]+)\s*([0-9]\.?[0-9]+)$


nb=10
x=1
var2=0 
while [ $x -le $nb ]
do
var=`./bin/atoms -v -O 1 -i 10 -n 1k |&  grep PERF | head -2 | tail -1 | cut -c 20-30 | tr -d ' '`
var2=$(( var2 + $var ))
#echo $var
#echo $var2
x=$(( $x + 1 ))
done
moy=$(( $var2 / $nb ))
echo -e "$moy"
