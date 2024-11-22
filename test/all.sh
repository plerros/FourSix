#! /bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

find $SCRIPT_DIR -type f -name '*.json' | sort -V > findjson.tmp

rm -f results.md

echo "| name | input obtuse | output obtuse | steiner | seconds |" >> results.md
echo "| --- | --- | --- | --- | --- |" >> results.md

while read filename; do
	base_name=$(basename ${filename})
	resultname=./"$base_name".result
	timeout 100s ./FourSix -i "$filename" > "$resultname"
	input_obtuse=$( cat $resultname | grep -o 'Input Obtuse: [[:digit:]]*' | grep -o '[[:digit:]]*' )
	obtuse=$( cat $resultname | grep -o 'obtuse [[:digit:]]*' | grep -o '[[:digit:]]*' )
	steiner=$( cat $resultname | grep -o 'steiner [[:digit:]]*' | grep -o '[[:digit:]]*' )
	ms=$( cat $resultname | grep -o 'ms [[:digit:]]*' | grep -o '[[:digit:]]*' )
	ms=$( echo "scale=2; $ms/1000" | bc -l )
	echo "| $base_name | $input_obtuse | $obtuse | $steiner | $ms |" >> results.md
	echo "| $base_name | $input_obtuse | $obtuse | $steiner | $ms |"
	rm "$resultname"
done <findjson.tmp

rm findjson.tmp