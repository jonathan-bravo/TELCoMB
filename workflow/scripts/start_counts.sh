input=$1
output=$2

echo $(zcat -c ${input} | wc -l)/4 | bc > ${output}
