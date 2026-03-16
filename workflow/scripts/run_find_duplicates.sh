outdir=$1
psl_indir=$2
threshold=$3

mkdir -p ${outdir};

for f in ${psl_indir}*;
do
    file="$(basename -- $f)";
    cluster=${file%.psl};
    python workflow/scripts/find_duplicates.py \
    -p $f \
    -s ${threshold} \
    -o ${outdir}${cluster}_dupes.txt;
done