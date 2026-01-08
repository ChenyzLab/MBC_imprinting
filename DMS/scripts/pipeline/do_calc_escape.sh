#!/bin/bash
set -e
IN=../outputs
INFO=../sampleInfo.csv
scripts=.
tmp=$IN/tmp
mkdir -p $IN/calc
cat $INFO | tail -n +2 | awk 'BEGIN{FS=","}{if($7 != "ref") print $0}' | while read x
do
    file=$(echo $x | cut -d',' -f 3)
    sample=$(echo $x | cut -d',' -f 1)
    lib=2mutBA5Tmerged
    antibody=$(echo $x | cut -d',' -f 2)
    ref=$(echo $x | cut -d',' -f 8)
    fraq=$(echo $x | cut -d',' -f 9)

    echo "$sample . Library: $lib . Antibody: $antibody"
    wt=SARS-CoV-2-BF7.fasta
    table=codon_variant_table_BF7.csv
    singlemut=expr_2mutBA5Tmerged.csv
    echo "${antibody}-${sample}_${lib}. Using reference $ref"

    cat << EOF > $tmp/sh/_02_${antibody}-${sample}_${lib}.tmp.sh
#!/bin/bash
python $scripts/02-call_escape_score.py \
           -r $ref \
           -e $IN/align/${antibody}-${sample}_${lib}/${sample}_variant_counts.csv \
           -exp FACS \
           -w $wt -t $table \
           -o $IN/calc \
           -p SARS-CoV-2-BF7 -S 330 \
           --mutbindexpr=$singlemut \
           --exprmin=-1 --bindmin=-999 --frac_escape=$fraq --cells_sorted=10000
EOF

    bash $tmp/sh/_02_${antibody}-${sample}_${lib}.tmp.sh > $tmp/outputs/_02_${antibody}-${sample}_${lib}.o.txt
done