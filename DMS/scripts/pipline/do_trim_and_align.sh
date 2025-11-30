#!/bin/bash
set -e

INFO=../sampleInfo.csv

cat $INFO | tail -n +2 | while read x
do
    raw_file=$(echo $x | cut -d',' -f 4)
    trimed_file=$(echo $x | cut -d',' -f 3)

    echo start processing $raw_file
    cutadapt -l 26 -o $trimed_file $raw_file
    echo has trimmed to $trimed_file

done

IN=../../fastq
OUT=../outputs
scripts=.

mkdir -p $OUT/align

tmp=$OUT/tmp

mkdir -p $tmp/sh
mkdir -p $tmp/outputs

cat $INFO | tail -n +2 | while read x
do
    sample=$(echo $x | cut -d',' -f 1)
    lib=2mutBA5Tmerged
    antibody=$(echo $x | cut -d',' -f 2)
    file=$(echo $x | cut -d',' -f 3)
    
    if [ ! -f $file ];then
        echo "Empty $sample"
        continue
    fi

    wt=SARS-CoV-2-BF7.fasta
    bclen=26
    table=codon_variant_table_BF7.csv

    alignout=$OUT/align/${antibody}-${sample}_${lib}
    mkdir -p $alignout

    echo "#!/bin/bash" > $tmp/sh/_01-align-${sample}.tmp.sh
    echo "python $scripts/01-align.py -i $file -o $alignout -t $table -w $wt -b $bclen" >> $tmp/sh/_01-align-${sample}.tmp.sh
    echo $sample,wt:$wt,bclen:$bclen,table:$table

    bash $tmp/sh/_01-align-${sample}.tmp.sh > $tmp/outputs/_01-align-${sample}.o.txt
done