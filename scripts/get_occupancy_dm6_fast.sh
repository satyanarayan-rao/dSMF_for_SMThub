chrom=$1
start=$2
end=$3
mid_point=`echo $2 $3 | awk '{print int (($1 + $2)/2)}'`
location="Distance from $4 ($1:${mid_point})"
annotation=$4
echo "$chrom $start $end" | awk '{print $1"\t"$2"\t"$3}' > ${annotation}.bed
#curl -JL "https://api.genome.ucsc.edu/getData/track?hubUrl=https://usegalaxy.org/api/datasets/f9cad7b01a472135e6ff84f7ae7d490c/display?to_ext=txt;genome=dm3;track=demo_hub;chrom=${chrom};start=${start};end=${end}" -o ${annotation}.json
curl -JL "https://api.genome.ucsc.edu/getData/track?hubUrl=https://usegalaxy.org/api/datasets/f9cad7b01a472135e6897af4f6ed17e1/display?to_ext=txt;genome=dm6;track=fp_and_mvec;chrom=${chrom};start=${start};end=${end}" -o ${annotation}.json
min_size=1000
actual_size=$(wc -c < ${annotation}.json)
if [ $actual_size -ge ${min_size} ]
then 
    cat ${annotation}.json | python scripts/parse_json_general.py fp_and_mvec > ${annotation}.all_fp.bed
    
    bedtools intersect -a ${annotation}.bed -b ${annotation}.all_fp.bed  -wa -wb -f 1  | \
    python scripts/expand_string.py | python scripts/add_edges.py  | \
    python scripts/assign_binding_states_first_check_percentage_methylation_around.py 30 30 \
    > ${annotation}.assigned_states.tsv
    cat ${annotation}.assigned_states.tsv | python scripts/extend_footprint_smthub.py 15 15 150 150 \
    ${annotation}.fp_extended.tsv ${annotation}.verbose.tsv
    sh scripts/order_footprints_single_binding_smthub.sh ${annotation}.verbose.tsv \
    ${annotation}.ordered_fp.tsv ${annotation}.ordered_m.tsv
    grep -v "#3" ${annotation}.ordered_fp.tsv | tac | Rscript \
    scripts/plot_matrix.R ${annotation}.pdf "$location"
    
    # cleanup
    #rm ${annotation}.ordered_fp.tsv ${annotation}.ordered_m.tsv 
    #rm ${annotation}.fp_extended.tsv ${annotation}.verbose.tsv 
    #rm ${annotation}.assigned_states.tsv ${annotation}.all_fp.bed
else
    echo -e "$chrom\t$start\t$end\tno hit"
fi
