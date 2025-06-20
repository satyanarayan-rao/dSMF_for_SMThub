rule bam2fragment_level_methylation_bedgz:
    input:
        bam_file = "suppressed_context/suppressed_merged_{cell_line}_to_{genome}.bam"
    params:
    output:
        read_level_methylation_status = "fragment_level_methylation/suppressed_merged_{cell_line}_to_{genome}_fragment_methylation_vec.bed.gz"
    shell:
        "samtools view {input.bam_file} | python scripts/prepare_methylation_string.py | gzip - > {output.read_level_methylation_status}" 
rule overlapping_or_adjacent_reads:
    input:
        read_level_methylation_status = "fragment_level_methylation/suppressed_merged_{cell_line}_to_{genome}_fragment_methylation_vec.bed.gz" 
    output:
        overlapping_or_adjacent = "overlapping_or_adjacent_fragments/suppressed_merged_{cell_line}_to_{genome}_overlapping_or_adjacent.bed.gz"
    shell:
        "sh scripts/overlapping_or_adjacent_reads.sh {input.read_level_methylation_status}"
        " {output.overlapping_or_adjacent}"

