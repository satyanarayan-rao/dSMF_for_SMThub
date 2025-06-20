rule making_compressed_string:
    input:
        footprints_on_fragments = "footprints_on_fragments/suppressed_merged_{cell_line}_to_{genome}_with_{setting}.bed.gz"
    output:
        compressed_cigar = "suppressed_merged/suppressed_merged_{cell_line}_to_{genome}_with_{setting}.cigar.bed.gz"
    shell:
        "zcat {input.footprints_on_fragments} | python scripts/compress_string.py | gzip -> {output.compressed_cigar}"

rule bedclip_cigar_strings:
    input:
        compressed_cigar = "suppressed_merged/suppressed_merged_{cell_line}_to_{genome}_with_{setting}.cigar.bed.gz",
        chrom_sizes = "metadata/{genome}.chrom.sizes" 
    output:
        bedclipped_cigar = "suppressed_merged/clipped_suppressed_merged_{cell_line}_to_{genome}_with_{setting}.cigar.bed"
    shell:
        "zcat {input.compressed_cigar} | scripts/bedClip stdin {input.chrom_sizes} {output.bedclipped_cigar}"
rule making_bigbed_file:
    input:
        bedclipped_cigar = "suppressed_merged/clipped_suppressed_merged_{cell_line}_to_{genome}_with_{setting}.cigar.bed",
        chrom_sizes = "metadata/{genome}.chrom.sizes"
    output:
        bigbed_file = "suppressed_merged/suppressed_merged_{cell_line}_to_{genome}_with_{setting}_and_mvec.bb"
    shell:
        "scripts/bedToBigBed -type=bed6+1 {input.bedclipped_cigar} {input.chrom_sizes} {output.bigbed_file}"
