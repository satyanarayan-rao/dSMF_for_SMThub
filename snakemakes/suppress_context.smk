rule suppress_context:
    input:
        bam_file = "bismark_mapped/merged_{cell_line}_to_{genome}.bam"
    output:
        out_bam = temp("suppressed_context/suppressed_merged_{cell_line}_to_{genome}.bam")
    shell:
        "sh scripts/suppress_context.sh {input.bam_file} {output.out_bam}"
