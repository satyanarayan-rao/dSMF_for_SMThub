rule prepare_bisulfite_genome:
    input:
        genome_fa = lambda wildcards: genomes_df.loc[wildcards.genome, "fasta_path"]
    params:
    output:
        ga_conversion = "ref_genome/{genome}/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa", 
        ct_conversion = "ref_genome/{genome}/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa"
    shell:
        "bismark_genome_preparation --bowtie2 --verbose ref_genome/{wildcards.genome}/"
 
rule bismark_align_pe:
    input:
        read1 = "trimmed/{sample}_val_1.fq.gz", 
        read2 = "trimmed/{sample}_val_2.fq.gz",
        ga_conversion = "ref_genome/{genome}/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa", 
        ct_conversion = "ref_genome/{genome}/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa"
        
    params:
        genome = lambda wildcards:  "ref_genome/{genome}/".format (genome = wildcards.genome), 
        out_dir = "bismark_mapped",
        extra = "--gzip --no_dovetail -p 4",
        basename = lambda wildcards: "--basename {b}_to_{genome}".format(b = wildcards.sample, genome = wildcards.genome)
    output:
        temp("bismark_mapped/{sample}_to_{genome}_pe.bam"), 
        "bismark_mapped/{sample}_to_{genome}_PE_report.txt"
    script:
        "wrapper/bismark/align_pe.py"

rule bismark2report: 
    input: 
        "bismark_mapped/{sample}_to_{genome}_PE_report.txt"
    output: 
        "bismark_align_report/{sample}_to_{genome}.html"
    shell:
        "bismark2report --output {output} --alignment_report {input}" 

rule sort_bismark_bam:
    input:
        bismark_mapped_bam = "bismark_mapped/{sample}_to_{genome}_pe.bam", 
    params:
    output:
        sorted_bam = temp("bismark_mapped/{sample}_to_{genome}_pe_sorted.bam")
    shell:
        "samtools sort -n {input.bismark_mapped_bam} -o {output.sorted_bam}" 
rule merge_sorted_bam: 
    input:
        bam_file_list = lambda wildcards: ["bismark_mapped/{sample}_to_{genome}_pe_sorted.bam".format(sample = s, genome = wildcards.genome) for s in config["bam_merge_config"][wildcards.cell_line]]
    params:
    output:
        merged_bam = "bismark_mapped/merged_{cell_line}_to_{genome}.bam"
    shell:
        "sh scripts/merge_bam.sh \"{input.bam_file_list}\" {output.merged_bam} {wildcards.cell_line}"
