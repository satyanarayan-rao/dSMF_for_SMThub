#!/bin/bash

zcat footprints_on_fragments/suppressed_merged_${1}_to_mm10_with_wobble_1_min_fp_10.bed.gz | python scripts/compress_string.py | gzip -> suppressed_merged/suppressed_merged_${1}_to_mm10_with_wobble_1_min_fp_10.cigar.bed.gz

../../SMThub/scripts/bedClip suppressed_merged/suppressed_merged_${1}_to_mm10_with_wobble_1_min_fp_10.cigar.bed.gz ../../SMThub/metadata/mm10.chrom.sizes suppressed_merged/suppressed_merged_clipped_${1}_to_mm10_with_wobble_1_min_fp_10.cigar.bed.gz

../../SMThub/scripts/bedToBigBed -type=bed6+1 suppressed_merged/suppressed_merged_clipped_${1}_to_mm10_with_wobble_1_min_fp_10.cigar.bed ../../SMThub/metadata/mm10.chrom.sizes suppressed_merged/suppressed_merged_${1}_to_mm10_with_wobble_1_min_fp_10.cigar_fp_and_mvec.bb
