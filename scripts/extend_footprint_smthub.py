import os
import sys
import pickle
import re
from scanf import scanf
# head -1 footprint.bed
#chr2L	11806729	11806730	chr2L:11806729-11806729^20	.	+	chr2L	11806460	11806747	SRR3133326.952549_952549/1_overlapping`99~147`99~147	.	........F...................FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF......................FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF.....FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF.....FFFFFFFFFFFFFFFFFFFFFFF.....FFF.....FFFFFFF......................................

# head -1 merged_S2_to_reordered_20clust_147_50PNE_open_v1_peak_cluster_7_lf_22_rf_19_methylation_matrix_for_site_peak_3484_nclust_3_sam_flag_83~163.tsv 
# chr3R:5416975-5416975^7`SRR3133326.727267_727267/1_overlapping`83~163`83~163#4	...................FFFFFFFFFFFFFF.....FFFF

#head -1 ../tmp/all_rows_with_same_columns.tsv 

#chr2L	480290	480320	chr2L	480021	480320	SRR3133329.41244179_41244179/1_adjacent`99~147	1.0	.	FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF.FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF.............................................................................FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF.FFFF	...................h.................................z..............................h..............h.........z....................Z....z...............................................H...X....H.Z........Z.H.....Z........X.........Z..................H.........Z..................x................H....	FFFFFFFFFFFFFFFFFFFFFFFFFF.FFFF	83.87	35	1	300

lflank = int(sys.argv[1])
rflank = int(sys.argv[2])
lextend = int(sys.argv[3])
rextend = int(sys.argv[4])

out_fp = open(sys.argv[5], "w")
out_verb_fp = open(sys.argv[6], "w")

for line in sys.stdin: 

    chrom, abs_roi_start, abs_roi_end, _, abs_read_start, abs_read_end, read_id, score, strand, complete_fp_str, complete_m_str, intersected_fp_str, max_percentage_overlap, length_of_max_contribution, label_on_read, complete_vec_len = scanf ("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%f\t%s\t%s\t%s\t%s\t%f\t%d\t%d\t%d\n", line)

    key_with_hash = chrom + ":" + str(abs_roi_start) + "-" + str(abs_roi_end) + "`" + read_id + "#" + str(label_on_read)
    existing_footprint = intersected_fp_str
    complete_footprint = complete_fp_str
    complete_mvec = complete_m_str
    complete_footprint_start = abs_read_start
    complete_footprint_end = abs_read_end
    peak_center = int((abs_roi_start + abs_roi_end)/2)
    m_vec_start = peak_center - complete_footprint_start -  lflank
    m_vec_stop = peak_center - complete_footprint_start  + rflank + 1
    ####
    # lf = -3, rf = 3 
    # m_vec_start = 26
    #  m_vec_stop = 33 
    # extend_start = m_vec_start - (lextend - lf) 
    # extend_stop = m_vec_stop + (rextend - rf)  
    # 
    # complete_footprint_start and complete_footprint_end are in closed form length of footprint = 288, `11806747 - 11806460 + 1` 
    # actual  2324252627282930313233343536 
    #          . . . . . F F F F . . . F F 
    #               -3-2-1 0 1 2 3
    # desired -5 7 
    #          . . . . . F F F F . . . F F
    #           -5-4-3-2-1 0 1 2 3 4 5 6 7 
    # 
    #### 
       
    desired_length = rextend - (0 - lextend) + 1 
    extend_start = peak_center - (lextend) # close interval
    extend_stop = peak_center + (rextend + 1 )  # open 
    extended_footprint = ""
    extended_mvec = ""
    #extended_bsseq = ""
    #print ([extend_start, extend_stop, complete_footprint_start, complete_footprint_end])
    #sys.exit(-1)
    if strand == ".":
        strand = "+"
    if (extend_start >= complete_footprint_start) and  (extend_stop <= complete_footprint_end + 1): 
        #print ([key_with_hash, "condition 1"])
        if strand == "+":
            extended_footprint = complete_footprint[extend_start - complete_footprint_start: extend_stop - complete_footprint_start]
            extended_mvec = complete_mvec[extend_start - complete_footprint_start: extend_stop - complete_footprint_start]
            #extended_bsseq = complete_bsseq[extend_start - complete_footprint_start: extend_stop - complete_footprint_start]
        elif strand == "-": 
            extended_footprint = complete_footprint[extend_start - complete_footprint_start: extend_stop - complete_footprint_start][::-1]
            extended_mvec = complete_mvec[extend_start - complete_footprint_start: extend_stop - complete_footprint_start][::-1]
            #extended_bsseq = complete_bsseq[extend_start - complete_footprint_start: extend_stop - complete_footprint_start][::-1]
            # recall that complete_footprint is on the watson strand always
            #print(["in", "in", key_without_hash,  extended_footprint, len(extended_footprint)])
    elif (extend_start >= complete_footprint_start) and  (extend_stop > complete_footprint_end + 1): 
        #print ([key_with_hash, "condition 2"])
        if strand == "+":
            extended_footprint = complete_footprint[extend_start - complete_footprint_start: complete_footprint_end - complete_footprint_start + 1]
            extended_mvec = complete_mvec[extend_start - complete_footprint_start: complete_footprint_end - complete_footprint_start + 1]
            #extended_bsseq = complete_bsseq[extend_start - complete_footprint_start: complete_footprint_end - complete_footprint_start + 1]
            # add `M` of how much short in the end 
            extended_footprint = extended_footprint + "M"*(desired_length - len(extended_footprint))
            extended_mvec = extended_mvec + "M"*(desired_length - len(extended_mvec))
            #extended_bsseq = extended_bsseq + "M"*(desired_length - len(extended_bsseq))
        elif strand == "-":
            extended_footprint = complete_footprint[extend_start - complete_footprint_start: complete_footprint_end - complete_footprint_start + 1][::-1] 
            extended_mvec = complete_mvec[extend_start - complete_footprint_start: complete_footprint_end - complete_footprint_start + 1][::-1] 
            #extended_bsseq = complete_bsseq[extend_start - complete_footprint_start: complete_footprint_end - complete_footprint_start + 1][::-1] 
            #print([peak_center, extend_start, complete_footprint_start, complete_footprint_end])
            extended_footprint = "M"*(desired_length - len(extended_footprint)) + extended_footprint
            extended_mvec = "M"*(desired_length - len(extended_mvec)) + extended_mvec 
            #extended_bsseq = "M"*(desired_length - len(extended_bsseq)) + extended_bsseq 
            #print(["in", "out", key_without_hash,  extended_footprint, len(extended_footprint)])
    elif (extend_start < complete_footprint_start) and (extend_stop <= complete_footprint_end + 1): 
        if strand == "+":
            extended_footprint = complete_footprint[0: extend_stop - complete_footprint_start]
            extended_mvec = complete_mvec[0: extend_stop - complete_footprint_start]
            #extended_bsseq = complete_bsseq[0: extend_stop - complete_footprint_start]
            extended_footprint = "M"*(desired_length - len(extended_footprint)) + extended_footprint 
            extended_mvec = "M"*(desired_length - len(extended_mvec)) + extended_mvec 
            #extended_bsseq = "M"*(desired_length - len(extended_bsseq)) + extended_bsseq
        #print(["out", "in", key_without_hash,  extended_footprint, len(extended_footprint)]) 
        elif strand == "-":
            extended_footprint = complete_footprint[0: extend_stop - complete_footprint_start][::-1]
            extended_mvec = complete_mvec[0: extend_stop - complete_footprint_start][::-1]
            #extended_bsseq = complete_bsseq[0: extend_stop - complete_footprint_start][::-1]
            extended_footprint = extended_footprint + "M"*(desired_length - len(extended_footprint)) 
            extended_mvec = extended_mvec + "M"*(desired_length - len(extended_mvec)) 
            #extended_bsseq = extended_bsseq + "M"*(desired_length - len(extended_bsseq)) 
            #print(["out", "in", key_without_hash,  extended_footprint, len(extended_footprint)]) 
        #print ([key_with_hash, "condition 3", extended_mvec])
        
    elif (extend_start < complete_footprint_start) and (extend_stop > complete_footprint_end + 1): 
        #print ([key_with_hash, "condition 4"])
        if strand == "+":
            to_add_in_left = "M"*(complete_footprint_start - extend_start)
            to_add_in_right = "M"*(extend_stop - (complete_footprint_end + 1) )
            extended_footprint = to_add_in_left + complete_footprint + to_add_in_right   
            extended_mvec = to_add_in_left + complete_mvec + to_add_in_right   
            #extended_bsseq = to_add_in_left + complete_bsseq + to_add_in_right   
        elif strand == "-":
            to_add_in_left = "M"*(complete_footprint_start - extend_start)
            to_add_in_right = "M"*(extend_stop - (complete_footprint_end + 1) )
            extended_footprint = to_add_in_right + complete_footprint[::-1] + to_add_in_left
            extended_mvec = to_add_in_right + complete_mvec[::-1] + to_add_in_left
            #extended_bsseq = to_add_in_right + complete_bsseq[::-1] + to_add_in_left
            #print(["out", "out", key_without_hash, extended_footprint, len(extended_footprint)])
    
    #print ([extend_start, extend_stop, complete_footprint_start, complete_footprint_end, extended_footprint])

    out_fp.write(key_with_hash + "\t" + extended_footprint + "\n")
    out_verb_fp.write(key_with_hash + "\t" + extended_footprint + "\n")
    out_verb_fp.write(key_with_hash + "\t" + extended_mvec + "\n")
    #out_verb_fp.write(key_with_hash + "\t" + extended_bsseq + "\n")

out_fp.close()
out_verb_fp.close()
    
#    prefix = [] 
#    for left in range(lextend - lflank): 
#        if m_vec_start - left - 1 >= complete_footprint_start:
#            prefix.append(m_vec_start - complete_footprint_start - left - 1)
#        else:
#            prefix.append("M") 
#    prefix_str = "".join(prefix)[::-1] 
#    for right in range() 
