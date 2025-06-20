import os
import sys
import re
import pickle
from collections import defaultdict
from collections import OrderedDict
import numpy as np
from length_and_loc_with_absolute import get_real_footprint_length_with_abs_start
import utils
from scanf import scanf

footprint_dict = OrderedDict()
footprint_abs_start_dict = defaultdict(list)
footprint_length_dict = defaultdict(list)
complete_footprint_length_dict = defaultdict(list)


label_dict = {
    "Naked-DNA" : "0",
    "TF" : "1",
    "Nuc":  "2",
    "discard" : "3"
}

left_extend = int(sys.argv[1])
right_extend = int(sys.argv[2])
for line in sys.stdin:

    #chr2L	480290	480320	chr2L	480021	480320	SRR3133329.41244179_41244179/1_adjacent`99~147	1	.	...................................................................................................................................FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF.............................................................................FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF.....	...................h.................................z..............................h..............h.........z....................Z....z...............................................H...X....H.Z........Z.H.....Z........X.........Z..................H.........Z..................x................H....
    #chr2L	480290	480320	chr2L	480028	480322	SRR3133326.43981750_43981750/1_overlapping`99~147	1	.	.......................................................................................................................................................................................................................................................................................................	............h.................................z..............................H..............H.........Z....................Z....z..............x.................h..............h.z.x....h.z........z.h.....z........x.........z..................h.........z..................x................h.....z
    #chr2L	480290	480320	chr2L	480028	480322	SRR3133329.39608881_39608881/1_overlapping`99~147	1	.	.................................................................................................................................FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF.......	............h.................................z..............................H..............H.........Z....................Z....Z..............x.................h..............h.z.x....h.z........z.h.....z........x.........z..................H.........z..................x................H.....z
    
    _, abs_roi_start, abs_roi_end, _, abs_read_start, abs_read_end, read_id, score, strand, complete_fp_str, complete_m_str = \
        scanf ("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%f\t%s\t%s\t%s\n", line)
    if strand == ".":
        strand = "+"
    rel_roi_start = abs_roi_start - abs_read_start
    rel_roi_end = abs_roi_end - abs_read_start + 1 
    intersected_fp_str = complete_fp_str[rel_roi_start: rel_roi_end]
    #print ([intersected_fp_str, rel_roi_start, rel_roi_end, complete_fp_str, line.strip()])
    intersected_m_str = complete_m_str[rel_roi_start: rel_roi_end]
    fp_length_list, _, _, abs_loc = get_real_footprint_length_with_abs_start(intersected_fp_str, rel_roi_start, rel_roi_end, complete_fp_str) 
    complete_fp_len = len(complete_fp_str)
    half = int ((rel_roi_end - rel_roi_start)/2)  
    #print ([complete_fp_str, complete_m_str]) 
    #print ([intersected_fp_str, complete_fp_str, abs_roi_start, abs_roi_end, abs_read_start, abs_read_end, left_extend, half, strand])
    
    left_neighbor = utils.extend_footprint_from_center (intersected_fp_str, complete_fp_str,
                                                  abs_roi_start, abs_roi_end,
                                                  abs_read_start, abs_read_end,
                                                  left_extend, half, strand) 
    right_neighbor = utils.extend_footprint_from_center (intersected_fp_str, complete_fp_str,
                                                  abs_roi_start, abs_roi_end,
                                                  abs_read_start, abs_read_end,
                                                  half, right_extend, strand) 
    
    
    per_capital_l, total_capital_l, total_small_l, total_l = utils.get_count_and_percentage_methylation (
                 left_neighbor)
    per_capital_r, total_capital_r, total_small_r, total_r = utils.get_count_and_percentage_methylation (
                 right_neighbor)
    # rightly assing percentange methylation level (inferred by the capital letters in the string)
    percentage_methylation = "NA"
    if (per_capital_l != "NA") and (per_capital_r != "NA"):                                                                
        if float(per_capital_l) <= float(per_capital_r):    
            percentage_methylation = per_capital_l # assigning minimum because, we have to check if the minimum is greater than a threshold
        else:
            percentage_methylation = per_capital_r
    elif (per_capital_l == "NA") and (per_capital_r != "NA"): 
        percentage_methylation = per_capital_r
    elif (per_capital_l != "NA") and (per_capital_r == "NA"):
        percentage_methylation = per_capital_l
    else: 
        percentage_methylation = "0" 
    complete_m_str_len, total_letters_on_complete_m_str, ratio_capital_letters, edges = \
        utils.capital_percentage_and_stretch_of_unbound(complete_m_str)  # see the output format in the definition
    max_of_edges = max(edges[0], edges[2])
    
    intersect_len = len(intersected_fp_str) 
    indiv_footprint_percent_overlap_with_intersect = utils.get_percent_for_each_footprint (intersected_fp_str)
    #print ([line.strip(), intersected_fp_str, intersect_len, intersected_fp_str])
    total_percent_footprint_in_intersect = round((intersected_fp_str.count('F')/intersect_len)*100, 3) # [(Total number of bases covered in footprint)/(intersect length)]*100.  
    
    # for each non-overlapping footprint in the intersect - label it with states, and finally choose the one with most contribution I
    
    labels_on_reads = defaultdict(list)
        
    
    lab = "NA"
    if len(indiv_footprint_percent_overlap_with_intersect) == 0:
        #print ([abs_start_loc_of_max, length_of_max_contribution, complete_fp_len])
        print ("\t".join([line.strip(), intersected_fp_str, "0.0", "0", label_dict["Naked-DNA"], str(complete_fp_len) ]))  
        lab = "Naked-DNA"
    else: # since we have realized there are fooptrints in the intersect with locus of interest - we should first check those lengths

        max_contribution_in_the_intersect = np.argmax (indiv_footprint_percent_overlap_with_intersect)
        max_percentage_overlap = indiv_footprint_percent_overlap_with_intersect[max_contribution_in_the_intersect] 
        length_of_max_contribution = fp_length_list[max_contribution_in_the_intersect]
        abs_start_loc_of_max = abs_loc [max_contribution_in_the_intersect]

        if total_percent_footprint_in_intersect <= 30 and float (percentage_methylation) > 25: 
            print ("\t".join([line.strip(), intersected_fp_str, str(max_percentage_overlap), str(length_of_max_contribution),label_dict["Naked-DNA"],  str(complete_fp_len)]))
            lab = "Naked-DNA"
        elif total_percent_footprint_in_intersect <= 30 and float(percentage_methylation) <= 25:
            print ("\t".join([line.strip(), intersected_fp_str, str(max_percentage_overlap), str(length_of_max_contribution),label_dict["Nuc"], str(complete_fp_len)]))
            lab = "Nuc"
        elif total_percent_footprint_in_intersect > 30 and length_of_max_contribution <= 50:
            if (abs_start_loc_of_max + length_of_max_contribution == complete_fp_len) or (abs_start_loc_of_max == 0):
                lab = "discard"
                print ("\t".join([line.strip(), intersected_fp_str, str(max_percentage_overlap), str(length_of_max_contribution), label_dict["discard"], str(complete_fp_len)]) )
            else:
                lab = "TF"
                print ("\t".join([line.strip(), intersected_fp_str, str(max_percentage_overlap), str(length_of_max_contribution), label_dict["TF"], str(complete_fp_len)]))
        elif length_of_max_contribution > 50:
            lab = "Nuc"
            print ("\t".join([line.strip(), intersected_fp_str, str(max_percentage_overlap), str(length_of_max_contribution),label_dict["Nuc"], str(complete_fp_len)]))
 
