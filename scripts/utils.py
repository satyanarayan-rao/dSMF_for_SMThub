def get_border_methylated_cytosines (m_vec):
    """
    take the methylation vector and reports the first methylated loci from both ends
    returns None when no methylated cytosines are found
    """
    len_vec = len(m_vec)
    left_border = None
    right_border = None
    cnt = 0
    for c in m_vec:
        if (c==".") or (c == c.lower()):
            cnt+=1
        elif c == c.upper():
            break
    left_border = cnt # point at the capital letter in python index scheme
    cnt = 0
    for c in m_vec[::-1]:
        if (c==".") or (c == c.lower()):
            cnt+=1
        elif c == c.upper():
            break
    right_border = len_vec - cnt - 1 # point at the capital letter in python index scheme

    if left_border == len_vec:
        left_border = None
        right_border = None
    return left_border, right_border

def get_count_and_percentage_methylation (m_vec):
    total = 0
    total_lower = 0
    total_upper = 0
    per_c = 0
    for c in m_vec:
        if c == ".":
            continue
        elif c == c.lower():
            total +=1
            total_lower +=1
        elif c == c.upper():
            total +=1
            total_upper +=1
    if total>0:
        per_c = str(round(total_upper/total, 3)*100)
    else:
        per_c = "NA"
    return per_c, total_upper, total_lower, total

def get_percent_for_each_footprint (fvec):
    flen_list = []
    cnt = 0
    f_switch = False
    fvec_len = len(fvec)
    for c in fvec:
        if c == 'F':
            cnt +=1
            f_switch = True
        else:
            if cnt!=0:
                flen_list.append(round((cnt/fvec_len)*100, 2))
                cnt = 0
    if cnt!= 0:
        flen_list.append(round((cnt/fvec_len)*100, 2))
    return flen_list



def extend_footprint_from_center (intersected_str, complete_str,
                                  abs_roi_start, abs_roi_end,
                                  abs_read_start, abs_read_end,
                                  lextend, rextend, strand): # pc : peak center 
    existing_footprint = intersected_str # existing footprint means that observed : extension can virtually be of any length
    complete_footprint = complete_str # just renaming 
    complete_footprint_start = abs_read_start
    complete_footprint_end = abs_read_end
    rel_roi_start = abs_roi_start - abs_read_start
    rel_roi_end = abs_roi_end - abs_read_start + 1
    peak_center = int ( (abs_roi_start + abs_roi_end )/2  )
    #print (peak_center)
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
    #lflank = int(sys.argv[3])
    #rflank = int(sys.argv[4])
    #lextend = int(sys.argv[5])
    #rextend = int(sys.argv[6])
    #    if strand == "-": 
    #        lexetend = int(sys.argv[6])
    #        rextend = int(sys.argv[5])
    #        lflank = int(sys.argv[4])
    #        rflank = int(sys.argv[3]) 
        
    desired_length = rextend - (0 - lextend) + 1 
    extend_start = peak_center - (lextend) # close interval
    extend_stop = peak_center + (rextend + 1 )  # open 
    extended_footprint = ""
    extended_mvec = ""
    extended_bsseq = ""
    #print ([extend_start, extend_stop, complete_footprint_start, complete_footprint_end])
    #sys.exit(-1)
    if (extend_start >= complete_footprint_start) and  (extend_stop <= complete_footprint_end + 1): 
        #print ([key_with_hash, "condition 1"])
        if strand == "+":
            extended_footprint = complete_footprint[extend_start - complete_footprint_start: extend_stop - complete_footprint_start]
        elif strand == "-": 
            extended_footprint = complete_footprint[extend_start - complete_footprint_start: extend_stop - complete_footprint_start][::-1]
            # recall that complete_footprint is on the watson strand always
            #print(["in", "in", key_without_hash,  extended_footprint, len(extended_footprint)])
    elif (extend_start >= complete_footprint_start) and  (extend_stop > complete_footprint_end + 1): 
        #print ([key_with_hash, "condition 2"])
        if strand == "+":
            extended_footprint = complete_footprint[extend_start - complete_footprint_start: complete_footprint_end - complete_footprint_start + 1]
            # add `M` of how much short in the end 
            extended_footprint = extended_footprint + "M"*(desired_length - len(extended_footprint))
        elif strand == "-":
            extended_footprint = complete_footprint[extend_start - complete_footprint_start: complete_footprint_end - complete_footprint_start + 1][::-1] 
            #print([peak_center, extend_start, complete_footprint_start, complete_footprint_end])
            extended_footprint = "M"*(desired_length - len(extended_footprint)) + extended_footprint
            #print(["in", "out", key_without_hash,  extended_footprint, len(extended_footprint)])
    elif (extend_start < complete_footprint_start) and (extend_stop <= complete_footprint_end + 1): 
        if strand == "+":
            extended_footprint = complete_footprint[0: extend_stop - complete_footprint_start]
            extended_footprint = "M"*(desired_length - len(extended_footprint)) + extended_footprint 
        #print(["out", "in", key_without_hash,  extended_footprint, len(extended_footprint)]) 
        elif strand == "-":
            extended_footprint = complete_footprint[0: extend_stop - complete_footprint_start][::-1]
            extended_footprint = extended_footprint + "M"*(desired_length - len(extended_footprint)) 
            #print(["out", "in", key_without_hash,  extended_footprint, len(extended_footprint)]) 
        #print ([key_with_hash, "condition 3", extended_mvec])
        
    elif (extend_start < complete_footprint_start) and (extend_stop > complete_footprint_end + 1): 
        #print ([key_with_hash, "condition 4"])
        if strand == "+":
            to_add_in_left = "M"*(complete_footprint_start - extend_start)
            to_add_in_right = "M"*(extend_stop - (complete_footprint_end + 1) )
            extended_footprint = to_add_in_left + complete_footprint + to_add_in_right   
        elif strand == "-":
            to_add_in_left = "M"*(complete_footprint_start - extend_start)
            to_add_in_right = "M"*(extend_stop - (complete_footprint_end + 1) )
            extended_footprint = to_add_in_right + complete_footprint[::-1] + to_add_in_left
            #print(["out", "out", key_without_hash, extended_footprint, len(extended_footprint)])
    
    #print ([extend_start, extend_stop, complete_footprint_start, complete_footprint_end, extended_footprint])
    return extended_footprint

def capital_percentage_and_stretch_of_unbound(mvec):
    """
    Input: Methylation string/vector : "...................h.................................z..............................h..............h.........z....................Z....z...............................................H...X....H.Z........Z.H.....Z........X.........Z..................H.........Z..................x................H...."
    Output: [String length, total number of nondot letter, ratio of capital letters, [occluded edge length left (until the first capital letter), total string length, occluded edge length right (until the first capital letter from right)]]  
          : [          300,                            20,                     0.65, [                                                       130,                 300,                                                                      4]]  
                               
    """
    len_vec = len(mvec)
    check_vec = "."*len_vec
    total_capital = 0
    total_small = 0
    left_extreme = 0
    right_extreme = len_vec
    first_found = False
    base_counter = -1
    last_capital = -1
    for c in mvec:
        base_counter +=1
        if c == ".":
            continue
        elif c == c.lower():
            total_small +=1
        elif c == c.upper():
             total_capital +=1
             if first_found == False:
                 left_extreme = base_counter
                 first_found = True
             else:
                 last_capital = base_counter
        else:
            continue
    naked_dna_stretches = []
    if first_found == True: # there was at least one methylation event
         if last_capital == -1: # there was only one capital
             naked_dna_stretches.append(left_extreme)
             naked_dna_stretches.append(0)
             naked_dna_stretches.append (right_extreme - left_extreme)
         else:
             naked_dna_stretches.append(left_extreme)
             naked_dna_stretches.append(last_capital - left_extreme - 1)
             naked_dna_stretches.append(right_extreme - last_capital - 1)
    else:
        naked_dna_stretches = [0, 0, 0]
    total = total_capital + total_small
    percentage_capital = "NA"
    if total > 0:
        percentage_capital = round(total_capital/total, 2)

    return len_vec, total, percentage_capital, naked_dna_stretches

