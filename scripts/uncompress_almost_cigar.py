def uncompress_almost_cigar(almost_cigar):
    from collections import defaultdict
    num_dict = defaultdict(lambda : False)
    for i in range(10):
        num_dict[str(i)] = True
    uncomp = []
    num_found = False
    rep_num = []
    for c in almost_cigar:
        if num_dict[c] == False:
            if num_found == True:
                rep = int("".join(rep_num))
                uncomp.append(uncomp[-1]*(rep - 1)) # rep - 1 because the character has been pushed once already 
                num_found = False
                rep_num = []
            uncomp.append(c)
        else:
            num_found = True
            rep_num.append(c)
    # handle the last instance - in case rep_num is not empty - then it has to be filled
    if len(rep_num) > 0:
        rep = int("".join(rep_num))
        uncomp.append(uncomp[-1]*(rep - 1))
        num_found = False
        return "".join(uncomp)
    else:
        return "".join(uncomp)    

