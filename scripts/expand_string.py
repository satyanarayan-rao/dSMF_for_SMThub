import sys
from uncompress_almost_cigar import uncompress_almost_cigar
for line in sys.stdin:
    #chr2L	480290	480320	chr2L	480021	480320	SRR3133329.41244179_41244179/1_adjacent`99~147|.131F52.77F35.5|.19h.33z.30h.14h.9z.20Z.4z.47H.3X.4H.Z.8Z.H.5Z.8X.9Z.18H.9Z.18x.16H.4	1	.
    #chr2L	480290	480320	chr2L	480028	480322	SRR3133326.43981750_43981750/1_overlapping`99~147|.295|.12h.33z.30H.14H.9Z.20Z.4z.14x.17h.14h.z.x.4h.z.8z.h.5z.8x.9z.18h.9z.18x.16h.5z	1	.
    #chr2L	480290	480320	chr2L	480028	480322	SRR3133329.39608881_39608881/1_overlapping`99~147|.129F159.7|.12h.33z.30H.14H.9Z.20Z.4Z.14x.17h.14h.z.x.4h.z.8z.h.5z.8x.9z.18H.9z.18x.16H.5z	1	.
    l_items = line.strip().split("\t")
    read_details, fp_comp, mvec_comp = l_items[6].split("|") 
    uncomp_fp = uncompress_almost_cigar(fp_comp)
    uncomp_mvec = uncompress_almost_cigar(mvec_comp)
    print ("\t".join(["\t".join(l_items[0:6]), 
                      read_details, l_items[7], l_items[8], uncomp_fp, uncomp_mvec]))

