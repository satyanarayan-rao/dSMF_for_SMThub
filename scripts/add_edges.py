import sys
from scanf import scanf
import utils
for line in sys.stdin:
    chr_a, abs_roi_start, abs_roi_end, chr_b, abs_read_start, abs_read_end, read_id, score, strand, complete_fp_str, complete_m_str = \
        scanf ("%s\t%d\t%d\t%s\t%d\t%d\t%s\t%f\t%s\t%s\t%s\n", line)
    left_border, right_border = utils.get_border_methylated_cytosines(complete_m_str)
    l_mvec = len (complete_m_str)
    new_vec = None
    if (left_border !=None) and (right_border!=None):
        new_vec = "F"*(left_border) + complete_fp_str[left_border: right_border+1] + 'F'*(l_mvec - right_border - 1)
    else:
        #No capital letters found - just assign the whole read as footprint
        new_vec =  'F'*(l_mvec)
    to_write = "\t".join( map(str, [chr_a, abs_roi_start, abs_roi_end, chr_b, abs_read_start, abs_read_end, read_id, score, strand, new_vec, complete_m_str]))
    print (to_write) 
