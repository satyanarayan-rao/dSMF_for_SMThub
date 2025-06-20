import sys
import json

intersect = json.load(sys.stdin)

for fp in intersect[sys.argv[1]]:
    chrom = fp["chrom"]
    st = fp["chromStart"]
    en = fp["chromEnd"]
    fp_info = fp["name"]
    score = fp["score"]
    strand = fp["strand"]
    m_vec = fp["field8"]
    print ("\t".join (map(str, [chrom, st, en, fp_info + "|" + m_vec, score, strand])))
