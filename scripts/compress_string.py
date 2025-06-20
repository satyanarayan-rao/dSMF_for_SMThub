import sys


def compress_string (s):
    final_string = s[0]
    current_counter = 1
    for c in s[1:]:
        if final_string[-1] == c:
            current_counter +=1
        else:
            if current_counter == 1:
                final_string += c 
                current_counter = 1
            else:
                final_string += str(current_counter)
                final_string += c
                current_counter = 1
    if current_counter > 1: 
        final_string += str(current_counter)
    return final_string

for line in sys.stdin:
    #chr2L	480021	480320	SRR3133329.41244179_41244179/1_adjacent`99~147	.	...................................................................................................................................FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF.............................................................................FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF.....	...................h.................................z..............................h..............h.........z....................Z....z...............................................H...X....H.Z........Z.H.....Z........X.........Z..................H.........Z..................x................H....	AGATTGTAATGTTGGTTGGTAATAAGTTAAAAGAGATTTAGGTGAATTTAAATTGTTTTTAGTAAATTTTTTAAAAGTGTATGGTTATTTTAAATGTTGTTATTTTTTTTGAGAAGATTTAAAAAATTACCGAAGTGATTTTTATTTGTGCAGGATATAGATATATATCTAGTAGGAGAGAGGCTAGCAGAGCACGTAATTAACGCATGGACGTTATTAGCTGTTTATTGCGGGGTTTTGTGTTTGGTGCAATAGAGAACGGGTTAGTAGAGGAGGAGTTGGGGTGGGTGTGGGGCTTTA
    l_items = line.strip().split("\t")
    footprint = l_items[5]
    methyl_vec = l_items[6]
    compressed_footprint = compress_string (footprint)
    compressed_vec = compress_string (methyl_vec)
    #print (compressed_footprint) 
    #print (compressed_vec)
    to_print = "\t".join ([l_items[0], l_items[1], l_items[2], 
                           l_items[3] + "|" + compressed_footprint, 
                           "1", ".", compressed_vec])
    print (to_print)
   
