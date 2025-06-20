import sys
import re
#max_windows = 'data_from_ReMap/sample_max_peaks_windows.fa'
#max_windows = 'data_from_ReMap/max_peaks_windows.fa'
#windows = open(max_windows)


header = ""

#for line in windows:
for line in sys.stdin:

    #>Max:Schneider-2@window1:12845::chr2L:12757-12787
    #GTATTTCTTTAGCAAGCTGCGCAGAAATTC
    #>Max:Schneider-2@window2:12845::chr2L:12787-12817
    #GGCGGGGCACGTGTGGTGGTGCATTGCCAC
    if line.startswith(">"):
       header = line.strip()
 
    else:
        sequence = line.strip() 
        if 'CACGTG' in sequence or 'cacgtg' in sequence or 'CATGTG' in sequence or 'catgtg' in sequence:
           all_cacgtg = [m.start() for m in re.finditer("cacgtg", sequence, re.IGNORECASE)]
           # we are interested in the first occurance only. 
           start = int(header.split(":")[-1].split ("-")[0])
           start = start + all_cacgtg[0] +  3 - 15  
           end = start + 30 
           l_items = header.split(":")
           
           header = ":".join (l_items[0: len(l_items) - 1]) + ":" + str(start) + "-" + str(end) + ":wMotif"
           
        else:
            header = header + ":woMotif"
 
        print(header)
        print(sequence)
#windows.close()
