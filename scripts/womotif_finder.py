#109349  chr2L   109514  109544  wMotif
#246920  chr2L   246911  246941  wMotif

wmotif_peaks = 'data_from_ReMap/wmotif_final_2.tsv'
w_peaks = []

file1 = open(wmotif_peaks)
for line in file1:
    line = line.strip().split("\t")
    w_peaks.append(line[0])
file1.close()

all_peaks = 'data_from_ReMap/womotif_final.tsv'
file2 = open(all_peaks)
for line in file2:
    l2 = line.strip().split("\t")
    for c in l2:
        if c in w_peaks:
            print(f"{l2[0]}\t{l2[1]}\t{l2[2]}\t{l2[3]}\t{l2[4]}")

        else:
            continue
file2.close()



