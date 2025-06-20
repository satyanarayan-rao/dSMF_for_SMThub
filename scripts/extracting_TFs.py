#chr2L   6464    6940    Sce:S2-DRSC     2       .       6697    6698    244,242,133
#chr2L   6470    6910    abd-A:S2-DRSC   1       .       6687    6688    16,200,243
#chr2L   6476    6912    Psc:third-instar,S2-DRSC        5       .       6711    6712    243,214,187

extracted_peaks = 'data_from_ReMap/extracted_s2_peaks.bed'
myTF = open(extracted_peaks)

for line in myTF:
    l1 = line.strip().split(":")[0].split("\t")
    print(l1[3])
