extracted_data = 'data_from_ReMap/sample.bed'
ext_peaks = open(extracted_data)


#chr2L  5086    6115    Ice1:Schneider-2        1       .       5690    5691    224,168,252
#chr2L  5088    6114    Sce:S2-DRSC     2       .       5818    5819    244,242,133
#chr2L  5103    6142    DnaJ-1:Schneider-2      1       .       5897    5898    62,241,1
for line in ext_peaks:
    l1 = line.strip("\t")
    print(l1)

