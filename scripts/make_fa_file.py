import sys

#fa_file = 'data_from_ReMap/max_annotated_windows.fa'
#file1 = open(fa_file)

for line in sys.stdin:
    l1 = line.strip()
    if 'cacgtg' in l1 or 'CACGTG' in l1 or 'wmotif' in l1:
        print(l1)
file1.close()
