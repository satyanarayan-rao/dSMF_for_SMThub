string_length = 'data_from_ReMap/dm6_windows.fa'
mystring = open(string_length)

for line in mystring:
    l1 = line.strip()
    print(len(l1))
