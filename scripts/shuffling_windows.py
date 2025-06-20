import sys

#>Max:Schneider-2@window1:12845::chr2L:12757-12787:woMotif
#GTATTTCTTTAGCAAGCTGCGCAGAAATTC
#>Max:Schneider-2@window2:12845::chr2L:12787-12817:wmotif
#GGCGGGGCACGTGTGGTGGTGCATTGCCAC

fa_file = 'data_from_ReMap/sample_max_annotated.fa'
file1 = open(fa_file)

for line in file1:
   if line.startswith('>'): 
      l1 = line.strip().split(':')
      l2 = l1[5].split('-')
      #print(l1)
      #print(l2)
   else:
      for line in file1:
         if 'cacgtg' in line or 'CACGTG' in line:
            l3 = line.strip().split()
            #l3 = list(line)
            substring = ['cacgtg'] or ['CACGTG']
            #index = l3.find(substring)
            #print(l3) 
            print(l3.index(substring))
