import sys
#we need to create the exact windows of 30 bp of chIP peaks (
# we can create different if conditions for when the no. of bases are not present in multiple of 30;
# i) we can either ignore the bases on the left side or on the right side
#ii) we can also ignore the equal no. of bases on the each side
#we need to annotate the different windows in ideal way, which could be w1, w2... or .1, .2... or windows1, windows2...

#extracted_data = 'data_from_ReMap/sample.bed'
#ext_peaks = open(extracted_data)


#chr2L	5086	6115	Ice1:Schneider-2	1	.	5690	5691	224,168,252
#chr2L	5088	6114	Sce:S2-DRSC	2	.	5818	5819	244,242,133
#chr2L	5103	6142	DnaJ-1:Schneider-2	1	.	5897	5898	62,241,1

#4       chr2L   5086    6115    Ice1:Schneider-2        1       .       5690    5691    224,168,252
#5       chr2L   5088    6114    Sce:S2-DRSC     2       .       5818    5819    244,242,133
#6       chr2L   5103    6142    DnaJ-1:Schneider-2      1       .       5897    5898    62,241,1


for line in sys.stdin:
    l1 = line.strip().split("\t")
    #print (l1)
    start = int(l1[1])
    end = int(l1[2])
    N = end-start 
    r = N%30
    a=0        
    if r==0:
        for i in range(start,end,30):
            a+=1    
            print(f"{l1[0]}\t{i}\t{i+30}\t{l1[3]}@window{a}")
        
        
    else:
       q = r%2
       p = int (r/2)
       if q==0:
           start = int(l1[1]) + p 
           end = int(l1[2]) - p
           for i in range(start,end,30):
               a+=1
               print(f"{l1[0]}\t{i}\t{i+30}\t{l1[3]}@window{a}") 

       elif q!=0:
           start = int(l1[1])+ p + 1 
           end = int(l1[2]) - p
           for i in range(start,end,30):
               a+=1
               print(f"{l1[0]}\t{i}\t{i+30}\t{l1[3]}@window{a}")
