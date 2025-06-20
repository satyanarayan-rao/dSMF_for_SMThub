#extract the lines that contain specific pattern

# Python can not read directly from file. You have to open the file and read line by line

# 1. You will read each line from s2_names - trim the "\n" and keep them in a list
# 2. for each line in the big file `remap2022_nr_macs2_dm6_v1_0.bed`, you will extract the column with system name and then extract the system. 
# 3. check for that system in the list that you prepared in 1. 
# 4. if the name is there in the list, you will print that line otherwise continue...

#for x in [/home/aashna_b.iitr/workplace/projects/SMThub/data_from_ReMap/remap2022_nr_macs2_dm6_v1_0.bed]
#	if x in [//home/aashna_b.iitr/workplace/projects/footprint_pipeline/metadata/s2_names] > extracted_from_remap_data.append [x]


#print ("done")


s2_names = ("data_from_ReMap/s2_names.txt") 
file_1 = open(s2_names)

s2_cell_lines_name = []

for line in file_1:
    line = line.strip().strip('"')
    s2_cell_lines_name.append(line)
file_1.close()


remap2022_nr_macs2_dm6_v1_0 = ("data_from_ReMap/remap2022_nr_macs2_dm6_v1_0.bed")
file_2 = open(remap2022_nr_macs2_dm6_v1_0)

for line in file_2:
    l_items = line.split("\t")
    cell_line =  l_items[3].split(":")[-1].split(",")
    for c in cell_line:
        if c in s2_cell_lines_name:
            print (line.strip())
            break
        else: 
            continue
 
file_2.close()
#remap2022_nr_macs2_dm6_v1_0 = ("data_from_Remap/remap2022_nr_macs2_dm6_v1_0.bed")
#file_2 = open(remap2022_nr_macs2_dm6_v1_0)
#
#spl_file = file_2.split("\t")
#spl_col = spl_file(":")
#
#for line in spl_col
#    if line == 
#	#if line in spl_coln 
##remap2022_nr_macs2_dm6_v1_0.bed

