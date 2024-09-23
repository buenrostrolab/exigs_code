# Zachary Chiang, Buenrostro Lab, Harvard University
# one of two index swapping filters, removes UMI-genomic location pairs based on patterns of PCR duplicates 

import sys

threshold = 0.25
current_location = ""

umi_dic = {}
umi_list = []

def consolidate(umi_list):

    total = sum([tup[1] for tup in umi_list])
    max_count = umi_list[0][1]
    for umi in umi_list:
        if umi[1] >= 0.5*total or (umi[1] > 0.2*total and total > 5) or umi[1] >= max_count:
            print "\t".join(umi_dic[(umi[0])])

for line in sys.stdin:

    column = line.rstrip().split()

    if column[0] == "multi_map" or column[0] == "unmap" or column[0] == "multi_umi":
        print line.rstrip()
        continue

    location = column[0] + "\t"  + column[1]
    umi = column[2]
    abundance = int(column[3])
    umi_id = column[4]

    if location == current_location:
	
	umi_list.append((umi, abundance))
        if umi not in umi_dic or (umi in umi_dic and abundance > int(umi_dic[umi][3])):
	    umi_dic[umi] = column

    else:

	if current_location != "":

	    umi_list = sorted(umi_list,key=lambda x:x[1],reverse=True)
	    consolidate(umi_list)

	current_location = location
	umi_list = [(umi, abundance)]
	umi_dic = {}
	umi_dic[umi] = column
