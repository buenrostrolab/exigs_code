# Zachary Chiang, Buenrostro Lab, Harvard University
# handles ex situ genomic reads that overlapped more than one repeat name, class, or family 

import sys

intersect_table = open(sys.argv[1],'r')
full_table = open(sys.argv[2],'r')

umi_name = {}
umi_class = {}
umi_family = {}

for line in intersect_table:

    column = line.rstrip().split("\t")
    
    umi = column[3]
    rep_name = column[8]
    rep_class = column[9]
    rep_family = column[10]

    if umi not in rep_name:
        umi_name[umi] = [rep_name]
        umi_class[umi] = [rep_class]
        umi_family[umi] = [rep_family]
    else:
        umi_name[umi].append(rep_name)
        umi_class[umi].append(rep_class)
        umi_family[umi].append(rep_family)

for line in full_table:

    column = line.rstrip().split("\t")

    umi = column[2]
    
    if umi in umi_name:
        column.append(";".join(set(umi_name[umi])))
        column.append(";".join(set(umi_class[umi])))
        column.append(";".join(set(umi_family[umi])))
        print "\t".join(column)
    else:
        column.append("N/A")
        column.append("N/A")
        column.append("N/A")
        print "\t".join(column)
