# Zachary Chiang, Buenrostro Lab, Harvard University
# adds multi-mapped and unmapped entries to ex situ UMI table

import sys

umi_file = open(sys.argv[1],'r')
multimap_file = open(sys.argv[2],'r')
unmap_file = open(sys.argv[3],'r')

max_index = 0

# uniquely mapped, get max UMI index

for line in umi_file:

	column = line.rstrip().split()
	print line.rstrip()

	if int(column[4])>max_index:
		max_index = int(column[4])

# multi-mapped, skip reads with N

for line in multimap_file:

        column = line.rstrip().split()
        if "N" in column[2]:
            continue
        max_index += 1
        column.append(str(max_index))
        print "\t".join(column)

# unmapped

for line in unmap_file:

        column = line.rstrip().split()
        max_index += 1
        column.append(str(max_index))
        print "\t".join(column)


	
