# Zachary Chiang, Buenrostro Lab, Harvard University
# adds genomic fragment size to table of ex situ UMIs 

import sys
import numpy as np
import math
import pysam

group_file = open(sys.argv[1],'r')
bam_file = pysam.AlignmentFile(sys.argv[2],'rb')

groups = {}
frag_sizes = {}

# put all valid UMI groups in dict

for line in group_file:

    column = line.rstrip().split("\t")
    groups[column[2]] = 1

reads = {}

# loop through .bam file

for read in bam_file.fetch():

    read_name = read.query_name

    # handle paired-end reads

    if read_name not in reads:
        reads[read_name] = read
        continue
    else:
        r1 = reads.pop(read_name)
        r2 = read

    # save frag len and UMI group

    frag_len = abs(r1.template_length) 

    umi_group = ""
    if r1.has_tag("BX"):
        umi_group = r1.get_tag("BX")
    if r2.has_tag("BX"):
        umi_group = r2.get_tag("BX")

    if umi_group in groups:
        frag_sizes[umi_group] = frag_len


group_file.seek(0,0)

# loop thorugh all valid UMI groups and append frag len

for line in group_file:

    column = line.rstrip().split()
    umi_group = column[2]

    frag_len = -1
    if umi_group in frag_sizes:
        frag_len = frag_sizes[umi_group]
            
    column.append(str(frag_len))


    print("\t".join(column))

