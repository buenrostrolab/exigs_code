# Zachary Chiang, Buenrostro Lab, Harvard University
# creates a .bed-esque file from uniquely and multi-mapped reads for overlap with repeat elements 

import sys

full_table = open(sys.argv[1],'r')
multimap_table = open(sys.argv[2],'r')

multimap_umis = {}

for line in full_table:

    column = line.rstrip().split('\t')

    chrom = column[0]
    pos = column[1]
    umi = column[2]

    if chrom == 'multi_map': 
        multimap_umis[umi] = 1
        continue
    elif chrom == 'multi_umi' or chrom == 'unmap':
        continue
    
    print "\t".join([chrom,pos,pos,umi])

for line in multimap_table:

    column = line.rstrip().split()

    umi = column[2]

    if umi in multimap_umis:

        chrom = column[0]
        pos = column[1]

        print "\t".join([chrom,pos,pos,umi])

