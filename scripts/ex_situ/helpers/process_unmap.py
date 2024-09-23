# Zachary Chiang, Buenrostro Lab, Harvard University
# handles different types of unmapped reads (w/ either one or both read pairs not uniquely mapping) 

import sys
import pysam

bam_file = pysam.AlignmentFile(sys.argv[1], "rb")
unmap_output = open(sys.argv[2],'w')
unmap_table = open(sys.argv[3],'w')
read_pairs = {}
umi_count = {}

for read in bam_file.fetch(until_eof=True):

    if read.query_name not in read_pairs:
        read_pairs[read.query_name] = read
        continue
    else:
        r1 = read
        r2 = read_pairs.pop(read.query_name)

    chrom_1 = r1.reference_name
    start_1 = r1.reference_start
    cigar_1 = r1.cigar
    seq_1 = r1.query

    chrom_2 = r2.reference_name
    start_2 = r2.reference_start
    seq_2 = r2.query

    umi = r1.get_tag("UM")

    if chrom_1 == None and start_1 == -1:
        chrom_1 = "*"
        start_1 = "*"
    if chrom_2 == None and start_2 == -1:
        chrom_2 = "*"
        start_2 = "*"

    print >> unmap_table, umi + "\t" + chrom_1 + "\t" + str(start_1) + "\t" + chrom_2 + "\t" + str(start_2) + "\t" + seq_1 + "\t" + seq_2
    
    if umi not in umi_count:
        umi_count[umi] = 1
    else:
        umi_count[umi] += 1

for umi in sorted(umi_count):

    print >> unmap_output, "unmap\t0\t" + "\t" + umi + "\t" + str(umi_count[umi])

sys.exit()

