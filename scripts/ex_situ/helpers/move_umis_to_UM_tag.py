# Zachary Chiang, Buenrostro Lab, Harvard University
# extracts UMI from split hairpin sequence
# indices depend on sequencing instrument, this is for Novaseq/Nextseq

import sys
import pysam

bam_file = pysam.AlignmentFile(sys.argv[1], "rb")
out_file = pysam.AlignmentFile(sys.argv[2], "wb", template=bam_file)

for read in bam_file.fetch(until_eof=True):

        header = read.query_name
        umi = header.split("+")[1]

	umi_ext = ""

	umi_ext = umi_ext + umi[98] # base 1
	umi_ext = umi_ext + umi[97] # base 2
	umi_ext = umi_ext + umi[96] # base 3
	umi_ext = umi_ext + umi[95] # base 4
	umi_ext = umi_ext + umi[94] # base 5
	umi_ext = umi_ext + umi[93] # base 6
	umi_ext = umi_ext + umi[74] # base 7
	umi_ext = umi_ext + umi[73] # base 8
	umi_ext = umi_ext + umi[72] # base 9
	umi_ext = umi_ext + umi[71] # base 10
	umi_ext = umi_ext + umi[70] # base 11
	umi_ext = umi_ext + umi[69] # base 12
	umi_ext = umi_ext + umi[52] # base 13
	umi_ext = umi_ext + umi[51] # base 14
	umi_ext = umi_ext + umi[50] # base 15
	umi_ext = umi_ext + umi[49] # base 16
	umi_ext = umi_ext + umi[48] # base 17
	umi_ext = umi_ext + umi[47] # base 18
	umi_ext = umi_ext + umi[30] # base 19
	umi_ext = umi_ext + umi[29] # base 20
	umi_ext = umi_ext + umi[28] # base 21
	umi_ext = umi_ext + umi[27] # base 22
    
        read.set_tag('UM',umi_ext)
        read.query_name = header.split("_")[0]
        out_file.write(read)
