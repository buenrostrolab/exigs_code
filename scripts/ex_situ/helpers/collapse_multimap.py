# Zachary Chiang, Buenrostro Lab, Harvard University
# handles ex situ UMIs with multi-mapped genomic fragments
# outputs a table of ex situ UMIs and a lookup table of all multi-mapped reads

import sys

reads = {}
frag_list = []
current_umi = ""

input_file = open(sys.argv[1],'r')
multimap_output = open(sys.argv[2],'w')
multimap_table = open(sys.argv[3],'w')

for line in input_file:

	column = line.rstrip().split()

        if column[0] == "read_id":
            continue

        read = ":".join(column[0].split(":")[-3:])
	chrom = column[1]
        pos = column[2]
	umi = column[6]
	abundance = column[7]
	umi_id = column[8]

	if umi == current_umi:

		frag_list.append((chrom,pos,umi,abundance,umi_id))
                if read not in reads:
                    reads[read] = 1

                if (chrom,pos,umi,abundance,umi_id) not in frag_reads:
                    frag_reads[(chrom,pos,umi,abundance,umi_id)] = [read]
                else:
                    frag_reads[(chrom,pos,umi,abundance,umi_id)].append(read)

	else:
		if current_umi != "":
                    uniq_frags = set(frag_list)
		    sort_frags = sorted(uniq_frags,key=lambda x:x[4])

                    for frag in sort_frags:
                        frag_read = ";".join(frag_reads[frag])
                        print >> multimap_table, "\t".join(frag) + "\t" + frag_read
		    
                    print >> multimap_output, "multi_map\t0\t" + current_umi + "\t" + str(len(reads))

		current_umi = umi
                reads = {}
                frag_reads = {}
                frag_list = []
		
                frag_list.append((chrom,pos,umi,abundance,umi_id))

                if read not in reads:
                    reads[read] = 1

                if (chrom,pos,umi,abundance,umi_id) not in frag_reads:
                    frag_reads[(chrom,pos,umi,abundance,umi_id)] = [read]
                else:
                    frag_reads[(chrom,pos,umi,abundance,umi_id)].append(read)

