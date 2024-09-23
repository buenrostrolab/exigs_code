# Zachary Chiang, Buenrostro Lab, Harvard University
# basic logic here is that if a low quality, primary alignment is grouped with a high quality alignment
# (i.e. it has the same UMI + genomic position) it is not a true low quality/multi-mapping read and should
# be removed otherwise the UMIs will collide when high/low/unmapped data are combined
# other considerations: usually primary alignment is the best but not always, so technically should check all
# secondary alignments, but by that logic we should have more than 4 secondary alignments too so *shrug*
# oh also we should probably increment the count in the high quality, primary alignment group but that's a
# pain in the ass, we can't use the all_pri file in it's place because we don't want legit low qual reads in there

import sys

low_pri_group = open(sys.argv[1],'r')
pri_group = open(sys.argv[2],'r')
low_all_group = open(sys.argv[3],'r')

low_pri_lib = {}
remove = {}

for line in low_pri_group:

    if "read_id" in line:
        continue

    column = line.rstrip().split()
    umi = column[6]
    count = int(column[7])
    
    if umi not in low_pri_lib:
        low_pri_lib[umi] = count

for line in pri_group:

    if "read_id" in line:
        continue

    column = line.rstrip().split()
    umi = column[6]
    count = int(column[7])
    
    if umi in low_pri_lib and count > low_pri_lib[umi]:
        remove[umi] = 1

for line in low_all_group:
    
    if "read_id" in line:
        continue

    column = line.rstrip().split()
    umi = column[6]

    if umi in remove:
        continue

    print(line.rstrip())
