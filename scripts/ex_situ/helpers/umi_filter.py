# Zachary Chiang, Buenrostro Lab, Harvard University
# one of two index swapping filters, removes UMI-genomic location pairs based on patterns of PCR duplicates 

import sys

threshold = 0.25
distance = 500
current_umi = ""

location_dic = {}
abundance_dic = {}
location_list = []

singleton_count = 0
non_singletons = []

input_file = open(sys.argv[1],'r')
multiumi_output = open(sys.argv[2],'w')
multiumi_table = open(sys.argv[3],'w')

def consolidate(location_list):

	location_list2 = location_list
	
        #if len(location_list) > 1:
        #    print location_list

	for i in range(len(location_list)):
		for j in range(len(location_list2)):
					
			chr_i = location_list[i][0].split()[0]
			pos_i = int(location_list[i][0].split()[1])
			num_i = int(location_list[i][1])
			uid_i = int(location_list[i][2])
			chr_j = location_list[j][0].split()[0]
			pos_j = int(location_list[j][0].split()[1])
			num_j = int(location_list[j][1])
			uid_j = int(location_list[j][2])
			#print i,chr_i,pos_j,num_i,uid_i,j,chr_j, pos_j,num_j,uid_j
			if i != j and chr_i == chr_j and abs(pos_i - pos_j) < distance and (num_i > num_j or (num_i == num_j and uid_i < uid_j)) and num_i > 0 and num_j > 0:
				location_list[i] = (location_list[i][0], num_i+num_j, uid_i)
				location_list[j] = (location_list[j][0], 0, uid_j)
				
				#print "hi", location_list
        
        #if len(location_list) > 1:
        #    print location_list

        location_list_filt = []

	total = sum([tup[1] for tup in location_list])
        max_count = max([tup[1] for tup in location_list])
        total_filt = 0
        first_index = location_list[0][2]
	for location in location_list:
		if location[1] >= 0.5*total or (location[1] > 0.2*total and total > 5 or location[1] == max_count):
			#print location
			location_dic[(location[0],str(location[2]))][3] = str(location[1])
			#print "\t".join(location_dic[(location[0],str(location[2]))])
                        location_list_filt.append("\t".join(location_dic[(location[0],str(location[2]))]))
                        total_filt += location[1]
        
        if len(location_list_filt) > 1:
            for location in location_list_filt:
                    print >> multiumi_table, location
            print >> multiumi_output, "multi_umi\t0\t" + current_umi + "\t" + str(total_filt) + "\t" + str(first_index)
        else:
            print >> multiumi_output, location_list_filt[0]

for line in input_file:

	column = line.rstrip().split()
	location = column[0] + "\t"  + column[1]
	umi = column[2]
	abundance = int(column[3])
	umi_id = column[4]

	if umi == current_umi:
	
		location_list.append((location, abundance, umi_id))
		location_dic[(location,umi_id)] = column

	else:

		if current_umi != "":

			#print sorted(location_list,key=lambda x:x[1],reverse=True)
			#if current_umi == "TAGCCGATGTACCTATGTTT":
				#print location_list
			location_list = sorted(location_list,key=lambda x:x[1],reverse=True)
			consolidate(location_list)


		current_umi = umi
		location_list = [(location, abundance, umi_id)]
		location_dic = {}
		location_dic[(location,umi_id)] = column
