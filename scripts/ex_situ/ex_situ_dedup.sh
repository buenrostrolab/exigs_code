run=$1
bam1=$2
bam2=$3
step=$4

mkdir -p $run

# 1. Merge BAMs from different sequencing runs

if  [ $step -le 1 ]; then
	
	echo "$run - Merge BAMs from different sequencing runs (1)"
	
	samtools merge $run/$run.merge_all.bam $bam1/$bam1.merge_all.bam $bam2/$bam2.merge_all.bam
	samtools index $run/$run.merge_all.bam

fi

# 2. Split merged BAM into high quality, low quality, and unmapped

if  [ $step -le 2 ]; then
	
	echo "$run - Split BAM into high and low quality (2)"

	samtools view $run/$run.merge_all.bam -b -q 30 -U $run/$run.low_all.bam > $run/$run.high_all.bam
	samtools index $run/$run.low_all.bam
	
	samtools1.16 view $run/$run.merge_all.bam -f 12 -F 256 -b > $run/$run.unmap.bam # forget why it needs a different version but it does
	samtools1.16 index $run/$run.unmap.bam

fi

# 3. Select primary alignment reads 

if  [ $step -le 3 ]; then
	
	echo "$run - Select primary alignment reads (3)"

	samtools view $run/$run.merge_all.bam -b -F 256 > $run/$run.pri.bam
	samtools index $run/$run.pri.bam
	
	samtools view $run/$run.high_all.bam -b -F 256 > $run/$run.high_pri.bam
	samtools index $run/$run.high_pri.bam

	samtools view $run/$run.low_all.bam -b -F 256 > $run/$run.low_pri.bam
	samtools index $run/$run.low_pri.bam

fi

# 4. Run UMI-tools group

if  [ $step -le 4 ]; then

	echo "$run - Run UMI-tools group (4)"
	
	umi_tools group -I $run/$run.pri.bam --extract-umi-method=tag --umi-tag=UM:Z --method=cluster --edit-distance-threshold=2 \
		--group-out=$run/$run.pri.group.tsv --log=$run/$run.pri.group.log --paired  --output-bam -S $run/$run.pri.group.bam
	umi_tools group -I $run/$run.high_pri.bam --extract-umi-method=tag --umi-tag=UM:Z --method=cluster --edit-distance-threshold=2 \
		--group-out=$run/$run.high_pri.group.tsv --log=$run/$run.high_pri.group.log --paired  --output-bam -S $run/$run.high_pri.group.bam
	umi_tools group -I $run/$run.low_pri.bam --extract-umi-method=tag --umi-tag=UM:Z --method=cluster --edit-distance-threshold=2 \
		--group-out=$run/$run.low_pri.group.tsv --log=$run/$run.low_pri.group.log --paired  --output-bam -S $run/$run.low_pri.group.bam
	umi_tools group -I $run/$run.low_all.bam --extract-umi-method=tag --umi-tag=UM:Z --method=cluster --edit-distance-threshold=2 \
		--group-out=$run/$run.low_all.group.tsv --log=$run/$run.low_all.group.log --paired  --output-bam -S $run/$run.low_all.group.bam

fi

# 5. Collapse reads by group

if  [ $step -le 5 ]; then

	echo "$run - Collapse reads by group (5)"

	cat $run/$run.high_pri.group.tsv | grep -v "contig" | cut -f2,3,7,8,9 | sort | uniq | sort -k5,5 -n | grep -P -v "N" > $run/$run.high.txt

	python bin/filter_low_qual.py $run/$run.low_pri.group.tsv $run/$run.pri.group.tsv $run/$run.low_all.group.tsv > $run/$run.low_filt.group.tsv
	cat $run/$run.low_filt.group.tsv | grep -v "contig" | sort -k7,7 >  $run/$run.low_filt.sort.tsv
	python2 bin/collapse_multimap.py $run/$run.low_filt.sort.tsv $run/$run.multimap.txt $run/$run.multimap_table.txt

	python2 bin/process_unmap.py $run/$run.unmap.bam $run/$run.unmap.txt $run/$run.unmap_table.txt
	python2 bin/add_multi_and_unmapped.py $run/$run.high.txt $run/$run.multimap.txt $run/$run.unmap.txt > $run/$run.reads_all.txt

fi

# 5. Index swapping filters

if  [ $step -le 6 ]; then

	echo "$run - Index swapping filters (6)"

	python2 bin/pos_filter.py < $run/$run.reads_all.txt | sort -k3,3 > $run/$run.reads_filt_pos_all.txt
	python2 bin/umi_filter.py $run/$run.reads_filt_pos_all.txt $run/$run.reads_filt_umi_all.txt $run/$run.reads_multi_umi_all.txt
	cat $run/$run.reads_filt_umi_all.txt | sort -n -k5,5 | grep -v "N" > $run/$run.reads_filt_both_all.txt

fi

# 7. Convert UMI to colospace

if  [ $step -le 7 ]; then
	
	echo "$run - Convert UMI to colorspace (7)"

	python2 bin/seq2umi_sbs_nextseq.py < $run/$run.reads_filt_both_all.txt > $run/$run.umis_all.txt

fi

# 8. Add fragment lengths

if  [ $step -le 8 ]; then

	echo "$run - Add fragment lengths (8)"

	samtools index $run/$run.high_pri.group.bam
	python2 bin/add_frag_size.py $run/$run.umis_all.txt $run/$run.high_pri.group.bam > $run/$run.umis_all.frag_len.txt

fi

# 8. Add repeat annotations

if  [ $step -le 9 ]; then

	echo "$run - Add repeat annotations (9)"

	python2 bin/get_all_genomic_locs.py $run/$run.umis_all.frag_len.txt $run/$run.multimap_table.txt > $run/$run.all_mapped.bed
	bedtools intersect -a $run/$run.all_mapped.bed -b hg38_repeat_masker.bed -wa -wb > $run/$run.all_intersect.bed 
	python2 bin/collapse_repeat_annotations.py $run/$run.all_intersect.bed $run/$run.umis_all.frag_len.txt > $run/$run.umis_all.repeats.txt

fi
exit


