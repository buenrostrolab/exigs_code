
## 0. Set up environment

home_dir='/mnt/users/zack/projects/xgs_alignment'
fastq_name=$1
run=$2
data_dir=$4
genome='/mnt/genomes/bowtie2/hg38/hg38'

step=$3
num_cores=4
split_size=20000000 # 5M reads, use for NovaSeq (less files to check)

bin='/mnt/users/zack/projects/xgs_alignment/bin'
tmp=$home_dir/tmp/

cd $home_dir
mkdir -p $run

r1_path=`ls $data_dir/$fastq_name"_"*R1*.fastq.gz`
r2_path=`echo $r1_path | sed 's/R1/R2/g'`
r3_path=`echo $r1_path | sed 's/R1/R3/g'`

r1=$run"_R1"
r2=$run"_R2"
r3=$run"_R3"

## 1. Split FASTQ

mkdir -p $tmp/split/

if  [ $step -le 1 ]; then

	echo "$run - Spliting FASTQs (1)"
	bin/fastp -i $r1_path -o $tmp/split/$r1.fastq.gz -S $split_size --thread 1 -d 4 -A -G -L -Q 2>logs/split_R1.log & 
	bin/fastp -i $r2_path -o $tmp/split/$r2.fastq.gz -S $split_size --thread 1 -d 4 -A -G -L -Q 2>logs/split_R2.log &
	bin/fastp -i $r3_path -o $tmp/split/$r3.fastq.gz -S $split_size --thread 1 -d 4 -A -G -L -Q 2>logs/split_R3.log &
	wait

fi

ls $tmp/split/ | grep $run | grep "R1" | grep -P -o "^[0-9]{4}" > $run/split_list.txt
num_split=`tail $run/split_list.txt`

## 2. Extract barcode (.ext)

# not needed in this version

## 3. Move R2 barcode to header (.umi)

mkdir -p $tmp/umi/

if  [ $step -le 3 ]; then
	
	echo "$run - Moving R2 barcode to header (3)"
	parallel --will-cite --jobs $num_cores --colsep '\t' \
		python2 $bin/add_umis_from_R2.py \
		$tmp/split/{1}.$r1.fastq.gz $tmp/split/{1}.$r2.fastq.gz $tmp/split/{1}.$r3.fastq.gz \
		$tmp/umi/{1}.$r1.umi.fastq.gz $tmp/umi/{1}.$r2.umi.fastq.gz \
		2>logs/umi.log :::: $run/split_list.txt

fi

## 4. Trim adapters (.trim)

mkdir -p $tmp/trim/

if  [ $step -le 4 ]; then

	echo "$run - Trimming adapters (4)"
	parallel --will-cite --jobs $num_cores --colsep '\t' \
		python2 $bin/pyadapter_trim.py \
		-a $tmp/umi/{1}.$r1.umi.fastq.gz -b $tmp/umi/{1}.$r2.umi.fastq.gz \
		-c $tmp/trim/{1}.$r1.trim.fastq.gz -d $tmp/trim/{1}.$r2.trim.fastq.gz \
		1>logs/trim.log :::: $run/split_list.txt
fi

## 5. Alignment (.aln)

mkdir -p $tmp/aln/
mkdir -p logs/align/

if  [ $step -le 5 ]; then

	echo "$run - Aligning (5)"
	parallel --will-cite --jobs $num_cores --colsep '\t' \
		bowtie2 -X2000 -p1 --rg-id $run \
		-x $genome \
		--very-sensitive \
		-k5 \
		-1 $tmp/trim/{1}.$r1.trim.fastq.gz \
		-2 $tmp/trim/{1}.$r2.trim.fastq.gz '|' \
		samtools view -bS - -o $tmp/aln/{1}.$run.aln_all.bam \
		2>logs/align_all/$run.align_all.log :::: $run/split_list.txt

fi

## 6. Sort (.sort)

mkdir -p $tmp/sort/

if  [ $step -le 6 ]; then


	echo "$run - Sorting (6)"
	parallel --will-cite --jobs $num_cores --colsep '\t' \
		samtools sort $tmp/aln/{1}.$run.aln_all.bam -o $tmp/sort/{1}.$run.sort_all.bam \
		2>logs/sort.log :::: $run/split_list.txt
	parallel --will-cite --jobs $num_cores --colsep '\t' \
		samtools index $tmp/sort/{1}.$run.sort_all.bam \
		2>logs/index.log :::: $run/split_list.txt

fi

## 7. Filter (.flt)

# not needed in this version

## 8. Add barcode to UM tag (.tag)

mkdir -p $tmp/tag/

if  [ $step -le 8 ]; then

	echo "$run - Adding barcode to UM tag (8)"
	parallel --will-cite --jobs $num_cores --colsep '\t' \
    		python2 $bin/move_umis_to_UM_tag_nextseq.py \
		$tmp/sort/{1}.$run.sort_all.bam \
		$tmp/tag/{1}.$run.tag_all.bam \
		2>logs/tag.log :::: $run/split_list.txt
	parallel --will-cite --jobs $num_cores --colsep '\t' \
		samtools index $tmp/tag/{1}.$run.tag_all.bam \
		2>logs/index.log :::: $run/split_list.txt

fi

## 9. Merge all (.merge)

if  [ $step -le 9 ]; then
	
	echo "$run - Merging (9)"
	ls $tmp/tag/*$run*bam > $run/merge_list.txt	
	samtools merge -f -b $run/merge_list.txt --threads $num_cores $run/$run.merge_all.bam
	samtools index $run/$run.merge_all.bam
fi

exit
