fastq_name='S83_1_Ad205_S1'
run='230530_sample83_1_novaseq_1'
fastq_name='S83_2_Ad206_S2'
run='230530_sample83_2_novaseq_1'
fastq_name='S83_3_Ad207_S3'
run='230530_sample83_3_novaseq_1'
fastq_name='S83_4_Ad208_S4'
run='230530_sample83_4_novaseq_1'

#fastq_name='S84_1_Ad213_S9'
#run='230601_sample84_1_novaseq_1'
#fastq_name='S84_2_Ad214_S10'
#run='230601_sample84_2_novaseq_1'
#fastq_name='S84_3_Ad215_S11'
#run='230601_sample84_3_novaseq_1'
#fastq_name='S84_4_Ad216_S12'
#run='230601_sample84_4_novaseq_1'

#fastq_name='S85_1_Ad217_S13'
#run='230601_sample85_1_novaseq_1'
#fastq_name='S85_2_Ad218_S14'
#run='230601_sample85_2_novaseq_1'
#fastq_name='S85_3_Ad219_S15'
#run='230601_sample85_3_novaseq_1'
#fastq_name='S85_4_Ad220_S16'
#run='230601_sample85_4_novaseq_1'

#fastq_name='S107_1_Ad233_S17'
#run='231004_sample107_1_novaseq_1'
#fastq_name='S107_2_Ad234_S18'
#run='231004_sample107_2_novaseq_1'
#fastq_name='S107_3_Ad235_S19'
#run='231004_sample107_3_novaseq_1'
#fastq_name='S107_4_Ad236_S20'
#run='231004_sample107_4_novaseq_1'

#fastq_name='S110_1_Ad237_S21'
#run='231002_sample110_1_novaseq_1'
#fastq_name='S110_2_Ad238_S22'
#run='231002_sample110_2_novaseq_1'

step=1
data_dir='/mnt/RawData/NovaSeq/20231122_IGS/fastqs/'
nohup bash ex_situ_alignment.sh $fastq_name $run $step $data_dir &

step=1
nohup bash ex_situ_dedup.sh $run $run '' $step &

