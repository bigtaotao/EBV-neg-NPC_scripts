#!/bin/sh
# *******************************************
# Script to perform TNscope variant calling
# using a matched paired Tumor+normal sample with fastq
# files named normal_1.fastq.gz, normal_2.fastq.gz
# tumor_1.fastq.gz, tumor_2.fastq.gz
# *******************************************


root="/BIGDATA1/sysucc_mszeng_3/lyt/wes_hknpc_last/batch_1"
for i in `cat fastq.list`;
#for i in 535253-A; #the test sample run
do
# Update with the fullpath location of your sample fastq
fastq_folder=$root"/clean"
tumor_fastq_1="${i}T_R1.fastq.gz"
tumor_fastq_2="${i}T_R2.fastq.gz" #If using Illumina paired data
tumor_sample="${i}T"
tumor_group="NPC_${i}T"
normal_fastq_1="${i}N_R1.fastq.gz"
normal_fastq_2="${i}N_R2.fastq.gz" #If using Illumina paired data
normal_sample="${i}N"
normal_group="NPC_${i}N"
platform="ILLUMINA"

# Update with the location of the reference data files
fasta="/BIGDATA1/sysucc_mszeng_3/Reference/hg19/hg19.fa"
dbsnp="/BIGDATA1/sysucc_mszeng_3/Reference/hg19/SNP_file/dbsnp_138.hg19.vcf"
known_Mills_indels="/BIGDATA1/sysucc_mszeng_3/Reference/hg19/SNP_file/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
known_1000G_indels="/BIGDATA1/sysucc_mszeng_3/Reference/hg19/SNP_file/1000G_phase1.indels.hg19.sites.vcf"

# Update with the location of the Sentieon software package and license file
#module load sentieon/201711
#release_dir=/BIGDATA1/app/sentieon/sentieon-genomics-201711.03
release_dir=$SENTIEON_INSTALL_DIR
# Other settings
nt=24 #number of threads to use in computation
workdir=$root"result/${i}" #Determine where the output files will be stored

# ******************************************
# 0. Setup
# ******************************************
mkdir -p $workdir
logfile=$workdir/run_tnscope.log
exec >$logfile 2>&1
cd $workdir

# ******************************************
# 1a. Mapping reads with BWA-MEM, sorting for tumor sample
# ******************************************
#The results of this call are dependent on the number of threads used. To have number of threads independent results, add chunk size option -K 10000000 
( $release_dir/bin/sentieon bwa mem -M -R "@RG\tID:$tumor_group\tSM:$tumor_sample\tPL:$platform" -t $nt -K 10000000 $fasta $fastq_folder/$tumor_fastq_1 $fastq_folder/$tumor_fastq_2 || echo -n 'error' ) | $release_dir/bin/sentieon util sort -o tumor_sorted.bam -t $nt --sam2bam -i -
# ******************************************
# 1b. Mapping reads with BWA-MEM, sorting for normal sample
# ******************************************
#The results of this call are dependent on the number of threads used. To have number of threads independent results, add chunk size option -K 10000000 
( $release_dir/bin/sentieon bwa mem -M -R "@RG\tID:$normal_group\tSM:$normal_sample\tPL:$platform" -t $nt -K 10000000 $fasta $fastq_folder/$normal_fastq_1 $fastq_folder/$normal_fastq_2 || echo -n 'error' ) | $release_dir/bin/sentieon util sort -o normal_sorted.bam -t $nt --sam2bam -i -

# ******************************************
# 2a. Metrics for tumor sample
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i tumor_sorted.bam --algo MeanQualityByCycle tumor_mq_metrics.txt --algo QualDistribution tumor_qd_metrics.txt --algo GCBias --summary tumor_gc_summary.txt tumor_gc_metrics.txt --algo AlignmentStat --adapter_seq '' tumor_aln_metrics.txt --algo InsertSizeMetricAlgo tumor_is_metrics.txt
$release_dir/bin/sentieon plot metrics -o tumor_metrics-report.pdf gc=tumor_gc_metrics.txt qd=tumor_qd_metrics.txt mq=tumor_mq_metrics.txt isize=tumor_is_metrics.txt
# ******************************************
# 2b. Metrics for normal sample
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i normal_sorted.bam --algo MeanQualityByCycle normal_mq_metrics.txt --algo QualDistribution normal_qd_metrics.txt --algo GCBias --summary normal_gc_summary.txt normal_gc_metrics.txt --algo AlignmentStat --adapter_seq '' normal_aln_metrics.txt --algo InsertSizeMetricAlgo normal_is_metrics.txt
$release_dir/bin/sentieon plot metrics -o normal_metrics-report.pdf gc=normal_gc_metrics.txt qd=normal_qd_metrics.txt mq=normal_mq_metrics.txt isize=normal_is_metrics.txt

# ******************************************
# 3a. Remove Duplicate Reads for tumor sample
# ******************************************
$release_dir/bin/sentieon driver -t $nt -i tumor_sorted.bam --algo LocusCollector --fun score_info tumor_score.txt
$release_dir/bin/sentieon driver -t $nt -i tumor_sorted.bam --algo Dedup --rmdup --score_info tumor_score.txt --metrics tumor_dedup_metrics.txt tumor_deduped.bam 
# ******************************************
# 3b. Remove Duplicate Reads for normal sample
# ******************************************
$release_dir/bin/sentieon driver -t $nt -i normal_sorted.bam --algo LocusCollector --fun score_info normal_score.txt
$release_dir/bin/sentieon driver -t $nt -i normal_sorted.bam --algo Dedup --rmdup --score_info normal_score.txt --metrics normal_dedup_metrics.txt normal_deduped.bam 

# ******************************************
# 4a. Indel realigner for tumor sample
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i tumor_deduped.bam --algo Realigner -k $known_Mills_indels -k $known_1000G_indels tumor_realigned.bam
# ******************************************
# 4b. Indel realigner for normal sample
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i normal_deduped.bam --algo Realigner -k $known_Mills_indels -k $known_1000G_indels normal_realigned.bam

# ******************************************
# 5a. Base recalibration for tumor sample
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i tumor_realigned.bam --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels tumor_recal_data.table
$release_dir/bin/sentieon driver -r $fasta -t $nt -i tumor_realigned.bam -q tumor_recal_data.table --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels tumor_recal_data.table.post --algo ReadWriter tumor_recalibrated.bam
$release_dir/bin/sentieon driver -t $nt --algo QualCal --plot --before tumor_recal_data.table --after tumor_recal_data.table.post tumor_recal.csv
$release_dir/bin/sentieon plot bqsr -o tumor_recal_plots.pdf tumor_recal.csv
# ******************************************
# 5b. Base recalibration for normal sample
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i normal_realigned.bam --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels normal_recal_data.table
$release_dir/bin/sentieon driver -r $fasta -t $nt -i normal_realigned.bam -q normal_recal_data.table --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels normal_recal_data.table.post --algo ReadWriter normal_recalibrated.bam
$release_dir/bin/sentieon driver -t $nt --algo QualCal --plot --before normal_recal_data.table --after normal_recal_data.table.post normal_recal.csv
$release_dir/bin/sentieon plot bqsr -o normal_recal_plots.pdf normal_recal.csv

# ******************************************
# 6. Corealignment of tumor and normal
# ******************************************
#$release_dir/bin/sentieon driver -r $fasta -t $nt -i tumor_realigned.bam -i normal_realigned.bam -q tumor_recal_data.table -q normal_recal_data.table --algo Realigner -k $known_Mills_indels -k $known_1000G_indels tn_corealigned.bam

# ******************************************
# 7. Somatic and Structural variant calling
# ******************************************
$release_dir/bin/sentieon driver -r $fasta -t $nt -i tumor_realigned.bam -i normal_realigned.bam -q tumor_recal_data.table -q normal_recal_data.table  --algo TNscope --tumor_sample $tumor_sample --normal_sample $normal_sample --dbsnp $dbsnp output_tnscope.vcf.gz
done
