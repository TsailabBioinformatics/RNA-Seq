#!/bin/sh
#SBATCH -J RNASeqDemoMain
#SBATCH --partition batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --array=1-2
#SBATCH --time=8:00:00
#SBATCH --mem=30gb
#SBATCH --mail-user=your_email@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
date

### File path and environment setup ###
masterfolder=$SLURM_SUBMIT_DIR
cd ${masterfolder}
mkdir clean
mkdir map
mkdir count
cleanfolder=${masterfolder}/clean
mapfolder=${masterfolder}/map
countfolder=${masterfolder}/count
readsfolder=${masterfolder}/fastq
echo $readsfolder
array_configure_file=file.list
sampleFile=`head -n $SLURM_ARRAY_TASK_ID ${array_configure_file} | tail -n1`
trimmo_path=/apps/eb/Trimmomatic/0.39-Java-11
rRNA_ref_path=/work/cjtlab/Database/rRNA_ref


genomefolder=${masterfolder}
### This is the clean and trim section ###
cd ${cleanfolder}
ml Trimmomatic
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads $SLURM_NTASKS_PER_NODE -phred33 \
${readsfolder}/${sampleFile}.R1.fastq ${readsfolder}/${sampleFile}.R2.fastq \
${sampleFile}_trimP_1.fq.gz ${sampleFile}_trimS_1.fq.gz ${sampleFile}_trimP_2.fq.gz ${sampleFile}_trimS_2.fq.gz \
ILLUMINACLIP:${trimmo_path}/adapters/TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:36 TOPHRED33 &>${sampleFile}_trim.log

# TODO investigate how important this part is
module load STAR
STAR --runThreadN $SLURM_NTASKS_PER_NODE --genomeDir ${rRNA_ref_path} \
--readFilesIn ${sampleFile}_trimP_1.fq.gz ${sampleFile}_trimP_2.fq.gz \
--readFilesCommand gunzip -c --outReadsUnmapped Fastx \
--outFileNamePrefix ${sampleFile}_STAR_ \
  --outFilterMismatchNmax 999 \
  --outFilterMismatchNoverLmax 0.4 \
  --outFilterScoreMinOverLread 0.4 \
  --outFilterMatchNminOverLread 0.4 \
--alignMatesGapMax 20 --alignIntronMax 20
STAR --runThreadN $SLURM_NTASKS_PER_NODE --genomeDir ${rRNA_ref_path} \
--readFilesIn ${sampleFile}_trimS_1.fq.gz \
--readFilesCommand gunzip -c --outReadsUnmapped Fastx \
--outFileNamePrefix ${sampleFile}_STAR_S1_ \
  --outFilterMismatchNmax 999 \
  --outFilterMismatchNoverLmax 0.4 \
  --outFilterScoreMinOverLread 0.4 \
  --outFilterMatchNminOverLread 0.4 \
--alignMatesGapMax 20 --alignIntronMax 20
STAR --runThreadN $SLURM_NTASKS_PER_NODE --genomeDir ${rRNA_ref_path} \
--readFilesIn ${sampleFile}_trimS_2.fq.gz \
--readFilesCommand gunzip -c --outReadsUnmapped Fastx \
--outFileNamePrefix ${sampleFile}_STAR_S2_ \
  --outFilterMismatchNmax 999 \
  --outFilterMismatchNoverLmax 0.4 \
  --outFilterScoreMinOverLread 0.4 \
  --outFilterMatchNminOverLread 0.4 \
--alignMatesGapMax 20 --alignIntronMax 20
printf "Trimming : " >>${sampleFile}_Final.log.txt
grep "^Input Read"  ${sampleFile}_trim.log>>${sampleFile}_Final.log.txt
printf "Mapping rRNA : " >>${sampleFile}_Final.log.txt
grep "Number of input read"  ${sampleFile}_STAR_Log.final.out  >>${sampleFile}_Final.log.txt
printf "Mapping rRNA : " >>${sampleFile}_Final.log.txt
grep "Uniquely mapped reads number"  ${sampleFile}_STAR_Log.final.out  >>${sampleFile}_Final.log.txt
printf "Mapping rRNA : " >>${sampleFile}_Final.log.txt
grep "Number of reads mapped to multiple loci"  ${sampleFile}_STAR_Log.final.out  >>${sampleFile}_Final.log.txt
printf "HPT : " >>${sampleFile}_Final.log.txt
grep "MARKER" ${sampleFile}_STAR_Aligned.out.sam |grep -v "^@"|grep  "MARKER_HPT"|cut -f 1|sort|uniq|wc -l>>${sampleFile}_Final.log.txt
printf "NPT : " >>${sampleFile}_Final.log.txt
grep "MARKER" ${sampleFile}_STAR_Aligned.out.sam |grep -v "^@"|grep  "MARKER_NPT"|cut -f 1|sort|uniq|wc -l>>${sampleFile}_Final.log.txt
printf "IRP : " >>${sampleFile}_Final.log.txt
grep "MARKER" ${sampleFile}_STAR_Aligned.out.sam |grep -v "^@"|grep  "MARKER_IRP"|cut -f 1|sort|uniq|wc -l>>${sampleFile}_Final.log.txt
sed '1~4 s|\s.\+$|/1|g' < ${sampleFile}_STAR_Unmapped.out.mate1 | paste - - - - |  sort -k1,1 -S 10G -T ./ -V | tr '\t' '\n' >${sampleFile}_clean_1.fq
sed '1~4 s|\s.\+$|/2|g' < ${sampleFile}_STAR_Unmapped.out.mate2 | paste - - - - |  sort -k1,1 -S 10G -T ./ -V | tr '\t' '\n' >${sampleFile}_clean_2.fq
sed -i '1~4 s|$|/1|g' ${sampleFile}_STAR_S1_Unmapped.out.mate1
sed -i '1~4 s|$|/1|g' ${sampleFile}_STAR_S2_Unmapped.out.mate1
mv ${sampleFile}_STAR_S1_Unmapped.out.mate1 ${sampleFile}_clean_S1.fq
sed '1~4 s|/1$|/2|g' < ${sampleFile}_STAR_S2_Unmapped.out.mate1  >${sampleFile}_clean_S2.fq
sed '1~4 s|/1$|/2|g' < ${sampleFile}_clean_S1.fq |  \
awk 'NR%4==1{print}NR%4==2{printf  "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n"}NR%4==3{print}NR%4==0{printf  "##############################\n"} '  > ${sampleFile}_clean_E2.fq
cat ${sampleFile}_STAR_S2_Unmapped.out.mate1 |  \
awk 'NR%4==1{print}NR%4==2{printf  "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n"}NR%4==3{print}NR%4==0{printf  "##############################\n"} '  > ${sampleFile}_clean_E1.fq
cat ${sampleFile}_clean_S1.fq >> ${sampleFile}_clean_1.fq
cat ${sampleFile}_clean_E2.fq >> ${sampleFile}_clean_2.fq
cat ${sampleFile}_clean_E1.fq >> ${sampleFile}_clean_1.fq
cat ${sampleFile}_clean_S2.fq >> ${sampleFile}_clean_2.fq
rm ${sampleFile}_STAR_SJ.out.tab
rm ${sampleFile}_STAR_Log.progress.out
rm ${sampleFile}_STAR_Log.out
rm ${sampleFile}_STAR_Log.final.out
rm ${sampleFile}_STAR_Aligned.out.sam
rm ${sampleFile}_STAR_S1_SJ.out.tab
rm ${sampleFile}_STAR_S1_Log.progress.out
rm ${sampleFile}_STAR_S1_Log.out
rm ${sampleFile}_STAR_S1_Log.final.out
rm ${sampleFile}_STAR_S1_Aligned.out.sam
rm ${sampleFile}_STAR_S2_SJ.out.tab
rm ${sampleFile}_STAR_S2_Log.progress.out
rm ${sampleFile}_STAR_S2_Log.out
rm ${sampleFile}_STAR_S2_Log.final.out
rm ${sampleFile}_STAR_S2_Aligned.out.sam
rm ${sampleFile}_clean_S1.fq
rm ${sampleFile}_clean_S2.fq
rm ${sampleFile}_clean_E1.fq
rm ${sampleFile}_clean_E2.fq
rm ${sampleFile}_STAR_S2_Unmapped.out.mate1
rm ${sampleFile}_STAR_Unmapped.out.mate1
rm ${sampleFile}_STAR_Unmapped.out.mate2
rm ${sampleFile}_trimP_1.fq.gz
rm ${sampleFile}_trimS_1.fq.gz
rm ${sampleFile}_trimP_2.fq.gz
rm ${sampleFile}_trimS_2.fq.gz
rm ${sampleFile}_trim.log
gzip -f ${sampleFile}_clean_1.fq
gzip -f ${sampleFile}_clean_2.fq


### Mapping with STAR ###
cd ${mapfolder}

### STAR mapping time
STAR --runThreadN $SLURM_NTASKS_PER_NODE \
--genomeDir ${genomefolder} --readFilesIn \
${cleanfolder}/${sampleFile}_clean_1.fq.gz ${cleanfolder}/${sampleFile}_clean_2.fq.gz --readFilesCommand gunzip -c \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix ${sampleFile}_ \
--alignMatesGapMax 20000 \
--alignIntronMax 10000 \
--outFilterScoreMinOverLread 0.1 \
--outFilterMatchNminOverLread 0.1
### Count by Subread ###
cd ${countfolder}
ml Subread
featureCounts -Q 2 -s 0 -T $SLURM_NTASKS_PER_NODE -p -C \
-a ${genomefolder}/gene.gtf \
-o ${sampleFile}_counts.txt ${mapfolder}/${sampleFile}_Aligned.sortedByCoord.out.bam


## get trim log
python ~/script/get_trim_sum.py ${masterfolder}/clean/
##Please notice you need to get the get_trim_sum.py from NGSclean github
## Get map log
grep "" ${masterfolder}/map/*Log.final.out > ${masterfolder}/all_mapping_logs.txt

date