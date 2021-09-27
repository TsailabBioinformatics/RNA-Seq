"""
This program creates a directory structure to
store and organize different species' genome data
"""
import os
import shutil
import pandas as pd
import numpy as np
import subprocess

# TODO - add resume option, so that no need to re-run the entire process, especially the copy part

# TODO - create directories for user, if they don't exist
root_directory = "./data/"
sub_directories = ["species", "clean", "count",
                   "fastq", "mapping", "miscelaneous", "reference"]

table = pd.read_csv("./input/sample_table.csv")
directory_ids = list(set(table["folder_id"].to_numpy()))
myID = table["myID"][0]
start_sh_list = []


def clean_csv_data():
    '''
    takes in a csv and cleans whitespace 
    '''
    # TODO
    pass


def make_directories():
    """
    Creates directory structure to store genome data

    @param species_name: the species to sequence
    """

    print("### creating directories ###")
    # TODO - table["species"] can be used to write readme files for human to read in each directory
    species = table["species"]


    os.mkdir("data")
    # TODO - skip if directory already exists
    for id in directory_ids:
        print(f"creating directory {id}")
        os.mkdir(root_directory + id)
        for dir in sub_directories:
            print(f"creating {root_directory + id + '/' + dir}")
            os.mkdir(root_directory + id + "/" + dir + "/")
    # we can skip this part for now
    # for i in range(len(table)):
    #     os.mkdir(root_directory +
    #              table.at[i, "folder_id"] + "/species/" + table.at[i, "sample_name"])


def extract_data():
    """
    As of now, extracts data from Phytozome and
    moves it to specific species directory

    @param species_name: the species to sequence
    """
    print("### extracting data ###")
    # copy genome reference, fasta and gff3 respectively
    for i in range(len(table)):
        print(f'copying {table.at[i, "genome_fa"]} to {"./data/" + table.at[i,"folder_id"] + "/reference/" + table.at[i, "folder_id"] + "_genome.fa"}')
        shutil.copyfile(table.at[i, "genome_fa"], "./data/" + table.at[i,
                                                                       "folder_id"] + "/reference/" + table.at[i, "folder_id"] + "_genome.fa")
    for i in range(len(table)):
        print(f'copying {table.at[i, "gene_gff3"]} to {"./data/" + table.at[i,"folder_id"] + "/reference/" + table.at[i, "folder_id"] + ".gff3"}')
        shutil.copyfile(table.at[i, "gene_gff3"], "./data/" + table.at[i,
                                                                       "folder_id"] + "/reference/" + table.at[i, "folder_id"] + ".gff3")
    # copy fastq files
    for i in range(len(table)):
        print(f'copying {table.at[i, "fastq_dir"] + table.at[i, "sample_name"] + ".R1.fastq"}')
        shutil.copyfile(table.at[i, "fastq_dir"] + table.at[i, "sample_name"] + ".R1.fastq",
                        "./data/" + table.at[i, "folder_id"] + "/fastq/" + table.at[i, "sample_name"] + ".R1.fastq")
        shutil.copyfile(table.at[i, "fastq_dir"] + table.at[i, "sample_name"] + ".R2.fastq",
                        "./data/" + table.at[i, "folder_id"] + "/fastq/" + table.at[i, "sample_name"] + ".R2.fastq")


def generate_scripts():
    # TODO - handle emails
    # chen's note, I add a column in the table called "myID" to handle this, so that different jobs can even be sent to different users ,if necessary
    
    print("### generating scripts ###")
    for element in directory_ids:
        """
        makes a script to prepare for array mapping
        """
        element = element.strip()
        start_sh_list.append("sbatch ./data/" + element + "/start.sh")
        with open("./data/" + element + "/start.sh", "wt") as script:
            script.write("#!/bin/sh\n")
            script.write("#SBATCH -J " + element + "_script\n")
            script.write("#SBATCH --partition batch\n")
            script.write("#SBATCH --nodes=1\n")
            script.write("#SBATCH --ntasks-per-node=12\n")
            script.write("#SBATCH --time=8:00:00\n")
            script.write("#SBATCH --mem=30gb\n")
            script.write(f"#SBATCH --mail-user={myID}@uga.edu\n")
            script.write("#SBATCH --mail-type=BEGIN,END,FAIL\n")
            script.write("date\n\n")

            script.write("module load gffread\n")
            script.write("gffread reference/" + element +
                         ".gff3 -T -o reference/" + element + "_gene.gtf\n")
            script.write("ml STAR\n")
            script.write("STAR \\ \n")
            script.write("--runThreadN 8 \\ \n")
            script.write("--runMode genomeGenerate \\ \n")
            script.write("--genomeSAindexNbases 13 \\ \n")
            script.write("--genomeDir . \\ \n")
            script.write("--genomeFastaFiles " + element + "_genome.fa \\ \n")
            script.write("--sjdbGTFfile reference/" + element + "_gene.gtf \n")
            script.write("sbatch " + "./data/" + element + "/map.sh\n")

        # handle array size and generate file.list
        sample_list = table[table["folder_id"]
                            == element].sample_name.to_list()
        requested_array_size = len(sample_list)
        with open("./data/" + element + "/file.list", "wt") as file_list:
            file_list.write("\n".join(sample_list))

        with open("./data/" + element + "/map.sh", "wt") as mapping_script:
            """
            makes array mappign script
            """
            # writing header
            mapping_script.write(f"""#!/bin/sh
#SBATCH -J {element}_array_mapping
#SBATCH --partition batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --array=1-{requested_array_size}
#SBATCH --time=8:00:00
#SBATCH --mem=30gb
#SBATCH --mail-user={myID}@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
date""")

            # TODO - reduce the unecessary reassignment of variables, on the python level or bash level?
            mapping_script.write("""
### File path and environment setup ###
masterfolder=$SLURM_SUBMIT_DIR
cd ${masterfolder}
cleanfolder=${masterfolder}/clean
mapfolder=${masterfolder}/map
countfolder=${masterfolder}/count
readsfolder=${masterfolder}/fastq
genomefolder=${masterfolder}/reference
echo $readsfolder
array_configure_file=file.list
sampleFile=`head -n $SLURM_ARRAY_TASK_ID ${array_configure_file} | tail -n1`
trimmo_path=/apps/eb/Trimmomatic/0.39-Java-11
rRNA_ref_path=/work/cjtlab/Database/rRNA_ref

### module load ###
# TODO - do we want to add version here?
module load STAR
module load trimmomatic
module load Subread
""")

# clean and trim reads
            mapping_script.write("""
### This is the clean and trim section ###
cd ${{cleanfolder}}
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads $SLURM_NTASKS_PER_NODE -phred33 \\ 
${{readsfolder}}/${{sampleFile}}.R1.fastq ${{readsfolder}}/${{sampleFile}}.R2.fastq \\
${{sampleFile}}_trimP_1.fq.gz ${{sampleFile}}_trimS_1.fq.gz ${{sampleFile}}_trimP_2.fq.gz ${{sampleFile}}_trimS_2.fq.gz \\ 
ILLUMINACLIP:${{trimmo_path}}/adapters/TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:36 TOPHRED33 &>${{sampleFile}}_trim.log
""")

# rRNA removal part
            mapping_script.write("""
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
""")

            mapping_script.write(f"""
### Mapping with STAR ###
cd ${{mapfolder}}

### STAR mapping ###
STAR --runThreadN $SLURM_NTASKS_PER_NODE \
--genomeDir ${{genomefolder}} --readFilesIn \
${{cleanfolder}}/${{sampleFile}}_clean_1.fq.gz ${{cleanfolder}}/${{sampleFile}}_clean_2.fq.gz --readFilesCommand gunzip -c \\
--outSAMtype BAM SortedByCoordinate \\
--outFileNamePrefix ${{sampleFile}}_ \\
--alignMatesGapMax 20000 \\
--alignIntronMax 10000 \\
--outFilterScoreMinOverLread 0.1 \\
--outFilterMatchNminOverLread 0.1

### Count by Subread ###
cd ${{countfolder}}
featureCounts -Q 2 -s 0 -T $SLURM_NTASKS_PER_NODE -p -C \\
-a ${{genomefolder}}{element}_gene.gtf \\
-o ${{sampleFile}}_counts.txt ${{mapfolder}}/${{sampleFile}}_Aligned.sortedByCoord.out.bam

## get trim log
python ./helper_scripts/get_trim_sum.py ${{masterfolder}}/clean/
## Get map log
grep "" ${{masterfolder}}/map/*Log.final.out > ${{masterfolder}}/all_mapping_logs.txt
            
date
""")

def submit_scripts():
    """
    Submit all scripts to the cluster
    """
    print("### Submitting scripts to the cluster ###")
    with open(f"submit_all_scripts.sh", "w") as submit_script:
        submit_script.write("\n".join(start_sh_list))
    bashCommand  = f"bash ./submit_all_scripts.sh"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    print(output)
    print(error)

# TODO - print out progress
make_directories()
extract_data()
generate_scripts()
submit_scripts()