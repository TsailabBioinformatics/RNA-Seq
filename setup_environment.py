"""
This program creates a directory structure to
store and organize different species' genome data
"""
import os
import shutil
import pandas as pd
import numpy as np

root_directory = "./data/"
sub_directories = ["species", "clean", "count", "fastq", "mapping", "miscelaneous", "reference"]

table = pd.read_csv("./input/sample_table.csv")
directory_ids = list(set(table["folder_id"].to_numpy()))

def clean_csv_data():
    '''
    takes in a csv and cleans whitespace 
    '''
    #TODO 
    pass

def make_directories():
    """
    Creates directory structure to store genome data

    @param species_name: the species to sequence
    """
    species = table["species"]

    for id in directory_ids:
        os.mkdir(root_directory + id)
        for dir in sub_directories:
            os.mkdir(root_directory + id + "/" + dir + "/")
    for i in range(len(table)):
        os.mkdir(root_directory + table.at[i, "folder_id"] + "/species/" + table.at[i, "sample_name"])


def extract_data():
    """
    As of now, extracts data from Phytozome and
    moves it to specific species directory

    @param species_name: the species to sequence
    """
    for i in range(len(table)):
        shutil.copyfile(table.at[i, "genome_fa"], "./data/" + table.at[i, "folder_id"] + "/reference/" + table.at[i, "folder_id"]  + "_genome.fa")
    for i in range(len(table)):
        shutil.copyfile(table.at[i, "gene_gff3"], "./data/" + table.at[i, "folder_id"] + "/reference/" + table.at[i, "folder_id"] + ".gff3")

def generate_scripts():
    #TODO - handle emails 


    for element in directory_ids:
        """
        makes a script to prepare for array mapping
        """
        with open("./data/" + element + "/start.sh", "wt") as script:
            script.write("#!/bin/sh\n")
            script.write("#SBATCH -J " + element + "_script\n")
            script.write("#SBATCH --partition batch\n")
            script.write("#SBATCH --nodes=1\n")
            script.write("#SBATCH --ntasks-per-node=12\n")
            script.write("#SBATCH --time=8:00:00\n")
            script.write("#SBATCH --mem=30gb\n")
            script.write("#SBATCH --mail-user=your_email@uga.edu\n")
            script.write("#SBATCH --mail-type=BEGIN,END,FAIL\n")
            script.write("date\n\n")

            script.write("module load gffread\n")
            script.write("gffread reference/" + element + ".gff3 -T -o reference/" + element + "_gene.gtf\n")
            script.write("ml STAR\n")
            script.write("STAR \\ \n")
            script.write("--runThreadN 8 \\ \n")
            script.write("--runMode genomeGenerate \\ \n")
            script.write("--genomeSAindexNbases 13 \\ \n")
            script.write("--genomeDir . \\ \n")
            script.write("--genomeFastaFiles " + element + "_genome.fa \\ \n")
            script.write("--sjdbGTFfile reference/" + element + "_gene.gtf \n")

            # user will provide path 
            #script.write("cp " + table.[path] + " + )
            script.write("tar xvfz eugra.tar.gz\n")
            script.write("cp -r ./eugra/* fastq/\n")
            script.write("sbatch array_mapping_script.sh\n")


        # handle array size 
        with open("./data/" + element + "/map.sh", "wt") as mapping_script:
            """
            makes array mappign script
            """
            mapping_script.write(f"""#!/bin/sh
#SBATCH -J " + element + "_script
#SBATCH --partition batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=8:00:00
#SBATCH --mem=30gb
#SBATCH --mail-user=your_email@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL
date

cd ${{cleanfolder}}
ml Trimmomatic
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads $SLURM_NTASKS_PER_NODE -phred33 \\ 
${{readsfolder}}/${{sampleFile}}.R1.fastq ${{readsfolder}}/${{sampleFile}}.R2.fastq \\
${{sampleFile}}_trimP_1.fq.gz ${{sampleFile}}_trimS_1.fq.gz ${{sampleFile}}_trimP_2.fq.gz ${{sampleFile}}_trimS_2.fq.gz \\ 
ILLUMINACLIP:${{trimmo_path}}/adapters/TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:36 TOPHRED33 &>${{sampleFile}}_trim.log

### Mapping with STAR ###
cd ${{mapfolder}}

### STAR mapping time
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
ml Subread
featureCounts -Q 2 -s 0 -T $SLURM_NTASKS_PER_NODE -p -C \\
-a ${{genomefolder}}/gene.gtf \\
-o ${{sampleFile}}_counts.txt ${{mapfolder}}/${{sampleFile}}_Aligned.sortedByCoord.out.bam

## get trim log
python ~/script/get_trim_sum.py ${{masterfolder}}/clean/
##Please notice you need to get the get_trim_sum.py from NGSclean github
## Get map log
grep "" ${{masterfolder}}/map/*Log.final.out > ${{masterfolder}}/all_mapping_logs.txt
            
date
""")


make_directories()
extract_data()
generate_scripts()
