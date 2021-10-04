# RNA-Seq analysis on sapelo2 at UGA

Tsai lab RNA-Seq script

## Description

RNA-Seq runs transcriptomic analysis. The input of this pipeline is a csv file. So to begin, prepare a table with the schema like `sample_table.csv`. This script will then run the RNA-Seq pipeline on the samples from different species in the table.

Here's how the pipeline flows:

![workflow](images/workflow.png)

It will create folder structure like:

```
RNA-Seq/data
|
 ----- {species_id}
        |
         ----- fastq (raw fastq files)
         ----- reference (genome reference from Phytozome or other sources)
               |
                ----- species_genome.fasta -> genome.fa
                ----- species_annotation.gff3 -> gene.gtf
        ----- clean (filtered reads)
        ----- map (alignment maps)
        ----- count (quantification results for further analysis)
        ----- miscelaneous
              |
               ----- *.txt
               ----- log....
        ----- file.list
        ----- start.sh
        ----- map.sh
        ----- get_trim_sum.py
```

## Instructions to Run
1. clone this repo
```bash
git clone https://github.com/TsailabBioinformatics/RNA-Seq.git;cp RNA-Seq/* .
```
- if cloning returned `fatal: Authentication failed`, then try again using a github personal access token. instructions to do so can be found here: <https://docs.github.com/en/enterprise-server@3.1/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token>
2. create an input csv file (like `sample_input`), name it `input.csv`, and add it to `RNA-Seq/input/`
- check to make sure `input.csv` exists within the folder `RNA-Seq/input`
3. change directories into `RNA-Seq`
5. load modules `python-utils` and `pandas` with commands:
```bash
ml python-utils
```
```bash
ml pandas
```
7. run the pipeline
```python3
python3 pipeline.py
```
