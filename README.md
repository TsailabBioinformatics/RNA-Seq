# RNA-Seq

Tsai lab RNA-Seq script

## demo using Eucalyptus data: instructions on sapelo2

we are gonna use interaction mode to process some data, then submit the main script to the cluster.

```bash
qlogin
```

1. clone this repo

```bash
git clone https://github.com/TsailabBioinformatics/RNA-Seq.git;cp RNA-Seq/* .
```

2. copy the genome fasta and annotatino gff3 file to the data folder

```bash
cp /work/cjtlab/Database/Egrandis/v2.0/annotation/Egrandis_297_v2.0.gene_exons.gff3 .
cp /work/cjtlab/Database/Egrandis/v2.0/assembly/Egrandis_297_v2.0.fa ./genome.fa

```

3. transfer the gff3 to `genome.gtf` using interactive mode

```bash
module load gffread
gffread Egrandis_297_v2.0.gene_exons.gff3 -T -o gene.gtf 
```

4. copy the fastq file to the data folder

```bash
cp /work/cjtlab/testing_data/eugra.tar.gz .
tar xvfz eugra.tar.gz
```

5. change the email address in the `main_script.sh`

6. submit the job

```bash
sbatch main_script.sh
```

# make it work for your own data (WIP)

1. change the file.list (see instruciton from [Ran's note](https://www.evernote.com/shard/s202/client/snv?noteGuid=070f6281-ef94-47c1-a4df-3dbb2083693c&noteKey=2e87d16e54db6d4b&sn=https%3A%2F%2Fwww.evernote.com%2Fshard%2Fs202%2Fsh%2F070f6281-ef94-47c1-a4df-3dbb2083693c%2F2e87d16e54db6d4b&title=RNAseq%2Bpipeline%2B%2528SLURM%2Bsystem%2B2020%2529))

2. copy your own fastq file to the data folder

3. copy your own genome fasta and genome.gtf to the data folder

4. change the array number: `#SBATCH --array=1-n`, `n` = number of lines in `file.list`

# TODO (WIP)

1. scale up for different species

2. make the script not restricted to `genome.fa` and `genome.gtf` this specific file name

3. turn the overall script into real pipeline using `Common Workflow Language (CWL)` or `Workflow Description Language (WDL)`

4. package stuff into container

5. automatic adjust the number of array requested

6. get rid of the unnecessary zipping and intermediate files removing commands
