cp /work/cjtlab/Database/Egrandis/v2.0/annotation/Egrandis_297_v2.0.gene_exons.gff3 .
cp /work/cjtlab/Database/Egrandis/v2.0/assembly/Egrandis_297_v2.0.fa ./genome.fa
module load gffread
gffread Egrandis_297_v2.0.gene_exons.gff3 -T -o gene.gtf 
STAR \
--runThreadN 8 \
--runMode genomeGenerate \
--genomeSAindexNbases 13 \
--genomeDir . \
--genomeFastaFiles genome.fa \
--sjdbGTFfile gene.gtf
cp /work/cjtlab/testing_data/eugra.tar.gz .
tar xvfz eugra.tar.gz
cp -r ./eugra/* .