#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("USAGE: fpkm_quick.R [COUNT] [SAMPLELIST] [OUTPUT]\nAt least one argument must be supplied (The count file)\n", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  args[3] = "out.txt"
}
library("DESeq2")
cts=read.table(args[1],T,row.names =1)
samples=read.table(args[2],F)
colnames(samples)=c("SAMPLES")
countData = as.data.frame(cts)
library("GenomicFeatures")
length <- nchar(getwd())
path <- paste("./reference", substring(getwd(), length - 5, length), "_gene.gtf", sep="")
txdb <- makeTxDbFromGFF(path, format = "gtf", circ_seqs = character())
ebg <- exonsBy(txdb, by="gene")

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = samples,
                              design = ~ SAMPLES)
rowRanges(dds) = ebg   
fpkm_out = as.data.frame(fpkm(dds))
nc=ncol(fpkm_out)
fpkm_out[,1:nc]=round(fpkm_out[,1:nc],3) #set digits to 3
write.table(fpkm_out,args[3],sep="\t",quote = FALSE)
sessionInfo() #print out the session information
