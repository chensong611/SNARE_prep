# SNARE_prep
dropseq-based **S**ingle **N**uclei genomic **A**ccessibility and **R**NA **E**xpression-seq (SNARE-seq) preprocessing pipeline

## Dependency
[Drop-seq_tools-1.13/Picard](https://github.com/broadinstitute/Drop-seq/releases)

[ATAC-seq pipeline](https://github.com/kundajelab/atac_dnase_pipelines#pipeline)

[Reference(mm10,hg38)](https://github.com/kundajelab/atac_dnase_pipelines#genome-data)

[.dict file](https://broadinstitute.github.io/picard/command-line-overview.html#CreateSequenceDictionary)

[pysam](https://github.com/pysam-developers/pysam)

[fastx_trimmer](http://hannonlab.cshl.edu/fastx_toolkit/index.html)

## Read(fastq) Structure 
Read1(30 bp): [12-bp barcode][9-bp UMI]TTTTTTTTT

Read2(75 bp): Nextera Read1

Read3(75 bp): Nextera Read2

##  Usage
`SNARE_Prep.sh -n $name -s $species -c $cells`
