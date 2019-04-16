#!/bin/bash

dropseq_dir=
picard_dir=
atac_dir=
ref_dir=

species='m'
cells=2000

#barcode end base position
bEnd=12

#UMI start and end base positions
uStart=$(($bEnd+1))
uEnd=21


while getopts ":n:s:c:" options; do
    case $options in
    	n ) name=$OPTARG;;
        s ) species=$OPTARG;;
	    c ) cells=$OPTARG;;
    esac
done
shift $(($OPTIND - 1))

if [ $species = 'h' ]; then
    ref_fasta=hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
    ref_dict=hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.dict
    species='hg38'
elif [ $species = 'm' ]; then
    ref_fasta=mm10/mm10_no_alt_analysis_set_ENCODE.fasta
    ref_dict=mm10/mm10_no_alt_analysis_set_ENCODE.dict
    species='mm10'
else
	echo "No reference found!"
fi

mkdir -p ATAC/$name
mkdir -p Reports
mkdir -p Tmp
mkdir -p DCA

#merge read1 and read2 together as a single read 
paste -d "" <(zcat $name'_R1_001.fastq.gz') <(zcat $name'_R2_001.fastq.gz' | sed '1~2s/.*//') > $name'_R12_001.fastq'

#make a pair-ended bam file, extract cell barcode and UMI information and store it in bam tags (XC/XM), and convert it back to pair-end fastq files
java -jar $picard_dir/picard.jar FastqToSam FASTQ=$name'_R12_001.fastq' FASTQ2=$name'_R3_001.fastq.gz' SAMPLE_NAME=$name OUTPUT=/dev/stdout | \
$dropseq_dir/TagBamWithReadSequenceExtended I=/dev/stdin O=/dev/stdout SUMMARY=Reports/$name'.cell_tag_report.txt' BASE_RANGE=1-$bEnd BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=false TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1 | \
$dropseq_dir/TagBamWithReadSequenceExtended I=/dev/stdin O=/dev/stdout SUMMARY=Reports/$name'.molecule_tag_report.txt' BASE_RANGE=$uStart-$uEnd BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=false TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1 | \
tee $name'.unaligned.tagged.bam' | \
java -jar $picard_dir/picard.jar SamToFastq I=/dev/stdin FASTQ=$name'_R12_002.fastq' SECOND_END_FASTQ=$name'_R3_001.fastq'

#remove read1 from merged single read
fastx_trimmer -f 31 -i $name'_R12_002.fastq' -o $name'_R2_001.fastq' -Q33

#call peaks from aggregated pair-end fastq files
bds $atac_dir/atac.bds -species $species -nth 16 -fastq1_1 $name'_R2_001.fastq' -fastq1_2 $name'_R3_001.fastq' -out_dir ./ATAC/$name -adapter1_1 CTGTCTCTTATACACATCT -adapter1_2 CTGTCTCTTATACACATCT -no_ataqc

#convert narrowPeak file to interval_list file
gunzip ATAC/$name/peak/macs2/overlap/conservative_set/$name'_rep1-pr.naive_overlap.filt.narrowPeak.gz'
java -jar $picard_dir/picard.jar BedToIntervalList I=ATAC/$name/peak/macs2/overlap/conservative_set/$name'_rep1-pr.naive_overlap.filt.narrowPeak' O=$name'.interval_list' SD=$ref_dir$ref_dict

#tag reads with peaks called, sort in query name order and retrive barcode/UMI information.
$dropseq_dir/TagReadWithInterval I=ATAC/$name/align/rep1/$name'_R2_001.trim.PE2SE.bam' O=/dev/stdout LOCI=$name'.interval_list' TAG=ZI | \
java -jar $picard_dir/picard.jar SortSam I=/dev/stdin O=$name'.sorted.interval.bam' SO=queryname TMP_DIR=Tmp
java -jar $picard_dir/picard.jar MergeBamAlignment REFERENCE_SEQUENCE=$ref_dir$ref_fasta UNMAPPED_BAM=$name'.unaligned.tagged.bam' ALIGNED_BAM=$name'.sorted.interval.bam' O=$name'.tagged.interval.bam' INCLUDE_SECONDARY_ALIGNMENTS=false PAIRED_RUN=true ATTRIBUTES_TO_RETAIN=ZI TMP_DIR=Tmp

#sort pair-end bam file in query name and convert it to single-end bam file labeled with peaks
java -jar $picard_dir/picard.jar SortSam I=$name'.tagged.interval.bam' O=$name'.paired.interval.bam' SO=queryname TMP_DIR=Tmp
SNARE_PairToSingle.py $name

#create digital count matrix of chromatin accessibility
$dropseq_dir/DigitalExpression I=$name'.single.interval.bam' O=DCA/$name'.counts.tsv' SUMMARY=DCA/$name'.count_summary.txt' NUM_CORE_BARCODES=$cells GENE_EXON_TAG=ZI EDIT_DISTANCE=1 USE_STRAND_INFO=false


