#!/usr/bin/env python

import sys
import pysam

sample_name=sys.argv[1]
pe_sam=pysam.AlignmentFile(sample_name+'.paired.interval.bam','rb')
se_sam=pysam.AlignmentFile(sample_name+'.single.interval.bam','wb', template=pe_sam)
r1_tag=""

#assume file is query name sorted
for read in pe_sam.fetch(until_eof=True):
    if read.is_read1:
        if read.has_tag("ZI"):
            r1_tag=read.get_tag("ZI")
        continue
    elif read.is_read2:
        s_read=pysam.AlignedSegment()
        s_read.query_name=read.query_name
        s_read.flag=read.flag
        s_read.reference_id=read.reference_id
        s_read.reference_start=read.reference_start
        s_read.mapping_quality=read.mapping_quality
        s_read.cigartuples=read.cigartuples
        s_read.next_reference_id=read.next_reference_id
        s_read.next_reference_start=read.next_reference_start
        s_read.template_length=read.template_length
        s_read.query_sequence=read.query_sequence
        s_read.query_qualities=read.query_qualities
        s_read.tags=read.tags

        if r1_tag:
            if read.has_tag("ZI")==False:
                s_read.set_tag('ZI', r1_tag, value_type='Z')
            elif read.get_tag("ZI")!=r1_tag:
                #print(read.query_name+" is skipped because of ambiguous peaks overlapped.")
                continue

        se_sam.write(s_read)
        r1_tag=""    
    else:
        print(read.query_name+" is a unpaired read/n.")

se_sam.close()
pe_sam.close()
