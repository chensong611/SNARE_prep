#!/usr/bin/env python

import os
import sys
import pysam

#open pair-end bam file, extract records by matching cellular barcode tags with provided barcode_list and store into sub_dir
tag="XC"
bam=pysam.AlignmentFile(sys.argv[1],'rb')
barcode_list=open(sys.argv[2],'r')
sub_dir=sys.argv[3]
os.mkdir(sub_dir)

d={}
for barcode in barcode_list:
    barcode=barcode.rstrip()
    f=pysam.AlignmentFile(sub_dir+'/'+barcode+'.bam','wb', template=bam)
    d[barcode]=f

for read in bam.fetch(until_eof=True):
    if read.is_read1:
        read1=read
    elif read.is_read2:
        f=d.get(read.get_tag(tag))
        if f!=None:
            read2=read
            f.write(read1)
            f.write(read2)
    else:
        print(read.query_name+" is a unpaired read/n.")

barcode_list.seek(0)
for barcode in barcode_list:
    barcode=barcode.rstrip()
    d[barcode].close()

barcode_list.close()
bam.close()
