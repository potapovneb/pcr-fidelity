#!/bin/bash

#$ -S /bin/bash
#$ -N inversion
#$ -o /dev/null
#$ -j yes

module load smrtlink
module load samtools

### run directory
mkdir -p $rundir
cd $rundir

### strand-specific analysis
for strand in fwd rev
do
    ### locate reads
    ccsreads="$sampledir/$sample/subreads_ccs.$strand.bam"
    
    ### align reads to reference
    blasr $ccsreads $reference -bestn 1 -bam -out aligned_reads.$strand.bam
    
    ### extract ZMW stats
    $scriptdir/bam2csv.pl aligned_reads.$strand.bam aligned_reads.$strand.csv
    
    ### combine info
    $scriptdir/extract_aln.pl aligned_reads.$strand.bam $reference aligned_reads.$strand.csv > summary.$strand.csv
done

### merge data
head --lines 1 summary.fwd.csv > summary.merged.csv
cat summary.fwd.csv summary.rev.csv | egrep -v "^Movie" >> summary.merged.csv
