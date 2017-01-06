#!/bin/bash

#$ -S /bin/bash
#$ -P longrun
#$ -pe smp 8

cat $JOB_SCRIPT >> $SGE_STDOUT_PATH

echo ""
echo "Started on $(date) on $(uname -n)"

module load bwa
module load p7zip
module load samtools

reference="$root/references/$amplicon/sequence/$amplicon.fasta"

### top run directory
rundir=`printf "%s/samples/%05i/chimeric" $root $sampleId`

### strand-specific analysis
echo ""
for strand in fwd rev
do
    ### strand-specific run directory
    echo ""
    mkdir -p $rundir/$strand
    cd $rundir/$strand
    echo "Task 1 ($strand) completed at $(date)"
    
    ### copy and index reference fasta
    echo ""
    cp "$reference" reference.fasta
    bwa index reference.fasta
    echo "Task 2 ($strand) completed at $(date)"
    
    ### convert CCS reads to FASTQ format
    echo ""
    ccsreads=`printf "%s/samples/%05i/subreads_ccs.$strand.bam" $root $sampleId`
    $root/bin/bam2fq.pl $ccsreads subreads_ccs.$strand.fastq
    echo "Task 3 ($strand) completed at $(date)"

    ### map reads
    echo ""
    bwa mem -t 8 reference.fasta subreads_ccs.$strand.fastq > aligned_ccs.$strand.sam
    echo "Task 4 ($strand) completed at $(date)"
    
    ### extract and save names of chimeric reads according to BWA
    echo ""
    samtools view -f 2048 -q60 aligned_ccs.$strand.sam | cut -f1 | cut -d'/' -f1,2 --output-delimiter=, | sort -u | sort -t, -k2 -n > chimeric.csv
    echo "Task 5 ($strand) completed at $(date)"

    ### compress output file
    echo ""
    7za a -t7z -m0=lzma -mx=9 -mfb=64 -md=32m -ms=on chimeric.csv.7z chimeric.csv
    echo "Task 6 ($strand) completed at $(date)"
    
    ### cleanup
    echo ""
    rm -f reference.*
    rm -f aligned_ccs.*
    rm -f subreads_ccs.*
    rm -f chimeric.csv
    echo "Task 7 ($strand) completed at $(date)"
done

echo ""
echo "Finished on $(date)"
