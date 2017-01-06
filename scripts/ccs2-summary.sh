#!/bin/bash

#$ -S /bin/bash
#$ -P longrun
#$ -pe smp 2

cat $JOB_SCRIPT >> $SGE_STDOUT_PATH

echo ""
echo "Started on $(date) on $(uname -n)"

module load p7zip

reference="$root/references/$amplicon/sequence/$amplicon.fasta"

### run directory
rundir=`printf "%s/samples/%05i/summary" $root $sampleId`
mkdir -p $rundir
cd $rundir

### tally up mutations
echo ""
sampledir=`printf "%s/samples/%05i" $root $sampleId`
$root/bin/ccs2-summary.pl \
    --np 15 \
    --qv 93 \
    --lb 40 \
    --ub 40 \
    --rlp 0.80 \
    $reference $sampledir >summary.csv 2>summary.log
echo "Task 1 completed at $(date)"

echo ""
echo "Finished on $(date)"
