#!/bin/bash

module load smrtlink
module load samtools

PROJDIR=

sampledir=$PROJDIR/samples
resultdir=$PROJDIR/extra/pcr_mediated_recombination
scriptdir=$PROJDIR/extra/pcr_mediated_recombination/scripts

while read line
do
    if [[ ! "$line" =~ "SampleID" ]]
    then
	### extract sample info
	sampleid=`echo $line | cut -d, -f1`
	enzyme=`echo $line | cut -d, -f2`
	amplicon=`echo $line | cut -d, -f3`
	input=`echo $line | cut -d, -f4`
	yield=`echo $line | cut -d, -f5`

	### define markers
	if [ "$amplicon" == "DNA-1" ] ; then
	    markers="--cue=77,T,C --cue=153,T,A --cue=279,G,T --cue=416,C,A --cue=543,G,A --cue=672,C,G --cue=798,A,T --cue=926,G,C --cue=1003,A,C"
	    amplicon_length=1076
	fi

	if [ "$amplicon" == "DNA-2" ] ; then
    	    markers="--cue=59,A,C --cue=121,G,A --cue=221,C,A --cue=317,C,G --cue=423,A,T --cue=521,T,G --cue=625,T,A --cue=729,G,C --cue=821,G,A --cue=916,T,C --cue=1020,A,G"
	    amplicon_length=1078
	fi

	### path to reference sequence
	reference=$PROJDIR/references/$amplicon/sequence/$amplicon.fasta

	echo ""
	echo "Processing sample $enzyme ($amplicon)..."

    	### define rundir
    	rundir=`printf "%s/results/%05i" $resultdir $sampleid`
    	mkdir -p $rundir
    	cd $rundir
	
     	for strand in fwd rev
    	do
    	    ### strand-specific consensus reads
    	    ccsreads=`printf "%s/%05i/subreads_ccs.%s.bam" $sampledir $sampleid $strand`
	    
    	    ### extract ZMW stats
    	    $scriptdir/bam2csv.pl $ccsreads subreads_ccs.$strand.csv
	    
    	    ### map reads
    	    blasr $ccsreads $reference -bam -out aligned_reads.$strand.bam -bestn 1 -nproc 2 2>aligned_reads.$strand.log

    	    ### extract bases at marker sites
    	    $scriptdir/markers.pl --np 3 --match 2 $markers aligned_reads.$strand.bam $reference subreads_ccs.$strand.csv > aligned_reads.$strand.csv
    	done
	
	### count recombination events
    	$scriptdir/recombinations.pl $markers aligned_reads.{fwd,rev}.csv > summary.csv
	
	### build tables
	Rscript $scriptdir/calculte.R summary.csv $amplicon_length $input $yield
    fi
done < samples.csv
