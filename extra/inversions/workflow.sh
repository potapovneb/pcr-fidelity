#!/bin/bash

### print usage instructions
if [ "$#" -ne 1 ] || ( [ "$1" != "execute" ] && [ "$1" != "extract" ] )
then
    echo "usage: $0 execute" >&2
    echo "       $0 extract" >&2
    exit 1
fi

### define top project directory
PROJDIR=

### define subdirectories
sampledir=$PROJDIR/samples
scriptdir=$PROJDIR/extra/inversions/scripts
resultdir=$PROJDIR/extra/inversions/results

### analyse lacZ inversions
cd $PROJDIR/extra/inversions

mkdir -p $resultdir

### generate reference sequences
export PATH=$scriptdir:$PATH
reference=$PROJDIR/references/LacZ-2/sequence/LacZ-2.fasta

reference.sh $reference   45  183 > $resultdir/reference.0045_0183.fasta
reference.sh $reference 1285 1305 > $resultdir/reference.1285_1305.fasta

for region in 0045_0183 1285_1305
do
    echo ""
    echo $region

    if [ "$1" == "extract" ] ; then
	echo "Enzyme,Correct,Double template switching,Single template switching (forward),Single template switching (reverse)"
    fi

    while read line
    do
	if [[ ! $line =~ "SampleID" ]]
	then
	    ### extract sampleId
	    sampleId=`echo $line | cut -d, -f1`
	    enzyme=`echo $line | cut -d, -f2`
	    
	    ### format
	    sample=`printf %05i $sampleId`
	    
	    ### define output directory
	    rundir=$resultdir/$region/$sample
	    
	    ### set of references to use
	    reference=$resultdir/reference.$region.fasta
	    
	    if [ "$1" == "execute" ] ; then
		### map reads
		qsub -v sample="$sample",reference="$reference",sampledir="$sampledir",scriptdir="$scriptdir",rundir="$rundir" $scriptdir/mapreads.sh
	    else
		### extract stats
		n1=`Rscript $scriptdir/summary.R $rundir/summary.merged.csv correct`
		n2=`Rscript $scriptdir/summary.R $rundir/summary.merged.csv double_switch`
		n3=`Rscript $scriptdir/summary.R $rundir/summary.merged.csv single_switch_forward`
		n4=`Rscript $scriptdir/summary.R $rundir/summary.merged.csv single_switch_reverse`
		echo $enzyme,$n1,$n2,$n3,$n4
	    fi
	fi
    done < samples.csv
done
