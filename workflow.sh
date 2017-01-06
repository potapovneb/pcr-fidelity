#!/bin/bash

if [ "$#" -ne 2 ] || ( [ "$2" != "view" ] && [ "$2" != "consensus" ] && [ "$2" != "chimeric" ] && [ "$2" != "mutation" ] && [ "$2" != "summary" ] && [ "$2" != "tabulate" ] )
then
    echo "usage: $0 samples.csv view"      >&2
    echo "       $0 samples.csv consensus" >&2
    echo "       $0 samples.csv chimeric"  >&2
    echo "       $0 samples.csv mutation"  >&2
    echo "       $0 samples.csv summary"   >&2
    echo "       $0 samples.csv tabulate"  >&2
    exit 1
fi

samples=$1
command=$2

### define top project directory
PROJDIR=

cd $PROJDIR

if [ "$command" == "tabulate" ] ; then
    echo "SampleID,Enzyme,Amplicon,Input,Yield,collectionPathUri,GroupID,AA,AC,AT,AG,CA,CC,CT,CG,TA,TC,TT,TG,GA,GC,GT,GG,Deletion,Insertion"
fi

while read line
do
    if [[ ! $line =~ "SampleID" ]]
    then
	### extract SampleID, Amplicon, and data path for each sample
	sampleId=`echo $line | cut -d, -f1`
	amplicon=`echo $line | cut -d, -f3`
	collectionPathUri=`echo $line | cut -d, -f6`
	
	if [ "$command" == "view" ]
	then
	    ### preview run info
	    echo ""
	    echo "sampleId=$sampleId"
	    echo "amplicon=$amplicon"
	    echo "collectionPathUri=$collectionPathUri"
	else
	    if [ "$command" != "tabulate" ] ; then
		logdir=`printf "%s/samples/%05i" $PROJDIR $sampleId`
		mkdir -p "$logdir"
    		qsub -v root="$PROJDIR",amplicon="$amplicon",sampleId="$sampleId",collectionPathUri="$PROJDIR/$collectionPathUri" -N "ccs$sampleId" -o $logdir/$command.log -j yes $PROJDIR/scripts/ccs2-$command.sh
	    else
		rundir=`printf "%s/samples/%05i/summary" $PROJDIR $sampleId`
		echo $line | tr '\n' ","
		tail --lines 1 $rundir/summary.csv
	    fi
	fi
    fi
done < $samples
