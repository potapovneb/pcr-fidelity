#!/bin/bash

reference=$1
start=$2
end=$3

w=30

echo ">correct"
extract_seq.pl $reference $((start-w)) $((end+w))
echo ""

echo ">double_switch"
extract_seq.pl $reference $((start-w)) $((start-1))
extract_seq.pl --reverse --complement $reference $start $end
extract_seq.pl $reference $((end+1)) $((end+w))
echo ""

echo ">single_switch_forward"
extract_seq.pl $reference $((start-w)) $((start-1))
extract_seq.pl --reverse --complement $reference $start $end
extract_seq.pl --reverse --complement $reference $((start-w)) $((start-1))
echo ""

echo ">single_switch_reverse"
extract_seq.pl --reverse --complement $reference $((end+1)) $((end+w))
extract_seq.pl $reference $start $end
extract_seq.pl $reference $((end+1)) $((end+w))
echo ""
