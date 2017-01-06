#!/usr/bin/perl -w

################################################################################
# PCR Fidelity: Examining Sources of Error in PCR by Single-Molecule Sequencing
# Copyright (C) 2016 New England Biolabs, Inc.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

use strict;

### command-line arguments
if( @ARGV == 0 )
{
    print "usage: $0 input.bam reference.fasta zmws.csv\n\n";
    exit;
}

my $bamfile = shift @ARGV;
my $reffile = shift @ARGV;
my $zmwfile = shift @ARGV;

### load reference
my $refs = load_references($reffile);

### load ZMW stats
my $zmws = load_zmw_data($zmwfile);

# data structures for parsing CIGAR string
my %digit = (
    "0" => 1,
    "1" => 1,
    "2" => 1,
    "3" => 1,
    "4" => 1,
    "5" => 1,
    "6" => 1,
    "7" => 1,
    "8" => 1,
    "9" => 1 );

my %operator = (
    "M" => 1,
    "I" => 1,
    "D" => 1,
    "=" => 1,
    "X" => 1 );

print "Movie,ZMW,Reference,QueryCoverage,RefLength,Mismatches,NP,ReadLength\n";

my %data = ();

# go through BAM file
open(BAM,"samtools view $bamfile |") || die "Can't open '$bamfile'";

while( my $line = <BAM> )
{
    my $s1 = ""; # aligned reference
    my $s2 = ""; # aligned read
    my $bq = ""; # base quality

    ### parse fields
    my ($qname,
	$flag,
	$rname,
	$pos,
	$mapq,
	$cigar,
	$rnext,
	$pnext,
	$tlen,
	$read,
	$qual,
	@tag) = split(/\t/,$line);

    ### skip unmapped reads
    next if( $flag & 0x4 );

    ### skip secondary alignment
    next if( $flag & 0x100 );

    ### skip supplementary alignment
    next if( $flag & 0x800 );

    ### skip reads with low mapping quality
    next if( $mapq < 254 );

    # ----- process CIGAR ------------------------------------------------------

    my $rseq = $$refs{$rname};

    my $p1  = $pos - 1; # start of aligned reference
    my $p2  = 0;        # start of aligned read
    my $num = 0;        # length of the block to operate
    my $len = 0;        # read length
    
    my @cigar = split(//,$cigar);
    
    while( @cigar )
    {
	my $c = shift @cigar;
	
	if( exists $digit{$c} )
	{
	    $num .= $c;
	    next;
	}
	
	if( $c eq "M" || $c eq "=" || $c eq "X" )
	{
	    $s1 .= substr($rseq,$p1,$num);
	    $s2 .= substr($read,$p2,$num);
	    $bq .= substr($qual,$p2,$num);

	    $p1 += $num;
	    $p2 += $num;
	    
	    $len += $num;
	}
	elsif( $c eq "D" )
	{
	    $s1 .= substr($rseq,$p1,$num);
	    $s2 .= ("-" x $num);
	    $bq .= (" " x $num);
	    
	    $p1 += $num;
	    
	    $len += $num;
	}
	elsif( $c eq "I" )
	{
	    $s1 .= ("-" x $num);
	    $s2 .= substr($read,$p2,$num);
	    $bq .= substr($qual,$p2,$num);

	    $p2 += $num;
	}
	elsif( $c eq "S" )
	{
	    $p2 += $num;
	}
	else
	{
	    print STDERR "error: unknown operator in CIGAR string '$c'\n\n";
	    exit;
	}
	
	$num = "";
    }
    
    ### compute alignment stats
    my $query_coverage = 0;
    my $mismatches = 0;

    for( my $i = 0; $i < length($s1); $i++ )
    {
	if( substr($s1,$i,1) ne "-" )
	{
	    $query_coverage++;

	    if( substr($s1,$i,1) ne substr($s2,$i,1) )
	    {
		$mismatches++;
	    }
	}
    }

    my $ref_length = length($$refs{$rname});

    my ($movie,$zmw,@oth) = split(/\//,$qname);

    ### print alignment & ZMW stats
    print join( ",",
    		$movie,
    		$zmw,
    		$rname,
    		$query_coverage,
    		$ref_length,
    		$mismatches,
    		$$zmws{$movie}{$zmw}{"NP"},
    		$$zmws{$movie}{$zmw}{"ReadLength"} ), "\n";
}

sub load_references {
    my ($file) = @_;
    
    my %refs = ();
    my $name = "";
    
    open(FA,$file) || die "Can't open '$file': $!";
    
    while( my $line = <FA> )
    {
	chomp($line);
	
	if( substr($line,0,1) eq ">" )
	{
	    $name = substr($line,1);
	}
	else
	{
	    $refs{$name} .= $line;
	}
    }
    
    close(FA);
    
    return \%refs;
}

sub load_zmw_data {
    my ($file) = @_;
    
    my %zmws = ();
    my @head = ();

    open(IN,$file) || die "Can't open '$file': $!";
    
    while( my $line = <IN> )
    {
	chomp($line);
	
	### parse line
	my @tokens = split(/,/,$line);
	
	if( @head == 0 )
	{
	    ### extract column names
	    @head = @tokens;
	}
	else
	{
	    ### store data fields
	    my %entry = ();
	    
	    for( my $i = 0; $i < @head; $i++ )
	    {
		$entry{$head[$i]} = $tokens[$i];
	    }

	    $zmws{$entry{"Movie"}}{$entry{"ZMW"}} = \%entry;
	}
    }

    close(IN);

    return \%zmws;
}
