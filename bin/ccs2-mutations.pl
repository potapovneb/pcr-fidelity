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
    print "usage: $0 input.bam reference.fasta\n\n";
    exit;
}

my $bamfile = shift @ARGV;
my $reffile = shift @ARGV;

### load reference
my $refs = load_references($reffile);

### compute homopolymer regions
my $hpmap = build_hpmap($refs);

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

print "Movie,ZMW,Ref,Pos,Length,Type,QV,AtEnd,RefBP,AltBP,IndelType,InHP,Bases,HPLen,HPChar\n";

open(BAM,"samtools view $bamfile |") || die "Can't open '$bamfile'";

while( my $line = <BAM> )
{
    my $s1 = ""; # aligned reference
    my $s2 = ""; # aligned read
    my $bq = ""; # base quality

    ### parse fields
    my ( $qname,
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
	 @tag ) = split(/\t/,$line);

    ### fwd/rev strand
    my $strand = ( $flag & 0x10 ) ? 1 : 0;

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
    
    my $ref_start = $pos;
    my $ccs_start = 1;
    
    my ($movie,$zmw,$other) = split(/\//,$qname);
    
    mutations($rname,$movie,$zmw,$s1,$s2,$bq,$ref_start-1,$ccs_start-1,$strand);
}

################################################################################
#                             SUBROUTINES                                      #
################################################################################

sub mutations {
    my ($ref_name,$movie,$zmw,$rseq,$read,$qual,$ref_start,$ccs_start,$strand) = @_;
    
    my $len1 = length($rseq);
    my $len2 = length($read);
    
    if( $len1 != $len2 )
    {
	print "error :: sequence length does not match\n\n\n";
	exit;
    }
    
    ### convert strings to arrays
    my @rseq = split(//,uc($rseq));
    my @read = split(//,uc($read));
    
    ### convert quality values
    my @qual = map { ord($_) - 33 } split(//,$qual);
    
    ### init position in the reference sequence
    my $rp = $ref_start - 1;
    
    ### init position in the sequencing read
    my $sp = $ccs_start - 1;
    
    ### init indel position in the reference sequence
    my $pindel_ref = undef;
    
    ### init auxilary indel data
    my @bases = ();
    my @positions = ();
    
    for( my $i = 0; $i < $len1; $i++ )
    {
	my $rb = $rseq[$i];
	my $sb = $read[$i];
	
	### update reference position
	$rp++ if( $rb ne "-" );
	
	### update sequencing read position
	$sp++ if( $sb ne "-" );
	
	if( $rb ne "-" && $sb ne "-" )
	{
	    ### SNP
	    
	    ### reset indel data
	    @bases = ();
	    @positions = ();
	    $pindel_ref = undef;
	    
	    my $offset = 0;
	    
	    ### update position in homopolymer regions
	    if( $strand == 1 )
	    {
 		my $k = 1;
		
		### offset by number of gaps in homopolymer
		while( ($i+$k) < $len1 && $rseq[$i+$k] eq $rseq[$i] )
		{
		    $offset++ if( $read[$i+$k] eq "-" );
		    $k++;
		}
	    }
	    
	    ### print call info
	    printf( "%s,%s,%s,%i,%i,%s,%i,%s,%s,%s,%s,%s,%s,%s,%s\n",
		    $movie,
		    $zmw,
		    $ref_name,
		    $rp + $offset + 1,
		    1,
		    "SNP",
		    $qual[$sp],
		    ( $i == 0 || $i + 1 == $len1 ? "True" : "False" ),
		    $rb,
		    $sb,
		    "NA",
		    "NA",
		    "NA",
		    "NA",
		    "NA" );
	}
	else
	{
	    ### INDEL
	    if( $rb ne "-" && $sb eq "-" )
	    {
		### DELETION
		
		### indel position
		$pindel_ref = $rp-1 if( ! defined $pindel_ref );
		
		### indel bases
		push @bases, $rb;
		
		### indel positions
		push @positions, $sp;
		
		### [end of alignment] or [end of deletion]
		if( $i + 1 == $len1 || $read[$i+1] ne "-" )
		{
		    my $qv = undef;
		    
		    if( $strand == 0 )
		    {
			### first base after deletion (fwd read)
			$qv = $qual[$positions[0]+1];
		    }
		    else
		    {
			### first base before deletion (rev read)
			$qv = $qual[$positions[0]];
		    }
		    
		    ### update indel position in homopolymer regions
		    if( $strand == 1 )
		    {
			my $k = 1;
			
			### offset $pindel by number of gaps in homopolymer
			while( ($i-$k) >= 0 && $rseq[$i-$k] eq $rseq[$i] )
			{
			    $pindel_ref-- if( $read[$i-$k] ne "-" );
			    $k++;
			}
		    }

		    ### homopolymer info (not currently used in downstream analysis)
		    my $HPLen = $$hpmap{$ref_name}[$pindel_ref+1];
		    my $InHP = ( $HPLen > 1 ) ? "True" : "False";
		    my $HPChar = substr($$refs{$ref_name},$pindel_ref+1,1);
		    
		    ### print call info
		    printf( "%s,%s,%s,%i,%i,%s,%i,%s,%s,%s,%s,%s,%s,%s,%s\n",
			    $movie,
			    $zmw,
			    $ref_name,
			    $pindel_ref + 1,
			    scalar(@bases),
			    "INDEL",
			    $qv,
			    "False",
			    "NA",
			    "NA",
			    "Deletion",
			    $InHP,
			    join("",@bases),
			    $HPLen,
			    $HPChar );
		}
	    }
	    elsif( $rb eq "-" && $sb ne "-" )
	    {
		### INSERTION
		
		### indel position
		$pindel_ref = $rp if( ! defined $pindel_ref );
		
		### indel bases
		push @bases, $sb;
		
		### indel positions
		push @positions, $sp;
		
		### [end of alignment] or [end of insertion]
		if( $i + 1 == $len1 || $rseq[$i+1] ne "-" )
		{
		    ### check for a homopolymer insertion
		    my $hp = 1;
		    
		    for( my $k = 0; $k < @bases-1; $k++ )
		    {
			if( $bases[$k] ne $bases[$k+1] )
			{
			    $hp = 0;
			    last;
			}
		    }

		    ### assign QV
		    my $qv = undef;

		    if( $hp )
		    {
			if( $strand == 0 )
			{
			    ### first base in a homopolymer (fwd read)
			    $qv = $qual[$positions[0]];
			}
			else
			{
			    ### last base in a homopolymer (rev read)
			    $qv = $qual[$positions[-1]];
			}
		    }
		    else
		    {
			### otherwise, take the lowest QV value
			$qv = 93;
			
			for( my $k = 0; $k < @bases; $k++ )
			{
			    if( $qv > $qual[$positions[$k]] )
			    {
				$qv = $qual[$positions[$k]];
			    }
			}
		    }
		    
		    ### update indel position in homopolymer regions
		    if( $strand == 1 )
		    {
			my $k = 1;
			
			### offset $pindel by number of gaps in homopolymer
			while( ($i-$k) >= 0 && $read[$i-$k] eq $read[$i] )
			{
			    $pindel_ref-- if( $rseq[$i-$k] ne "-" );
			    $k++;
			}
		    }

		    ### homopolymer info (not currently used in downstream analysis)
		    my $HPLen = $$hpmap{$ref_name}[$pindel_ref+1];
		    my $InHP = ( $HPLen > 1 ) ? "True" : "False";
		    my $HPChar = substr($$refs{$ref_name},$pindel_ref+1,1);
		    
		    ### print call info
		    printf( "%s,%s,%s,%i,%i,%s,%i,%s,%s,%s,%s,%s,%s,%s,%s\n",
			    $movie,
			    $zmw,
			    $ref_name,
			    $pindel_ref + 1,
			    scalar(@bases),
			    "INDEL",
			    $qv,
			    "False",
			    "NA",
			    "NA",
			    "Insertion",
			    $InHP,
			    join("",@bases),
			    $HPLen,
			    $HPChar );
		}
	    }
	    else
	    {
		print STDERR "error :: rb='$rb' sb='$sb'\n\n";
		exit;
	    }
	}
    }
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

sub build_hpmap {
    my ($refs) = @_;
    
    my %hpmap = ();

    for my $name ( keys %$refs )
    {
	my $seq = $$refs{$name};
	
	for( my $i = 0; $i < length($seq); $i++ )
	{
	    my $len = 1;
	    my $k = $i + 1;
	    
	    while( $k < length($seq) && substr($seq,$i,1) eq substr($seq,$k,1) )
	    {
		$len++;
		$k++;
	    }

	    push @{$hpmap{$name}}, $len;
	}
    }

    return \%hpmap;
}
