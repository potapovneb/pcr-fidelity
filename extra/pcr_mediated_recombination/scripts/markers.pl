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
use Getopt::Long;

### command-line options
my $o_np = 0;
my $o_match = 0;
my @o_cues = ();
my %o_cues = ();
my $o_mapq = 254;

GetOptions(
    "np=f"    => \$o_np,
    "match=f" => \$o_match,
    "cue=s"   => \@o_cues,
    "mapq=f"  => \$o_mapq,
    );

for my $x ( split(/;/,join(";",@o_cues)) )
{
    my ($pos,$b0,$b1) = split(/,/,$x);

    $o_cues{$pos}{$b0} = 0;
    $o_cues{$pos}{$b1} = 1;
}

### command-line usage
if( @ARGV == 0 )
{
    print "usage: $0 ccsfile.bam reference.fasta\n";
    print "\n";
    print "options:\n";
    print "  --np\t\tminimum number of passes ($o_np)\n";
    print "  --match\tnumber of exact matches around marker ($o_match)\n";
    print "  --cues\tmarker sites (eg. --cue 77,T,C --cue 153,T,A --cue 279,G,T ...)\n";
    print "  --mapq\tminimum read mapping quality ($o_mapq)\n";
    print "\n";
    exit;
}

### command-line arguments
my $bamfile = shift @ARGV;
my $reffile = shift @ARGV;
my $zmwfile = shift @ARGV;

my $refs = load_references($reffile);
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

### CSV column names
printf( "%s,%s,%s,%s\n",
	"Movie",
	"ZMW",
	"NP",
	join(",",sort { $a <=> $b } keys %o_cues) );

open(BAM,"samtools view $bamfile |") || die "Can't open '$bamfile': $!\n";

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
	 @tag) = split(/\t/,$line);
    
    ### fwd/rev strand
    my $strand = ( $flag & 0x10 ) ? 1 : 0;
    
    ### skip unmapped reads
    next if( $flag & 0x4 );
    
    ### skip secondary alignment
    next if( $flag & 0x100 );
    
    ### skip supplementary alignment
    next if( $flag & 0x800 );
    
    ### skip reads with low mapping quality
    next if( $mapq < $o_mapq );

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
    
    ### ----- work with alignment ----------------------------------------------

    ### init base identities at marker sites
    my %data = map { $_ => "NA" } keys %o_cues;
    
    my $ref_pos = $pos - 1;

    ### extract base identities from aligned read
    for( my $i = 0; $i < length($s1); $i++ )
    {
	if( substr($s1,$i,1) ne "-" )
	{
	    ### update reference position
	    $ref_pos++;
	    
	    ### find cue
	    if( exists $o_cues{$ref_pos} )
	    {
		### make sure cue site has a correct base
		if( ! exists $o_cues{$ref_pos}{substr($s2,$i,1)} )
		{
		    last;
		}
		
		### require exact match at nearby positions
		for( my $k = $i - $o_match; $k <= $i + $o_match; $k++ )
		{
		    next if( $k == $i );
		    
		    if( substr($s1,$k,1) ne substr($s2,$k,1) || $k < 0 || $k >= length($s1) )
		    {
			last;
		    }
		}
		
		### record base identity at the specified position
		$data{$ref_pos} = substr($s2,$i,1);
	    }
	}
    }
    
    ### discard reads if any of the marker sites or nearby bases are incorrect
    my $skip = 0;

    for my $cue ( keys %data )
    {
	if( $data{$cue} eq "NA" )
	{
	    $skip = 1;
	    last;
	}
    }
    
    next if( $skip == 1 );
    
    my ($movie,$zmw,$other) = split(/\//,$qname);

    ### filter by number of passes
    next if( $$zmws{$movie}{$zmw}{"NP"} < $o_np );

    ### output all bases at all marker sites for each read
    printf( "%s,%s,%s,%s\n",
	    $movie,
	    $zmw,
	    $$zmws{$movie}{$zmw}{"NP"},
	    join(",", map { $data{$_} } (sort { $a <=> $b } keys %data)) );
}

close(BAM);

sub load_references {
    my ($file) = @_;
    
    my %refs = ();
    
    my $name = "";
    my $seq = "";
    
    open(FA,$file) || die "Can't open '$file': $!";
    
    while( my $line = <FA> )
    {
	chomp($line);
	
	if( substr($line,0,1) eq ">" )
	{
	    if( $name ne "" )
	    {
		$refs{$name} = uc($seq);
	    }
	    
	    $name = substr($line,1);
	    $seq = "";
	}
	else
	{
	    $seq .= $line;
	}
    }
    
    if( $name ne "" )
    {
	$refs{$name} = uc($seq);
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
