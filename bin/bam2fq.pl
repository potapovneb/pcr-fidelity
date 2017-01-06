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
my $o_qv = 0;
my $o_np = 1;
my $o_rl = 1;
my $o_zm = -1;
my $o_sna = 1e6;
my $o_snt = 1e6;
my $o_sng = 1e6;
my $o_snc = 1e6;
my $o_fasta = "";

GetOptions(
    "qv=f" => \$o_qv,
    "np=f" => \$o_np,
    "rl=f" => \$o_rl,
    "zm=f" => \$o_zm,
    "sna=f" => \$o_sna,
    "snc=f" => \$o_snc,
    "sng=f" => \$o_sng,
    "snt=f" => \$o_snt,
    "fasta" => \$o_fasta,
    );

### command-line arguments
if( @ARGV == 0 )
{
    print "usage: $0 input.bam output.fastq\n";
    print "\n";
    print "options:\n";
    print " --qv ($o_qv)\n";
    print " --np ($o_np)\n";
    print " --rl ($o_rl)\n";
    print " --zm ($o_zm)\n";
    print " --sna ($o_sna)\n";
    print " --snc ($o_snc)\n";
    print " --sng ($o_sng)\n";
    print " --snt ($o_snt)\n";
    print " --fasta ($o_fasta)\n";
    print "\n";
    exit;
}

my $bamfile = shift @ARGV;
my $outfile = shift @ARGV;

### open BAM file
open(BAM,"samtools view $bamfile |") || die "Can't open '$bamfile': $!";

open(FQ,">",$outfile) || die "Can't write '$outfile': $!";

while( my $line = <BAM> )
{
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
    
    my $rl = length($read);
    my $np = -1;
    my $rq = -1;
    my $zm = -1;
    my @sn = (-1,-1,-1,-1);
    
    for my $tag ( @tag )
    {
	my ($name,$type,$value) = split(/:/,$tag);
	
	if( $name eq "np" )
	{
	    $np = $value;
	}
	elsif( $name eq "rq" )
	{
	    $rq = $value;
	}
	elsif( $name eq "zm" )
	{
	    $zm = $value;
	}
	elsif( $name eq "sn" )
	{
	    my @val = split(/,/,$value);
	    shift @val;
	    @sn = @val;
	}
    }
    
    ### filtering (optional)
    next if( $np < $o_np );
    next if( $rl < $o_rl );
    next if( $zm != $o_zm && $o_zm != -1 );
    next if( $sn[0] > $o_sna );
    next if( $sn[1] > $o_snc );
    next if( $sn[2] > $o_sng );
    next if( $sn[3] > $o_snt );
    
    if( $o_fasta )
    {
	print FQ ">", $qname, "\n";
	print FQ $read, "\n";
    }
    else
    {
	print FQ "\@", $qname, "\n";
	print FQ $read, "\n";
	print FQ "+\n";
	print FQ $qual, "\n";
    }
}

close(FQ);

close(BAM);
