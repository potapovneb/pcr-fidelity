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
    print "usage: $0 input.bam output.csv\n\n";
    exit;
}

my $bamfile = shift @ARGV;
my $outfile = shift @ARGV;

### generate CSV headers
open(CSV,">",$outfile) || die "Can't write '$outfile': $!";
print CSV "Movie,ZMW,SnrA,SnrC,SnrG,SnrT,RQ,AvgZScore,ReadLength,NP\n";

open(BAM,"samtools view $bamfile |") || die "Can't open '$bamfile': $!";

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
    my $za = -1;
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
	elsif( $name eq "za" )
	{
	    $za = $value;
	}
	elsif( $name eq "sn" )
	{
	    my @val = split(/,/,$value);
	    shift @val;
	    @sn = @val;
	}
    }
    
    my ($movie,$zmw,@other) = split(/\//,$line);

    printf( CSV "%s,%s,%f,%f,%f,%f,%f,%f,%i,%i\n",
	    $movie,
	    $zmw,
	    @sn,
	    $rq,
	    $za,
	    $rl,
	    $np );
}

close(BAM);

close(CSV);
