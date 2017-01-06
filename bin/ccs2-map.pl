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

if( @ARGV == 0 )
{
    print "usage: $0 input.bam\n\n";
    exit;
}

my $bamfile = shift @ARGV;

### read BAM file
open(BAM,"samtools view $bamfile |") || die "Can't open '$bamfile': $!\n";

while( my $line = <BAM> )
{
    chomp($line);
    
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
	 $seq,
	 $qual,
	 @tag ) = split(/\t/,$line);
    
    ### skip unmapped reads
    next if( $flag & 0x4 );
    
    ### mapping direction
    if( $flag & 0x10 )
    {
	printf( "%s,%i\n", $qname, 1 );
    }
    else
    {
	printf( "%s,%i\n", $qname, 0 );
    }
}

close(BAM);
