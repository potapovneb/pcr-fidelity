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
my $o_fwd = "subreads.fwd.sam";
my $o_rev = "subreads.rev.sam";

### command-line arguments
if( @ARGV == 0 )
{
    print "usage: $0 [options] input.bam clusters.csv.bz2\n";
    print "\n";
    print "options:\n";
    print "  --fwd\t\toutput file for forward reads ($o_fwd)\n";
    print "  --rev\t\toutput file for reverse reads ($o_rev)\n";
    print "\n";
    exit;
}

my $bamfile = shift @ARGV;
my $csvfile = shift @ARGV;

### read mapping information
my %csv = ();

open(CSV,"bzip2 -cd $csvfile |") || die "Can't open '$csvfile': $!";

while( my $line = <CSV> )
{
    chomp($line);
    
    my ($qname,$flag) = split(/,/,$line);

    ### save mapping direction
    $csv{$qname} = $flag;
}

close(CSV);

### split BAM file
open(FWD,">",$o_fwd) || die "Can't write '$o_fwd': $!";
open(REV,">",$o_rev) || die "Can't write '$o_rev': $!";

open(BAM,"samtools view -h $bamfile |") || die "Can't open '$bamfile': $!";

while( my $line = <BAM> )
{
    if( substr($line,0,1) eq "@" )
    {
	### copy BAM headers to both files
	print FWD $line;
	print REV $line;
    }
    else
    {
	my ($qname,@other) = split(/\t/,$line);
	
	if( exists $csv{$qname} )
	{
	    if( $csv{$qname} == 0 )
	    {
		print FWD $line;
	    }
	    else
	    {
		print REV $line;
	    }
	}
    }
}

close(BAM);

close(FWD);
close(REV);
