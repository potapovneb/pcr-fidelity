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
my @o_cues = ();
my %o_cues = ();

GetOptions( "cue=s" => \@o_cues );

for my $x ( split(/;/,join(";",@o_cues)) )
{
    my ($pos,$b0,$b1) = split(/,/,$x);

    $o_cues{$pos}{$b0} = 0;
    $o_cues{$pos}{$b1} = 1;
}

### command-line usage
if( @ARGV == 0 )
{
    print "usage: $0 markers1.csv markers2.csv ...\n";
    print "\n";
    print "options:\n";
    print "  --cues\tmarker sites (eg. --cue 77,T,C --cue 153,T,A --cue 279,G,T ...)\n";
    print "\n";
    exit;
}

### command-line arguments
my @files = @ARGV;

my %stat = ();

my @sorted = sort { $a <=> $b } keys %o_cues;

for my $file ( @files )
{
    ### load precomputed markers
    my $data = load_zmw_data($file);
    
    for my $movie ( keys %$data )
    {
	for my $zmw ( sort { $a <=> $b } keys %{$$data{$movie}} )
	{
	    ### {0,1} vector
	    my @cross = map { $o_cues{$_}{$$data{$movie}{$zmw}{$_}} } @sorted;
	    
	    ### count chimera switch points
	    my $cnum = 0;
	    
	    for( my $i = 1; $i < @cross; $i++ )
	    {
		if( $cross[$i-1] ne $cross[$i] )
		{
		    $cnum++;
		}
	    }
	    
	    ### cnum statistics
	    $stat{$cnum}++;
	}
    }
}

my $total_reads = 0;

for my $cnum ( keys %stat )
{
    $total_reads += $stat{$cnum};
}

print "Recombinations,Reads,Fraction\n";
for my $cnum ( sort { $a <=> $b } keys %stat )
{
    printf( "%i,%i,%.7f\n",
	    $cnum,
	    $stat{$cnum},
	    $stat{$cnum}/$total_reads );
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
