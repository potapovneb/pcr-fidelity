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
my $o_reverse = "";
my $o_complement = "";

GetOptions(
    "reverse" => \$o_reverse,
    "complement" => \$o_complement,
    );

### command-line arguments
my $file  = shift @ARGV;
my $start = shift @ARGV;
my $end   = shift @ARGV;

### obtain sequence
my $sequence = "";

open(IN,$file) || die "Can't open '$file': $!";

while( my $line = <IN> )
{
    chomp($line);

    next if( index($line,">") == 0 );

    $sequence .= $line;
}

close(IN);

### extract fragment
my $length = $end - $start + 1;
my $fragment = substr($sequence,$start-1,$length);

### transform (optinal)
if( $o_complement )
{
    $fragment = complement($fragment);
}

if( $o_reverse )
{
    $fragment = scalar reverse($fragment);
}

### print resulting string
print $fragment;

sub complement {
    my ($seq) = @_;

    my %bases = (
	"A" => "T",
	"C" => "G",
	"G" => "C",
	"T" => "A",
	"a" => "t",
	"c" => "g",
	"g" => "c",
	"t" => "a",
	"-" => "-",
	);

    my $complement = join("",map { $bases{$_} } split(//,$seq));

    return $complement;
}
