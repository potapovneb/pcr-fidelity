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
my $o_qv  = 0;
my $o_np  = 1;
my $o_rla = 0;
my $o_rlp = 0.0;
my $o_zm  = -1;
my $o_sna = 1e6;
my $o_snc = 1e6;
my $o_sng = 1e6;
my $o_snt = 20.0;
my $o_lb  = -1;
my $o_ub  = -1;
my $o_log = "";

GetOptions(
    "qv=f"  => \$o_qv,
    "np=f"  => \$o_np,
    "rla=f" => \$o_rla,
    "rlp=f" => \$o_rlp,
    "zm=f"  => \$o_zm,
    "sna=f" => \$o_sna,
    "snc=f" => \$o_snc,
    "sng=f" => \$o_sng,
    "snt=f" => \$o_snt,
    "lb=f"  => \$o_lb,
    "ub=f"  => \$o_ub,
    "log"   => \$o_log,
    );

### command-line arguments
if( @ARGV == 0 )
{
    print "usage: $0 [options] reference.fasta sampleDir1 SampleDir2 ...\n";
    print "\n";
    print "options:\n";
    print " --qv\t\tminimum quality value ($o_qv)\n";
    print " --np\t\tminimum number of passes ($o_np)\n";
    print " --rlp\t\tminimum read length ($o_rlp% of expected length)\n";
    print " --rla\t\tminimum read length ($o_rla nucleotides)\n";
    print " --zm\t\tanalyse a specific ZMW ($o_zm)\n";
    print " --sna\t\tmaximum SNR-A ($o_sna)\n";
    print " --snc\t\tmaximum SNR-C ($o_snc)\n";
    print " --sng\t\tmaximum SNR-G ($o_sng)\n";
    print " --snt\t\tmaximum SNR-T ($o_snt)\n";
    print " --lb\t\tskip first N bases ($o_lb)\n";
    print " --ub\t\tskip last N bases ($o_ub)\n";
    print " --log\t\tprint log info ($o_log)\n";
    print "\n";
    exit;
}

my $reference = shift @ARGV;
my @sampleDirs = @ARGV;

### load reference data
my $refs = load_references($reference);

### init counts
my @type = qw(AA AC AT AG CA CC CT CG TA TC TT TG GA GC GT GG Deletion Insertion);
my %count = map { $_ => 0 } @type;

### process samples one-by-one
for my $sampleDir ( @sampleDirs )
{
    for my $strand ( qw/fwd rev/ )
    {
	my $variants_file = sprintf( "%s/mutation/%s/variants.csv.7z", $sampleDir, $strand );
	my $zmws_file     = sprintf( "%s/mutation/%s/zmws.csv.7z",     $sampleDir, $strand );
	my $chimeric_file = sprintf( "%s/chimeric/%s/chimeric.csv.7z", $sampleDir, $strand );
	
	if( $o_log )
	{
	    print STDERR $variants_file, "\n";
	    print STDERR $chimeric_file, "\n";
	    print STDERR $zmws_file, "\n";
	}

	### load pacbio read stats
	my $zmws = load_zmw_data($zmws_file);
	
	### load chimeric reads
	my $blacklist = load_blacklist($chimeric_file);
	
	### filter variants
	filter($variants_file, $zmws, $blacklist, $strand);
    }
}

### print column headers
print join(",", @type) ,"\n";
print join(",", map { $count{$_} } @type), "\n";

sub filter {
    my ($variants_file,$zmws,$blacklist,$strand) = @_;

    my @head = ();
    
    open(VARIANTS,"7za e -so $variants_file |") || die "Can't open '$variants_file': $!";
    
    while( my $line = <VARIANTS> )
    {
	chomp($line);
	
	### parse line
	my @tokens = split(/,/,$line);

	my %entry = ();
	
	if( @head == 0 )
	{
	    @head = @tokens;
	    next;
	}
	else
	{
	    for( my $i = 0; $i < @tokens; $i++ )
	    {
		$entry{$head[$i]} = $tokens[$i];
	    }
	}
	
	if( ! exists $$refs{$entry{"Ref"}} )
	{
	    print STDERR "error :: missing reference '$entry{'Ref'}'\n\n";
	    exit;
	}
	
	### filter
	my $refLength = length($$refs{$entry{"Ref"}});
	
	### frequently used vars
	my $movie = $entry{"Movie"};
	my $zmw = $entry{"ZMW"};

	### filter by number of passes
	next if( $$zmws{$movie}{$zmw}{"NP"} < $o_np );
	
	### filter by read length
	next if( $$zmws{$movie}{$zmw}{"ReadLength"} < $o_rla );
	next if( $$zmws{$movie}{$zmw}{"ReadLength"} < $o_rlp * $refLength );
	
	### filter by SNR
	next if( $$zmws{$movie}{$zmw}{"SnrA"} > $o_sna );
	next if( $$zmws{$movie}{$zmw}{"SnrC"} > $o_snc );
	next if( $$zmws{$movie}{$zmw}{"SnrG"} > $o_sng );
	next if( $$zmws{$movie}{$zmw}{"SnrT"} > $o_snt );
	
	### filter by QV
	next if( $entry{"QV"} < $o_qv );
	
	### skip bases at the begining and at the end of reference (primer sites)
	next if( $o_lb != -1 && $entry{"Pos"} < $o_lb );
	next if( $o_ub != -1 && $entry{"Pos"} > $refLength - $o_ub );
	
	### exclude blacklisted reads
	next if( exists $$blacklist{$movie}{$zmw} );
	
	### count substitutions, deletions, insertions
	if( $entry{"Type"} eq "SNP" )
	{
	    my $rb = uc($entry{"RefBP"});
	    my $cb = uc($entry{"AltBP"});
	    
	    if( $strand eq "fwd" )
	    {
		$rb = complement($rb);
		$cb = complement($cb);
	    }
	    
	    $count{$rb.$cb}++;
	}
	elsif( $entry{"Type"} eq "INDEL" )
	{
	    $count{$entry{"IndelType"}} += $entry{"Length"};	    
	}
	
	### print out data lines (optional)
	if( $o_log )
	{
	    print STDERR $line, "\n";
	}
    }
    
    close(VARIANTS);
}

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

sub load_blacklist {
    my ($file) = @_;

    my %data = ();

    open(SR,"7za e -so $file |") || die "Can't open '$file': $!";

    while( my $line = <SR> )
    {
	chomp($line);

	my ($movie,$zmw) = split(/,/,$line);

	$data{$movie}{$zmw} = 1;
    }

    close(SR);

    return \%data;
}

sub load_zmw_data {
    my ($file) = @_;
    
    my %zmws = ();
    my @head = ();

    open(IN,"7za x -so $file |") || die "Can't open '$file': $!";
    
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
