#!/usr/bin/perl
# May 2017 Kuan-lin Huang

use strict;
use warnings;

use IO::File;
use FileHandle;

my $usage = 'perl filter_PDB.pl <pairwise file> <pass PDB file>
';

die $usage , unless @ARGV == 2;

my $pairwise = shift;
my $pass_PDB = shift;
my $pass_pairwise_out = $pairwise.".pass";

# pass PDB structures
my $pass_PDB_structures = parse_pass_PDB($pass_PDB);
# valid phosphosites
my %phosphorylated_aa = ('S' => 1,'T' => 1,'Y' => 1,'D' => 1,'H' => 1); 

open ( IN , "<$pairwise" ) or die "Cannot open $pairwise: $!";
my $OUT = FileHandle->new( $pass_pairwise_out , "w" );
if ( not defined $OUT ) { die "ADSERROR: Could not open/write $pass_pairwise_out\n"; }

my $head = <IN>;
chomp $head; 
print "$head\n";

while ( <IN> )
{
	chomp;
	my @line = split "\t" , $_;
	# filter out structures; include only passing structures
	my @distanceInfo = split " ", $line[18];
	my $PDB_structure = $distanceInfo[1];
	if (! exists($pass_PDB_structures->{$PDB_structure})){
		print STDERR "Discard pairs due to not-passed PDB structure: $PDB_structure\n";
		print STDERR "$_\n";
		next;
	}

	# filter out phosphosites that is not S, T, Y, D, or H (known phosphorylated aa in human)
	my ($site1, $feature1, $site2, $feature2) = ($line[4],$line[7],$line[12],$line[15]);
	my $aa1 = substr($site1, 2, 1); my $aa2 = substr($site2, 2, 1); 
	if ($feature1 =~ /Phospho/ && !exists($phosphorylated_aa{$aa1})){
		print STDERR "Discard phosphosite due to non-phosphorylated amino acid: $aa1\n";
		print STDERR "$_\n";
		next;
	}
	if ($feature2 =~ /Phospho/ && !exists($phosphorylated_aa{$aa2})){
		print STDERR "Discard phosphosite due to non-phosphorylated amino acid: $aa1\n";
		print STDERR "$_\n";
		next;
	}
	
	#print "$_\n";
	print $OUT "$_\n";
}

close(IN);
$OUT->close();

sub parse_pass_PDB {
	my $pass_PDB =shift;
	my $pass_PDB_structures ={};
		open(FILE, $pass_PDB ) or die "Unable to open file $pass_PDB due to $!";
		while(<FILE>) {
			chomp;
			my @line = split "\t" , $_ ;
			my $structure = $line[2];
			$pass_PDB_structures->{$structure}=1;
		}
		close FILE;
		return $pass_PDB_structures;
}