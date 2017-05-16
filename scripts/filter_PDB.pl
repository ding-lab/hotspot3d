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
#print "$head\n";

my $hasBadSites = {};
while ( <IN> )
{
	chomp;
	my @line = split "\t" , $_;
	# filter out structures; include only passing structures
	my $structureFail = 0;
	foreach my $info ( split /\|/ , $line[-1] ) {
		my $PDB_structure = ( split( /\ / , $info ) )[1];
		if ( not exists($pass_PDB_structures->{$PDB_structure} ) ) {
			print STDERR "Discard pairs due to not-passed PDB structure: $PDB_structure\n";
			print STDERR "$_\n";
			$structureFail = 1;
			last;
		}
	}
	next if ( $structureFail );

	# filter out phosphosites that is not S, T, Y, D, or H (known phosphorylated aa in human)
	# musites columns 12,15; sites columns 3,6,11,14
	if ( $pairwise =~ /\.musites$/ ) {
		my ( $site2 , $feature2 ) = @line[12,15];
		my $ok = &filterSite( $site2 , $feature2 , $_ , $hasBadSites );
		if ( not $ok ) { next; }
	} elsif ( $pairwise =~ /\.sites$/ ) {
		my ( $site1 , $feature1 , $site2 , $feature2 ) = @line[3,6,11,14];
		my $ok = &filterSite( $site1 , $feature1 , $_ , $hasBadSites );
		if ( not $ok ) { next; }
		$ok = &filterSite( $site2 , $feature2 , $_ , $hasBadSites );
		if ( not $ok ) { next; }
	}
	
	#print "$_\n";
	print $OUT "$_\n";
}

close(IN);
$OUT->close();

print STDOUT "# These structures had bad sites but were not filtered by PDB ID:\n";
foreach my $pdb ( sort keys %{$hasBadSites} ) {
	print STDOUT $pdb."\n";
}

sub filterSite {
	my ( $site , $feature , $line , $hasBadSites ) = @_;
	my $aa = substr( $site , 2 , 1 );
	if ( $feature =~ /Phospho/ && !exists( $phosphorylated_aa{$aa} ) ) {
		my $distInfo = ( split( /\t/ , $line ) )[-1];
		foreach my $info ( split( /\|/ , $distInfo ) ) {
			print STDERR $info."\n";
			my $pdb = ( split( /\ / , $info ) )[1];
			$hasBadSites->{$pdb} = 1;
		}
		print STDERR "Discard phosphosite due to non-phosphorylated amino acid: $aa\n";
		print STDERR "$line\n";
		return 0;
	}
	return 1;
}

sub parse_pass_PDB {
	my $pass_PDB =shift;
	my $pass_PDB_structures ={};
	open(FILE, $pass_PDB ) or die "Unable to open file $pass_PDB due to $!";
	while(<FILE>) {
		chomp;
		my @line = split "\t" , $_ ;
		next if ( scalar @line < 4 );
		my $structure = $line[2];
		$pass_PDB_structures->{$structure}=1;
	}
	close FILE;
	return $pass_PDB_structures;
}
