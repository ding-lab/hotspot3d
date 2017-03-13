#!/usr/bin/perl
#24 June 2016 - Adam Scott - 

use strict;
use warnings;

use IO::File;
use FileHandle;

my $usage = 'perl clusterPDBPresence.drug.pl <drugport_pairs> <clusters> <output_prefix> 
';

die $usage , unless @ARGV == 3;
my ( $drugportFile , $clustersFile , $output ) = @ARGV;

my $IN1 = FileHandle->new( "$drugportFile" , "r" );
if ( not defined $IN1 ) { die "ADSERROR: Could not open/read $drugportFile\n"; }
my $IN2 = FileHandle->new( "$clustersFile" , "r" );
if ( not defined $IN2 ) { die "ADSERROR: Could not open/read $clustersFile\n"; }
my $OUT1 = FileHandle->new( "$output.chains" , "w" );
if ( not defined $OUT1 ) { die "ADSERROR: Could not open/write $output.chains\n"; }
my $OUT2 = FileHandle->new( "$output.xmer" , "w" );
if ( not defined $OUT2 ) { die "ADSERROR: Could not open/write $output.xmer\n"; }

my %structures;
my %represent;
my $minMax = {};
while ( my $line = <$IN1> ) {
	chomp( $line );
	next if ( $line =~ /Drug/ );
	#1	Drug
	#3	PDB_ID
	#4	Chain
	#5	Compound_Location
	#7	Gene
	#11	Amino_Acid_Change
	#12	Chain
	#13	Mutation_Location_In_PDB
	#14	Res_Name
	my ( $drug , $pdb , $chain1 , $drugPosition , $gene2 , $mutation2 , $chain2 , $residue2 ) = (split( "\t" , $line ))[0,2,3,4,6,10,11,12];
	$chain1 =~ s/\[(.*)\]/$1/;
	$chain2 =~ s/\[(.*)\]/$1/;
	if ( 0 ) {
		print( $gene2."\t" );
		print( $drug."\t" );
		print( $pdb."\t" );
		print( $chain1."\t" );
		print( $drugPosition."\t" );
		print( $chain2."\t" );
		print( $residue2."\n" );
	}
	$structures{$gene2}{$drug}{$pdb}{$chain1} = $drugPosition;
	$structures{$gene2}{$mutation2}{$pdb}{$chain2} = $residue2;
}
$IN1->close();

while ( my $line = <$IN2> ) {
	if ( $line !~ /Cluster/ ) {
		chomp( $line );
		my ( $cluster , $cgene , $mutation , $recurrence ) = (split( "\t" , $line ))[0,1,2,6];
		if ( $mutation !~ /^p\./ ) {
			my $temp = $cgene;
			$cgene = $mutation;
			$mutation = $temp;
		}
		foreach my $pdb ( keys %{$structures{$cgene}{$mutation}} ) {
			foreach my $chain ( keys %{$structures{$cgene}{$mutation}{$pdb}} ) {
				$represent{$pdb}{$cgene}{$chain}{$cluster}{$mutation.":".$structures{$cgene}{$mutation}{$pdb}{$chain}} = $recurrence;
			}
		}
	}
}
$IN2->close();

$OUT1->print( "PDB_ID\tGene\tChain\tCluster\tnMutations\tnResidues\tTotalRecurrence\tMutations|Position\n" );
my %complex;
foreach my $pdb ( sort keys %represent ) {
	foreach my $gene ( sort keys %{$represent{$pdb}} ) {
		foreach my $chain ( sort keys %{$represent{$pdb}{$gene}} ) {
			foreach my $cluster ( sort keys %{$represent{$pdb}{$gene}{$chain}} ) {
				my ( $mutres , $mutation , $position );
				my @mutations;
				my %residues;
				my $recurrence = 0;
				foreach my $mutpos ( sort keys %{$represent{$pdb}{$gene}{$chain}{$cluster}} ) {
					( $mutation , $position ) = split( ":" , $mutpos );
					if ( $mutation =~ m/p\./ ) {
						$recurrence += $represent{$pdb}{$gene}{$chain}{$cluster}{$mutpos};
					}
					$mutres = join( "|" , ( $mutation , $position ) );
					push @mutations , $mutres;
					$residues{$position} = 1;
				}
				my @logline = ( $pdb , $gene , $chain , $cluster , $mutation , $position , $recurrence );
				#print join( "\t" , @logline )."\n";
				my @outline = ( $pdb , $gene , $chain , $cluster , scalar @mutations , scalar keys %residues , $recurrence , join( ";" , @mutations ) );
				$complex{$cluster}{$pdb}{$gene}{$chain} = \@outline;
				$OUT1->print( join( "\t" , @outline )."\n" );
			}
		}
	} #foreach pdb in represent => cluster
} #foreach cluster in represent
$OUT1->close();

$OUT2->print( "Cluster\tPDB_ID\tGene\tChain\tnMutationsDrugs\tnPositions\tTotalRecurrence\tMutationsDrugs|Position\n" );
foreach my $cluster ( sort keys %complex ) {
	foreach my $pdb ( sort keys %{$complex{$cluster}} ) {
		my ( @mutations , @geneChains );
		my ( $mutations , $residues , $recurrence );
		foreach my $gene ( sort keys %{$complex{$cluster}{$pdb}} ) {
			my @chains;
			foreach my $chain ( sort keys %{$complex{$cluster}{$pdb}{$gene}} ) {
				my @outline = @{$complex{$cluster}{$pdb}{$gene}{$chain}};
				$mutations += $outline[4];
				$residues += $outline[5];
				$recurrence += $outline[6];
				push @chains , $chain;
				push @mutations , $chains[-1]."\\";
				$mutations[-1] .= $outline[-1];
			} #foreach chain
			push @geneChains , $gene."|".join( "/" , @chains );
		} #foreach gene
		my @complexLine = ( $cluster , $pdb , join( ";" , @geneChains ) , $mutations , $residues , $recurrence , join( ";" , @mutations ) );
		$OUT2->print( join( "\t" , @complexLine )."\n" );
	} #foreach pdb
} #foreach cluster
$OUT2->close();
