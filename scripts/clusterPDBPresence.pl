#!/usr/bin/perl
#26 April 2016 - Adam Scott - 

use strict;
use warnings;

use IO::File;
use FileHandle;

my $usage = 'perl clusterPDBPresence.pl <pairwise> <clusters> <output_prefix> 
';

die $usage , unless @ARGV == 3;
my ( $pairwiseFile , $clustersFile , $output ) = @ARGV;

my $IN1 = FileHandle->new( "$pairwiseFile" , "r" );
if ( not defined $IN1 ) { die "ADSERROR: Could not open/read $pairwiseFile\n"; }
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
	my ( $gene1 , $mutation1 , $chain1 , $residue1 , $gene2 , $mutation2 , $chain2 , $residue2 , $pdbInfos ) = (split( "\t" , $line ))[0,4,5,6,9,13,14,15,19];
	$chain1 =~ s/\[(.*)\]/$1/;
	$chain2 =~ s/\[(.*)\]/$1/;
	my @pdbInfos = split( '\|' , $pdbInfos );
	foreach my $pdbInfo ( @pdbInfos ) {
		my $pdb = (split( ' ' , $pdbInfo ))[1];
		$structures{$gene1}{$mutation1}{$pdb}{$chain1} = $residue1;
		$structures{$gene2}{$mutation2}{$pdb}{$chain2} = $residue2;
		&checkMin( $minMax , $pdb , $gene1 , $chain1 , $residue1 );
		&checkMin( $minMax , $pdb , $gene2 , $chain2 , $residue2 );
		&checkMax( $minMax , $pdb , $gene1 , $chain1 , $residue1 );
		&checkMax( $minMax , $pdb , $gene2 , $chain2 , $residue2 );
	} #foreach pdbInfo in pdbInfos
}
$IN1->close();

sub checkMin {
	my ( $minMax , $pdb , $gene , $chain , $residue ) = @_;
	return unless ( $residue =~ /^\d+$/ ); 
	if ( exists $minMax->{$pdb}->{$gene}->{$chain}->{'min'} ) {
		if ( $minMax->{$pdb}->{$gene}->{$chain}->{'min'} <= $residue ) {
			return;
		}
	}
	$minMax->{$pdb}->{$gene}->{$chain}->{'min'} = $residue;
}

sub checkMax {
	my ( $minMax , $pdb , $gene , $chain , $residue ) = @_;
	return unless ( $residue =~ /^\d+$/ ); 
	if ( exists $minMax->{$pdb}->{$gene}->{$chain}->{'max'} ) {
		if ( $minMax->{$pdb}->{$gene}->{$chain}->{'max'} >= $residue ) {
			return;
		}
	}
	$minMax->{$pdb}->{$gene}->{$chain}->{'max'} = $residue;
}

while ( my $line = <$IN2> ) {
	if ( $line !~ /Cluster/ ) {
		chomp( $line );
		my ( $cluster , $cgene , $mutation , $recurrence ) = (split( "\t" , $line ))[0,1,2,6];
		foreach my $pdb ( keys %{$structures{$cgene}{$mutation}} ) {
			foreach my $chain ( keys %{$structures{$cgene}{$mutation}{$pdb}} ) {
				$represent{$pdb}{$cgene}{$chain}{$cluster}{$mutation.":".$structures{$cgene}{$mutation}{$pdb}{$chain}} = $recurrence;
			}
		}
	}
}
$IN2->close();

$OUT1->print( "PDB_ID\tGene\tChain\tCluster\tMinResidue\tMaxResidue\tnMutations\tnResidues\tTotalRecurrence\tMutations|Position\n" );
my %complex;
foreach my $pdb ( sort keys %represent ) {
	foreach my $gene ( sort keys %{$represent{$pdb}} ) {
		foreach my $chain ( sort keys %{$represent{$pdb}{$gene}} ) {
			$chain =~ m/\[(.*)\]/;
			foreach my $cluster ( sort keys %{$represent{$pdb}{$gene}{$chain}} ) {
				my ( $mutres , $mutation , $position );
				my @mutations;
				my %residues;
				my $recurrence = 0;
				foreach my $mutpos ( sort keys %{$represent{$pdb}{$gene}{$chain}{$cluster}} ) {
					( $mutation , $position ) = split( ":" , $mutpos );
					$recurrence += $represent{$pdb}{$gene}{$chain}{$cluster}{$mutpos};
					$mutres = join( "|" , ( $mutation , $position ) );
					push @mutations , $mutres;
					$residues{$position} = 1;
				}
				my @logline = ( $pdb , $gene , $chain , $cluster , $mutation , $position , $recurrence );
				print join( "\t" , @logline )."\n";
				if ( exists $minMax->{$pdb}->{$gene}->{$chain} ) {
					my $min = $minMax->{$pdb}->{$gene}->{$chain}->{'min'};
					my $max = $minMax->{$pdb}->{$gene}->{$chain}->{'max'};
					my @outline = ( $pdb , $gene , $chain , $cluster , $min , $max , scalar @mutations , scalar keys %residues , $recurrence , join( ";" , @mutations ) );
					$complex{$cluster}{$pdb}{$gene}{$chain} = \@outline;
					$OUT1->print( join( "\t" , @outline )."\n" );
				}
			}
		}
	} #foreach pdb in represent => cluster
} #foreach cluster in represent
$OUT1->close();

$OUT2->print( "Cluster\tPDB_ID\tGene\tChain\tnMutations\tnResidues\tTotalRecurrence\tMutations|Position\n" );
foreach my $cluster ( sort keys %complex ) {
	foreach my $pdb ( sort keys %{$complex{$cluster}} ) {
		my ( @mutations , @geneChains );
		my ( $mutations , $residues , $recurrence );
		foreach my $gene ( sort keys %{$complex{$cluster}{$pdb}} ) {
			my @chains;
			foreach my $chain ( sort keys %{$complex{$cluster}{$pdb}{$gene}} ) {
				my @outline = @{$complex{$cluster}{$pdb}{$gene}{$chain}};
				$mutations += $outline[6];
				$residues += $outline[7];
				$recurrence += $outline[8];
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
