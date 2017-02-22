#!/usr/bin/perl
#08 February 2017 - Adam D Scott - 

use strict;
use warnings;

use IO::File;
use FileHandle;

my $usage = 'perl nStructures.pl <hupt> <output> 
';

die $usage , unless @ARGV == 2;
my ( $hupt , $output ) = @ARGV;

my $IN1 = FileHandle->new( "$hupt" , "r" );
if ( not defined $IN1 ) { die "ADSERROR: Could not open/read $hupt\n"; }

my $OUT = FileHandle->new( "$output" , "w" );
if ( not defined $OUT ) { die "ADSERROR: Could not open/write $output\n"; }

my $nStructures = 0;
my $nGenesWithStructure = 0;
my $nGenes = 0;
while ( my $line = <$IN1> ) {
	chomp( $line );
	my @line = split( /\t/ , $line );
	next if ( $line[2] eq "N/A" );
	my @pdbs = split( /\ / , $line[2] );
	my $nStructureForGene = scalar @pdbs;
	if ( $nStructureForGene > 0 ) {
		$nGenesWithStructure += 1;
		$nStructures += $nStructureForGene;
		$OUT->print( join( "\t" , ( $line[0] , $nStructureForGene ) )."\n" );
	}
}
$IN1->close();
$OUT->close();
print join( "\t" , ( $nStructures , $nGenesWithStructure , ( $nStructures / $nGenesWithStructure ) ) )."\n";
