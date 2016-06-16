#!/usr/bin/perl
#15 June 2016 - Adam D Scott - get .pdb files from a gene list via HGNC download file
#  make sure your HGNC download includes UniProt ID and the style of gene name you work with
#  HGNC downloads can be found here: http://www.genenames.org/cgi-bin/download

use strict;
use warnings;

use IO::File;
use FileHandle;

use LWP::Simple;

my $usage = 'perl getPDBs.pl <geneList> <HGNC> <output> [ (\d=0) geneColumn , (\d=2) uniprotColumn , (bool=0) writePDB , (bool=1) compressPDB ]

ABOUT: get .pdb files from a gene list via HGNC download file

TIPS:  1) Many .pdb files take a lot of space, so please BE CAREFUL with this.
	   2) Run this script wherever you want your .pdb files.
       3) Make sure your HGNC download includes UniProt ID and the style of gene name you work with
          HGNC downloads can be found here: http://www.genenames.org/cgi-bin/download
';

die $usage , unless @ARGV >= 3;
my ( $geneList , $HGNC , $output , $geneColumn , $uniprotColumn , $writePDB , $compress ) = @ARGV;
if ( not defined $geneColumn ) {
	$geneColumn = 0;
	print STDERR "ADSWarning: gene column not provided for HGNC file. Assuming it is column 0 (0-base).\n";
}
if ( not defined $uniprotColumn ) {
	$uniprotColumn = 2;
	print STDERR "ADSWarning: uniprot column not provided for HGNC file. Assuming it is column 2 (0-base).\n";
}
if ( not defined $writePDB ) {
	$writePDB = 0;
	print STDERR "ADSWarning: write PDB status not provided. Assuming not to write.\n";
}
if ( not defined $compress ) {
	$compress = 1;
	print STDERR "ADSWarning: compress PDB status not provided. Assuming to compress.\n";
}

my $IN1 = FileHandle->new( "$geneList" , "r" );
if ( not defined $IN1 ) { die "ADSERROR: Could not open/read $geneList\n"; }

my $IN2 = FileHandle->new( "$HGNC" , "r" );
if ( not defined $IN2 ) { die "ADSERROR: Could not open/write $HGNC\n"; }

my $OUT = FileHandle->new( "$output" , "w" );
if ( not defined $OUT ) { die "ADSERROR: Could not open/write $output\n"; }

my @genes;
while ( my $line = <$IN1> ) {
	chomp( $line );
	push @genes , $line;
}
$IN1->close();

my %map;
while ( my $line = <$IN2> ) {
	chomp( $line );
	my ( $hugo , $uniprot ) = (split( "\t" , $line ))[$geneColumn,$uniprotColumn];
	if ( grep{ $_ eq $hugo } @genes ) {
		$map{$hugo} = $uniprot;
	}
}
$IN2->close();

my $uniprotBase = "http://www.uniprot.org/uniprot/";
$OUT->print( "Gene\tUniProt\tPDB_ID\tMethod\tResolution\tChain_Start\tChain_End\tChains\n" );
my ( $pdbURL , $pdbPage , $PDB );
my $pdbBase = "http://www.rcsb.org/pdb/files/";
foreach my $hugo ( keys %map ) {
	my ( $url , $page , $pdb , $method , $resolution , $positions , $chains , $residues , $begin , $end );
	my $uniprot = $map{$hugo};
	print STDOUT $hugo."\t";
	if ( defined $uniprot ) {
		$url = $uniprotBase.$uniprot.".txt";
		print STDOUT $url."\t";
		$page = get( $url );
		if ( $page ) {
			print STDOUT "success\t";
			foreach my $line ( split( "\n" , $page ) ) {
				if ( $line =~ m/^DR/ ) {
					if ( $line =~ m/PDB;/ ) { #PDBsum lines also have PDB ID
						chomp( $line );
						my @line = split( /\s*;/ , $line );
						( $pdb , $method , $resolution , $positions ) = @line[1,2,3,4];
						$pdb =~ s/\s*(\w\w\w\w)/$1/;
						print STDOUT $pdb."\t";
						( $chains , $residues ) = split( "=" , $positions );
						if ( $residues =~ /-/ ) {
							( $begin , $end ) = split( "-" , $residues );
						} else {
							$begin = $residues;
							$end = $residues;
						}
						$end =~ s/(\d*)\./$1/;
						$OUT->print( join( "\t" , ( $hugo , $uniprot , $pdb , $method , $resolution , $begin , $end , $chains ) )."\n" );

						if ( $writePDB ) {
							$pdbURL = $pdbBase.$pdb.".pdb";
							if ( $compress ) {
								$pdbURL .= ".gz";
							} else {
								$PDB = FileHandle->new( $pdb.".pdb" , "w" );
								if ( not defined $PDB ) { die "ADSERROR: Could not open/write ".$pdb.".pdb\n"; }
							}
							print STDOUT $pdbURL."\t";
							$pdbPage = get( $pdbURL );
							if ( not $compress and $pdbPage ) {
								print STDOUT "writing\t";
								foreach my $line ( split( "\n" , $pdbPage ) ) {
									$PDB->print( $line."\n" );
								}
								$PDB->close();
							} else {
								print STDOUT "downloaded\t";
							}
						}
						print STDOUT "\n";
					}
				} #if uniprot line starts with DR
			} #foreach line in uniprot page
		} #if uniprot page
	} #if uniprot id
	print STDOUT "\n";
} #foreach hugo
$OUT->close();
