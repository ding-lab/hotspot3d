#!/bin/perl
#3 April 2015 - Adam D Scott - 
 
use strict;
use warnings;
 
my $usage = 'perl annotate_clusters_domains.pl <pairwise> <drugClean> <clusters>
';
 
die $usage , unless @ARGV == 3;
my ( $pairwise , $drugClean , $clusters ) = @ARGV;
 
my %domains;
my %variants;
open ( IN , "<$pairwise" ) or die "Cannot open $pairwise: $!";
while( <IN> ) {
	chomp; if ( /Gene/ ) { next; }
	my @line = split( "\t" , $_ );
	my $gene1 = $line[0];
	my $AA1 = $line[4];
	my $domain1 = $line[7];
	$domain1 = &filter_eco( $domain1 );
	my $gene2 = $line[9];
	my $AA2 = $line[13];
	my $domain2 = $line[16];
	$domain2 = &filter_eco( $domain2 );

	$domains{$gene1}{$AA1}{$domain1} = 1;
	$domains{$gene2}{$AA2}{$domain2} = 1;
}
close IN;

open ( IN , "<$drugClean" ) or die "Cannot open $drugClean: $!";
while( <IN> ) {
	chomp; if ( /Gene/ ) { next; }
	my @line = split( "\t" , $_ );
	my $gene1 = $line[5];
	my $AA1 = $line[6];
	my $domain1 = $line[7];
	$domain1 = &filter_eco( $domain1 );

	$domains{$gene1}{$AA1}{$domain1} = 1;
}
close IN;

open ( OUT , ">domains.$clusters" ) or die "Cannot open domains.$clusters: $!";
open ( IN , "<$clusters" ) or die "Cannot open $clusters: $!";
while ( <IN> ) {
	chomp;	
	if ( /Cluster/ ) { 
		print OUT $_."\tProtein_Domain\n";
		next; 
	}
	my @line = split( "\t" , $_ );
	my $genedrug = $line[1];
	my $AAgene = $line[2];

	if ( exists $domains{$genedrug} ) {
		my @domains = keys %{$domains{$genedrug}{$AAgene}};
		print OUT join( "\t" , ( @line , join( "|" , sort @domains ) ) )."\n";
	} else {
		print OUT join( "\t" , @line )."\tNULL\n";
	}
}
close IN;
close OUT;

sub filter_eco {
	my ( $domain ) = @_;

	$domain =~ s/ {ECO.*//;
	$domain =~ s/^{ECO.*//;
	if ( length( $domain ) == 0 || $domain =~ /N\/A/ ) {
		$domain = "NULL";
	}

	return $domain;
}
