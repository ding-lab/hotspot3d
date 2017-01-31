#!/bin/perl
#15 March 2015 - Adam D Scott - 
 
use strict;
use warnings;
 
my $usage = 'perl annotate_families.pl <hgnc> <protein_kinases> <clusters> <output>  
';
 
die $usage , unless @ARGV == 4;
my ( $hgnc , $prokin , $clusters , $output ) = @ARGV;
 
my $proteinkinase = "Protein Kinase";
my $unclassified = "Unclassified";
my $unclassifiedkinase = $unclassified.", Protein Kinase";
my $notgene = "NA";

my %families;
my %fams;
open ( IN , "<$hgnc" ) or die "Cannot open $hgnc: $!";
while ( <IN> ) {
	chomp;
	my @line = split "\t" , $_;
 
	my $gene = $line[1];
	my $family = $line[10];
	$families{$gene}{$family} = 1;
	$fams{$family} = 0;
}
close IN;

open ( IN , "<$prokin" ) or die "Cannot open $prokin: $!";
while ( <IN> ) {
	chomp;
	my @line = split "\t" , $_;
 
	my $gene = $line[0];
	$families{$gene}{$proteinkinase} = 1;
	$fams{$proteinkinase} = 0;
}
close IN;

my %counted;
my $total = 0;
my %lines;
open ( OUT , ">$output" ) or die "Cannot open $output: $!";
open ( IN2 , "<$clusters" ) or die "Cannot open $clusters: $!";
while ( <IN2> ) {
	chomp;
	if ( /Cluster/ ) {
		print OUT $_."\tGene_Families\n";
		next;
	}
	my @line = split "\t" , $_;
	
	my $id = $line[0];
	my $gene = $line[1];
	my $aachange = $line[2];
	my $mutations = $line[6];
	my @families;
	if ( $aachange =~ /p\./ ) {
		if ( exists $families{$gene} ) {
			@families = sort keys %{$families{$gene}};
			foreach my $family ( keys %{$families{$gene}} ) {
				if ( exists $fams{$family} ) {
					$fams{$family} += $mutations;
				} else {
					$fams{$family} = $mutations;
				}
			}
		} else {
			@families = ( $unclassified );
			if ( exists $fams{$unclassified} ) {
				$fams{$unclassified} += $mutations;
			} else {
				$fams{$unclassified} = $mutations;
			}
		}
		if ( not exists $counted{$gene} ) {
			$total++;
		}
	} else { 
		@families = ( $notgene );
	}
	print OUT join( "\t" , ( @line , join( "|" , @families ) ) )."\n";
}
close IN2;
close OUT;

print "Gene_Family\tMutations\tPercentage\n";
foreach my $family ( keys %fams ) {
	my $n = $fams{$family};
	if ( $n > 0 ) {
		my $percent = 100*$n/$total;
		print $family."\t".$n."\t";
		printf "%.3f\n" , $percent;
	}
}
