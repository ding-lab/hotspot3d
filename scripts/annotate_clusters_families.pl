#!/bin/perl
#15 March 2015 - Adam D Scott - 
 
use strict;
use warnings;
 
my $usage = 'perl annotate_families.pl <families> <clusters> <output>  
';
 
die $usage , unless @ARGV == 3;
my ( $families , $clusters , $output ) = @ARGV;
 
my %families;
my %fams;
open ( IN , "<$families" ) or die "Cannot open $families: $!";
while ( <IN> ) {
	chomp;
	my @line = split "\t" , $_;
 
	my $gene = $line[1];
	my $family = $line[10];
	$families{$gene}{$family} = 1;
	$fams{$family} = 0;
}
close IN;

my %counted;
my $unclassified = "Unclassified";
my $notgene = "NA";
my $total = 0;
my %lines;
open ( OUT , ">$output" ) or die "Cannot open $output: $!";
open ( IN2 , "<$clusters" ) or die "Cannot open $clusters: $!";
while ( <IN2> ) {
	chomp;
	if ( /Cluster/ ) {
		print OUT $_."\tHGNC_Gene_Families\n";
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
