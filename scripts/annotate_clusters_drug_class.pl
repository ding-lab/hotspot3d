#!/bin/perl
#26 February 2015 - Adam D Scott - 
 
use strict;
use warnings;
 
my $usage = 'perl drug_class_annotate_clusters.pl <clusters> <drugclasses> <output>  
';
 
die $usage , unless @ARGV == 3;
my ( $clusters , $drugclasses , $output ) = @ARGV;
 
my %clusters;
my %classes;

open ( IN , "<$drugclasses" ) or die "Cannot open $drugclasses: $!";
while ( <IN> ) {
	chomp;
	my @line = split( "\t" , $_ );
 	if ( /Total/ ) { next; }
	##from two-column association list: drug \t class
	my $drug = $line[0];
	my $class = $line[1];
	$class =~ s/\"//g;
	$class =~ s/\;/\|/g;
	if ( $class ) {
		$classes{$drug}{$class} = 1;
	} else {
		$classes{$drug}{"Unclassified"} = 1;
	}
}
close IN;

my $unclassified = "Unclassified";
my $notdrug = "NA";
my %listclasses;
open ( IN , "<$clusters" ) or die "Cannot open $clusters: $!";
while ( <IN> ) {
	chomp; if ( /Cluster/ ) { next; }
	my @line = split( "\t" , $_ );
 
 #print join( "\t" , @line )."\n";
 	my $id = $line[0];
 	my $drug = $line[1];
	my $gene = $line[1];
	my $AA = $line[2];
	if ( exists $classes{$drug} ) {
		my @class = sort keys %{$classes{$drug}};
		$clusters{$id}{$drug} = join( "\t" , ( @line , join( "|" , @class ) ) );
	} else {
		if ( $AA =~ /^p\./ ) {
			$clusters{$id}{$gene.$AA} = join( "\t" , ( @line , $notdrug ) );
		} else {
			$clusters{$id}{$gene.$AA} = join( "\t" , ( @line , $unclassified ) );
		}
	}
}
close IN;

open ( OUT , ">$output" );
foreach my $id ( keys %clusters ) {
	foreach my $spec ( keys %{$clusters{$id}} ) {
		print OUT $clusters{$id}{$spec}."\n";
	}
}
close OUT;

#open ( OUT , ">$output.table" );
#print OUT "Drug\tClass\n";
#foreach my $drug ( keys %classes ) {
#	foreach my $class ( keys %{$classes{$drug}} ) {
#		print OUT join( "\t" , ( $drug , $class ) )."\n";
#	}
#}
#close OUT;
