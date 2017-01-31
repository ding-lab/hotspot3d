#!/bin/perl
#10 February 2015 - Adam D Scott - 
 
use strict;
use warnings;
 
my $usage = 'perl determine_transcript_lengths.pl <gtf>  
';
 
die $usage , unless @ARGV == 1;
my ( $gtf ) = @ARGV;
 
my %transcripts;

open ( IN , "<$gtf" ) or die "Cannot open $gtf: $!";
while ( <IN> )
{
	chomp;
	my @line = split( "\t" , $_ );
 
 	if ( $line[2] eq "CDS" )
	{
		my $length = $line[4] - $line[3];
		my $gene = $line[-1];
		$gene =~ s/.*gene_name \"([\w\\\$.-]+)\"; gene_biotype.*/$1/;
		my $trans = $line[-1];
		$trans =~ s/.*transcript_id \"(ENST\d+)\"; exon_number.*/$1/;
		my $exon = $line[-1];
		$exon =~ s/.*exon_number \"(\d+)\"\; gene_name.*/$1/;
		#print $gene."\t".$trans."\t".$exon."\t".$length."\n";
		$transcripts{$gene}{$trans}{$exon} = $length+1;
	}
}
close IN;

foreach my $gene ( sort keys %transcripts )
{
	foreach my $trans ( sort keys %{$transcripts{$gene}} )
	{
		my $length = 0;
		foreach my $exon ( keys %{$transcripts{$gene}{$trans}} )
		{
			$length += $transcripts{$gene}{$trans}{$exon};
		}
		print $gene."\t".$trans."\t".$length."\n";
	}
}
