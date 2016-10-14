#!/bin/perl
#16 January 2015 - Adam D Scott - 
 
use strict;
use warnings;
 
my $usage = 'perl collapse_transcripts.pl <MAF> <clusters> <gtl> <output>
';
die $usage , unless @ARGV == 4;
my ( $maf , $clusters , $gtl , $output ) = @ARGV;
my $barcode = 16;
my $chrcol = 33;
my $startcol = 34;
my $refcol = 36;
my $varcol = 37;
my $genecol = 39;
my $transcol = 40;
my $typecol = 46;
my $cposcol = 47;
my $AAcol = 48;

my %mutations;
my %alternates;
my %transcripts;

open ( IN0 , "<$gtl" ) or die "Cannot open $gtl: $!";
while ( <IN0> ) {
	chomp;
	my @line = split( "\t" , $_ );
	my $gene = $line[0];
	my $trans = $line[1];
	my $length = $line[2];

	$transcripts{$gene}{$trans} = $length;
}
close IN0;

my %usetrans;
foreach my $gene ( keys %transcripts ) {
	my @ts = keys %{$transcripts{$gene}}; my $long = $ts[0]; $usetrans{$gene}{$long} = $transcripts{$gene}{$long};
	foreach my $trans ( keys %{$transcripts{$gene}} ) {
		if ( $usetrans{$gene}{$long} < $transcripts{$gene}{$trans} ) {
			$usetrans{$gene}{$trans} = $transcripts{$gene}{$trans};
			delete $usetrans{$gene}{$long};
			$long = $trans;
		}
	}
}

open ( IN , "<$maf" ) or die "Cannot open $maf: $!";
while ( <IN> ) {
	chomp; if ( /Hugo/ ) { next; }
	my @line = split( "\t" , $_ );

	my $chr = $line[$chrcol-1];
	my $start = $line[$startcol-1];
	my $stop = $line[$startcol];
	my $ref = $line[$refcol-1];
	my $var = $line[$varcol-1];
	my $gene = $line[$genecol-1];
	my $trans = $line[$transcol-1];
	my $type = $line[$typecol-1];
	my $cpos = $line[$cposcol-1];
	my $AA = $line[$AAcol-1];
	if ( $type =~ /missense|in_frame/ ) {
		my $cv1 = $line[-2];
		my $cv2 = $line[-1];
		if ( $cv2 !~ /ClinVar/ ) {
			$cv1 = "NULL";
			$cv2 = "NULL";
		}

		my $vari = $chr."\t".$start."\t".$stop."\t".$ref."\t".$var;
		$mutations{$gene}{$AA}{$vari}{$trans}{$cpos} = $cv1."\t".$cv2;
		$alternates{$gene}{$vari}{$trans}{$cpos} = $AA;
	}
}
close IN;

my %clusterlines;
open ( OUT , ">$output" );
print OUT "Cluster_ID\tGene\tAAchange\tDegree\tCloseness_Centrality\tGeodesic\tFrequency\tTranscript\tc_position\tChromosome\tStart\tStop\tReference\tVariant\tClinVarAnnotation\tClinVarCitation\tLongest_Transcript\tLongest_AAchange\tLongest_c_position\n";
open ( IN2 , "<$clusters" ) or die "Cannot open $clusters: $!";
while ( <IN2> ) {
	chomp; if ( /Cluster/ ) { next; }
	my $line = $_;
	my @line = split( "\t" , $_ );

	my ( $id , $gene , $AA , $deg , $Cc , $geo , $freq ) = @line;

	my $vari = "NULL";
	my $trans = "NULL";
	my $cpos = "NULL";
	my $others = "";
	my $val = "NULL";
	my @l;
	my $long = "NULL";
	my @othercpos;
	my $othercpos = "NULL";
	my $otherAA = "NULL";
	if ( exists $mutations{$gene}{$AA} ) {
		foreach $vari ( keys %{$mutations{$gene}{$AA}} ) {
			@l = keys %{$usetrans{$gene}}; 
			$long = $l[0];
			foreach $trans ( sort keys %{$mutations{$gene}{$AA}{$vari}} ) {
				foreach $cpos ( sort keys %{$mutations{$gene}{$AA}{$vari}{$trans}} ) {
					@othercpos = keys %{$alternates{$gene}{$vari}{$long}}; 
					$othercpos = $othercpos[0];
					$otherAA = $alternates{$gene}{$vari}{$long}{$othercpos};
					#$line[2] = $otherAA;
					$val = $mutations{$gene}{$AA}{$vari}{$trans}{$cpos};
					$clusterlines{$id}{$gene}{$vari}{$long} = join( "\t" , ( @line , $trans , $cpos , $vari , $val , $long , $otherAA , $othercpos ) );#$transcripts{$gene}{$trans} ) );
				}
			}
		}
	} else {
		print "This did not match a MAF mutation: ".join( "\t" , @line ). "\n";
	}
	print OUT join( "\t" , ( @line , $trans , $cpos , $vari , $val , $long , $otherAA , $othercpos ) );
}
close IN2;

#foreach my $id ( keys %clusterlines )
#{
#	foreach my $gene ( sort keys %{$clusterlines{$id}} )
#	{
#		foreach my $AA ( sort keys %{$clusterlines{$id}{$gene}} )
#		{
#			foreach my $trans ( sort keys %{$usetrans{$gene}} )
#			{
#				print OUT $clusterlines{$id}{$gene}{$AA}{$trans}."\n";
#			}
#		}
#	}
#}
close OUT;
