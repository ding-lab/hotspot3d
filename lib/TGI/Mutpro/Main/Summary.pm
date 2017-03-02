package TGI::Mutpro::Main::Summary;
#
#----------------------------------
# $Authors: Adam Scott
# $Date: 2015-05-21
# $Revision:  $
# $URL: $
# $Doc: $ summarize clusters
# 
#----------------------------------
#
use strict;
use warnings;

use Carp;
use Cwd;
use Getopt::Long;

use List::MoreUtils qw( uniq );
use List::Util qw( min max );

use IO::File;
use FileHandle;

sub new {
    my $class = shift;
    my $this = {};
    $this->{'clusters_file'} = '3D_Proximity.pairwise.clusters';
    $this->{'output_prefix'} = undef;
    $this->{mutationmass} = {};
    $this->{recurrencemass} = {};
    $this->{drugmass} = {};
    $this->{degrees} = {};
    $this->{centralities} = {};
    $this->{geodesics} = {};
    $this->{centroids} = {};
    $this->{genes} = {};
    $this->{aas} = {};
    $this->{transcripts} = {};
    bless $this, $class;
    $this->process();
    return $this;
}

sub process {
    my $this = shift;
	$this->setOptions();
	$this->readClustersFile();
	$this->writeSummary();
}

sub setOptions {
	my ( $this ) = shift;
    my ( $help, $options );
    unless( @ARGV ) { die $this->help_text(); }
    $options = GetOptions (
        'output-prefix=s' => \$this->{'output_prefix'},
        'clusters-file=s' => \$this->{'clusters_file'},
        'help' => \$help,
    );
    if ( $help ) { print STDERR help_text(); exit 0; }
    unless( $options ) { die $this->help_text(); }
    unless( $this->{'clusters_file'} ) { warn 'You must provide a clusters file! ', "\n"; die $this->help_text(); }
    unless( -e $this->{'clusters_file'} ) { warn "The input clusters file (".$this->{'clusters_file'}.") does not exist! ", "\n"; die $this->help_text(); }
	return;
}

sub generateOutputFileName {
	my ( $this ) = @_;
	my $outFilename = "";
	if ( defined $this->{'output_prefix'} ) {
		$outFilename = $this->{'output_prefix'};
	} else {
		$outFilename = $this->{'clusters_file'};
	}
	$outFilename .= ".summary";
	return $outFilename;
}

sub readClustersFile {
	my ( $this ) = @_;
	my $infh = new FileHandle;
	unless( $infh->open( $this->{'clusters_file'} , "r" ) ) { die "Could not open clusters file $! \n" };
	my @cols;
	while ( my $line = <$infh> ) {
		chomp( $line );
		if ( $line =~ /^Cluster/ ) {
			my $i = 0;
			my %cols = map{ ( $_ , $i++ ) } split( /\t/ , $line );
			unless( defined( $cols{"Cluster"} )
				and defined( $cols{"Gene/Drug"} )
				and defined( $cols{"Mutation/Gene"} )
				and defined( $cols{"Degree_Connectivity"} )
				and defined( $cols{"Closeness_Centrality"} )
				and defined( $cols{"Geodesic_From_Centroid"} )
				and defined( $cols{"Recurrence"} || $cols{"Weight"} ) 
				and defined( $cols{"Chromosome"} ) 
				and defined( $cols{"Start"} ) 
				and defined( $cols{"Reference"} ) 
				and defined( $cols{"Alternate"} ) 
				and defined( $cols{"Transcript"} ) ) {
				die "Not a valid clusters file!\n";
			}
			@cols = ( 	$cols{"Cluster"} , 
						$cols{"Gene/Drug"} , 
						$cols{"Mutation/Gene"} , 
						$cols{"Degree_Connectivity"} , 
						$cols{"Closeness_Centrality"} , 
						$cols{"Geodesic_From_Centroid"} , 
						($cols{"Recurrence"} || $cols{"Weight"}) ,
						$cols{"Chromosome"} , 
						$cols{"Start"} ,
						$cols{"Stop"} ,
						$cols{"Reference"} ,
						$cols{"Alternate"} ,
						$cols{"Transcript"} 
						 );
		} else {
			my ( $id , $genedrug , $aagene , $degree , $centrality , $geodesic , $recurrence, $chromosome, $start, $stop, $ref, $alt, $trans, $alt_trans ) = (split( "\t" , $line ))[@cols];
			$this->sum( 'degrees' , $id , $degree );
			$this->sum( 'centralities' , $id , $centrality );
			$this->sum( 'geodesics' , $id , $geodesic );
			$this->{genes}->{$id}->{$genedrug} = 1;
			$this->{transcripts}->{$id}->{$trans} = 1;
                        if (! exists $this->{aas}->{$id}){
				$this->{aas}->{$id} = ();
                        } 
			my @list;
			if ( $aagene =~ /p\./ ) {
				push @{$this->{aas}->{$id}}, $aagene; 
				if ( $geodesic == 0 ) { $this->{centroids}->{$id} = $genedrug.":".$aagene; }
				$this->sum( 'mutationmass' , $id , 1 );
				$this->sum( 'recurrencemass' , $id , $recurrence );
				#if ( $families ) { &sumlist( $family , $id , $families , $recurrence ); }
				#if ( $doms ) { &sumlist( $domains , $id , $doms , $recurrence ); }
				#$this->list( 'family' , $id , $families );
				#$this->list( 'domains' , $id , $doms );
				#if ( $dmkb !~ /NULL/ ) {
				#	&sum( $DMKB , $id , $recurrence );
				#	$DMKBs->{$id}->{$genedrug.":".$aagene} = 1;
				#}
			} else {
				if ( $geodesic == 0 ) { $this->{centroids}->{$id} = $genedrug; }
				$this->sum( 'drugmass' , $id , 1 );
				#if ( $classes_DrugBank ) { &sumlist( $DrugBank , $id , $classes_DrugBank , $recurrence ); }
				#if ( $classes_NIH ) { &sumlist( $NIH , $id , $classes_NIH , $recurrence ); }
				#&list( $DrugBank , $id , $classes_DrugBank );
				#&list( $NIH , $id , $classes_NIH );
			}
		}
	}
	$infh->close();
	return;
}

sub writeSummary {
	my ( $this ) = @_;
	my $outFilename = $this->generateOutputFileName();
	my $fh = new FileHandle;
	unless( $fh->open( $outFilename , "w" ) ) { die "Could not open $outFilename $! \n"; }
	my $fill = "%.3f"."\t";
	$fh->print( "Cluster_ID\tCentroid\tAvg_Degree\tCentrality\tAvg_Centrality\tAvg_Geodesic\tRecurrence_Mass\tAvg_Recurrence" );
	$fh->print( "\tMutations_(Unique_AAchanges)" );
	#$fh->print( "\tKnown_Mutations_(Unique_Known)" );
	$fh->print( "\tTotal_Drugs\tGenes_Drugs\tAA_Mutations\tTranscripts");
	#$fh->print( "\tHGNC_Gene_Families\tPfam_Domains\tDrugBank_Classes\tNIH_Classes" );
	$fh->print( "\n" );
	foreach my $id ( sort { $a <=> $b } keys %{$this->{mutationmass}} ) {
		$fh->print( $id."\t" ); #Cluster_ID
		if ( exists $this->{centroids}->{$id} ) {
            $fh->print( $this->{centroids}->{$id}."\t" ); #Centroid
        } else {
            print STDERR $id." has no centroid\n";
			$fh->print( "NULL\t" );
        }
		$fh->printf( $fill , $this->avg( 'degrees' , $id , 'mutationmass' ) ); #AVG_Degree (pairs)
		$fh->printf( $fill , $this->{centralities}->{$id} ); #Centrality (cluster closeness)
		$fh->printf( $fill , $this->avg( 'centralities' , $id , 'recurrencemass' ) ); #Avg_Frequency (average recurrence)
		$fh->printf( $fill , $this->avg( 'geodesics' , $id , 'mutationmass' ) ); #Avg_Geodesic (average geodesic from centroid)
		$fh->printf( $fill , $this->{recurrencemass}->{$id} ); #Recurrence_Mass (sum recurrence in cluster)
		$fh->printf( $fill , $this->avg( 'recurrencemass' , $id , 'mutationmass' ) ); #Avg_Frequency (average recurrence)
		$fh->print( $this->{recurrencemass}->{$id}." (".$this->{mutationmass}->{$id}.")\t" ); #Mutations_(Unique_AAchanges)
		#$fh->print( $DMKB->{$id}." (".(scalar keys %{$DMKBs->{$id}}).")\t" ); #Known_Mutations_(Unique_Known)#known druggable
		if ( exists $this->{drugmass}->{$id} ) {
			$fh->print( $this->{drugmass}->{$id}."\t" );
		} else {
			$fh->print( "NA\t" );
		}
		if ( exists $this->{genes}->{$id} && $this->{genes}->{$id} ne "" ) {
			$fh->print( join( "; " , sort keys %{$this->{genes}->{$id}} )."\t" );
		} else {
			$fh->print( "NA\t" );
		} 
		if ( exists $this->{aas}->{$id} && $this->{aas}->{$id} ne ""){
			$fh->print( join("," , uniq(@{$this->{aas}->{$id}}))."\t");
		} else {
			$fh->print( "NA\t" );
		}
		if ( exists $this->{transcripts}->{$id} && $this->{transcripts}->{$id} ne ""){
			$fh->print( join( "; " , sort keys %{$this->{transcripts}->{$id}} )."\t" );
		} else {
			$fh->print( "NA\t" );
		} 

		 #list of genes & drugs
		#if ( exists $family->{$id} && $family->{$id} ne "" ) {
		#	$fh->print( join( "; " , sort keys %{$family->{$id}} )."\t" );
		#} else {
		#	$fh->print( "NA\t" );
		#} #list of gene families
		#if ( exists $domains->{$id} && $domains->{$id} ne "" ) {
		#	$fh->print( join( "; " , sort keys %{$domains->{$id}} )."\t" );
		#} else {
		#	$fh->print( "NA\t" );
		#} #list of gene families
		#if ( exists $DrugBank->{$id} && $DrugBank->{$id} ne "" ) {
		#	$fh->print( join( "; " , sort keys %{$DrugBank->{$id}} )."\t" );
		#} else {
		#	$fh->print( "NA\t" );
		#} #list of drug classes
		#if ( exists $NIH->{$id} && $NIH->{$id} ne "" ) {
		#	$fh->print( join( "; " , sort keys %{$NIH->{$id}} )."\n" );
		#} else {
		#	$fh->print( "NA\n" );
		#} #list of drug classes
		$fh->print( "\n" );
	}
	$fh->close();
	return;
}

sub sum {
	my ( $this , $measure , $id , $sample ) = @_;
	if ( exists $this->{$measure}->{$id} ) {
		$this->{$measure}->{$id} += $sample;
	} else {
		$this->{$measure}->{$id} = $sample;
	}
	return 1;
}

sub sum2 {
	my ( $this , $measure , $id , $key , $sample ) = @_;
	if ( exists $this->{$measure}->{$id}->{$key} ) {
		$this->{$measure}->{$id}->{$key} += $sample;
	} else {
		$this->{$measure}->{$id}->{$key} = $sample;
	}
	return 1;
}

sub list {
	my ( $this , $measure , $id , $thing ) = @_;
	my @list;
	if ( $thing =~ /\|/ ) {
		@list = split( /\|/ , $thing );
	} else {
		@list = ( $thing );
	}
	foreach my $type ( @list ) {
		$this->{$measure}->{$id}->{$type} = 1;
	}
	return 1;
}

sub avg {
	my ( $this , $measure , $id , $N ) = @_;
	if ( exists $this->{$N}->{$id} && exists $this->{$measure}->{$id} ) {
		my $quo = $this->{$measure}->{$id}/$this->{$N}->{$id};
		return $quo;
	}
	return 0;
}

sub help_text{
	my $this = shift;
	return <<HELP

Usage: hotspot3d summary [options]

                             REQUIRED
--clusters-file              Clusters file

                             OPTIONAL
--output-prefix              Output prefix

--help                       this message

HELP

}

1;
