package TGI::Mutpro::Main::Summary;
#
#----------------------------------
# $Authors: Adam Scott and Matthew Bailey
# $Date: 2015-05-21
# $Revision: 0.3 $
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

#use List::MoreUtils qw( uniq );
use List::Util qw( min max );

use IO::File;
use FileHandle;

my $PTM = "ptm";
my $NULL = "NA";
my $fill = "%.3f";
my $count = "%i";

sub new {
    my $class = shift;
    my $this = {};
    $this->{'clusters_file'} = undef;
    $this->{isrecurrence} = 0;
    $this->{'output_prefix'} = undef;
    $this->{mutationmass} = {};
    $this->{recurrencemass} = {};
    $this->{drugmass} = {};
    $this->{vertexmass} = {};
    $this->{degrees} = {};
    $this->{centralities} = {};
    $this->{geodesics} = {};
    $this->{centroids} = {};
    $this->{genes} = {};
    $this->{genomicmutations} = {};
    $this->{aas} = {};
    $this->{sites} = {};
    $this->{residues} = {};
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
	unless( $infh->open( $this->{'clusters_file'} , "r" ) ) {
		die "Could not open clusters file $! \n";
	}
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
			if ( $line =~ /Recurrence/ ) {
				$this->{isrecurrence} = 1;
			}
		} else {
			my @line = split( "\t" , $line );
			my ( $id , $genedrug , $aagene , $degree , $centrality , 
				 $geodesic , $recurrence, $chromosome, $start, $stop, 
				 $ref, $alt, $trans, $alt_trans 
			) = @line[@cols];

			$this->sum( 'degrees' , $id , $degree );
			$this->sum( 'centralities' , $id , $centrality );
			$this->sum( 'geodesics' , $id , $geodesic );
			$this->{transcripts}->{$id}->{$trans} += 1;
			my @list;
			if ( $aagene =~ /p\./ ) {
				$this->{genes}->{$id}->{$genedrug} += 1;
				my $mutationKey = $this->makeKey( $genedrug , 
						$chromosome , $start , $stop , $ref , $alt );
				my $proteinKey = $this->makeKey( $genedrug , $aagene );
				if ( $mutationKey =~ m/$PTM/ ) {
					$this->{sites}->{$id}->{$proteinKey} += 1;
				} else {
					if ( $this->{isrecurrence} ) {
						$this->{genomicmutations}->{$id}->{$mutationKey} += $recurrence;
						$this->{aas}->{$id}->{$proteinKey} += $recurrence;
						$this->sum( 'recurrencemass' , $id , $recurrence );
					} else {
						$this->{genomicmutations}->{$id}->{$mutationKey} += 1;
						$this->{aas}->{$id}->{$proteinKey} += 1;
						$this->sum( 'recurrencemass' , $id , 1 );
					}
					$this->sum( 'mutationmass' , $id , 1 );
				}
				my $position = $1 if $aagene =~ m/p\.\D*(\d+).*/;
				$this->{residues}->{$id}->{$position} += 1;
				if ( $geodesic == 0 ) { $this->{centroids}->{$id} = $genedrug.":".$aagene; }
				if ( $this->{isrecurrence} ) {
					$this->sum( 'vertexmass' , $id , $recurrence );
				} else {
					$this->sum( 'vertexmass' , $id , 1 );
				}
			} else {
				$this->{drugs}->{$id}->{$genedrug} += 1;
				if ( $geodesic == 0 ) { $this->{centroids}->{$id} = $genedrug; }
				$this->sum( 'drugmass' , $id , 1 );
				$this->sum( 'vertexmass' , $id , 1 );
			}
		}
	}
	$infh->close();
	return;
}

sub makeKey {
	my $this = shift;
	return &combine( @_ );
}

sub combine {
	return join( ":" , @_ );
}

sub writeSummary {
	my ( $this ) = @_;
	my $outFilename = $this->generateOutputFileName();
	my $fh = new FileHandle;
	unless( $fh->open( $outFilename , "w" ) ) { die "Could not open $outFilename $! \n"; }
	$this->printHeader( $fh );
	foreach my $id ( sort { $a cmp $b } keys %{$this->{vertexmass}} ) {
		$this->printClusterID( $fh , $id );
		$this->printCentroid( $fh , $id );
		$this->printCentrality( $fh , $id );
		$this->printRecurrenceMass( $fh , $id );
		$this->printAvgCentrality( $fh , $id );
		$this->printAvgRecurrenceMass( $fh , $id );
		$this->printAvgDegree( $fh , $id );
		$this->printAvgGeodesic( $fh , $id );
		$this->printNGenes( $fh , $id );
		$this->printNVertices( $fh , $id );
		$this->printNGenomicMutations( $fh , $id );
		$this->printNProteinMutations( $fh , $id );
		$this->printNProteinSites( $fh , $id );
		$this->printNProteinPositions( $fh , $id );
		$this->printNDrugs( $fh , $id );
		$this->printGenes( $fh , $id );
		$this->printTranscripts( $fh , $id );
		$this->printGenomicMutations( $fh , $id );
		$this->printProteinMutations( $fh , $id );
		$this->printProteinSites( $fh , $id );
		$this->printProteinPositions( $fh , $id );
		$this->printDrugs( $fh , $id );

		$fh->print( "\n" );
	}
	$fh->close();
	return;
}

sub printHeader {
	my ( $this , $fh ) = @_;
	$fh->print( join( "\t" , ( 	"Cluster_ID" , 
								"Centroid" , 
								"Cluster_Closeness" , 
								"Recurrence_Mass" , 
								"Avg_Centrality" , 
								"Avg_Recurrence" , 
								"Avg_Degree" , 
								"Avg_Geodesic" , 
								"N_Genes" , 
								"N_Vertices" , 
								"N_Genomic_Mutations" , 
								"N_Protein_Mutations" , 
								"N_Protein_Sites" , 
								"N_Protein_Positions" , 
								"N_Drugs" , 
								"Genes" , 
								"Transcripts" ,
								"Genomic_Mutations" , 
								"Protein_Mutations" , 
								"Protein_Sites" , 
								"Protein_Positions" , 
								"Drugs" 
							)
					)
	);
	$fh->print( "\n" );
	return;
}


sub printClusterID {
	my ( $this , $fh , $id ) = @_;
	$fh->print( $id ); #Cluster_ID
	$fh->print( "\t" );
	return;
}


sub printCentroid {
	my ( $this , $fh , $id ) = @_;
	if ( exists $this->{centroids}->{$id} ) {
		$fh->print( $this->{centroids}->{$id} ); #Centroid
	} else {
		print STDERR $id." has no centroid\n";
		$fh->print( $NULL );
	}
	$fh->print( "\t" );
	return;
}


sub printAvgDegree {
	my ( $this , $fh , $id ) = @_;
	if ( exists $this->{degrees}->{$id} ) {
		$fh->printf( $fill , $this->avg( 'degrees' , $id , 'vertexmass' ) ); #AVG_Degree (pairs)
	} else {
		$fh->print( $NULL );
	}
	$fh->print( "\t" );
	return;
}

sub printCentrality {
	my ( $this , $fh , $id ) = @_;
	if ( exists $this->{centralities}->{$id} ) {
		$fh->printf( $fill , $this->{centralities}->{$id} ); #Centrality (cluster closeness)
	} else {
		$fh->print( $NULL );
	}
	$fh->print( "\t" );
	return;
}

sub printAvgCentrality {
	my ( $this , $fh , $id ) = @_;
	if ( exists $this->{centralities}->{$id} ) {
		$fh->printf( $fill , $this->avg( 'centralities' , $id , 'vertexmass' ) ); #Avg_Frequency (average recurrence)
	} else {
		$fh->print( $NULL );
	}
	$fh->print( "\t" );
	return;
}

sub printAvgGeodesic {
	my ( $this , $fh , $id ) = @_;
	if ( exists $this->{geodesics}->{$id} ) {
		$fh->printf( $fill , $this->avg( 'geodesics' , $id , 'vertexmass' ) ); #Avg_Geodesic (average geodesic from centroid)
	} else {
		$fh->print( $NULL );
	}
	$fh->print( "\t" );
	return;
}

sub printRecurrenceMass {
	my ( $this , $fh , $id ) = @_;
	if ( exists $this->{recurrencemass}->{$id} ) {
		$fh->printf( $count , $this->{recurrencemass}->{$id} ); #Recurrence_Mass (sum recurrence in cluster)
	} else {
		$fh->print( $NULL );
	}
	$fh->print( "\t" );
	return;
}

sub printAvgRecurrenceMass {
	my ( $this , $fh , $id ) = @_;
	if ( exists $this->{recurrencemass}->{$id} ) {
		$fh->printf( $fill , $this->avg( 'recurrencemass' , $id , 'mutationmass' ) ); #Avg_Frequency (average recurrence)
	} else {
		$fh->print( $NULL );
	}
	$fh->print( "\t" );
	return;
}

sub printNGenomicMutations {
	my ( $this , $fh , $id ) = @_;
	if ( exists $this->{genomicmutations}->{$id} ) {
		$fh->printf( $count , ( scalar keys %{$this->{genomicmutations}->{$id}} ) ); #Total_Mutations
	} else {
		$fh->print( $NULL );
	}
	$fh->print( "\t" );
	return;
}

sub printNProteinMutations {
	my ( $this , $fh , $id ) = @_;
	if ( exists $this->{aas}->{$id} ) {
		$fh->printf( $count , ( scalar keys %{$this->{aas}->{$id}} ) ); #Total_Mutations
	} else {
		$fh->print( $NULL );
	}
	$fh->print( "\t" );
	return;
}

sub printNProteinSites {
	my ( $this , $fh , $id ) = @_;
	if ( exists $this->{sites}->{$id} ) {
		$fh->printf( $count , ( scalar keys %{$this->{sites}->{$id}} ) ); #Total_Mutations
	} else {
		$fh->print( $NULL );
	}
	$fh->print( "\t" );
	return;
}

sub printNProteinPositions {
	my ( $this , $fh , $id ) = @_;
	if ( exists $this->{residues}->{$id} ) {
		$fh->printf( $count , ( scalar keys %{$this->{residues}->{$id}} ) ); #Total_Mutations
	} else {
		$fh->print( $NULL );
	}
	$fh->print( "\t" );
	return;
}

sub printNGenes {
	my ( $this , $fh , $id ) = @_;
	if ( exists $this->{genes}->{$id} ) {
		$fh->print( ( scalar keys %{$this->{genes}->{$id}} ) );
	} else {
		$fh->print( $NULL );
	}
	$fh->print( "\t" );
	return;
}

sub printNDrugs {
	my ( $this , $fh , $id ) = @_;
	if ( exists $this->{drugs}->{$id} ) {
		$fh->print( ( scalar keys %{$this->{drugs}->{$id}} ) );
	} else {
		$fh->print( $NULL );
	}
	$fh->print( "\t" );
	return;
}

sub printNVertices {
	my ( $this , $fh , $id ) = @_;
	my ( $nDrugs , $nGenomicMutations , $nSites ) = (0) x 3;
	#if ( exists $this->{drugs}->{$id} ) {
	#	$nDrugs = scalar keys %{$this->{drugs}->{$id}};
	#}
	#if ( exists $this->{genomicmutations}->{$id} ) {
	#	$nGenomicMutations = scalar @{$this->{genomicmutations}->{$id}};
	#}
	#if ( exists $this->{sites}->{$id} ) {
	#	$nSites = scalar @{$this->{sites}->{$id}};
	#}
	#my $nVertices = $nDrugs + $nGenomicMutations + $nSites;
	if ( exists $this->{vertexmass}->{$id} ) {
		my $nVertices = $this->{vertexmass}->{$id};
		$fh->print( $nVertices );
	} else {
		$fh->print( $NULL );
	}
	$fh->print( "\t" );
	return;
}

sub printGenomicMutations {
	my ( $this , $fh , $id ) = @_;

	if ( exists $this->{genomicmutations}->{$id} ) {
		$fh->print( join( "," , sort keys %{$this->{genomicmutations}->{$id}} ) );
	} else {
		$fh->print( $NULL );
	}
	$fh->print( "\t" );
	return;
}

sub printProteinMutations {
	my ( $this , $fh , $id ) = @_;
	if ( exists $this->{aas}->{$id} ) {
		$fh->print( join( "," , sort keys %{$this->{aas}->{$id}} ) );
	} else {
		$fh->print( $NULL );
	}
	$fh->print( "\t" );
	return;
}

sub printProteinPositions {
	my ( $this , $fh , $id ) = @_;
	if ( exists $this->{residues}->{$id} ) {
		$fh->print( join( "," , sort keys %{$this->{residues}->{$id}} ) );
	} else {
		$fh->print( $NULL );
	}
	$fh->print( "\t" );
	return;
}

sub printProteinSites {
	my ( $this , $fh , $id ) = @_;
	if ( exists $this->{sites}->{$id} ) {
		$fh->print( join( "," , sort keys %{$this->{sites}->{$id}} ) );
	} else {
		$fh->print( $NULL );
	}
	$fh->print( "\t" );
	return;
}

sub printGenes {
	my ( $this , $fh , $id ) = @_;
	if ( exists $this->{genes}->{$id} && $this->{genes}->{$id} ne "" ) {
		my $out = "";
		foreach my $gene ( sort keys %{$this->{genes}->{$id}} ) {
			$out .= $gene."(".$this->{genes}->{$id}->{$gene}."),";
		}
		$out = substr( $out , 0 , -1 );
		$fh->print( $out );
	} else {
		$fh->print( $NULL );
	} 
	$fh->print( "\t" );
	return;
}

sub printTranscripts {
	my ( $this , $fh , $id ) = @_;
	if ( exists $this->{transcripts}->{$id} && $this->{transcripts}->{$id} ne ""){
		my $out = "";
		foreach my $transcript ( sort keys %{$this->{transcripts}->{$id}} ) {
			$out .= $transcript."(".$this->{transcripts}->{$id}->{$transcript}."),";
		}
		$out = substr( $out , 0 , -1 );
		$fh->print( $out );
	} else {
		$fh->print( $NULL );
	} 
	$fh->print( "\t" );
	return;
}

sub printDrugs {
	my ( $this , $fh , $id ) = @_;
	if ( exists $this->{drugs}->{$id} && $this->{drugs}->{$id} ne ""){
		$fh->print( join( "," , sort keys %{$this->{drugs}->{$id}} ) );
	} else {
		$fh->print( $NULL );
	}
	$fh->print( "\t" );
	return;
}


sub count {
	my ( $this , $type , $id ) = @_;
	return scalar @{$this->{$type}->{$id}};
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
