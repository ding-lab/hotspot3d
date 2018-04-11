package TGI::Mutpro::Main::Network;
#
#----------------------------------
# $Authors: Adam Scott, Sohini Sengupta, & Amila Weerasinghe
# $Date: 2014-01-14 14:34:50 -0500 (Tue Jan 14 14:34:50 CST 2014) $
# $Revision: 8.0 $
# $URL: $
# $Doc: $ determine mutation clusters from HotSpot3D inter, intra, and druggable data
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
use File::Basename;
use File::Temp;
use Parallel::ForkManager;

use Data::Dumper;

use TGI::Variant;
use TGI::ProteinVariant;
use TGI::Mutpro::Main::Cluster;

our @ISA = qw( TGI::Mutpro::Main::Cluster );

my $WEIGHT = "weight";
my $RECURRENCE = "recurrence";
my $UNIQUE = "unique";
my $SITE = "site";
my $PVALUEDEFAULT = 0.05;
my $DISTANCEDEFAULT = 10;
my $MAXDISTANCE = 10000;
my $AVERAGEDISTANCE = "average";
my $SHORTESTDISTANCE = "shortest";
my $NETWORK = "network";
my $DENSITY = "density";
my $INDEPENDENT = 0; #"independent";
my $DEPENDENT = 1; #"dependent";
my $ANY = "any";
my $MINCLUSTERSIZE = 2;
my $NULL = "-";

my $CENTRALITY = "centrality";
my $EXPONENTIALS = "exponentials";
my $MAXWEIGHT = 20;
my $MAXGEODESIC = 10;

my $MULTIMER = "multimer";
my $MONOMER = "monomer";
my $HOMOMER = "homomer";
my $HETEROMER = "heteromer";
my $UNSPECIFIED = "unspecified";
my $INTRA = "intra";
my $INTER = "inter";

my $BSUB = "bsub";
my $LOCAL = "local";
my $NONE = "none";

sub new {
    my $class = shift;
    my $this = shift;
    bless $this, $class;
    $this->process();
    return $this;
}

sub process {
    my $this = shift;
	my $temp_distance_matrix = {};
 	my $temp_mutations = {};
	my $distance_matrix = {};
 	my $mutations = {};
 	my $siteVertexMap = {}; # a hash to store refAlt and pKey info when "SITE" specific clustering is done.

	$this->readMAF( $temp_mutations );
	#$this->readSite( $temp_mutations );
	$this->getDrugMutationPairs( $temp_mutations , $temp_distance_matrix );
	$this->getMutationSitePairs( $temp_mutations , $temp_distance_matrix );
	$this->getSiteSitePairs( $temp_mutations , $temp_distance_matrix );
	$this->getMutationMutationPairs( $temp_distance_matrix );
	$this->vertexFilter( $temp_mutations , $temp_distance_matrix , $mutations , $distance_matrix, $siteVertexMap );
	#$this->initializeSameSiteDistancesToZero( $distance_matrix );
	$this->{"siteVertexMap"} = $siteVertexMap; # store the reference to siteVertexMap
	$this->networkClustering( $mutations , $distance_matrix );
	#print Dumper $distance_matrix;

    return 1;
}

#####
#	functions
#####
sub networkClustering {
	#$this->finalize( $clusterings , $mutations , $distance_matrix );
	my ( $this , $mutations , $distance_matrix ) = @_;
	print STDOUT "HotSpot3D::Cluster::networkClustering\n";
	my $clusterings = {};
	$this->link( $clusterings , $distance_matrix );
	if ( $this->{'parallel'} eq $LOCAL ) {
		$this->localParallelNetworkClustering( $clusterings , $mutations , $distance_matrix );
	#} elsif ( $this->{'parallel'} eq $BSUB ) {
	#	$this->bsubParallelNetworkClustering( $clusterings , $mutations , $distance_matrix );
	} else {
		$this->noParallelNetworkClustering( $clusterings , $mutations , $distance_matrix );
	}
	return;
}

sub localParallelNetworkClustering {
	my ( $this , $clusterings , $mutations , $distance_matrix ) = @_;
	print STDOUT "HotSpot3D::Cluster::localParallelNetworkClustering\n";
	print STDOUT "\tparallel clustering over structures with up to ".$this->{'max_processes'}." processes\n";
	my $fh = $this->writeHeader();
	my $tempDir = File::Temp->newdir( TEMPLATE => 'hs3dXXXXX' );
	my $pm = Parallel::ForkManager->new( $this->{'max_processes'} , $tempDir );
	my %finalLines;
	$pm->run_on_finish( 
		sub {
			my ( $pid , $exit_code , $ident , $exit_signal , $core_dump , $data ) = @_;
			if ( defined $data ) {
				if ( exists $clusterings->{$data->[0]} ) {
					foreach my $superClusterID ( sort keys %{$data->[1]} ) {
						print STDOUT $superClusterID;
						$this->writeClustersFile( $fh , $data->[1]->{$superClusterID} );
					}
				}
			} else {
				print qq|No data from child $pid!\n|;
			}
		}
	);
	DATA_LOOP:
	foreach my $structure ( sort keys %{$distance_matrix} ) {
		my $pid = $pm->start and next DATA_LOOP;
		#my $lines = "";
		my $superLines = {};
		foreach my $superClusterID ( sort keys %{$clusterings->{$structure}} ) {
			my $lines = {};
			$this->determineStructureClusters( $clusterings , $mutations , $distance_matrix , 
					$structure , $superClusterID , $lines );
			$superLines->{$superClusterID} = $lines;
		}
		print STDOUT "children: ".(scalar keys %{$superLines})."\n";
		my @pair = ( $structure , $superLines );
		#my @pair = ( $structure , $lines );
		$pm->finish( 0 , \@pair );
	}
	$pm->wait_all_children;
	$fh->close();
	return;
}

#TODO finish bsub parallelization
#sub bsubParallelNetworkClustering {
#	my ( $this , $clusterings , $mutations , $distance_matrix ) = @_;
#	print STDOUT "bsub parallel clustering over structures\n";
#	foreach my $structure ( sort keys %{$distance_matrix} ) {
#		$this->{'results'}->{$structure} = "";
#		foreach my $superClusterID ( sort keys %{$clusterings->{$structure}} ) {
#			my $subClusterID = 0;
#			$this->determineStructureClusters( $clusterings , $mutations , $distance_matrix , 
#					$structure , $superClusterID , $subClusterID );
#		}
#	}
#	return;
#}
#
sub writeHeader {
	my ( $this ) = shift;
	print STDOUT "HotSpot3D::Cluster::writeHeader\n";
    my $outFilename = $this->generateFilename();
    print STDOUT "\tCreating cluster output file: ".$outFilename."\n";
    my $fh = new FileHandle;
    die "Could not create clustering output file\n" unless( $fh->open( $outFilename , "w" ) );
	my $score = "Closeness_Centrality";
	if ( $this->{'vertex_score'} eq $EXPONENTIALS ) { $score = "Exponentials_Score"; }
    $fh->print( join( "\t" , (  "Cluster" , "Gene/Drug" , "Mutation/Gene" ,
                                "Degree_Connectivity" , $score ,
                                "Geodesic_From_Centroid" , "Weight" ,
                                "Chromosome" , "Start" , "Stop" ,
                                "Reference" , "Alternate" ,
                                "Transcript" , "Alternative_Transcripts"
                             )
                    )."\n"
              );
	return $fh
}

sub noParallelNetworkClustering {
	my ( $this , $clusterings , $mutations , $distance_matrix ) = @_;
	print STDOUT "serially clustering over structures\n";
	my $fh = $this->writeHeader();
	foreach my $structure ( sort keys %{$distance_matrix} ) {
		foreach my $superClusterID ( sort keys %{$clusterings->{$structure}} ) {
			my $lines = {};
			$this->determineStructureClusters( $clusterings , $mutations , $distance_matrix , 
					$structure , $superClusterID , $lines );
			$this->writeClustersFile( $fh , $lines );
		}
	}
	$fh->close();
	return;
}

## NETWORK CLUSTERING - AGGLOMERATIVE HIERARCHICAL CLUSTERING
sub link {
	#$this->link( $clusterings , $distance_matrix );
	my ( $this, $clusterings , $distance_matrix , $mutations ) = @_;
	print STDOUT "HotSpot3D::Cluster::link\n";
	foreach my $structure ( sort keys %{$distance_matrix} ) {
		print "\t".$structure."\n";
		foreach my $mutationKey1 ( sort keys %{$distance_matrix->{$structure}} ) {
			foreach my $mutationKey2 ( sort keys %{$distance_matrix->{$structure}} ) { #->{$mutationKey1}} ) {
				my $distance;
				if ( $this->isSameProteinPosition( $mutations , $mutationKey1 , $mutationKey2 ) == 1 ) {
					$distance = 0;
				} else {
					$distance = $this->getElementByKeys( $distance_matrix , $structure , $mutationKey1 , $mutationKey2 );
					if ( $distance == $MAXDISTANCE ) {
						next;
					}
				}
				#print join( "\t" , ( $mutationKey1 , $mutationKey2 , $distance ) )."\n";
				my @mutations = ( $mutationKey1 , $mutationKey2 );
				my ( $combine1 , $combine2 , $id ); 
				my @combine;
				foreach $id ( keys %{$clusterings->{$structure}} ) { #each cluster
					if ( exists $clusterings->{$structure}->{$id}->{$mutationKey1} ) {
						push @combine , $id;
					}
					if ( exists $clusterings->{$structure}->{$id}->{$mutationKey2} ) {
						push @combine , $id;
					}
				}
				TGI::Mutpro::Main::Cluster::numSort( \@combine );
				if ( scalar @combine > 0 ) { #collapse clusters into one
					my $collapse_to = $combine[0]; #cluster type
					#print "collapsing to (".$collapse_to.") ".join( ", " , @combine )."\n";;
					foreach my $otherClusters ( @combine ) {
						if ( $otherClusters != $collapse_to ) {
							foreach my $mutationKey ( keys %{$clusterings->{$structure}->{$otherClusters}} ) {
								$clusterings->{$structure}->{$collapse_to}->{$mutationKey} = 1;
								delete $clusterings->{$structure}->{$otherClusters}->{$mutationKey};
							}
							delete $clusterings->{$structure}->{$otherClusters};
						}
					}
					$clusterings->{$structure}->{$collapse_to}->{$mutationKey1} = 1;
					$clusterings->{$structure}->{$collapse_to}->{$mutationKey2} = 1;
				} else { #new cluster
					#print "\tnew cluster\n";
					my @ids = keys %{$clusterings->{$structure}};
					if ( scalar @ids > 0 ) {
						TGI::Mutpro::Main::Cluster::numSort( \@ids );
						$id = $ids[-1] + 1;
					} else { $id = 0; }
					$clusterings->{$structure}->{$id}->{$mutationKey1} = 1;
					$clusterings->{$structure}->{$id}->{$mutationKey2} = 1;
				}
			} #foreach mutation2
		} #foreach mutation1
#		my $nsuper = scalar keys %{$clusterings->{$structure}};
#		print "there are ".$nsuper." superclusters in ".$structure."\n";
	} #foreach structure
	#foreach my $s ( sort keys %{$clusterings} ) {
	#	foreach my $id ( sort keys %{$clusterings->{$s}} ) {
	#		foreach my $m ( sort keys %{$clusterings->{$s}->{$id}} ) {
	#			print "clusterings: ".join( "\t" , ( $s , $id , $m ) )."\n";
	#		}
	#	}
	#}
    return;
}

sub initializeGeodesics {
	my ( $this , $clusterings , $superClusterID , $structure , $distance_matrix , $mutations ) = @_;
	#print STDOUT "HotSpot3D::Cluster::initializeGeodesics\n";
	my $geodesics = {};
	my $nInitialized = 0;
	my $nMutations = scalar keys %{$clusterings->{$structure}->{$superClusterID}};
	#print $nMutations." mutations to initialize geodesics\n";
	foreach my $mutationKey1 ( sort keys %{$clusterings->{$structure}->{$superClusterID}} ) { #initialize geodesics
		next if ( $this->hasBeenProcessed( $structure , $mutationKey1 ) );
	#foreach my $mutationKey1 ( sort keys %{$distance_matrix->{$structure}} ) { #initialize geodesics
		#print "need to process mutationKey1 = ".$mutationKey1.": \n";
		foreach my $mutationKey2 ( sort keys %{$clusterings->{$structure}->{$superClusterID}} ) {
			next if ( $this->hasBeenProcessed( $structure , $mutationKey2 ) );
		#foreach my $mutationKey2 ( sort keys %{$distance_matrix->{$structure}} ) {
			#next if ( exists $dist{$mutationKey1}{$mutationKey2} );
			#print "need to process mutationKey2 = ".$mutationKey2.": \n";
			if ( $this->isSameProteinPosition( $mutations , $mutationKey1 , $mutationKey2 ) == 1 ) {
				#print "same site: ".$mutationKey1."\t".$mutationKey2."\n";
				$geodesics->{$structure}->{$mutationKey1}->{$mutationKey2} = 0;
				$geodesics->{$structure}->{$mutationKey2}->{$mutationKey1} = 0;
				#$distance_matrix->{$structure}->{$mutationKey1}->{$mutationKey2} = 0;
				#$distance_matrix->{$structure}->{$mutationKey2}->{$mutationKey1} = 0;
				$nInitialized += 1;
			} else {
				my $distance = $this->getElementByKeys( $distance_matrix , 
										$structure , $mutationKey1 , 
										$mutationKey2 );
				if ( $distance != $MAXDISTANCE ) { #NOTE if amino acids are neighboring, then a distance is given
					$nInitialized += 1;
					$geodesics->{$structure}->{$mutationKey1}->{$mutationKey2} = $distance;
					$geodesics->{$structure}->{$mutationKey2}->{$mutationKey1} = $distance;
					#$distance_matrix->{$structure}->{$mutationKey1}->{$mutationKey2} = $distance;
					#$distance_matrix->{$structure}->{$mutationKey2}->{$mutationKey1} = $distance;
					#print join( "\t" , ( "known" , $mutationKey1 , $mutationKey2 , $distance ) )."\n";
				} #if distance not known
			}
		} #foreach mutationKey2
	} #foreach mutationKey1

#	print "GEODESICS\n";
#	foreach my $mu1 ( sort keys %{$geodesics->{$structure}} ) {
#		foreach my $mu2 ( sort keys %{$geodesics->{$structure}} ) {
#			print "\t".join( "  " , ( $mu1 , $mu2 , $geodesics->{$structure}->{$mu1}->{$mu2} ) )."\n";
#		}
#	}

	#print "\tnInitialized = ".$nInitialized."\n";
	return $geodesics;
	#return $distance_matrix;
}

sub isRadiusOkay {
	my ( $this , $geodesics , $structure , $mutationKey1 , $mutationKey2 ) = @_;
	#print "isRadiusOkay of d(".$mutationKey1.",".$mutationKey2.") = ".$geodesics->{$structure}->{$mutationKey1}->{$mutationKey2}.": ";
	if ( $geodesics->{$structure}->{$mutationKey1}->{$mutationKey2} <= $this->{'max_radius'} ) {
		#print "OKAY\n";
		return 1;
	}
	#print "TOO LONG\n";
	return 0;
}

sub determineCentroid {
	my ( $this , $mutationKey , $newScore , 
		 $currentCentroid , $currentScore ) = @_;
	if ( not $currentScore ) {
		return ( $mutationKey , $newScore );
	}
	if ( abs( $newScore ) > abs( $currentScore ) ) {
		return ( $mutationKey , $newScore );
	}
	return ( $currentCentroid , $currentScore );
}

sub determineSubClusters {
	my ( $this , $scores , $geodesics , $structure ) = @_;
	my $subClusters = {};
	my $moreToFind = 1;
	my $subClusterID = -1;
	while ( $moreToFind == 1 ) {
		my $score = 0;
		my $centroidScore;
		my $centroid;
		$subClusterID++;
		print "Finding subcluster: ".$subClusterID."\n";
		my $checked = 0;
		foreach my $mutationKey ( keys %{$scores} ) {
			next if $this->hasBeenProcessed( $structure , $mutationKey );
			$checked++;
			foreach my $refAlt ( keys %{$scores->{$mutationKey}} ) {
				$score = $scores->{$mutationKey}->{$refAlt};
				( $centroid , $centroidScore ) = $this->determineCentroid( 
						$mutationKey , $score ,	$centroid , $centroidScore );
			}
		}
		my $vertices = 0;
		if ( defined $centroid ) {
			foreach my $mutationKey ( keys %{$geodesics->{$structure}->{$centroid}} ) {
				next if $this->hasBeenProcessed( $structure , $mutationKey );
				if ( $mutationKey eq $centroid ) {
					#print $subClusterID." - CENTROID: ".$mutationKey."\n";
					$this->setProcessStatus( $structure , $mutationKey , 1 );
					$subClusters->{$structure}->{$subClusterID}->{'centroid'} = $mutationKey;
					$vertices++;
					next;
				}
				if ( $this->isRadiusOkay( $geodesics , $structure , 
							$centroid , $mutationKey ) ) {
					#print $subClusterID." - NEIGHBOR: ".$mutationKey."\n";
					$this->setProcessStatus( $structure , $mutationKey , 1 );
					$subClusters->{$structure}->{$subClusterID}->{'neighbors'}->{$mutationKey} = 1;
					$vertices++;
				}
			}
		}
		if ( $vertices < $MINCLUSTERSIZE ) {
			if ( not $this->anyFiniteGeodesicsRemaining( $geodesics , $structure ) ) {
				#print "Checked ".$checked.", but found no new cluster\n";
				$moreToFind = 0;
			} else {
				#print "Found ".$vertices." from ".$checked." checked in ".$subClusterID." with centroid ".$centroid."\n";
				print( "WARNING: vertices remain, but no centroids identified\n" );
				$moreToFind = 0;
			}
		} else {
			print "Found ".$vertices." from ".$checked." checked in ".$subClusterID." with centroid ".$centroid."\n";
		}
	}
	return $subClusters;
}

sub calculateVertexScore {
	my ( $this , $mutations , $geodesics , $structure ) = @_;
	#my ( $this , $mutations , $geodesics , $structure ) = @_;
	my $scores = {};
	#print STDOUT "HotSpot3D::Cluster::calculateVertexScore\n";#.$x." by ";
	my $max=0;
	my $centroid = "";
	my ( $mutationKey1 , $mutationKey2 , $weight );
	my $x = scalar keys %{$geodesics->{$structure}};
	foreach $mutationKey1 ( keys %{$geodesics->{$structure}} ) {
		my $y = scalar keys %{$geodesics->{$structure}->{$mutationKey1}};
		#print $y."\n";
		#print "mutationKey1 = ".$mutationKey1."\t";
#TODO if alternative transcripts have different proref & proalt, then double counted
		foreach my $refAlt1 ( sort keys %{$mutations->{$mutationKey1}} ) {
			#print $refAlt1."\n";
			my $C = 0;
			my @proteinKeys1 = sort keys %{$mutations->{$mutationKey1}->{$refAlt1}};
			my $proteinKey1 = shift @proteinKeys1;
			#print $mutationKey1."|".$proteinKey1."\n";
			foreach $mutationKey2 ( keys %{$geodesics->{$structure}->{$mutationKey1}}) {
#TODO take only contributions within radius of centroid
				my $geodesic = $geodesics->{$structure}->{$mutationKey1}->{$mutationKey2};
				my $term = 0;
				#next if ( not $this->isRadiusOkay( $geodesics , $structure , $mutationKey1 , $mutationKey2 ) );
				foreach my $refAlt2 ( sort keys %{$mutations->{$mutationKey2}} ) {
					#print $refAlt2."\n";
					my @proteinKeys2 = sort keys %{$mutations->{$mutationKey2}->{$refAlt2}};
					my $proteinKey2 = shift @proteinKeys2;
					#print "\t".$mutationKey2."|".$proteinKey2."\t";
					$weight = 1;
					$term = 0;
					if ( $this->{'vertex_type'} ne $UNIQUE ) {
						if ( exists $mutations->{$mutationKey2} ) {
							$weight = $mutations->{$mutationKey2}->{$refAlt2}->{$proteinKey2};
						}
					}
					#print join( "\t" , ( $weight , $geodesics->{$structure}->{$mutationKey1}->{$mutationKey2} ) )."\t";
					if ( $mutationKey1 ne $mutationKey2 ) { #geodesic is non-zero
						if ( $this->{'vertex_score'} eq $EXPONENTIALS ) {
							$term = $this->exponential( $weight , $geodesic );
							$C += $term;
							#print( "weight: ".$weight."\t\tgeodesic: ".$geodesic."\t\tterm: ".$term."\t\tnewC: ".$C."\n" );
						} else {
							$C += $this->closenessCentrality( $weight , $geodesic );
						}
					} else { #mutationKey1 is same as mutationKey2, geodesic is zer
						if ( $this->{'vertex_score'} eq $EXPONENTIALS ) {
							$term = $this->exponential( $weight , $geodesic );
							$C += $term;
							#print( "weight: ".$weight."\t\tgeodesic: ".$geodesic."\t\tterm: ".$term."\t\tnewC: ".$C."\n" );
						} else {
							if ( $this->{'vertex_type'} eq $WEIGHT ) { 
								$C += $weight;
							} else {
								if ( $refAlt1 ne $refAlt2 ) {
									$C += $weight;
								} else {
									$C += $weight - 1;
								}
							}
						}
					}
					$scores->{$mutationKey1}->{$refAlt1} = $C;
					#print join( "\t" , ( "cid=".$superClusterID.".".$subClusterID , "cent=".$centroid , "maxCc=".$max , "Cc=".$C ) )."\n";
				} #foreach refAlt2
			} #foreach mutationKey2
			#print "Score of ".$mutationKey1.":".$refAlt1." = ".$C."\n";
		} #foreach refAlt1
	} #foreach mutationKey1
	#print join( "\t" , ( "result from calculation: " , $centroid , $max ) )."\n";
#	foreach my $super ( sort keys %{$centrality} ) {
#		foreach my $sub ( sort keys %{$centrality->{$super}} ) {
#			foreach my $mu ( sort keys %{$centrality->{$super}->{$sub}} ) {
#				print "c = ".join( "  " , ( $super.".".$sub , $mu , $centrality->{$super}->{$sub}->{$mu} ) )."\n";
#			}
#		}
#	}

	return $scores;
}

sub closenessCentrality {
	my ( $this , $weight , $geodesic ) = @_;
	return $weight / ( 2**$geodesic );
}

sub exponential {
	my ( $this , $weight , $geodesic ) = @_;
	my $weightExp = abs( $weight / $this->{'weight_scale'} );
	my $geodesicExp = $geodesic / $this->{'length_scale'};
	my $score = $weight * exp( $weightExp - $geodesicExp );
	return $score;
}

sub anyFiniteGeodesicsRemaining {
	my ( $this , $geodesics , $structure ) = @_;
	foreach my $mutationKey1 ( keys %{$geodesics->{$structure}} ) {
		next if ( $this->hasBeenProcessed( $structure , $mutationKey1 ) );
		foreach my $mutationKey2 ( keys %{$geodesics->{$structure}->{$mutationKey1}} ) {
			next if ( $this->hasBeenProcessed( $structure , $mutationKey2 ) 
				and not $this->isRadiusOkay( $geodesics , $structure , $mutationKey1 , $mutationKey2 ) );
			return 1;
		}
	}
	return 0;
}

sub determineStructureClusters {
	my ( $this , $clusterings , $mutations , $distance_matrix ,
		 $structure , $superClusterID , $linesToWrite ) = @_;
	print STDOUT "HotSpot3D::Cluster::determineStructureClusters\n";
	print "\t".$structure."\t".$superClusterID."\n";
	my $geodesics = $this->initializeGeodesics( $clusterings , $superClusterID ,
							$structure , $distance_matrix , $mutations );
	$this->floydWarshall( $geodesics , $structure );
	my $scores = $this->calculateVertexScore( $mutations , $geodesics , $structure );
	my $subClusters = $this->determineSubClusters( $scores , $geodesics , $structure );
	$this->collectOutputLines( $mutations , $geodesics , $structure , 
			$superClusterID , $scores , $subClusters , $linesToWrite );
	return;
}

#TODO use this method to recalculate closeness centralities of acceptable region 
#		or else closeness centralities must be described as measured for the remaining
#		super cluster nodes within range
#sub carveOutSubCluster {
#	my ( $this , $mutations , $geodesics , $structure , $centroid ,
#			$centrality , $superClusterID , $subClusterID ) = @_;
#	#	$this->carveOutSubCluster( $mutations , $geodesics , $structure , $centroid ,
#	#			$centrality , $superClusterID , $subClusterID );
#	my $geods = {};
#	foreach my $mutationKey ( keys %{$geodesics->{$structure}->{$centroid}} ) {
#		$geods->{$structure}->{$centroid}->{$mutationKey} = $geodesics->{$structure}->{$centroid}->{$mutationKey};
#	}
#}

sub collectOutputLines {
	#$this->collectOutputLines( $fh , $mutations , $geodesics , $structure , 
	#		$superClusterID , $subClusterID , $centroid , $scores );
	my ( $this , $mutations , $geodesics , $structure ,
		 $superClusterID , $scores , $subClusters , $linesToWrite ) = @_;
	my $clusterID = $superClusterID;
	foreach my $subClusterID ( sort {$a<=>$b} keys %{$subClusters->{$structure}} ) {
		#print STDOUT "HotSpot3D::Cluster::collectOutputLines\n"; 
		my $centroid = $subClusters->{$structure}->{$subClusterID}->{'centroid'};
		if ( $this->{'structure_dependence'} eq $DEPENDENT 
			 or $this->{'subunit_dependence'} eq $DEPENDENT ) {
			$clusterID = join( "." , ( $superClusterID , $subClusterID , $structure ) );
		} else {
			$clusterID = join( "." , ( $superClusterID , $subClusterID ) );
		}
		print "\tcluster ID ".$clusterID." has centroid ".$centroid."\n";
		my $geodesic = 0;
		my $degrees = $this->getDegrees( $geodesics , $structure , $centroid );
		my ( $gene , $chromosome , $start , $stop ) = @{$this->splitMutationKey( $centroid )};
		my @alternateAnnotations;
		my $proteinChanges = {};
		my ( $reportedTranscript , $reportedAAChange );
		my $weight;
		foreach my $refAlt ( sort keys %{$mutations->{$centroid}} ) {
			my $score = $scores->{$centroid}->{$refAlt};
#TODO make sure this works for in_frame_ins
			my ( $reference , $alternate ) = @{ TGI::Mutpro::Main::Cluster::uncombine( $refAlt ) };
			@alternateAnnotations = sort keys %{$mutations->{$centroid}->{$refAlt}};
			my $reported = shift @alternateAnnotations;
			$weight = $mutations->{$centroid}->{$refAlt}->{$reported};
			( $reportedTranscript , $reportedAAChange ) = @{$this->splitProteinKey( $reported )};
			my $alternateAnnotations = join( "|" , @alternateAnnotations );
			if ( not defined $alternateAnnotations ) {
				$alternateAnnotations = $NULL;
			}
			my $out = join( "\t" , ( $clusterID , $gene , $reportedAAChange , 
									   $degrees , $score , 
									   $geodesic , $weight ,
									   $chromosome , $start , $stop ,
									   $reference , $alternate ,
									   $reportedTranscript , $alternateAnnotations
									 )
							);
			$linesToWrite->{$out} = 1;
		} #foreach refAlt
		$this->setProcessStatus( $structure , $centroid , 2 );
		foreach my $mutationKey2 ( sort keys %{$subClusters->{$structure}->{$subClusterID}->{'neighbors'}} ) {
			$geodesic = $geodesics->{$structure}->{$centroid}->{$mutationKey2};
			#print $centroid." geodesic to ".$mutationKey2."\t".$geodesic."\n";
			$degrees = $this->getDegrees( $geodesics , $structure , $mutationKey2 );
			( $gene , $chromosome , $start , $stop ) = @{$this->splitMutationKey( $mutationKey2 )};
			$proteinChanges = {};
			foreach my $refAlt ( sort keys %{$mutations->{$mutationKey2}} ) {
#TODO make sure this works for in_frame_ins
				my $score = $scores->{$mutationKey2}->{$refAlt};
				my ( $reference , $alternate ) = @{ TGI::Mutpro::Main::Cluster::uncombine( $refAlt ) };
				@alternateAnnotations = sort keys %{$mutations->{$mutationKey2}->{$refAlt}};
				my $reported = shift @alternateAnnotations;
				$weight = $mutations->{$mutationKey2}->{$refAlt}->{$reported};
				( $reportedTranscript , $reportedAAChange ) = @{$this->splitProteinKey( $reported )};
				my $alternateAnnotations = $NULL;
				if ( scalar @alternateAnnotations > 0 ) {
					$alternateAnnotations = join( "|" , @alternateAnnotations );
				}
				my $out = join( "\t" , ( $clusterID , $gene , $reportedAAChange , 

										   $degrees , $score , 
										   $geodesic , $weight ,
										   $chromosome , $start , $stop ,
										   $reference , $alternate ,
										   $reportedTranscript , $alternateAnnotations
										 )
								);
				$linesToWrite->{$out} = 1;
			} #foreach refAlt
			$this->setProcessStatus( $structure , $mutationKey2 , 2 );
		} #foreach other vertex in network
	} #foreach subClusterID
	return;
}

sub getDegrees {
	my ( $this , $geodesics , $structure , $mutationKey ) = @_;
	my $degrees = 0;
	foreach my $neighbor ( keys %{$geodesics->{$structure}->{$mutationKey}} ) {
		my $geodesic = $geodesics->{$structure}->{$mutationKey}->{$neighbor};
		if ( $geodesic < $this->{'3d_distance_cutoff'} ) {
			$degrees++;
		}
	}
	return $degrees;
}

sub limitSubClusterSize {
	my ( $this , $results ) = @_;
	my $outLines = scalar @{$results};
	if ( $outLines < $MINCLUSTERSIZE ) {
		return 0;
	}
	return $outLines;
}

sub floydWarshall {
	my ( $this , $geodesics , $structure ) = @_;
	#print STDOUT "HotSpot3D::Cluster::floydWarshall\n";
	foreach my $mu_k ( keys %{$geodesics->{$structure}} ) {
		#print "\t".$mu_k."\n";
		foreach my $mu_i ( keys %{$geodesics->{$structure}} ) {
			#print "\t\t".$mu_i."\n";
			my ( $dist_ik , $dist_ij , $dist_kj );
			$dist_ik = $this->getElementByKeys( $geodesics , $structure , $mu_i , $mu_k );
			next if ( $dist_ik == $MAXDISTANCE );
			foreach my $mu_j ( keys %{$geodesics->{$structure}} ) {
				if ( exists $geodesics->{$structure}->{$mu_i}->{$mu_j} ) {
					$dist_ij = $geodesics->{$structure}->{$mu_i}->{$mu_j};
				} else {
					$dist_ij = $MAXDISTANCE;
				}
				next if ( $dist_ij == 0 );
				$dist_kj = $this->getElementByKeys( $geodesics , $structure , $mu_k , $mu_j );
				next if ( $dist_ik == $MAXDISTANCE );
				if ( $dist_ij > $dist_ik + $dist_kj ) {
					$geodesics->{$structure}->{$mu_i}->{$mu_j} = $dist_ik + $dist_kj;
					$geodesics->{$structure}->{$mu_j}->{$mu_i} = $dist_ik + $dist_kj;
#					print join( " -- " , ( $mu_i , $mu_j , 
#							$dist_ij , $dist_ik , $dist_kj , 
#							$geodesics->{$structure}->{$mu_i}->{$mu_j} 
#							) )."\n";
				}
			}
		}
	}
	return;
}

1;
