package TGI::Mutpro::Main::Cluster;
#
#----------------------------------
# $Authors: Adam Scott & Sohini Sengupta
# $Date: 2014-01-14 14:34:50 -0500 (Tue Jan 14 14:34:50 CST 2014) $
# $Revision:  $
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

use Data::Dumper;

use TGI::Mutpro::Main::Density;

my $WEIGHT = "weight";
my $RECURRENCE = "recurrence";
my $UNIQUE = "unique";
my $PVALUEDEFAULT = 0.05;
my $DISTANCEDEFAULT = 10;
my $MAXDISTANCE = 100;
my $AVERAGEDISTANCE = "average";
my $SHORTESTDISTANCE = "shortest";
my $NETWORK = "network";
my $DENSITY = "density";
my $INDEPENDENT = "independent";
my $DEPENDENT = "dependent";
my $ANY = "any";

sub new {
    my $class = shift;
    my $this = {};
    $this->{'pairwise_file'} = '3D_Proximity.pairwise';
    $this->{'collapsed_file'} = '3D_Proximity.pairwise.singleprotein.collapsed';
    $this->{'drug_clean_file'} = undef;
    $this->{'output_prefix'} = undef;
    $this->{'p_value_cutoff'} = undef;
    $this->{'3d_distance_cutoff'} = undef;
    $this->{'linear_cutoff'} = 0;
	$this->{'max_radius'} = 10;
	$this->{'vertex_type'} = $RECURRENCE;
	$this->{'distance_measure'} = $AVERAGEDISTANCE;
    $this->{'maf_file'} = undef;
    $this->{'amino_acid_header'} = "amino_acid_change";
    $this->{'transcript_id_header'} = "transcript_name";
    $this->{'weight_header'} = $WEIGHT;
    $this->{'clustering'} = undef;
    $this->{'structure_dependence'} = undef;
    #$this->{'pairwise_file'} = '3D_Proximity.pairwise';
    $this->{'Epsilon'} = undef;
    $this->{'MinPts'} = undef;
    $this->{'number_of_runs'} = undef;
    $this->{'probability_cut_off'} = undef;
    bless $this, $class;
    $this->process();
    return $this;
}

sub process {
    my $this = shift;
	$this->setOptions();
    ## processing procedure
	my $clusterings = {};
	my $distance_matrix = {};
	my $pdb_loc = {};
	my $aa_map = {};
	my $master = {};
	my $mut_chrpos = {};
	my $locations = {};
	my $Variants = {};
	my $transcripts = {};
	my $variants_from_pairs = {};
 	my $mutations = {};
	my $WEIGHT = "weight";
	$this->getDrugMutationPairs( $pdb_loc , $aa_map , $locations , 
								 $mut_chrpos , $master );
	$this->getMutationMutationPairs( $pdb_loc , $aa_map , 
									 $locations , $mut_chrpos , 
									 $master , $clusterings , 
									 $distance_matrix );
	$this->cleanUpClusters( $clusterings );
	$this->getVariantWeight( $variants_from_pairs , $mut_chrpos , 
							 $clusterings , $Variants , 
							 $mutations , $transcripts );
	$this->writeOutput( $clusterings , $Variants , 
						$distance_matrix , $transcripts );
    return 1;
}
#####
#	sub functions
#####
sub setOptions {
	my $this = shift;
    my ( $help, $options );
    unless( @ARGV ) { die $this->help_text(); }
    $options = GetOptions (
        'output-prefix=s' => \$this->{'output_prefix'},
        'pairwise-file=s' => \$this->{'pairwise_file'},
        'collapsed-file=s' => \$this->{'collapsed_file'},
        'drug-clean-file=s' => \$this->{'drug_clean_file'},
        'p-value-cutoff=f' => \$this->{'p_value_cutoff'},
        '3d-distance-cutoff=f' => \$this->{'3d_distance_cutoff'},
        'linear-cutoff=f' => \$this->{'linear_cutoff'},
        'max-radius=f' => \$this->{'max_radius'},
        'vertex-type=s' => \$this->{'vertex_type'},
        'distance-measure=s' => \$this->{'distance_measure'},
        'maf-file=s' => \$this->{'maf_file'},
        'amino-acid-header=s' => \$this->{'amino_acid_header'},
        'transcript-id-header=s' => \$this->{'transcript_id_header'},		
        'weight-header=s' => \$this->{'weight_header'},
        'clustering=s' => \$this->{'clustering'},
        'structure-dependence=s' => \$this->{'structure_dependence'},		
        'help' => \$help,

        'Epsilon=f' => \$this->{'Epsilon'},
        'MinPts=f' => \$this->{'MinPts'},
        'number-of-runs=f' => \$this->{'number_of_runs'},
        'probability-cut-off=f' => \$this->{'probability_cut_off'},
    );
    unless( $options ) { die $this->help_text(); }
	if ( not defined $this->{'clustering'} ) {
		$this->{'clustering'} = $NETWORK;
		warn "HotSpot3D::Cluster warning: no clustering option given, setting to default network\n";
	}
	if ( not defined $this->{'structure_dependence'} ) {
		$this->{'structure_dependence'} = $INDEPENDENT;
		warn "HotSpot3D::Cluster warning: no structure-dependence option given, setting to default independent\n";
	}
	if ( $this->{'clustering'} eq $DENSITY ) {
		if ( $help ) { print STDERR density_help_text(); exit 0; }
		else{
			TGI::Mutpro::Main::Density->new($this);
			exit;
		}
	}
	if ( $help ) { print STDERR help_text(); exit 0; }
	if ( not defined $this->{'p_value_cutoff'} ) {
		if ( not defined $this->{'3d_distance_cutoff'} ) {
			warn "HotSpot3D::Cluster warning: no pair distance limit given, setting to default p-value cutoff = 0.05\n";
			$this->{'p_value_cutoff'} = $PVALUEDEFAULT;
			$this->{'3d_distance_cutoff'} = $MAXDISTANCE;
		} else {
			$this->{'p_value_cutoff'} = 1;
		}
	} else {
		if ( not defined $this->{'3d_distance_cutoff'} ) {
			$this->{'3d_distance_cutoff'} = $MAXDISTANCE;
		}
	}
	print STDOUT "p-value-cutoff = ".$this->{'p_value_cutoff'};
	print STDOUT " & 3d-distance-cutoff = ".$this->{'3d_distance_cutoff'}."\n";
    if ( ( not defined $this->{'collapsed_file'} )
		 and ( not defined $this->{'drug_clean_file'} ) ) {
		warn 'You must provide a collapsed pairs file or drug pairs file! ', "\n";
		die $this->help_text();
	}
	if ( not defined $this->{'drug_clean_file'} ) {
		if ( not -e $this->{'collapsed_file'} ) { 
			warn "The input collapsed pairs file (".$this->{'collapsed_file'}.") does not exist! ", "\n";
			die $this->help_text();
		}
	} elsif ( not defined $this->{'collapsed_file'} ) {
		if ( not -e $this->{'drug_clean_file'} ) { 
			warn "The input drug pairs file (".$this->{'drug_clean_file'}.") does not exist! ", "\n";
			die $this->help_text();
		}
	}
    unless( $this->{'pairwise_file'} ) {
		warn 'You must provide a pairwise file! ', "\n";
		die $this->help_text();
	}
    unless( -e $this->{'pairwise_file'} ) {
		warn "The input pairwise file (".$this->{'pairwise_file'}.") does not exist! ", "\n";
		die $this->help_text();
	}
	if ( $this->{'vertex_type'} ne $RECURRENCE
		 and $this->{'vertex_type'} ne $UNIQUE
		 and $this->{'vertex_type'} ne $WEIGHT ) {
		warn "vertex-type option not recognized as \'recurrence\', \'unique\', or \'weight\'\n";
		warn "Using default vertex-type = \'recurrence\'\n";
		$this->{'vertex_type'} = $RECURRENCE;
	}
	if ( $this->{'distance_measure'} ne $AVERAGEDISTANCE
		 and $this->{'distance_measure'} ne $SHORTESTDISTANCE ) {
		warn "distance-measure option not recognized as \'shortest\' or \'average\'\n";
		warn "Using default distance-measure = \'average\'\n";
		$this->{'distance_measure'} = $AVERAGEDISTANCE;
	}
	if ( $this->{'vertex_type'} ne $UNIQUE ) {
		unless( $this->{'maf_file'} ) {
			warn 'You must provide a .maf file if not using unique vertex type! ', "\n";
			die $this->help_text();
		}
		unless( -e $this->{'maf_file'} ) {
			warn "The input .maf file )".$this->{'maf_file'}.") does not exist! ", "\n";
			die $this->help_text();
		}
	}
	return;
}

sub getMutationMutationPairs {
	#$this->getMutationMutationPairs( $pdb_loc , $aa_map , 
	#								 $locations , $mut_chrpos , 
	#								 $master , $clusterings , 
	#								 $distance_matrix );
	my ( $this , $pdb_loc , $aa_map , 
		 $locations , $mut_chrpos , 
		 $master , $clusterings , 
		 $distance_matrix ) = @_;
	print STDOUT "HotSpot3D::Cluster getting pairwise data\n";
	$this->readPairwise( $pdb_loc , $aa_map , $locations , 
						 $mut_chrpos , $master , $distance_matrix );
	#Pick longest transcript representation of unique
	$this->longestRepresentation( $pdb_loc , $aa_map , $locations , 
						 		  $mut_chrpos , $master );
#####
#	collapsed pairwise data
#####
	$this->readCollapsed( $aa_map , $clusterings , $distance_matrix );
	return;
}

sub readPairwise {
	#$this->readPairwise( $pdb_loc , $aa_map , $locations , 
	#					 $mut_chrpos , $master , $distance_matrix );
	my ( $this , $pdb_loc , $aa_map , $locations , 
		 $mut_chrpos , $master , $distance_matrix ) = @_;
	print STDOUT "\nReading in pairwise data ... \n";
	my $fh = new FileHandle;
	unless( $fh->open( $this->{'pairwise_file'} , "r" ) ) { die "Could not open pairwise file $! \n" };
	map {
		my ( $gene1 , $chr1 , $start1 , $stop1 , $aa_1 , $loc_1 , 
			 $gene2 , $chr2 , $start2 , $stop2 , $aa_2 , $loc_2 , 
			 $infos ) = (split /\t/)[0..4,6,9..13,15,19];
		my $gmu1 = $gene1.":".$aa_1;
		my $gmu2 = $gene2.":".$aa_2;
		$mut_chrpos->{$gmu1}->{$chr1."_".$start1."_".$stop1} = 1;
		$mut_chrpos->{$gmu2}->{$chr2."_".$start2."_".$stop2} = 1;
		$this->redundant( $pdb_loc , $aa_map , $locations , $gene1, $aa_1, $loc_1);
		$this->redundant( $pdb_loc , $aa_map , $locations , $gene2, $aa_2, $loc_2);
		if ( $this->{'distance_measure'} eq $AVERAGEDISTANCE ) {
			$this->getAverageDistance( $gmu1 , $gmu2 , $distance_matrix , $infos );
		}
	} $fh->getlines;
	$fh->close();
	return;
}

sub longestRepresentation {
	#$this->longestRepresentation( $pdb_loc , $aa_map , $locations , 
	#					 		  $mut_chrpos , $master );
	my ( $this , $pdb_loc , $aa_map , $locations , 
		 $mut_chrpos , $master ) = @_;
	foreach my $gene ( keys %{$aa_map} ) {
		foreach my $aa ( @{$aa_map->{$gene}} ) {
			if ( $aa ne 'NA' ) {
				if ( exists $locations->{$gene}->{$aa} ) {
					my @locs = @{$locations->{$gene}->{$aa}};
					my $orig_loc = $aa;
					my $orig_letters = $aa;
					$orig_letters =~ s/[^A-Z]//g;
					$orig_loc =~ s/\D//g;
					my @pdb_locations=@{$pdb_loc->{$gene}};
					foreach my $current_loc ( @locs ) {
						if ( grep{ $_ eq $current_loc } @pdb_locations ) {
							my @idx = grep{ $pdb_locations[$_] eq $current_loc } 0..$#pdb_locations;
							foreach my $current_idx ( @idx ) {
								my $current = $aa_map->{$gene}->[$current_idx];
								my $current_int = $current;
								my $current_letters = $current;
								$current_int =~ s/\D//g;
								$current_letters =~ s/[^A-Z]//g;
								if ( $orig_loc > $current_int && $orig_letters eq $current_letters ) {
									$aa_map->{$gene}->[$current_idx] = $aa;
								} #foreach good AA, keep track
							}
						}
					} #foreach location
				}
			} #if not na
		} #foreach aa
	} #foreach gene
	return;
}

sub readCollapsed{
	#$this->readCollapsed( $aa_map , $clusterings , $distance_matrix );
	my ( $this , $aa_map , $clusterings , $distance_matrix ) = @_;
	print STDOUT "HotSpot3D::Cluster reading/getting collapsed data\n";
	if ( $this->{'collapsed_file'} ) {
		my $fh = new FileHandle;
		unless( $fh->open( $this->{'collapsed_file'} , "r" ) ) {
			die "Could not open collapsed file $! \n"
		};
		map {
			chomp;
			my ( $gene1 , $m1 , $gene2 , $m2 , $lindists , $dist , 
				 $pdb , $pval ) = (split /\t/)[0,1,5,6,11..14];
			my @mus1 = split /,/, $m1;
			my @mus2 = split /,/, $m2;
			my @gm1 = (); my @gm2 = ();
			my $structure = $this->structureDependence( $pdb );
			print "RC - ".$structure."\n";
			for ( my $i = 0 ; $i < @mus1 ; $i++ ) { 
				my $mut1 = $mus1[$i];
				my $mut2 = $mus2[$i];
				if ( ( grep{ $_ eq $mut1 } @{$aa_map->{$gene1}} )
					 && ( grep{ $_ eq $mut2 } @{$aa_map->{$gene2}} ) ) {
					my $first = $gene1.":".$mut1;
					my $second = $gene2.":".$mut2;
					if ( $this->checkPair( $dist , $pval ) == 1 ) {
						my @lindists = split /\,/ , $lindists;
						my $lin_dist = $lindists[0];
						if ( $lin_dist eq "N/A" ) {
							push @gm1 , ( $first ); #unique mutation
							push @gm2 , ( $second ); #unique mutation
							$distance_matrix->{$structure}->{$first}->{$second} = $dist;
							$distance_matrix->{$structure}->{$second}->{$first} = $dist;
						} elsif ( $lin_dist > $this->{'linear_cutoff'} ) {
							push @gm1 , ( $first ); #unique mutation
							push @gm2 , ( $second ); #unique mutation
							$distance_matrix->{$structure}->{$first}->{$second} = $dist;
							$distance_matrix->{$structure}->{$second}->{$first} = $dist;
						} #if linear distance okay
					} #if spatial significance okay
				} #if okay transcript representations
			} #foreach transcript representation of mutations
			my @mutations = @gm1;
			push @mutations , @gm2;
			$this->AHC( $pval , $dist , $clusterings , \@mutations , $structure );
		} $fh->getlines;
		$fh->close();
	} #if using collapsed pairs file
	return;
}

sub cleanUpClusters {
	#$this->cleanUpClusters( $clusterings );
	my ( $this , $clusterings ) = @_;
	print STDOUT "HotSpot3D::Cluster clean up clusters\n";
    my $i = 0;
    ## REASSIGN CLUSTER IDS BY +1 INCREMENTS (REMOVE GAPS IN ID LIST)
	foreach my $structure ( sort keys %{$clusterings} ) {
		print "CUC - ".$structure."\n";
		foreach my $id ( sort {$a<=>$b} keys %{$clusterings->{$structure}} ) {
			my $size = scalar @{$clusterings->{$structure}->{$id}};
			if ( scalar @{$clusterings->{$structure}->{$id}} < 2 ) {
				delete $clusterings->{$structure}->{$id};
			}
		} #assure empty or singleton clusters (if any) are removed
		my @keys = sort { $a <=> $b } keys %{$clusterings->{$structure}};
		while ( $i < scalar keys %{$clusterings->{$structure}} ) {
			if ( $keys[$i] != $i ) {
				$clusterings->{$structure}->{$i} = $clusterings->{$structure}->{$keys[$i]};
				delete $clusterings->{$structure}->{$keys[$i]};
			}
			$i++;
		}
		#@keys = sort { $a <=> $b } keys %{$clusterings->{$structure}};
	}
	return;
}

sub getVariantWeight {
	#$this->getVariantWeight( $variants_from_pairs , $mut_chrpos , 
	#						 $clusterings , $Variants , 
	#						 $mutations , $transcripts );
	my ( $this , $variants_from_pairs , $mut_chrpos , 
		 $clusterings , $Variants , 
		 $mutations , $transcripts ) = @_;
    ## SHAVE CLUSTERS TO CORE
    #use distance_matrix matrix
	print STDOUT "HotSpot3D::Cluster get variant weighting\n";
	if ( $this->{'vertex_type'} ne $UNIQUE ) {
		print STDOUT "Reading .maf to get ";
		if ( $this->{'vertex_type'} eq $RECURRENCE ) { print "recurrence\n"; }
		elsif ( $this->{'vertex_type'} eq $WEIGHT ) { print "weight\n"; }
		foreach my $structure ( keys %{$clusterings} ) {
			print "GVW - ".$structure."\n";
			foreach my $id ( keys %{$clusterings->{$structure}} ) {
				foreach my $gene_mut ( @{$clusterings->{$structure}->{$id}} ) {
					my ( $gene , $mutation ) = split /\:/ , $gene_mut;
					if ( exists $mut_chrpos->{$gene_mut} ) {
						foreach my $chr_start_stop ( keys %{$mut_chrpos->{$gene_mut}} ) {
							my $variant = join( "_" , ( $gene , $mutation , $chr_start_stop ) );
							$variants_from_pairs->{$variant} = 1;
						} #foreach chr_start_stop in mut_chrpos->gene_mut
					} #if gene_mut in mut_chrpos
				} #foreach gene_mut in clusterings->structure->id
			} #foreach id in clusterings->structure
		} #foreach structure in clusterings
		##Mutation recurrence or weight from MAF
		$this->readMAF( $variants_from_pairs , $Variants , 
						$mutations , $transcripts );
	}
	return;
}

sub readMAF {
	#$this->readMAF( $variants_from_pairs , $Variants , 
	#				 $mutations , $transcripts );
	my ( $this , $variants_from_pairs , $Variants , 
		 $mutations , $transcripts ) = @_;
	my $fh = new FileHandle;
	die "Could not open .maf file\n" unless( $fh->open( $this->{'maf_file'} , "r" ) );
	print STDOUT "\nReading in .maf ...\n";
	my $mafi = 0;
	my $headline = $fh->getline(); chomp( $headline );
	my %mafcols = map{ ( $_ , $mafi++ ) } split( /\t/ , $headline );
	unless(		defined( $mafcols{"Hugo_Symbol"} )
			and defined( $mafcols{"Chromosome"} )
			and defined( $mafcols{"Start_Position"} )
			and defined( $mafcols{"End_Position"} )
			and defined( $mafcols{"Reference_Allele"} )
			and defined( $mafcols{"Tumor_Seq_Allele2"} )
			and defined( $mafcols{"Tumor_Sample_Barcode"} )
			and defined( $mafcols{$this->{"transcript_id_header"}} )
			and defined( $mafcols{$this->{"amino_acid_header"}} ) ) {
		die "not a valid .maf file! Check transcript and amino acid change headers.\n";
	}
	my @mafcols = ( $mafcols{"Hugo_Symbol"},
					$mafcols{"Chromosome"},
					$mafcols{"Start_Position"},
					$mafcols{"End_Position"},
					$mafcols{"Reference_Allele"},
					$mafcols{"Tumor_Seq_Allele2"},
					$mafcols{"Tumor_Sample_Barcode"},
					$mafcols{$this->{"transcript_id_header"}},
					$mafcols{$this->{"amino_acid_header"}} );
	if ( $this->{'vertex_type'} eq $WEIGHT ) {
		unless( defined( $mafcols{$this->{"weight_header"}} ) ) { die "HotSpot3D::Cluster error: weight vertex-type chosen, but weight-header not recocgnized\n"; };
		push @mafcols , $mafcols{$this->{"weight_header"}};
	}
	map {
		chomp;
		my @line = split /\t/;
		if ( $#line >= $mafcols[-1] && $#line >= $mafcols[-2] ) { #makes sure custom maf cols are in range
			my ( $gene , $chr , $start , $stop , $reference , $tumorAllele , $barID , $transcript_name , $aachange );
			my $weight = 1;
			if ( $this->{'vertex_type'} eq $WEIGHT ) {
				( $gene , $chr , $start , $stop , $reference , $tumorAllele , $barID , $transcript_name , $aachange , $weight ) = @line[@mafcols];
			} else {
				( $gene , $chr , $start , $stop , $reference , $tumorAllele , $barID , $transcript_name , $aachange ) = @line[@mafcols];
			}
			my $variant = join( "_" , ( $gene , $aachange , $chr , $start , $stop ) );
			if ( exists $variants_from_pairs->{$variant} ) {
				my $gene_aachange = $gene.":".$aachange;
				$transcripts->{$gene_aachange}->{$transcript_name.":".$aachange} = $chr.":".$start.":".$stop;
				if ( exists $Variants->{$gene_aachange} ) {
					if ( $this->{'vertex_type'} ne $WEIGHT ) {
						if ( not exists $mutations->{$variant}->{$barID} ) {
							$Variants->{$gene_aachange} += $weight;
							$mutations->{$variant}->{$barID} += $weight;
						}
					}
				} else {
					$Variants->{$gene_aachange} = $weight;
					$mutations->{$variant}->{$barID} = $weight;
				}
			}
		}
	} $fh->getlines;
	$fh->close();
	return;
}

sub writeOutput {
	#$this->writeOutput( $clusterings , $Variants , 
	#					$distance_matrix , $transcripts );
	my ( $this , $clusterings , $Variants , $distance_matrix , $transcripts ) = @_;
	print STDOUT "HotSpot3D::Cluster write cluster output\n";
	my $outFilename = $this->generateFilename();
	print STDOUT "Creating cluster output file: ".$outFilename."\n";
	my $fh = new FileHandle;
	die "Could not create clustering output file\n" unless( $fh->open( $outFilename , "w" ) );
	$fh->print( join( "\t" , ( 	"Cluster" , "Gene/Drug" , "Mutation/Gene" , 
								"Degree_Connectivity" , "Closeness_Centrality" , 
								"Geodesic_From_Centroid" , "Recurrence" , 
								"Chromosome" , "Start" , "Stop" , 
								"Transcript" , "Alternative_Transcripts"
							 )
					)."\n"
			  );
	print STDOUT "Getting Cluster ID's & Centroids\n";
	my $numstructures = scalar keys %{$clusterings};
	my $numclusters = 0;
	foreach my $structure ( keys %{$clusterings} ) {
		$numclusters += scalar keys %{$clusterings->{$structure}};
		foreach my $clus_num ( keys %{$clusterings->{$structure}} ) {
			my @clus_mut = @{$clusterings->{$structure}->{$clus_num}};
			print $structure."\t".$clus_num."\n";
			$this->centroid( $Variants , $distance_matrix , $clus_num , 
							 \@clus_mut , $fh , 0 , 1 , $transcripts , $structure );
		} #foreach cluster
	} #foreach structure
	if ( $this->{'structure_dependence'} eq $DEPENDENT ) {
		print STDOUT "Found ".$numclusters." super-clusters on ".$numstructures." structures\n";
	} else {
		print STDOUT "Found ".$numclusters." super-clusters\n";
	}
	return;
}

sub generateFilename {
	my $this = shift;
	my @outFilename;
	if ( $this->{'output_prefix'} ) {
		push @outFilename , $this->{'output_prefix'};
	} else {
		if ( $this->{'maf_file'} ) {
			my $maf = basename( $this->{'maf_file'} );
			push @outFilename , $maf;
		}
		if ( $this->{'collapsed_file'} ) {
			my $collapsed = basename( $this->{'collapsed_file'} );
			push @outFilename , $collapsed;
		}
		if ( defined $this->{'drug_clean_file'} ) {
			my $clean = basename( $this->{'drug_clean_file'} );
			if ( $clean ne '' and scalar @outFilename > 1 ) {
				push @outFilename , $clean;
			} elsif ( $clean ne '' ) {
				push @outFilename , $clean;
			}
		}
		push @outFilename , "l".$this->{'linear_cutoff'};
		my $m = "a";
		if ( $this->{'distance_measure'} eq $SHORTESTDISTANCE ) { $m = "s"; }
		if ( $this->{'3d_distance_cutoff'} != $MAXDISTANCE ) {
            if ( $this->{'p_value_cutoff'} != 1 ) {
                push @outFilename , "p".$this->{'p_value_cutoff'};
                push @outFilename , $m."d".$this->{'3d_distance_cutoff'};
            } else {
                push @outFilename , $m."d".$this->{'3d_distance_cutoff'};
            }
        } else {
            if ( $this->{'p_value_cutoff'} != 1 ) {
                push @outFilename , "p".$this->{'p_value_cutoff'};
            }
        }
		push @outFilename , "r".$this->{'max_radius'};
	}
	push @outFilename , "clusters";
	return join( "." , @outFilename );
}

sub getDrugMutationPairs {
	#$this->getDrugMutationPairs( $pdb_loc , $aa_map , $locations , 
	#							 $mut_chrpos , $master );
	my ( $this , $pdb_loc , $aa_map , $locations , $mut_chrpos , $master ) = shift;
	print STDOUT "HotSpot3D::Cluster get drug mutation pairs\n";
	$this->readDrugClean( $pdb_loc , $aa_map , $locations , $mut_chrpos , $master );
	#Pick longest transcript representation of unique
	$this->getBestRepresentation( $aa_map , $locations );
	#cluster drug-mutation pairs and build distance matrix
	$this->createDistanceMatrix( $master , $aa_map );
	return;
}

sub readDrugClean {
	#$this->readDrugClean( $pdb_loc , $aa_map , $locations , $mut_chrpos , $master );
	my ( $this , $pdb_loc , $aa_map , $locations , $mut_chrpos , $master ) = shift;
	print STDOUT "HotSpot3D::Cluster read drug clean\n";
    if ( $this->{'drug_clean_file'} ) { #if drug pairs included
		my $fh = new FileHandle;
		unless( $fh->open( $this->{'drug_clean_file'} , "r" ) ) {
			die "Could not open drug pairs data file $! \n"
		};
		my $rxi = 0;
		my $headline = $fh->getline(); chomp( $headline );
		my %drugcols = map{ ( $_ , $rxi++ ) } split( /\t/ , $headline );
		my @required = ( "Drug" , "PDB_ID" , "Gene" , "Chromosome" , "Start" , 
						 "Stop" , "Amino_Acid_Change" , 
						 "Mutation_Location_In_PDB" ,
						 "3D_Distance_Information" );
		unless(	defined( $drugcols{"Drug"} )										#0
			and defined( $drugcols{"PDB_ID"} )										#2
			and defined( $drugcols{"Gene"} )										#6
			and defined( $drugcols{"Chromosome"} )									#7
			and defined( $drugcols{"Start"} )										#8
			and defined( $drugcols{"Stop"} )										#9
			and defined( $drugcols{"Amino_Acid_Change"} )							#10
			and defined( $drugcols{"Mutation_Location_In_PDB"} )					#12
			and defined( $drugcols{"3D_Distance_Information"} ) ) {					#17
			die "not a valid drug-clean file\n";
		}
		my @wantrxcols = (	$drugcols{"Drug"} ,
							$drugcols{"PDB_ID"} ,
							$drugcols{"Gene"} ,
							$drugcols{"Chromosome"} ,
							$drugcols{"Start"} ,
							$drugcols{"Stop"} ,
							$drugcols{"Amino_Acid_Change"} ,
							$drugcols{"Mutation_Location_In_PDB"} ,
							$drugcols{"3D_Distance_Information"} );
		map { 
				chomp;
				my @line = split /\t/; 
				map{ $_ =~ s/"//g } @line;
				my ( $drug, $pdb, $gene2, $chr , $start , $stop , $m2, $loc, $infor_3d ) = @line[@wantrxcols];
				my ( $dist, $pdb2, $pval ) = split / /, $infor_3d;
				$mut_chrpos->{$gene2.":".$m2}->{$chr."_".$start."_".$stop} = 1;
				$pval =~ s/"//g;
				my $first = $drug.":".$gene2;
				my $second = $gene2.":".$m2;
				my $info = $dist.":".$pval;
				my $structure = &checkStructureDependence( $pdb );
				$master->{$structure}->{$first}->{$second} = $info; #store necessary pair info
				$this->redundant($pdb_loc, $aa_map, $locations , $gene2, $m2, $loc); #filter transcripts
		} $fh->getlines; 
		$fh->close();
	} #if drug pairs included
	return;
}

sub getBestRepresentation {
	#$this->getBestRepresentation( $aa_map , $locations );
	my ( $this , $aa_map , $locations , $pdb_loc ) = @_;
	print STDOUT "HotSpot3D::Cluster::getBestRepresentation\n";
	foreach my $gene ( keys %{$aa_map} ) {
		foreach my $aa ( @{$aa_map->{$gene}} ) {
			if ( $aa ne 'NA' ) {
				if ( exists $locations->{$gene}->{$aa} ) {
					my @locs = @{$locations->{$gene}->{$aa}};
					my $orig_loc = $aa;
					my $orig_letters = $aa;
					$orig_letters =~ s/[^A-Z]//g;
					$orig_loc =~ s/\D//g;
					my @pdb_locations=@{$pdb_loc->{$gene}};
					foreach my $current_loc ( @locs ) {
						if ( grep{ $_ eq $current_loc } @pdb_locations ) {
							my @idx = grep{ $pdb_locations[$_] eq $current_loc } 0..$#pdb_locations;
							foreach my $current_idx ( @idx ) {
								my $current = $aa_map->{$gene}->[$current_idx];
								my $current_int = $current;
								my $current_letters = $current;
								$current_int =~ s/\D//g;
								$current_letters =~ s/[^A-Z]//g;
								if ( $orig_loc > $current_int && $orig_letters eq $current_letters ) {
									$aa_map->{$gene}->[$current_idx] = $aa;
								} #foreach good AA, keep track
							}
						}
					} #foreach location
				} #if location exists
			} #if not na
		} #foreach aa
	} #foreach gene
	return ( $aa_map );
}

sub createDistanceMatrix {
	#$this->createDistanceMatrix( $master , $aa_map );
	my ( $this , $master , $aa_map , $clusterings ) = @_;
	print STDOUT "HotSpot3D::Cluster::createDistanceMatrix\n";
	my $distance_matrix = {};
	foreach my $structure ( keys %{$master} ) {
		foreach my $first ( keys %{$master->{$structure}} ) {
			foreach my $second ( keys %{$master->{$structure}->{$first}} ) {
				my ( $gene2 , $m2 ) = split ":" , $second;
				if ( grep{ $_ eq $m2 } @{$aa_map->{$gene2}} ) { 
					my @mutations = ();
					push @mutations , $first;
					if ( $this->{'vertex_type'} eq $UNIQUE ) { $m2 =~ s/\D+(\d+)\D+/$1/g; }
					$second = $gene2.":".$m2;
					push @mutations , $second; #@mus2;
					my ( $dist , $pval ) = split ":" , $master->{$structure}->{$first}->{$second};
					$this->AHC( $pval , $dist , $clusterings->{$structure} , 
								\@mutations , $structure );
					if ( $this->checkPair( $dist , $pval ) ) {
						$distance_matrix->{$structure}->{$first}->{$second} = $dist;
						$distance_matrix->{$structure}->{$second}->{$first} = $dist;
					}
				}
			}
		}
	} #foreach in master
	return;
}

sub centroid{
	my ( $this , $Variants , $distance_matrix , $clus_num , $clus_mut , 
		 $fh , $recluster , $counter , $transcripts , $structure )=@_;
	my %dist;
	foreach my $mut1 ( @{$clus_mut} ) { #initialize geodesics
		my @mu1 = split( ":" , $mut1 );
		foreach my $mut2 ( @{$clus_mut} ) {
			my @mu2 = split( ":" , $mut2 );
			if ( $mu1[1] =~ /p\./ ) { $mu1[1] =~ s/\D//g; }
			if ( $mu2[1] =~ /p\./ ) { $mu2[1] =~ s/\D//g; }
			if ( exists $distance_matrix->{$structure}->{$mut1}->{$mut2} ) {
				print join( "\t" , ( $structure , $mut1 , $mut2 , $distance_matrix->{$structure}->{$mut1}->{$mut2} ) )."\n";
				$dist{$mut1}{$mut2} = $distance_matrix->{$structure}->{$mut1}->{$mut2};
			} elsif ( ( $mu1[0] eq $mu2[0] ) && ( $mu1[1] eq $mu2[1] ) ) {
				$dist{$mut1}{$mut2} = 0;
			} else {
				$dist{$mut1}{$mut2} = 1000000;
			}
		}
	}
	&floydwarshall( \%dist, \@{$clus_mut}); #get geodesics
	my %centrality;
	my $max=0;
	my $centroid='';
	my $count=0;
	my ( $other , $weight );
	foreach my $current ( keys %dist ) {
		my $C = 0;
		foreach $other ( keys %{$dist{$current}}) {
			$weight = 1; #stays as 1 if vertex_type eq 'unique'
			if ( $this->{'vertex_type'} ne $UNIQUE ) {
				if ( exists $Variants->{$other} ) {
					$weight = $Variants->{$other};
				}
			}
			if ( $current ne $other ) {
				if ( $dist{$current}{$other} <= $this->{'max_radius'} ) {
					$C += $weight/( 2**$dist{$current}{$other} );
					$count+=1;
				} else{
					$recluster=1;
				}
			} else { #current is same as other
				if ( $this->{'vertex_type'} eq $WEIGHT ) { 
					$C += $weight;
				} else {
					$C += $weight -1;
				}
			}
			$centrality{$clus_num}{$current} = $C;
			if ( $C > $max ) {
				$max = $C;
				$centroid = $current;
			}
		} #foreach other
	} #foreach current

	my $cluster_size=scalar @{$clus_mut};
	$count=0;
	if ( exists $dist{$centroid} ) {
		foreach $other ( keys %{$dist{$centroid}} ) {
			my $geodesic = $dist{$centroid}{$other};
			my $degrees = scalar keys %{$dist{$other}};
			my $closenesscentrality = $centrality{$clus_num}{$other};
			my ( $gene , $mutation ) = split /\:/ , $other;
			my ( $reportedTranscript , $altTranscript , 
				 $chromosome , $start , $stop ) = $this->getTranscriptInfo( $transcripts , $other );
			$weight = 1;

			if ( $geodesic <= $this->{'max_radius'} ) {
				if ( exists $Variants->{$other} ) {
					$weight = $Variants->{$other};
				}
				$count+=1;
				my $clusterID = $clus_num;
				if ($recluster==1){
					if ( $this->{'structure_dependence'} eq $DEPENDENT ) {
						$clusterID = join( "." , ( $clus_num , $counter , $structure ) );
					} else {
						$clusterID = join( "." , ( $clus_num , $counter ) );
					}
					$fh->print( join( "\t" , ( $clusterID , $gene , $mutation , 
											   $degrees , $closenesscentrality , 
											   $geodesic , $weight ,
											   $chromosome , $start , $stop  ,
											   $reportedTranscript , $altTranscript 
											 )
									)."\n"
							  );
					my $index=0;
					$index++ until $clus_mut->[$index] eq $other;
					splice(@{$clus_mut}, $index, 1);
				} else {
					if ( $this->{'structure_dependence'} eq $DEPENDENT ) {
						$clusterID = join( "." , ( $clus_num , $structure ) );
					}
					$fh->print( join( "\t" , ( 	$clusterID , $gene , $mutation , $degrees , 
												$closenesscentrality , $geodesic , $weight , 
												$chromosome , $start , $stop  ,
												$reportedTranscript , $altTranscript 
											 )
									)."\n"
							  );
				}
			}
		} #foreach other vertex in network
	} #if dist for centroid
	if ($count<2) { return; }
	if ($recluster==1) {
		$counter+=1;
		$this->centroid( $Variants , $distance_matrix , $clus_num , 
						 $clus_mut , $fh , $recluster , $counter , 
						 $transcripts , $structure );
	}
}

## CLUSTERING FUNCTION - AGGLOMERATIVE HIERARCHICAL CLUSTERING (AHC)
sub AHC {
	#$this->AHC( $pval , $dist , $clusterings , \@mutations );
    my ( $this, $pval , $dist , $clusterings , $mutations , $structure ) = @_;
    if ( $this->checkPair( $dist , $pval ) ) { #meets desired significance
        my ( @temp, @found, @combine ); 
        my ( @uniq, $c );
        foreach $c ( keys %{$clusterings->{$structure}} ) { #each cluster
            my @mus_in_cluster = @{$clusterings->{$structure}->{$c}};
            foreach my $mu ( @{$mutations} ) {
				foreach ( @mus_in_cluster ) {
					if ( $mu eq $_ ) { push @combine , $c; }
				}
			}
        }
        my @uniqcombo = uniq @combine; #cluster types
        if ( scalar @uniqcombo > 0 ) { #collapse clusters into one
            my $collapse_to = min @uniqcombo; #cluster type
            my $j = 0; #iterator type
            while ( $j < scalar @uniqcombo ) {
                push @{$clusterings->{$structure}->{$collapse_to}} , @{$clusterings->{$structure}->{$uniqcombo[$j]}}; #mutation types
                push @{$clusterings->{$structure}->{$collapse_to}} , @{$mutations}; #
                if ( $collapse_to ne $uniqcombo[$j] ) { delete $clusterings->{$structure}->{$uniqcombo[$j]}; }
                $j++;
            }
            @{$clusterings->{$structure}->{$collapse_to}} = uniq @{$clusterings->{$structure}->{$collapse_to}};
        } else { #new cluster
            if ( scalar keys %{$clusterings->{$structure}} > 0 ) {
                $c = ( max keys %{$clusterings->{$structure}} ) + 1;
            } else { $c = 0; }
            push @temp , @{$mutations};
            @uniq = uniq @temp; #mutation types
            $clusterings->{$structure}->{$c} = \@uniq;
        }
        $c = scalar keys %{$clusterings->{$structure}};
    } #if pval significant
    return;
}

sub redundant {
    my ( $this, $pdb_loc, $aa_map, $locations , $gene, $aa, $loc ) = @_;
    my $aa_orig = $aa;
	if ( not grep{ $_ eq $loc} @{$locations->{$gene}->{$aa}} ) {
		push @{$locations->{$gene}->{$aa}} , $loc;
	}
    if ( exists $pdb_loc->{$gene} ) { #if pdb_loc has gene
        if ( grep{$_ eq $loc} @{$pdb_loc->{$gene}} ) { #if primary location list has this mapping location
            my @array = @{$pdb_loc->{$gene}}; #each gene has list of mutations
            my @idx = grep{$array[$_] eq $loc} 0..$#array; #get mutation locations matching mapping location
            my $pos = $idx[0]; #
            my $prev_aa = $aa_map->{$gene}[$pos];
            $prev_aa =~ s/\D//g;
            $aa =~ s/\D//g;
            $loc =~ s/\D//g;
            if ( abs($prev_aa-$loc)>abs($aa-$loc) ) {
                foreach( @idx ) {
                    $pdb_loc->{$gene}[$_] = "NA";
                    $aa_map->{$gene}[$_] = "NA";
                }
                push(@{$pdb_loc->{$gene}}, $loc);
                push(@{$aa_map->{$gene}}, $aa_orig);
            } elsif ( abs($prev_aa-$loc)==abs($aa-$loc) && !(grep(/^\Q$aa_orig\E/,@{$aa_map->{$gene}})) ) {
                push(@{$pdb_loc->{$gene}}, $loc);
                push(@{$aa_map->{$gene}}, $aa_orig);
            }
        } else {
            push( @{$pdb_loc->{$gene}}, $loc );
            push( @{$aa_map->{$gene}}, $aa_orig );
        } #if loc among 
    } else { #don't have gene yet
        push( @{$pdb_loc->{$gene}}, $loc );
        push( @{$aa_map->{$gene}}, $aa_orig );
    } #if have gene

    return 1;
}

sub floydwarshall {
	my ( $dist , $clus_mut ) = @_;
	foreach my $mu_k ( @{$clus_mut} ) {
		foreach my $mu_i ( @{$clus_mut} ) {
			my $dist_ik = $dist->{$mu_i}->{$mu_k};
			foreach my $mu_j ( @{$clus_mut} ) {
				my $dist_ij = $dist->{$mu_i}->{$mu_j};
				my $dist_kj = $dist->{$mu_k}->{$mu_j};
				if ( $dist_ij > $dist_ik + $dist_kj ) {
					$dist->{$mu_i}->{$mu_j} = $dist_ik + $dist_kj;
				}
			}
		}
	}
}

sub getTranscriptInfo {
	my ( $this , $transcripts , $other ) = @_;
	my ( $reportedTranscript , $altTranscript , $transcript );
	my ( $chromosome , $start , $stop );
	my ( $gene , $mu ) = split /\:/ , $other;
	foreach my $tranmu ( sort keys %{$transcripts->{$other}} ) {
		my ( $transcript , $mutation ) = split /\:/ , $tranmu;
		my $css = $transcripts->{$other}->{$tranmu};
		( $chromosome , $start , $stop ) = split /\:/ , $css;
		if ( $mu eq $mutation ) {
			if ( not $reportedTranscript ) {
				$reportedTranscript = $transcript;
			} else {
				$reportedTranscript .= "|".$transcript;
			}
		} else {
			if ( not $altTranscript ) {
				$altTranscript = $transcript.":".$mutation;
			} else {
				$altTranscript .= "|".$transcript.":".$mutation;
			}
		}
	} #foreach tranmu

	if ( not $reportedTranscript ) { $reportedTranscript = "NA"; }
	if ( not $altTranscript ) { $altTranscript = "NA"; }
	return ( $reportedTranscript , $altTranscript , $chromosome , $start , $stop );
}

sub checkPair {
	my ( $this , $dist , $pval ) = @_;
	if ( $this->{'3d_distance_cutoff'} == $MAXDISTANCE ) { #3d-dist undef & p-val def
		if ( $pval < $this->{'p_value_cutoff'} ) {
			return 1;
		}
	} elsif ( $this->{'p_value_cutoff'} == 1 ) { #3d-dist def & p-val undef
		if ( $dist < $this->{'3d_distance_cutoff'} ) {
			return 1;
		}
	} else { #3d-dist def & p-val def
		if ( $dist < $this->{'3d_distance_cutoff'} and $pval < $this->{'p_value_cutoff'} ) {
			return 1;
		}
	}
	return 0;
}

sub getAverageDistance {
	my ( $this , $gmu1 , $gmu2 , $distance_matrix , $infos ) = @_;
	my @infos = split /\|/ , $infos;
	my $avgDistance = {};
	foreach my $info ( @infos ) {
		chomp( $info );
		next unless ( $info );
		my ( $distance , $pdbID , $pvalue ) = split / / , $info;
		my $structure = $this->structureDependence( $pdbID );
		$avgDistance->{$structure} += $distance;
	}
	foreach my $structure ( keys %{$avgDistance} ) {
		$avgDistance->{$structure} = $avgDistance->{$structure}  / ( scalar @infos );
		$distance_matrix->{$structure}->{$gmu1}->{$gmu2} = $avgDistance->{$structure};
		$distance_matrix->{$structure}->{$gmu2}->{$gmu1} = $avgDistance->{$structure};
	}
	return;
}

sub structureDependence {
	my $this = shift;
	my $structure = shift;
	if ( $this->{'structure_dependence'} eq $DEPENDENT ) {
		return $structure;
	}
	return $ANY;
}

sub density_help_text{
    my $this = shift;
        return <<HELP

Usage: hotspot3d density [options]

                             REQUIRED
--pairwise-file              3D pairwise data file

                             OPTIONAL
--Epsilon                    Epsilon value, default: 10
--MinPts                     MinPts, default: 4
--number-of-runs             Number of density clustering runs to perform before the cluster membership probability being calculated, default: 10
--probability-cut-off        Clusters will be formed with variants having at least this probability, default: 100 

--help                       this message

HELP
}

sub help_text{
    my $this = shift;
        return <<HELP

Usage: hotspot3d cluster [options]

                             REQUIRED
--pairwise-file              3D pairwise data file

                             AT LEAST ONE
--collapsed-file             Either (or concatenated) pairwise.singleprotein.collapsed & pairwise.complex.collapsed data
--drug-clean-file            Either (or concatenated) drugs.target.clean & drugs.nontarget.clean data

                             OPTIONAL
--output-prefix              Output prefix, default: 3D_Proximity
--p-value-cutoff             P_value cutoff (<), default: 0.05 (if 3d-distance-cutoff also not set)
--3d-distance-cutoff         3D distance cutoff (<), default: 100 (if p-value-cutoff also not set)
--linear-cutoff              Linear distance cutoff (> peptides), default: 0
--max-radius                 Maximum cluster radius (max network geodesic from centroid, <= Angstroms), default: 10
--clustering                 Cluster using network or density-based methods (network or density), default: network
--vertex-type                Graph vertex type for network-based clustering (recurrence, unique, or weight), default: recurrence
--distance-measure           Pair distance to use (shortest or average), default: average
--structure-dependence       Clusters for each structure or across all structures (dependent or independent), default: independent
--maf-file                   .maf file used in proximity search step (used if vertex-type = recurrence)
--transcript-id-header       .maf file column header for transcript id's, default: transcript_name
--amino-acid-header          .maf file column header for amino acid changes, default: amino_acid_change 
--weight-header              .maf file column header for mutation weight, default: weight (used if vertex-type = weight)

--help                       this message

HELP

}

1;
