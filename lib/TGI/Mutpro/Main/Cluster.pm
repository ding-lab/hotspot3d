package TGI::Mutpro::Main::Cluster;
#
#----------------------------------
# $Authors: Adam Scott & Sohini Sengupta
# $Date: 2014-01-14 14:34:50 -0500 (Tue Jan 14 14:34:50 CST 2014) $
# $Revision: 4 $
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
use TGI::Mutpro::Main::Density;

my $WEIGHT = "weight";
my $RECURRENCE = "recurrence";
my $UNIQUE = "unique";
my $PVALUEDEFAULT = 0.05;
my $DISTANCEDEFAULT = 10;
my $MAXDISTANCE = 10000;
my $AVERAGEDISTANCE = "average";
my $SHORTESTDISTANCE = "shortest";
my $NETWORK = "network";
my $DENSITY = "density";
my $INDEPENDENT = "independent";
my $DEPENDENT = "dependent";
my $ANY = "any";

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
    my $this = {};
    $this->{'pairwise_file'} = '3D_Proximity.pairwise';
    $this->{'maf_file'} = undef;
    $this->{'drug_clean_file'} = undef;
    $this->{'output_prefix'} = undef;
    $this->{'p_value_cutoff'} = undef;
    $this->{'3d_distance_cutoff'} = undef;
    $this->{'linear_cutoff'} = 0;
	$this->{'max_radius'} = 10;
	$this->{'vertex_type'} = $RECURRENCE;
	$this->{'distance_measure'} = $AVERAGEDISTANCE;
    $this->{'amino_acid_header'} = "amino_acid_change";
    $this->{'transcript_id_header'} = "transcript_name";
    $this->{'weight_header'} = $WEIGHT;
    $this->{'clustering'} = undef;

	$this->{'parallel'} = undef;
	$this->{'max_processes'} = undef;

    $this->{'structure_dependence'} = undef;
    $this->{'subunit_dependence'} = undef;
#TODO add meric based on protein or gene id, if protein, need hugo.uniprot.pdb.transcript.csv file
    $this->{'meric_type'} = undef;
	
	$this->{'processed'} = undef;
	$this->{'distance_matrix'} = undef;
	$this->{'mutations'} = undef;
	$this->{'results'} = undef;

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
	my $distance_matrix = {};
 	my $mutations = {};
	my $WEIGHT = "weight";

	$this->readMAF( $mutations );
#	foreach my $mk ( sort keys %{$mutations} ) {
#		foreach my $ra ( sort keys %{$mutations->{$mk}} ) {
#			foreach my $pk ( sort keys %{$mutations->{$mk}->{$ra}} ) {
#				print join( " -- " , ( $mk , $ra , $pk , $mutations->{$mk}->{$ra}->{$pk} ) )."\n";
#			}
#		}
#	}
	$this->getDrugMutationPairs( $distance_matrix );
	$this->getMutationMutationPairs( $distance_matrix );
	#$this->initializeSameSiteDistancesToZero( $distance_matrix );
	$this->networkClustering( $mutations , $distance_matrix );

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
        'drug-clean-file=s' => \$this->{'drug_clean_file'},
        'p-value-cutoff=f' => \$this->{'p_value_cutoff'},
        '3d-distance-cutoff=f' => \$this->{'3d_distance_cutoff'},
        'linear-cutoff=f' => \$this->{'linear_cutoff'},
        'max-radius=f' => \$this->{'max_radius'},
        'vertex-type=s' => \$this->{'vertex_type'},
        'number-of-runs=f' => \$this->{'number_of_runs'},
        'probability-cut-off=f' => \$this->{'probability_cut_off'},
        'distance-measure=s' => \$this->{'distance_measure'},
        'maf-file=s' => \$this->{'maf_file'},
        'amino-acid-header=s' => \$this->{'amino_acid_header'},
        'transcript-id-header=s' => \$this->{'transcript_id_header'},
        'weight-header=s' => \$this->{'weight_header'},
        'clustering=s' => \$this->{'clustering'},
        'structure-dependence=s' => \$this->{'structure_dependence'},
        'subunit-dependence=s' => \$this->{'subunit_dependence'},
        'parallel=s' => \$this->{'parallel'},
        'max-processes=i' => \$this->{'max_processes'},
        'meric-type=s' => \$this->{'meric_type'},

        'help' => \$help,
    );
    unless( $options ) { die $this->help_text(); }
	if ( not defined $this->{'clustering'} ) {
		$this->{'clustering'} = $NETWORK;
		warn "HotSpot3D::Cluster::setOptions warning: no clustering option given, setting to default network\n";
	}
	if ( $help ) { print STDERR help_text(); exit 0; }
	my $acceptableParallels = { $LOCAL => 1 , $NONE => 1 }; #, $BSUB => 1 };
	if ( not defined $this->{'parallel'} ) { 
		$this->{'parallel'} = $NONE;
		warn "HotSpot3D::Cluster::setOptions warning: no parallel option given, setting to default none\n";
	} elsif ( not exists $acceptableParallels->{$this->{'parallel'}} ) {
		warn "HotSpot3D::Cluster::setOptions error: unrecognized parallel option given, please use either none, local, or bsub\n";
		die $this->help_text();
	} else {
		if ( $this->{'parallel'} eq $LOCAL and ( not defined $this->{'max_processes'} ) ) {
			warn "HotSpot3D::Cluster::setOptions error: local parallelization specified, but max-processes was not set\n";
			die $this->help_text();
		}
	}
	if ( not defined $this->{'structure_dependence'} ) {
		$this->{'structure_dependence'} = $INDEPENDENT;
		warn "HotSpot3D::Cluster::setOptions warning: no structure-dependence option given, setting to default independent\n";
	}
	if ( not defined $this->{'subunit_dependence'} ) {
		$this->{'subunit_dependence'} = $INDEPENDENT;
		warn "HotSpot3D::Cluster::setOptions warning: no subunit-dependence option given, setting to default independent\n";
	}
	if ( not defined $this->{'meric_type'} ) {
		$this->{'meric_type'} = $UNSPECIFIED;
		warn "HotSpot3D::Cluster::setOptions warning: no meric-type option given, setting to default unspecified\n";
	}
	if ( not defined $this->{'p_value_cutoff'} ) {
		if ( not defined $this->{'3d_distance_cutoff'} ) {
			warn "HotSpot3D::Cluster::setOptions warning: no pair distance limit given, setting to default p-value cutoff = 0.05\n";
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
	if ( defined $this->{'drug_clean_file'} ) {
		if ( not -e $this->{'drug_clean_file'} ) { 
			warn "The input drug pairs file (".$this->{'drug_clean_file'}.") does not exist! ", "\n";
			die $this->help_text();
		}
	} else {
		warn "HotSpot3D::Cluster::setOptions warning: no drug-clean-file included (cannot produce drug-mutation clusters)!\n";
	}
    if ( defined $this->{'pairwise_file'} ) {
		if ( not -e $this->{'pairwise_file'} ) { 
			warn "HotSpot3D::Cluster::setOptions error: the input pairwise file (".$this->{'pairwise_file'}.") does not exist!\n";
			die $this->help_text();
		}
	} else {
		warn "HotSpot3D::Cluster::setOptions error: must provide a pairwise-file!\n";
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
	unless( $this->{'maf_file'} ) {
		warn 'You must provide a .maf file! ', "\n";
		die $this->help_text();
	}
	unless( -e $this->{'maf_file'} ) {
		warn "The input .maf file )".$this->{'maf_file'}.") does not exist! ", "\n";
		die $this->help_text();
	}
	my $tempMericHash = {}; ### check for correct meric option 
	$tempMericHash->{$MULTIMER} = 0;
	$tempMericHash->{$MONOMER} = 0;
	$tempMericHash->{$HOMOMER} = 0;
	$tempMericHash->{$HETEROMER} = 0;
	$tempMericHash->{$UNSPECIFIED} = 0;
	$tempMericHash->{$INTRA} = 0;
	$tempMericHash->{$INTER} = 0;
	if ( not exists $tempMericHash->{$this->{'meric_type'}} ) {
		die "Error: meric-type should be one of the following: intra, monomer, homomer, inter, heteromer, multimer, unspecified\n";
	}
	print STDOUT "=====Parameters=====\n";
	print STDOUT " p-value-cutoff        = ".$this->{'p_value_cutoff'}."\n";
	print STDOUT " 3d-distance-cutoff    = ".$this->{'3d_distance_cutoff'}."\n";
	print STDOUT " max-radius            = ".$this->{'max_radius'}."\n";
	print STDOUT " vertex-type           = ".$this->{'vertex_type'}."\n";
	print STDOUT " distance-measure      = ".$this->{'distance_measure'}."\n";
	print STDOUT " structure-dependence  = ".$this->{'structure_dependence'}."\n";
	print STDOUT " subunit-dependence    = ".$this->{'subunit_dependence'}."\n";
	print STDOUT " meric-type            = ".$this->{'meric_type'}."\n";
	print STDOUT " parallel              = ".$this->{'parallel'}."\n";
	if ( defined $this->{'max_processes'} ) {
		print STDOUT " max-processes         = ".$this->{'max_processes'}."\n";
	}
	print STDOUT "====================\n";

	if ( $this->{'clustering'} eq $DENSITY ) {
		if ( $help ) {
			die density_help_text();
		}
		else{
			TGI::Mutpro::Main::Density->new($this);
			exit;
		}
	}
	return;
}

sub getMutationMutationPairs {
	#$this->getMutationMutationPairs( $distance_matrix );
	my ( $this , $distance_matrix ) = @_;
	print STDOUT "HotSpot3D::Cluster::getMutationMutationPairs\n";
	$this->readPairwise( $distance_matrix );
	return;
}

sub readPairwise { # shared
	#$this->readPairwise( $distance_matrix );
	my ( $this , $distance_matrix ) = @_;
	print STDOUT "HotSpot3D::Cluster::readPairwise\n";
	my $fh = new FileHandle;
	unless( $fh->open( $this->{'pairwise_file'} , "r" ) ) { die "Could not open pairwise file $! \n" };
	my $pdbCount;
	map {
		my ( $gene1 , $chromosome1 , $start1 , $stop1 , $aa_1 , $chain1 , $loc_1 , $domain1 , $cosmic1 , 
			 $gene2 , $chromosome2 , $start2 , $stop2 , $aa_2 , $chain2 , $loc_2 , $domain2 , $cosmic2 , 
			 $linearDistance , $infos ) = split /\t/ , $_;

		if ( $this->checkMeric( $gene1 , $gene2 , $chain1 , $chain2 ) ) {
			#print $_."\n";
			$chain1 =~ s/\[(\w)\]/$1/g;
			$chain2 =~ s/\[(\w)\]/$1/g;
			my $proteinMutation = TGI::ProteinVariant->new();
			my $mutation1 = TGI::Variant->new();
			$mutation1->gene( $gene1 );
			$mutation1->chromosome( $chromosome1 );
			$mutation1->start( $start1 );
			$mutation1->stop( $stop1 );
			$proteinMutation->aminoAcidChange( $aa_1 );
			$mutation1->addProteinVariant( $proteinMutation );

			my $mutation2 = TGI::Variant->new();
			$mutation2->gene( $gene2 );
			$mutation2->chromosome( $chromosome2 );
			$mutation2->start( $start2 );
			$mutation2->stop( $stop2 );
			$proteinMutation->aminoAcidChange( $aa_2 );
			$mutation2->addProteinVariant( $proteinMutation );
#			print "pair structures ".join( "  " , ( $infos , $chain1 , $chain2 ) )."\n";
			$this->setDistance( $distance_matrix , $mutation1 , $mutation2 , 
								$chain1 , $chain2 , $infos , $pdbCount );
		}
	} $fh->getlines;
	$fh->close();
#	print "DISTANCE MATRIX\n";
#	foreach my $s ( sort keys %{$distance_matrix} ) {
#		foreach my $mu1 ( sort keys %{$distance_matrix->{$s}} ) {
#			foreach my $mu2 ( sort keys %{$distance_matrix->{$s}->{$mu1}} ) {
#				print "\t".join( "  " , ( $s , $mu1 , $mu2 , $this->getElementByKeys( $distance_matrix , $s , $mu1 , $mu2 ) ) )."\n";
#			}
#		}
#	}
	return;
}

sub checkMeric {
	my ( $this , $gene1 , $gene2 , $chain1 , $chain2 ) = @_;
#	print join( "  " , ( $gene1 , $chain1 , $gene2 , $chain2 , "" ) );
	if ( $this->{'meric_type'} eq $UNSPECIFIED ) {
#		print "unspec okay\n";
		return 1;
	} elsif ( $this->{'meric_type'} eq $INTRA ) {
		if ( $gene1 eq $gene2 ) {
#			print "intra okay\n";
			return 1;
		}
	} elsif ( $this->{'meric_type'} eq $INTER ) {
		if ( $chain1 ne $chain2 and $gene1 ne $gene2 ) {
#			print "inter okay\n";
			return 1;
		}
	} elsif ( $this->{'meric_type'} eq $MONOMER ) {
		if ( $chain1 eq $chain2 and $gene1 eq $gene2 ) {
#			print "mono okay\n";
			return 1;
		}
	} elsif ( $this->{'meric_type'} eq $MULTIMER ) {
		if ( $chain1 ne $chain2 ) {
#			print "multi okay\n";
			return 1;
		}
	} elsif ( $this->{'meric_type'} eq $HOMOMER ) {
		if ( $chain1 ne $chain2 and $gene1 eq $gene2 ) {
#			print "homo okay\n";
			return 1;
		}
	} elsif ( $this->{'meric_type'} eq $HETEROMER ) {
		if ( $chain1 ne $chain2 and $gene1 ne $gene2 ) {
#			print "heter okay\n";
			return 1;
		}
	} else {
#		print "not okay\n";
		return 0;
	}
#	print "backup not okay\n";
	return 0;
}

sub readMAF{ # shared
	#$this->readMAF( $mutations );
	my ( $this , $mutations ) = @_;
	print STDOUT "HotSpot3D::Cluster::readMAF\n";
	my $fh = new FileHandle;
	die "Could not open .maf file\n" unless( $fh->open( $this->{'maf_file'} , "r" ) );
	my $headline = $fh->getline(); chomp( $headline );
	my $mafi = 0;
	my %mafcols = map{ ( $_ , $mafi++ ) } split( /\t/ , $headline );
	unless( defined( $mafcols{"Hugo_Symbol"} )
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
					$mafcols{"Variant_Classification"},
					$mafcols{"Reference_Allele"},
					$mafcols{"Tumor_Seq_Allele2"},
					$mafcols{"Tumor_Sample_Barcode"},
					$mafcols{$this->{"transcript_id_header"}},
					$mafcols{$this->{"amino_acid_header"}} );
	if ( $this->{'vertex_type'} eq $WEIGHT ) {
		unless( defined( $mafcols{$this->{"weight_header"}} ) ) {
			die "HotSpot3D::Cluster::readMAF error: weight vertex-type chosen, but weight-header not recocgnized\n";
		};
		push @mafcols , $mafcols{$this->{"weight_header"}};
	}
	map {
		chomp;
		my @line = split /\t/;
		#print $_."\n";
		if ( $#line >= $mafcols[-1] && $#line >= $mafcols[-2] ) { #makes sure custom maf cols are in range
			my ( $gene , $chromosome , $start , $stop , $classification , $reference ,
				 $alternate , $barID , $transcript_name , $aachange );
			my $weight = 1;
			if ( $this->{'vertex_type'} eq $WEIGHT ) {
				( $gene , $chromosome , $start , $stop , $classification , $reference ,
				  $alternate , $barID , $transcript_name , $aachange , $weight
				) = @line[@mafcols];
			} else {
				( $gene , $chromosome , $start , $stop , $classification , $reference ,
				  $alternate , $barID , $transcript_name , $aachange 
				) = @line[@mafcols];
			}
			if ( $classification =~ /Missense/ 
				or $classification =~ /In_Frame/ ) {
				my $mutation = TGI::Variant->new();
				$mutation->gene( $gene );
				$mutation->chromosome( $chromosome );
				$mutation->start( $start );
				$mutation->stop( $stop );
				$mutation->reference( $reference );
				$mutation->alternate( $alternate );
				my $proteinMutation = TGI::ProteinVariant->new();
				$proteinMutation->transcript( $transcript_name );
				$proteinMutation->aminoAcidChange( $aachange );
				$mutation->addProteinVariant( $proteinMutation );
				$this->setMutation( $mutations , $mutation , $barID , $weight );
			} #if mutation is missense or in frame
		} #if columns in range
	} $fh->getlines; #map
	$fh->close();
	return;
}

sub setMutation {
	my ( $this , $mutations , $mutation , $barID , $weight ) = @_;

	my $mutationKey = $this->makeMutationKey( $mutation , "" );
	my $refAlt = &combine( $mutation->reference() , $mutation->alternate() );
	my $proteinKey = $this->makeProteinKey( $mutation );
	#print "setMutation: ".$refAlt."\n";
	#print join( "\t" , ( $mutationKey , $proteinKey , $barID , $weight ) )."\n";

	if ( exists $mutations->{$mutationKey}->{$refAlt}->{$proteinKey} ) {
		#print "existing\t";
		if ( $this->{'vertex_type'} ne $WEIGHT ) {
			$mutations->{$mutationKey}->{$refAlt}->{$proteinKey} += 1;
		} else {
			$mutations->{$mutationKey}->{$refAlt}->{$proteinKey} = $weight;
		}
	} else {
		#print "new\t";
		if ( $this->{'vertex_type'} eq $WEIGHT ) {
			$mutations->{$mutationKey}->{$refAlt}->{$proteinKey} = $weight;
		} else {
			$mutations->{$mutationKey}->{$refAlt}->{$proteinKey} += 1;
		}
	} #if mutation exists
	#my $w = $mutations->{$mutationKey}->{$refAlt}->{$proteinKey};
	#print $w."\t".$n."\n";
	return;
}

sub getMutationInfo {
	my ( $this , $mutations , $mutationKey ) = @_;
	my $weights = {};
	foreach my $refAlt ( sort keys %{$mutations->{$mutationKey}} ) {
		foreach my $proteinKey ( sort keys %{$mutations->{$mutationKey}->{$refAlt}} ) {
			my ( $transcript , $aaChange ) = @{$this->splitProteinKey( $proteinKey )};
			$aaChange =~ m/p\.(\D*)(\d+)(\D*)/;
			$weights->{$refAlt}->{$proteinKey} = $mutations->{$mutationKey}->{$refAlt}->{$proteinKey};
		}
	}
	return $weights;
}

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
	$this->writeClustersFile();
	return;
}

sub localParallelNetworkClustering {
	my ( $this , $clusterings , $mutations , $distance_matrix ) = @_;
	print STDOUT "HotSpot3D::Cluster::localParallelNetworkClustering\n";
	print STDOUT "\tparallel clustering over structures with up to ".$this->{'max_processes'}." processes\n";
	my $tempDir = File::Temp->newdir( TEMPLATE => 'hs3dXXXXX' );
	my $pm = Parallel::ForkManager->new( $this->{'max_processes'} , $tempDir );
	my %finalLines;
	$pm->run_on_finish( 
		sub {
			my ( $pid , $exit_code , $ident , $exit_signal , $core_dump , $data ) = @_;
			if ( defined $data ) {
				if ( exists $clusterings->{$data->[0]} ) {
					$finalLines{$data->[0]} = $data->[1];
				}
			} else {
				print qq|No data from child $pid!\n|;
			}
		}
	);
	DATA_LOOP:
	foreach my $structure ( sort keys %{$distance_matrix} ) {
		my $pid = $pm->start and next DATA_LOOP;
		my %lines;
		#my $lines = "";
		foreach my $superClusterID ( sort keys %{$clusterings->{$structure}} ) {
			my @lines;
			my $subClusterID = 0;
			#$lines .= $this->determineStructureClusters( $clusterings , $mutations , $distance_matrix , 
					#$structure , $superClusterID , $subClusterID , () );
			$lines{$superClusterID} = $this->determineStructureClusters( $clusterings , $mutations , $distance_matrix , 
					$structure , $superClusterID , $subClusterID , \@lines );
		}
		my @finalLines;
		foreach my $superClusterID ( keys %lines ) {
			push @finalLines , @{$lines{$superClusterID}};
		}
		my @pair = ( $structure , \@finalLines );
		#my @pair = ( $structure , $lines );
		$pm->finish( 0 , \@pair );
	}
	$pm->wait_all_children;
	foreach my $structure ( keys %finalLines ) {
		#if ( exists $this->{'distance_matrix'}->{$structure} ) {
			#$this->{'results'}->{$structure} = join( "\n" , @{$lines{$structure}} )."\n";
			$this->{'results'}->{$structure} = $finalLines{$structure};
		#}
	}
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

sub noParallelNetworkClustering {
	my ( $this , $clusterings , $mutations , $distance_matrix ) = @_;
	print STDOUT "serially clustering over structures\n";
	foreach my $structure ( sort keys %{$distance_matrix} ) {
		my %lines;
		#my $lines = "";
		foreach my $superClusterID ( sort keys %{$clusterings->{$structure}} ) {
			my @lines;
			my $subClusterID = 0;
			$lines{$superClusterID} = $this->determineStructureClusters( $clusterings , $mutations , $distance_matrix , 
					$structure , $superClusterID , $subClusterID , \@lines );
			# $lines .= $this->determineStructureClusters( $clusterings , $mutations , $distance_matrix , 
			#		$structure , $superClusterID , $subClusterID , 0 );
		}
		#$this->{'results'}->{$structure} = join( "\n" , @lines );
		#$this->{'results'}->{$structure}  = $lines;
		foreach my $superClusterID ( keys %lines ) {
			push @{$this->{'results'}->{$structure}} , @{$lines{$superClusterID}};
		}
	}
	return;
}

sub writeClustersFile {
	my ( $this ) = @_;
	print STDOUT "HotSpot3D::Cluster::writeClustersFile\n";
	my $outFilename = $this->generateFilename();
	print STDOUT "\tCreating cluster output file: ".$outFilename."\n";
	my $fh = new FileHandle;
	die "Could not create clustering output file\n" unless( $fh->open( $outFilename , "w" ) );
	$fh->print( join( "\t" , ( 	"Cluster" , "Gene/Drug" , "Mutation/Gene" , 
								"Degree_Connectivity" , "Closeness_Centrality" , 
								"Geodesic_From_Centroid" , "Weight" , 
								"Chromosome" , "Start" , "Stop" , 
								"Reference" , "Alternate" ,
								"Transcript" , "Alternative_Transcripts"
							 )
					)."\n"
			  );
	foreach my $structure ( %{$this->{'results'}} ) {
		if ( exists $this->{'results'}->{$structure} ) {
			my @mutations = uniq( @{$this->{'results'}->{$structure}} );
			my $nMutations = scalar( uniq( @mutations ) );
			if ( $nMutations > 0 ) {
				print STDOUT "\twriting for structure: ".$structure;
				print STDOUT "\thas ".$nMutations." mutations\n";
				foreach my $line ( @mutations ) {
					$fh->print( $line."\n" );
				}
			} else {
				print STDOUT "\tno mutations to write for structure: ".$structure."\n";
			}
		}
	}
	$fh->close();
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
	#$this->getDrugMutationPairs( $distance_matrix );
	my ( $this , $distance_matrix ) = shift;
	print STDOUT "HotSpot3D::Cluster::getDrugMutationPairs\n";
	$this->readDrugClean( $distance_matrix );
	return;
}

sub readDrugClean {
	#$this->readDrugClean( $distance_matrix );
	my ( $this , $distance_matrix ) = shift;
	print STDOUT "HotSpot3D::Cluster::readDrugClean\n";
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
		unless(	defined( $drugcols{"Drug"} )						#0
			and defined( $drugcols{"PDB_ID"} )						#2
			and defined( $drugcols{"Chain1"} )						#2
			and defined( $drugcols{"Gene"} )						#6
			and defined( $drugcols{"Chromosome"} )					#7
			and defined( $drugcols{"Start"} )						#8
			and defined( $drugcols{"Stop"} )						#9
			and defined( $drugcols{"Amino_Acid_Change"} )			#10
			and defined( $drugcols{"Chain2"} )						#2
			and defined( $drugcols{"Mutation_Location_In_PDB"} )	#12
			and defined( $drugcols{"3D_Distance_Information"} ) ) {	#17
			die "not a valid drug-clean file\n";
		}
		my @wantrxcols = (	$drugcols{"Drug"} ,
							$drugcols{"PDB_ID"} ,
							$drugcols{"Chain1"} ,
							$drugcols{"Gene"} ,
							$drugcols{"Chromosome"} ,
							$drugcols{"Start"} ,
							$drugcols{"Stop"} ,
							$drugcols{"Amino_Acid_Change"} ,
							$drugcols{"Chain2"} ,
							$drugcols{"Mutation_Location_In_PDB"} ,
							$drugcols{"3D_Distance_Information"} );
		my $pdbCount = {};
		map { 
				chomp;
				my @line = split /\t/; 
				map{ $_ =~ s/"//g } @line;
				my ( $drug, $pdb , $chain1 , $gene, 
					 $chromosome , $start , $stop , $aaChange , $chain2 , 
					 $loc, $infos ) = @line[@wantrxcols];
				#my ( $dist, $pdb2, $pval ) = split / /, $infos;
				$infos =~ s/"//g;

				my $structure = $this->makeStructureKey( $pdb , $chain1 , $chain2 );
				my $proteinMutation = TGI::ProteinVariant->new();
				my $mutation = TGI::Variant->new();
				$mutation->gene( $gene );
				$mutation->chromosome( $chromosome );
				$mutation->start( $start );
				$mutation->stop( $stop );
				$proteinMutation->aminoAcidChange( $aaChange );
				$mutation->addProteinVariant( $proteinMutation );
				my $soCalledMutation = TGI::Variant->new(); #for the purpose of making keys, will look like a mutation
				$soCalledMutation->gene( $gene );
				$soCalledMutation->chromosome( $drug );
				$soCalledMutation->start( $structure );
				$soCalledMutation->stop( $structure );
				$proteinMutation->aminoAcidChange( $aaChange );
				$soCalledMutation->addProteinVariant( $proteinMutation );
				$this->setDistance( $distance_matrix , $soCalledMutation , 
									$mutation , $chain1 , $chain2 , 
									$infos , $pdbCount );
		} $fh->getlines; 
		$fh->close();
	} #if drug pairs included
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
				&numSort( \@combine );
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
						&numSort( \@ids );
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

sub calculateClosenessCentrality {
	my ( $this , $mutations , $geodesics , $structure , $superClusterID , $subClusterID ) = @_;
	#my ( $this , $mutations , $geodesics , $structure ) = @_;
	my $centrality = {};
	#print STDOUT "HotSpot3D::Cluster::calculateClosenessCentrality\n";#.$x." by ";
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
				#next if ( not $this->isRadiusOkay( $geodesics , $structure , $mutationKey1 , $mutationKey2 ) );
				foreach my $refAlt2 ( sort keys %{$mutations->{$mutationKey2}} ) {
					#print $refAlt2."\n";
					my @proteinKeys2 = sort keys %{$mutations->{$mutationKey2}->{$refAlt2}};
					my $proteinKey2 = shift @proteinKeys2;
					#print "\t".$mutationKey2."|".$proteinKey2."\t";
					$weight = 1;
					if ( $this->{'vertex_type'} ne $UNIQUE ) {
						if ( exists $mutations->{$mutationKey2} ) {
							$weight = $mutations->{$mutationKey2}->{$refAlt2}->{$proteinKey2};
						}
					}
					#print join( "\t" , ( $weight , $geodesics->{$structure}->{$mutationKey1}->{$mutationKey2} ) )."\t";
					if ( $mutationKey1 ne $mutationKey2 ) {
						$C += $weight/( 2**$geodesics->{$structure}->{$mutationKey1}->{$mutationKey2} );
					} else { #mutationKey1 is same as mutationKey2
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
					$centrality->{$superClusterID}->{$subClusterID}->{$mutationKey1} = $C;
					#$centrality->{$mutationKey1} = $C;
					if ( $C > $max ) {
						$max = $C;
						$centroid = $mutationKey1;
					}
					#print join( "\t" , ( "cid=".$superClusterID.".".$subClusterID , "cent=".$centroid , "maxCc=".$max , "Cc=".$C ) )."\n";
				} #foreach refAlt2
			} #foreach mutationKey2
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

	return ( $centroid , $centrality );
}

#sub determineCentroid {
#	my ( $this , $centrality , $superClusterID ) = @_;
#	my $max = 0;
#	my $centroid = "";
#	foreach my $mutationKey ( keys %{$centrality->{$superClusterID}} ) {
#		my $C = $centrality->{$superClusterID}->{$mutationKey};
#		if ( $C > $max ) {
#			$max = $C;
#			$centroid = $mutationKey1;
#		}
#	}
#	return ( $centroid , $max );
#}

sub checkProcessedDistances {
	my ( $this , $distance_matrix , $structure ) = @_;
	my $count = 0;
	#print "CHECK PROCESSED: \n";
	foreach my $mutationKey1 ( keys %{$distance_matrix->{$structure}} ) {
		#print "\t".$mutationKey1;
		if ( not $this->hasBeenProcessed( $structure , $mutationKey1 ) ) {
			#print "no\n";
			$count += 1;
		#} else {
			#print "yes\n";
		}
	}
	return $count;
}

sub mutationWeights {
	my ( $this , $mutations , $mutationKey ) = @_;
	my $weight = [];
	foreach my $refAlt ( sort keys %{$mutations->{$mutationKey}} ) {
		my @proteinKeys = sort keys %{$mutations->{$mutationKey}->{$refAlt}};
		my $proteinKey = shift @proteinKeys;
		push @{$weight} , $mutations->{$mutationKey}->{$refAlt}->{$proteinKey};
	}
	return $weight;
}

sub sum {
	my $array = shift;
	my $sum = 0;
	foreach my $value ( @{$array} ) {
		$sum += $value;
	}
	return $sum;
}

sub anyFiniteGeodesicsRemaining {
	my ( $this , $mutations , $geodesics , $structure ) = @_;
	foreach my $mutationKey1 ( keys %{$geodesics->{$structure}} ) {
		next if ( $this->hasBeenProcessed( $structure , $mutationKey1 ) );
		foreach my $mutationKey2 ( keys %{$geodesics->{$structure}->{$mutationKey1}} ) {
			next if ( $this->hasBeenProcessed( $structure , $mutationKey2 ) 
				and not $this->isRadiusOkay( $geodesics , $structure , $mutationKey1 , $mutationKey2 ) );
			#my $sumWeights1 = &sum( $this->mutationWeights( $mutations , $mutationKey1 ) );
			#my $sumWeights2 = &sum( $this->mutationWeights( $mutations , $mutationKey2 ) );
			#next if ( $sumWeights1 == 1 and $sumWeights2 == 1 );
			#print "FINITE: ".join( "..." , ( $mutationKey1 , $mutationKey2 , $geodesics->{$structure}->{$mutationKey1}->{$mutationKey2} ) )."\n";
			return 1;
		}
	}
	return 0;
}

sub determineStructureClusters {
	my ( $this , $clusterings , $mutations , $distance_matrix ,
		 $structure , $superClusterID , $subClusterID , $linesToWrite ) = @_;
	#print STDOUT "HotSpot3D::Cluster::determineStructureClusters\n";
	#print "\t".$structure."\t".$superClusterID.".".$subClusterID."\n";
	my $geodesics = $this->initializeGeodesics( $clusterings , $superClusterID ,
							$structure , $distance_matrix , $mutations );
	if ( $this->anyFiniteGeodesicsRemaining( $mutations , $geodesics , $structure ) ) {
		$this->floydWarshall( $geodesics , $structure );
		my ( $centroid , $centrality ) = $this->calculateClosenessCentrality( 
												$mutations , $geodesics , 
												$structure , $superClusterID , 
												$subClusterID );
		#my ( $centroid , $centrality ) = $this->calculateClosenessCentrality( 
		#										$mutations , $geodesics , 
		#										$structure );
		#my $writtenLines = $this->writeCluster( $mutations , $geodesics , $structure , 
		my ( $writtenLines , $lines ) = $this->writeCluster( $mutations , $geodesics , $structure , 
				$superClusterID , $subClusterID , $centroid , $centrality );
		
		my $count = $this->checkProcessedDistances( $geodesics , $structure );
		#print join( "\t" , ( "recluster?" , $count , 
		#		$superClusterID.".".$subClusterID , $structure ) );
		if ( $count >= 2 and $writtenLines ) {
			push @{$linesToWrite} , @{$lines};
			#$linesToWrite .= join( "\n" , @{$lines} );
			$subClusterID += 1;
			#print " yes\n";
			#$linesToWrite .= $this->determineStructureClusters( $clusterings , 
								#$mutations , $geodesics , $structure , 
								#$superClusterID , $subClusterID , $linesToWrite )};
			push @{$linesToWrite} , @{$this->determineStructureClusters( $clusterings , 
								$mutations , $geodesics , $structure , 
								$superClusterID , $subClusterID , $linesToWrite )};
		#} else {
			#print " no\n";
			#return $linesToWrite;
		}
	#	my $numstructures = scalar keys %{$distance_matrix->{$structure}};
	#	if ( $this->{'structure_dependence'} eq $DEPENDENT ) {
	#		print STDOUT "Found ".$numclusters." super-clusters on ".$numstructures." structures\n";
	#	} else {
	#		print STDOUT "Found ".$numclusters." super-clusters\n";
	#	}
	} 
	return $linesToWrite;
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

sub setProcessStatus { # shared
	my ( $this , $structure , $mutationKey , $status ) = @_;
	$this->{'processed'}->{$structure}->{$mutationKey} = $status;
	return $this->{'processed'}->{$structure}->{$mutationKey};
}

sub hasBeenProcessed { # shared
	my ( $this , $structure , $mutationKey ) = @_;
	if ( $this->{'processed'}->{$structure}->{$mutationKey} ) {
		return 1;
	}
	return 0;
}

sub resetProcessed {
    my $this = shift;
	my $structure = shift;
    my $SetOfNodes = {};
    if ( @_ ) {
        $SetOfNodes = shift;
    } else {
        $SetOfNodes = $this->{'processed'}->{$structure};
    }

    foreach my $mutationKey ( keys %{$SetOfNodes} ) {
		#print join( "\t" , ( $mutationKey , $SetOfNodes->{$mutationKey} , $this->{'processed'}->{$mutationKey} ) )."\t";
        $this->setProcessStatus( $structure , $mutationKey , 0 );
		#print $this->{'processed'}->{$mutationKey}."\n";
    }

    return 1;
}

sub writeCluster {
	#$this->writeCluster( $fh , $mutations , $geodesics , $structure , 
	#		$superClusterID , $subClusterID , $centroid , $centrality );
	my ( $this , $mutations , $geodesics , $structure ,
		 $superClusterID , $subClusterID , $centroid , $centrality ) = @_;
	my $writtenLines = 0;
	my @linesToWrite;
	my $clusterID = $superClusterID;
	#print STDOUT "HotSpot3D::Cluster::writeCluster\n"; 
	if ( $this->{'structure_dependence'} eq $DEPENDENT 
		 or $this->{'subunit_dependence'} eq $DEPENDENT ) {
		$clusterID = join( "." , ( $superClusterID , $subClusterID , $structure ) );
	} else {
		$clusterID = join( "." , ( $superClusterID , $subClusterID ) );
	}
	#print "\tcluster ID ".$clusterID."\n";
	my $geodesic = 0;
	my $degrees = scalar keys %{$geodesics->{$structure}->{$centroid}}; #TODO update to only count subcluster nodes
	my $closenessCentrality = $centrality->{$superClusterID}->{$subClusterID}->{$centroid};
	#my $closenessCentrality = $centrality->{$centroid};
	my ( $gene , $chromosome , $start , $stop ) = @{$this->splitMutationKey( $centroid )};
	my @alternateAnnotations;
	my $proteinChanges = {};
	my ( $reportedTranscript , $reportedAAChange );
	my $weight; # = $weights->{$proteinKey};
	foreach my $refAlt ( sort keys %{$mutations->{$centroid}} ) {
#TODO make sure this works for in_frame_ins
		my ( $reference , $alternate ) = @{&uncombine( $refAlt )};
		@alternateAnnotations = sort keys %{$mutations->{$centroid}->{$refAlt}};
		#print join( "," , ( @alternateAnnotations ) )."\n";
		my $reported = shift @alternateAnnotations;
		#print "\t".$reported."\n";
		$weight = $mutations->{$centroid}->{$refAlt}->{$reported};
		( $reportedTranscript , $reportedAAChange ) = @{$this->splitProteinKey( $reported )};
		my $alternateAnnotations = join( "|" , @alternateAnnotations );
		#$fh->print( join( "\t" , ( $clusterID , $gene , $reportedAAChange , 
		push @linesToWrite , join( "\t" , ( $clusterID , $gene , $reportedAAChange , 
								   $degrees , $closenessCentrality , 
								   $geodesic , $weight ,
								   $chromosome , $start , $stop ,
								   $reference , $alternate ,
								   $reportedTranscript , $alternateAnnotations
								 )
						);
				  #);
		#$writtenLines += 1;
	} #foreach refAlt
	$this->setProcessStatus( $structure , $centroid , 1 );
	foreach my $mutationKey2 ( sort keys %{$geodesics->{$structure}->{$centroid}} ) {
		next if ( $this->hasBeenProcessed( $structure , $mutationKey2 ) );
		$geodesic = $geodesics->{$structure}->{$centroid}->{$mutationKey2};
		next if ( $geodesic > $this->{'max_radius'} ); 
		#print $centroid." geodesic to ".$mutationKey2."\t".$geodesic."\n";
		$degrees = scalar keys %{$geodesics->{$structure}->{$mutationKey2}}; #TODO update to only count subcluster nodes
		$closenessCentrality = $centrality->{$superClusterID}->{$subClusterID}->{$mutationKey2};
		#$closenessCentrality = $centrality->{$mutationKey2};
		( $gene , $chromosome , $start , $stop ) = @{$this->splitMutationKey( $mutationKey2 )};
		$proteinChanges = {};
		foreach my $refAlt ( sort keys %{$mutations->{$mutationKey2}} ) {
#TODO make sure this works for in_frame_ins
			my ( $reference , $alternate ) = @{&uncombine( $refAlt )};
			@alternateAnnotations = sort keys %{$mutations->{$mutationKey2}->{$refAlt}};
			#print join( "," , ( @alternateAnnotations ) )."\n";
			my $reported = shift @alternateAnnotations;
			#print "\t".$reported."\n";
			$weight = $mutations->{$mutationKey2}->{$refAlt}->{$reported};
			( $reportedTranscript , $reportedAAChange ) = @{$this->splitProteinKey( $reported )};
			my $alternateAnnotations = join( "|" , @alternateAnnotations );
			#$fh->print( join( "\t" , ( $clusterID , $gene , $reportedAAChange , 
			push @linesToWrite , join( "\t" , ( $clusterID , $gene , $reportedAAChange , 

									   $degrees , $closenessCentrality , 
									   $geodesic , $weight ,
									   $chromosome , $start , $stop ,
									   $reference , $alternate ,
									   $reportedTranscript , $alternateAnnotations
									 )
							);
					  #);
			#$writtenLines += 1;
		} #foreach refAlt
		#print "deleting: ".$mutationKey2." and distances with centroid ".$centroid."\n";
		$this->setProcessStatus( $structure , $mutationKey2 , 1 );
	} #foreach other vertex in network
#	print "LINES TO WRITE = ".(scalar @linesToWrite)."\n";
	if ( scalar @linesToWrite > 1 ) {
		#$this->{'results'}->{$structure} .= join( "\n" , @linesToWrite );
		#$fh->print( $line );
		$writtenLines = scalar @linesToWrite;
	}
	return ( $writtenLines , \@linesToWrite );
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

## MUTATIONS
#sub printForTranscript {
#	my ( $this , $mutations , $mutationKey , $transcript ) = @_;
#	my $str = $mutationKey."=";
#	foreach my $ra ( sort keys %{$mutations->{$mutationKey}} ) {
#		foreach my $pk ( sort keys %{$mutations->{$mutationKey}->{$ra}} ) {
#			if ( $pk =~ /$transcript/ ) {
#				$str .= $ra.":".join( "|" , ( sort keys %{$mutations->{$mutationKey}->{$ra}} ) );
#			}
#		}
#		$str .= ";";
#	}
#	print $str;
#}

sub isSameProteinPosition { # shared
	my ( $this , $mutations , $mutationKey1 , $mutationKey2 ) = @_;
	#print join( "\t" , ( "begin" , $mutationKey1 , $mutationKey2 ) )."\n";
	if ( $mutationKey1 eq $mutationKey2 ) { return 1; }
	foreach my $refAlt1 ( sort keys %{$mutations->{$mutationKey1}} ) {
		foreach my $proteinKey1 ( sort keys %{$mutations->{$mutationKey1}->{$refAlt1}} ) {
			my ( $transcript1 , $aaChange1 ) = @{$this->splitProteinKey( $proteinKey1 )};
			my ( $aaReference1 , $aaPosition1 , $aaAlternate1 );
			#print "\tproteinKey1: ".$proteinKey1;
			if ( $aaChange1 =~ m/p\.\D\D*(\d+)\D*/ ) {
				$aaPosition1 =  $1;
			} else {
				print "...next, no match aaChange1\n";
				next;
			}
			foreach my $refAlt2 ( sort keys %{$mutations->{$mutationKey2}} ) {
#TODO make sure this works for in_frame_ins
				foreach my $proteinKey2 ( sort keys %{$mutations->{$mutationKey2}->{$refAlt2}} ) {
					my ( undef , $aaChange2 ) = @{$this->splitProteinKey( $proteinKey2 )};
					my ( $aaReference2 , $aaPosition2 , $aaAlternate2 );
					#print "\tproteinKey2: ".$proteinKey2."\t";
					if ( $aaChange2 =~ m/p\.\D\D*(\d+)\D*/ ) {
						$aaPosition2 =  $1;
					} else {
						print "...next, no match aaChange2\n";
						next;
					}
					next if ( $proteinKey2 !~ /$transcript1/ );
					#if ( $this->checkProteinPosition( $proteinKey1 , $proteinKey2 ) ) {
					#	return 1;
					#}
					if ( $aaPosition1 eq $aaPosition2 ) {
						#$this->printForTranscript( $mutations , $mutationKey1 , $transcript1 );
						#print "\n";
						#$this->printForTranscript( $mutations , $mutationKey2 , $transcript1 );
						#print join( "  " , ( $mutationKey1 , $refAlt1 , $proteinKey1 , $mutationKey2 , $refAlt2 , $proteinKey2 ) );
						#print "\n^--same aaPosition\n";#diff = ".$diff."\n";
						return 1;
					} else {
						next;
					}
					if ( $this->checkGenomicPositionNearby( $mutationKey1 , $mutationKey2 ) ) {
#TODO may fail at splice sites
						return 1;
					}
				} #foreach proteinKey2
			} #foreach refAlt2
			#print "\n";
		} #foreach proteinKey1
	} #foreach refAlt1
	#print "not the same site : ".join( "   " , ( $mutationKey1 , $mutationKey2 ) )."\n";
	return 0;
}

sub checkProteinPosition {
	my ( $this , $proteinKey1 , $proteinKey2 ) = @_;
	#print join( "\t" , ( $aaPosition1 , $aaReference1 , $aaPosition2 , $aaReference2 ) )."\t";
	return 0;
}

sub checkGenomicPositionNearby {
	my ( $this , $mutationKey1 , $mutationKey2 ) = @_;
	my ( undef , $chromosome1 , $start1 , undef ) = @{$this->splitMutationKey( $mutationKey1 )};
	next if ( !$start1 );
	my ( undef , $chromosome2 , $start2 , undef ) = @{$this->splitMutationKey( $mutationKey2 )};
	next if ( !$start2 );
	my $diff = 3;
	if ( $start1 <= $start2 ) {
		#print "ASDF21 ".$start2."\t".$start1."\n";
		$diff = $start2 - $start1;
	} else {
		#print "ASDF12 ".$start1."\t".$start2."\n";
		$diff = $start1 - $start2;
	}
	if ( $chromosome1 eq $chromosome2 
		 and $diff <= 2 ) {
#TODO make sure that the transcript details assure this works
		#print join( "  " , ( $mutationKey1 , $mutationKey2 ) );
		#print "<--same protein ref/alt ".$start2." - ".$start1." = ".$diff."\n";
		return 1;
	}
	return 0;
}

sub makeMutationKey {
	my ( $this , $mutation , $chain ) = @_;
	#my $proteinMutation = $mutation->proteinVariant( 0 );
	#my $mutationKey = &combine( $mutation->gene() , $proteinMutation->transcript() );
	#$mutationKey = &combine( $mutationKey , $proteinMutation->aminoAcidChange() );
	#$mutationKey = &combine( $mutationKey , $mutation->chromosome() );
	my $mutationKey = &combine( $mutation->gene() , $mutation->chromosome() );
	$mutationKey = &combine( $mutationKey , $mutation->start() );
	$mutationKey = &combine( $mutationKey , $mutation->stop() );
	#if ( $mutation->reference() ) {
	#	if ( $mutation->alternate() ) {
	#		$mutationKey = &combine( $mutationKey , $mutation->reference() );
	#		$mutationKey = &combine( $mutationKey , $mutation->alternate() );
	#	}
	#}
	return $mutationKey;
}

sub makeRefAltKey {
	my ( $this , $mutation ) = @_;
	return &combine( $mutation->reference() , $mutation->alternate() );
}

sub makePairKey {
	my ( $this , $mutation1 , $mutation2 , $chain1 , $chain2 ) = @_;
	my $mutationKey1 = $this->makeMutationKey( $mutation1 , $chain1 );
	my $mutationKey2 = $this->makeMutationKey( $mutation2 , $chain2 );
	return $mutationKey1."_".$mutationKey2;
}

sub makeProteinKey {
	my ( $this , $mutation ) = @_;
	#print "makeProteinKey: ";
	#print $this->makeMutationKey( $mutation ).": ";
	my $proteinVariant = $mutation->proteinVariant();
	my $proteinKey = &combine( $proteinVariant->transcript() , $proteinVariant->aminoAcidChange() );
	#print $proteinKey."\n";
	return $proteinKey;
}

sub splitMutationKey { # shared
	my ( $this , $mutationKey ) = @_;
	return &uncombine( $mutationKey );
}

sub splitRefAltKey { # shared
	my ( $this , $refAlt ) = @_;
	return &uncombine( $refAlt );
}

sub splitPairKey { 
	my ( $this , $pairKey ) = @_;
	return (split /_/ , $pairKey);
}

sub splitProteinKey { # shared
	my ( $this , $proteinKey ) = @_;
	my @split = @{&uncombine( $proteinKey )};
	return \@split;
}

sub getPairKeys {
	my ( $this , $mutation1 , $mutation2 ) = @_;
	my $pairKeyA = $this->makePairKey( $mutation1 , $mutation2 );
	my $pairKeyB = $this->makePairKey( $mutation2 , $mutation1 );
	return ( $pairKeyA , $pairKeyB );
}

## DISTANCE MATRIX
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

sub setShortestDistance {
	my ( $this , $distance_matrix , $mutation1 , $mutation2 , $chain1 , $chain2 , $infos ) = @_;

	my @infos = split /\|/ , $infos;
	my $nStructures = scalar @infos;
	print "HotSpot3D::Cluster::setShortestDistance\n";
	if ( $this->{'structure_dependence'} eq $DEPENDENT ) {
		foreach my $info ( @infos ) {
			chomp( $info );
			next if ( $info eq "" ); 
			my ( $distance , $pdbID , $pvalue ) = split / / , $info;
			#print "info ".$info." == ( ".join( " " , ( $distance , $pdbID , $pvalue ) )." )\n";
#TODO corresponding to update in setAverageDistance
			if ( $this->checkPair( $distance , $pvalue ) ) {
				my $structure = $this->makeStructureKey( $pdbID , $chain1 , $chain2 );
				my $oldDistance = $this->getElement( $distance_matrix , $structure , $mutation1 , $mutation2 , $chain1 , $chain2 );
				if ( $distance < $oldDistance ) {
					$this->setElement( $distance_matrix , $structure , $mutation1 , $mutation2 , $distance );
				}
			}
		}
	} else {
		my ( $distance , $pdbID , $pvalue ) = split / / , $infos[0];
#TODO corresponding to update in setAverageDistance
		if ( $this->checkPair( $distance , $pvalue ) ) {
			my $structure = $this->makeStructureKey( $pdbID , $chain1 , $chain2 );
			my $oldDistance = $this->getElement( $distance_matrix , $structure , $mutation1 , $mutation2 , $chain1 , $chain2 );
			if ( $distance < $oldDistance ) {
				$this->setElement( $distance_matrix , $structure , $mutation1 , $mutation2 , $distance , $chain1 , $chain2 );
			}
		}
	}
	return;
}

#TODO make sure set for average works as well as for shortest
sub setAverageDistance {
	my ( $this , $distance_matrix , $mutation1 , $mutation2 , $chain1 , $chain2 , $infos , $pdbCount ) = @_;
	my @infos = split /\|/ , $infos;
	my $sumDistances = {};
	my $nStructures = {};
	foreach my $info ( @infos ) {
		chomp( $info );
		next unless ( $info );
		my ( $distance , $pdbID , $pvalue ) = split / / , $info;
#TODO average over ALL pairs or just the ones satisfying conditions (leave checkPair here if latter)
		if ( $this->checkPair( $distance , $pvalue ) ) {
			my $structure = $this->makeStructureKey( $pdbID , $chain1 , $chain2 );
			$sumDistances->{$structure} += $distance;
			$nStructures->{$structure} += 1;
		}
	}
	$this->calculateAverageDistance( $distance_matrix , $mutation1 , $mutation2 , $pdbCount , $sumDistances , $nStructures );
	return;
}

sub calculateAverageDistance {
	my ( $this , $distance_matrix , $mutation1 , $mutation2 , $pdbCount , $sumDistances , $nStructures , $chain1 , $chain2 ) = @_;
	my $mutationKey1 = $this->makeMutationKey( $mutation1 , $chain1 );
	my $mutationKey2 = $this->makeMutationKey( $mutation2 , $chain2 );
	foreach my $structure ( keys %{$sumDistances} ) {
		if ( exists $pdbCount->{$structure}->{$mutationKey1} ) {
			if ( exists $pdbCount->{$structure}->{$mutationKey1}->{$mutationKey2} ) { #have seen the pair on a prior line
				my $count = $pdbCount->{$structure}->{$mutationKey1}->{$mutationKey2};
				my $oldDistance = $this->getElement( $distance_matrix , $structure , $mutation1 , $mutation2 , $chain1 , $chain2 );
				$sumDistances->{$structure} += $count*$oldDistance;
				$nStructures->{$structure} += $count;
			} else {
				$pdbCount->{$structure}->{$mutationKey1}->{$mutationKey2} = $nStructures;
			}
		} else {
			$pdbCount->{$structure}->{$mutationKey1}->{$mutationKey2} = $nStructures;
		}
		my $finalDistance = $sumDistances->{$structure} / $nStructures->{$structure};

		$this->setElement( $distance_matrix , $structure , $mutation1 , $mutation2 , $finalDistance , $chain1 , $chain2 );
	}
	return;
}

sub setDistance {
	my ( $this , $distance_matrix , $mutation1 , $mutation2 , 
		 $chain1 , $chain2 , $infos , $pdbCount ) = @_;
	if ( $this->{'distance_measure'} eq $AVERAGEDISTANCE ) {
		#print "get average distance\n";
		$this->setAverageDistance( $distance_matrix , $mutation1 , $mutation2 , 
								   $chain1 , $chain2 , $infos , $pdbCount );
	} else {
		#print "get shortest distance\n";
		$this->setShortestDistance( $distance_matrix , $mutation1 , 
									$mutation2 , $chain1 , $chain2 , $infos );
	}
	return;
}

sub setElement {
	my ( $this , $distance_matrix , $structure , $mutation1 , $mutation2 , $distance , $chain1 , $chain2 ) = @_;

#TODO corresponding to update in setAverageDistance
	#if ( $this->checkPair( $distance , $pvalue ) ) { #meets desired significance
		my $mutationKey1 = $this->makeMutationKey( $mutation1 );
		my $mutationKey2 = $this->makeMutationKey( $mutation2 );

		$distance_matrix->{$structure}->{$mutationKey1}->{$mutationKey2} = $distance;
		$distance_matrix->{$structure}->{$mutationKey2}->{$mutationKey1} = $distance;

		$this->setProcessStatus( $structure , $mutationKey1 , 0 );
		$this->setProcessStatus( $structure , $mutationKey2 , 0 );
	#}
	return;
}

sub initializeSameSiteDistancesToZero {
	my ( $this , $distance_matrix , $mutations ) = @_;
	foreach my $structure ( keys %{$distance_matrix} ) {
		foreach my $mutationKey1 ( keys %{$distance_matrix->{$structure}} ) {
			foreach my $mutationKey2 ( keys %{$distance_matrix->{$structure}} ) {
#TODO if mutationKeys are equal, but protein change is different
				next if ( $mutationKey1 eq $mutationKey2 );
				if ( $this->isSameProteinPosition( $mutationKey1 , $mutationKey2 ) ) {
					$distance_matrix->{$structure}->{$mutationKey1}->{$mutationKey2} = 0;
					$distance_matrix->{$structure}->{$mutationKey2}->{$mutationKey1} = 0;
				}
			}
		}
	}
}

sub getElement {
	my ( $this , $distance_matrix , $structure , $mutation1 , $mutation2 , $chain1 , $chain2 ) = @_;
	my $mutationKey1 = $this->makeMutationKey( $mutation1 , $chain1 );
	my $mutationKey2 = $this->makeMutationKey( $mutation2 , $chain2 );

	return ( $this->getElementByKeys( $distance_matrix , $structure , $mutationKey1 , $mutationKey2 ) );
}

sub getElementByKeys {
	my ( $this , $distance_matrix , $structure , $mutationKey1 , $mutationKey2 ) = @_;
	my $distance = $MAXDISTANCE;
	if ( exists $distance_matrix->{$structure}->{$mutationKey1} ) {
		if ( exists $distance_matrix->{$structure}->{$mutationKey1}->{$mutationKey2} ) {
			$distance = $distance_matrix->{$structure}->{$mutationKey1}->{$mutationKey2};
		}
	} elsif ( exists $distance_matrix->{$structure}->{$mutationKey2} ) {
		warn "ASYMMETRIC DISTANCE MATRIX ".$structure." ".$mutationKey1." ".$mutationKey2."\n";
		if ( exists $distance_matrix->{$structure}->{$mutationKey2}->{$mutationKey1} ) {
			$distance = $distance_matrix->{$structure}->{$mutationKey2}->{$mutationKey1};
		}
	}
	return $distance;
}

sub makeStructureKey {
	my ( $this , $structure , $chain1 , $chain2 ) = @_;
	my @chains = sort( ( $chain1 , $chain2 ) );
	#print "makeStructureKey: ".join( "  " , ( $structure , $chain1 , $chain2 ) )."\n";
	if ( $this->{'structure_dependence'} eq $DEPENDENT ) {
		if ( $this->{'subunit_dependence'} eq $DEPENDENT ) {
			#print "\t".&combine( $structure , &combine( $chains[0] , $chains[1] ) )."\n";
			return &combine( $structure , &combine( $chains[0] , $chains[1] ) );
		} else {
			#print "\t".$structure."\n";
			return $structure;
		}
	} else {
		if ( $this->{'subunit_dependence'} eq $DEPENDENT ) {
			#print "\t".&combine( $structure , &combine( $chains[0] , $chains[1] ) )."\n";
			return &combine( $structure , &combine( $chains[0] , $chains[1] ) );
		} else {
			#print "\t".$ANY."\n";
			return $ANY;
		}
	}
	#print "\tbackup ".$ANY."\n";
	return $ANY;
}

sub setSameSiteDistancesToZero { # used in density (this shuld be the same as initializeSameSiteDistancesToZero; but there was a problem when it was used)
	my ( $this , $distance_matrix , $mutations ) = @_;

	foreach my $structure ( sort keys %{$distance_matrix} ) {
		#print "\t".$structure."\n";
		foreach my $mutationKey1 ( sort keys %{$distance_matrix->{$structure}} ) {
			foreach my $mutationKey2 ( sort keys %{$distance_matrix->{$structure}} ) { #->{$mutationKey1}} ) {
				if ( $this->isSameProteinPosition( $mutations , $mutationKey1 , $mutationKey2 ) == 1 && $mutationKey1 ne $mutationKey2 ) {
					$distance_matrix->{$structure}->{$mutationKey1}->{$mutationKey2} = 0;
					$distance_matrix->{$structure}->{$mutationKey2}->{$mutationKey1} = 0;
				}
			}
		}
	}
}

## MISCELLANEOUS METHODS
sub numSort {
       my ( $list ) = @_;
       if ( scalar @{$list} > 0 ) {
               @{$list} = sort {$a <=> $b} @{$list};
       }
       return;
}

sub combine {
       my ( $a , $b ) = @_;
       return join( ":" , ( $a , $b ) );
}

sub uncombine {
       my $a = shift;
       my @split = split( /\:/ , $a );
       return \@split;
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
--distance-measure           Pair distance to use (shortest or average), default: average
--structure-dependence       Clusters for each structure or across all structures (dependent or independent), default: independent 

--help                       this message

HELP
}

sub help_text {
    my $this = shift;
        return <<HELP

Usage: hotspot3d cluster [options]

                             REQUIRED
--pairwise-file              3D pairwise data file
--maf-file                   .maf file used in proximity search step

                             OPTIONAL
--drug-clean-file            Either (or concatenated) drugs.target.clean & drugs.nontarget.clean data
--output-prefix              Output prefix, default: 3D_Proximity
--p-value-cutoff             P_value cutoff (<), default: 0.05 (if 3d-distance-cutoff also not set)
--3d-distance-cutoff         3D distance cutoff (<), default: 100 (if p-value-cutoff also not set)
--linear-cutoff              Linear distance cutoff (> peptides), default: 0
--max-radius                 Maximum cluster radius (max network geodesic from centroid, <= Angstroms), default: 10
--clustering                 Cluster using network or density-based methods (network or density), default: network
--vertex-type                Graph vertex type for network-based clustering (recurrence, unique, or weight), default: recurrence
--distance-measure           Pair distance to use (shortest or average), default: average
--structure-dependence       Clusters for each structure or across all structures (dependent or independent), default: independent
--subunit-dependence         Clusters for each subunit or across all subunits (dependent or independent), default: independent
--meric-type                 Clusters for each intra-molecular (both monomers and homomers), monomer, homomer, 
                                 inter-molecular (heteromers), heteromer, multimer (simultaneously homomer & heteromer), or any *mer 
                                 (intra, monomer, homomer, inter, heteromer, multimer, or unspecified), default: unspecified
--transcript-id-header       .maf file column header for transcript id's, default: transcript_name
--amino-acid-header          .maf file column header for amino acid changes, default: amino_acid_change 
--weight-header              .maf file column header for mutation weight, default: weight (used if vertex-type = weight)
--parallel                   Parallelization for structure and subunit dependent runs (none or local), default: none
--max-processes              Set if using parallel type local (CAUTION: make sure you know your max CPU processes)

--help                       this message

HELP

}

1;
