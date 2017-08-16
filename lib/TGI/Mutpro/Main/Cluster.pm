package TGI::Mutpro::Main::Cluster;
#
#----------------------------------
# $Authors: Adam Scott, Sohini Sengupta, & Amila Weerasinghe
# $Date: 2014-01-14 14:34:50 -0500 (Tue Jan 14 14:34:50 CST 2014) $
# $Revision: 6 $
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
use TGI::Mutpro::Main::Network;
use TGI::Mutpro::Main::Density;

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
my $NULL = "-";
my $PTM = "ptm";

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
    $this->{'maf_file'} = undef;
    $this->{'site_file'} = undef;
    $this->{'pairwise_file'} = undef;
    $this->{'drug_clean_file'} = undef;
    $this->{'sites_file'} = undef;
    $this->{'musites_file'} = undef;
    $this->{'output_prefix'} = undef;
    $this->{'p_value_cutoff'} = undef;
    $this->{'3d_distance_cutoff'} = undef;
    $this->{'linear_cutoff'} = 0;
	$this->{'max_radius'} = 10;
	$this->{'vertex_type'} = $SITE;
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

    $this->{'gene_list_file'} = undef;
    $this->{'structure_list_file'} = undef;
    $this->{'listOption'} = undef;
	
	$this->{'processed'} = undef;
	$this->{'distance_matrix'} = undef;
	$this->{'mutations'} = undef;
	$this->{'results'} = undef;

    $this->{'Epsilon'} = undef;
    $this->{'MinPts'} = undef;
    $this->{'number_of_runs'} = undef;
    $this->{'probability_cut_off'} = undef;

	$this->{'JSON_status'} = undef;
	$this->{'mutations_json_hash'} = undef;
	$this->{'distance_matrix_json_hash'} = undef;

    bless $this, $class;
    $this->process();
    return $this;
}

sub process {
    my $this = shift;
	$this->setOptions();
	$this->launchClustering(); 

##TODO bring the following commented out region back, and find a way to carry $distance_matrix and $mutations to clustering modules
	# my $temp_distance_matrix = {};
 # 	my $temp_mutations = {};
	# my $distance_matrix = {};
 # 	my $mutations = {};
	# my $WEIGHT = "weight";

	# $this->readMAF( $temp_mutations );
	# $this->getDrugMutationPairs( $temp_distance_matrix );
	# $this->getMutationMutationPairs( $temp_distance_matrix );
	# $this->vertexFilter( $temp_mutations , $temp_distance_matrix , $mutations , $distance_matrix );
	# $this->networkClustering( $mutations , $distance_matrix );

    return 1;
}

sub mafFile {
	my ( $this , $file ) = @_;
	return ( $this->getFile( 'maf_file' , $file ) );
}

sub geneListFile {
	my ( $this , $file ) = @_;
	return ( $this->getFile( 'gene_list_file' , $file ) );
}

sub structureListFile {
	my ( $this , $file ) = @_;
	return ( $this->getFile( 'structure_list_file' , $file ) );
}

sub siteListFile {
	my ( $this , $file ) = @_;
	return ( $this->getFile( 'site_file' , $file ) );
}

sub pairwiseFile {
	my ( $this , $file ) = @_;
	return ( $this->getFile( 'pairwise_file' , $file ) );
}

sub drugsCleanFile {
	my ( $this , $file ) = @_;
	return ( $this->getFile( 'drug_clean_file' , $file ) );
}

sub sitePairsFile {
	my ( $this , $file ) = @_;
	return ( $this->getFile( 'sites_file' , $file ) );
}

sub musitePairsFile {
	my ( $this , $file ) = @_;
	return ( $this->getFile( 'musites_file' , $file ) );
}

sub getFile {
	my ( $this , $name , $file ) = @_;
	if ( $file ) { $this->{$name} = $file; }
	if ( not defined $this->{$name} ) { 
		return undef; 
	}
	return $this->{$name};
}

sub printMutations {
	my ( $this , $mutations , $type ) = @_;
	foreach my $a ( sort keys %{$mutations} ) {
		foreach my $b ( sort keys %{$mutations->{$a}} ) {
			foreach my $c ( sort keys %{$mutations->{$a}->{$b}} ) {
				print join( "  " , ( $type , $a , $b , $c , $mutations->{$a}->{$b}->{$c} ) )."\n";
			}
		}
	}
	print scalar( keys %{$mutations})." of mu\n";
}

sub printDistanceMatrix {
	my ( $this , $distance_matrix , $type ) = @_;
	foreach my $a ( sort keys %{$distance_matrix} ) {
		foreach my $b ( sort keys %{$distance_matrix->{$a}} ) {
			foreach my $c ( sort keys %{$distance_matrix->{$a}->{$b}} ) {
				print join( "  " , ( $type , $a , $b , $c , $distance_matrix->{$a}->{$b}->{$c} ) )."\n";
			}
		}
		print scalar( keys %{$distance_matrix->{$a}})." of dm\n";
		last;
	}
	print scalar( keys %{$distance_matrix})." of dms\n";
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
        'sites-file=s' => \$this->{'sites_file'},
        'musites-file=s' => \$this->{'musites_file'},
        'p-value-cutoff=f' => \$this->{'p_value_cutoff'},
        '3d-distance-cutoff=f' => \$this->{'3d_distance_cutoff'},
        'linear-cutoff=f' => \$this->{'linear_cutoff'},
        'max-radius=f' => \$this->{'max_radius'},
        'vertex-type=s' => \$this->{'vertex_type'},
        'number-of-runs=f' => \$this->{'number_of_runs'},
        'probability-cut-off=f' => \$this->{'probability_cut_off'},
        'distance-measure=s' => \$this->{'distance_measure'},
        'maf-file=s' => \$this->{'maf_file'},
        'site-file=s' => \$this->{'site_file'} , 
        'amino-acid-header=s' => \$this->{'amino_acid_header'},
        'transcript-id-header=s' => \$this->{'transcript_id_header'},
        'weight-header=s' => \$this->{'weight_header'},
        'clustering=s' => \$this->{'clustering'},
        'structure-dependent' => \$this->{'structure_dependence'},
        'subunit-dependent' => \$this->{'subunit_dependence'},
        'parallel=s' => \$this->{'parallel'},
        'max-processes=i' => \$this->{'max_processes'},
        'meric-type=s' => \$this->{'meric_type'},
        'gene-list-file=s' => \$this->{'gene_list_file'},
        'structure-list-file=s' => \$this->{'structure_list_file'},
        'epsilon=f' => \$this->{'Epsilon'},
        'minPts=f' => \$this->{'MinPts'},
        'number-of-runs=f' => \$this->{'number_of_runs'},
        'probability-cut-off=f' => \$this->{'probability_cut_off'},
		'use-JSON' => \$this->{'JSON_status'},
		'mutations-hash-json-file=s' => \$this->{'mutations_json_hash'},
		'distance-matrix-json-file=s' => \$this->{'distance_matrix_json_hash'},

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
		warn "HotSpot3D::Cluster::setOptions warning: no structure-dependent option given, setting to default independent\n";
	} else {
		$this->{'structure_dependence'} = $DEPENDENT;
	}
	if ( not defined $this->{'subunit_dependence'} ) {
		$this->{'subunit_dependence'} = $INDEPENDENT;
		warn "HotSpot3D::Cluster::setOptions warning: no subunit-dependent option given, setting to default independent\n";
	} else {
		$this->{'subunit_dependence'} = $DEPENDENT;
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
	if ( defined $this->siteListFile() ) {
		if ( not -e $this->siteListFile() ) { 
			warn "The input site file (".$this->siteListFile().") does not exist! ", "\n";
			die $this->help_text();
		}
	} else {
		warn "HotSpot3D::Cluster::setOptions warning: no site-file included (cannot produce site-site clusters)!\n";
	}
	if ( defined $this->sitePairsFile() ) {
		#$this->requireSite();
		if ( not -e $this->sitePairsFile() ) { 
			warn "The input site-site pair file (".$this->sitePairsFile().") does not exist! ", "\n";
			die $this->help_text();
		}
	} else {
		warn "HotSpot3D::Cluster::setOptions warning: no sites-file included (cannot produce site-site clusters)!\n";
	}
	if ( defined $this->musitePairsFile() ) {
		#$this->requireSite();
		$this->requireMAF();
		if ( not -e $this->musitePairsFile() ) { 
			warn "The input mutation-site pair file (".$this->musitePairsFile().") does not exist! ", "\n";
			die $this->help_text();
		}
	} else {
		warn "HotSpot3D::Cluster::setOptions warning: no musites-file included (cannot produce mutation-site clusters)!\n";
	}
	if ( defined $this->drugsCleanFile() ) {
		$this->requireMAF();
		if ( not -e $this->drugsCleanFile() ) { 
			warn "The input drug pairs file (".$this->drugsCleanFile().") does not exist! ", "\n";
			die $this->help_text();
		}
	} else {
		warn "HotSpot3D::Cluster::setOptions warning: no drug-clean-file included (cannot produce drug-mutation clusters)!\n";
	}
    if ( defined $this->pairwiseFile() ) {
		$this->requireMAF();
		if ( not -e $this->pairwiseFile() ) { 
			warn "HotSpot3D::Cluster::setOptions error: the input pairwise file (".$this->pairwiseFile().") does not exist!\n";
			die $this->help_text();
		}
	} else {
		warn "HotSpot3D::Cluster::setOptions warning: no pairwise-file included (cannot produce mutation-mutation clusters)!\n";
	}
	if ( defined $this->{'JSON_status'} ) {
		$this->{'JSON_status'} = 1;
		warn "HotSpot3D::Cluster::setOptions warning: use-JSON flag used (will not look for pairwise data or maf file)\n";
		if ( not defined $this->{'mutations_json_hash'} or not defined $this->{'distance_matrix_json_hash'}  ) {
			die "HotSpot3D::Cluster::setOptions Error: use-JSON flag is used, but json file locations are not provided!\n";
		}
		elsif ( not -e $this->{'mutations_json_hash'} or not -e $this->{'mutations_json_hash'} ) {	
			die "HotSpot3D::Cluster::setOptions Error: use-JSON flag is used, but the provided JSON files do not exist!\n";
		}
	}
	else { $this->{'JSON_status'} = 0; }

	if (  not defined $this->pairwiseFile() and
		  not defined $this->sitePairsFile() and
		  not defined $this->musitePairsFile() and
		  not defined $this->drugsCleanFile() and not $this->{'JSON_status'} ) {
		warn "HotSpot3D::Cluster::setOptions error: no pair file provided. Need at least one of *.pairwise, *.clean, *.sites, *.musites.\n";
		die $this->help_text();
	}
	my $acceptableVertex = { $RECURRENCE => 1 , $UNIQUE => 1 , $SITE => 1 , $WEIGHT => 1 };
	if ( not defined $this->{'vertex_type'} ) {
		warn "vertex-type not specified. Using default vertex-type = \'recurrence\'\n";
		$this->{'vertex_type'} = $RECURRENCE;
	} elsif ( not exists $acceptableVertex->{$this->{'vertex_type'}} ) {
		warn "vertex-type option not recognized as \'recurrence\', \'unique\', \'site\', or \'weight\'\n";
		die $this->help_text();
	}
	if ( $this->{'distance_measure'} ne $AVERAGEDISTANCE
		 and $this->{'distance_measure'} ne $SHORTESTDISTANCE ) {
		warn "distance-measure option not recognized as \'shortest\' or \'average\'\n";
		warn "Using default distance-measure = \'average\'\n";
		$this->{'distance_measure'} = $AVERAGEDISTANCE;
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

	##### START gene-list and structure-list options
	if ( defined $this->geneListFile() ) {
		if ( not -e $this->geneListFile() ) { 
			warn "The input gene list file (".$this->geneListFile().") does not exist! ", "\n";
			die $this->help_text();
		}		
	}

	if ( defined $this->structureListFile() ) {
		if ( not -e $this->structureListFile() ) { 
			warn "The input structure list file (".$this->structureListFile().") does not exist! ", "\n";
			die $this->help_text();
		}
	}

	$this->initializeLists();
	##### END

	print STDOUT "=====Parameters=====\n";
	print STDOUT " linear-cutoff         = ".$this->{'linear_cutoff'}."\n";
	print STDOUT " p-value-cutoff        = ".$this->{'p_value_cutoff'}."\n";
	print STDOUT " 3d-distance-cutoff    = ".$this->{'3d_distance_cutoff'}."\n";
	print STDOUT " max-radius            = ".$this->{'max_radius'}."\n";
	print STDOUT " vertex-type           = ".$this->{'vertex_type'}."\n";
	print STDOUT " distance-measure      = ".$this->{'distance_measure'}."\n";
	print STDOUT " structure-dependent   = ".$this->{'structure_dependence'}."\n";
	print STDOUT " subunit-dependent     = ".$this->{'subunit_dependence'}."\n";
	print STDOUT " meric-type            = ".$this->{'meric_type'}."\n";
	print STDOUT " parallel              = ".$this->{'parallel'}."\n";
	if ( defined $this->{'max_processes'} ) {
		print STDOUT " max-processes         = ".$this->{'max_processes'}."\n";
	}
	print STDOUT "====================\n";
	print STDOUT "=====Data===========\n";
	print STDOUT " maf-file              = ".$this->mafFile()."\n";
	print STDOUT " transcript_id_header  = ".$this->{'transcript_id_header'}."\n";
	print STDOUT " amino-acid-header     = ".$this->{'amino_acid_header'}."\n";
	print STDOUT " weight-header         = ".$this->{'weight_header'}."\n";
	print STDOUT " site-file             = ".$this->siteListFile()."\n"; 
	print STDOUT " gene-list-file        = ".$this->geneListFile()."\n";
	print STDOUT " structure-list-file   = ".$this->structureListFile()."\n";
	print STDOUT " pairwise-file         = ".$this->pairwiseFile()."\n";
	print STDOUT " drug-clean-file       = ".$this->drugsCleanFile()."\n";
	print STDOUT " sites-file            = ".$this->sitePairsFile()."\n";
	print STDOUT " musites-file          = ".$this->musitePairsFile()."\n";
	print STDOUT " output-prefix         = ".$this->{'output_prefix'}."\n";
	print STDOUT "====================\n";

	return;
}

sub requireSite {
	my $this = shift;
	unless( $this->siteListFile() ) {
		warn "Provided sites-file requires a site-file, but no site-file was given\n";
		die $this->help_text();
	}
	unless ( -e $this->siteListFile() ) {
		warn "The site-file, ".$this->siteListFile().", does not exist!\n";
		die $this->help_text();
	}
	return;
}

sub requireMAF {
	my $this = shift;
	unless( $this->mafFile() ) {
		warn "You must provide a maf-file for the pair data provided.\n";
		die $this->help_text();
	}
	unless( -e $this->mafFile() ) {
		warn "The maf-file, ".$this->mafFile().", does not exist! ", "\n";
		die $this->help_text();
	}
	return;
}

sub initializeLists {
	my $this = shift;
	#print "at initializeLists\n";

	if ( defined $this->geneListFile() ) {
		my $fh = new FileHandle;
		unless( $fh->open( $this->geneListFile() , "r" ) ) { die "Could not open gene list file $! \n" };
		while ( my $gene = <$fh> ) {
    		chomp( $gene );
    		$this->{'providedGenes'}->{$gene} = 0; # hash with provided gene names as keys, will be used later in filterByProvidedLists
    	}
	}
	if ( defined $this->structureListFile() ) {
		my $fh = new FileHandle;
		unless( $fh->open( $this->structureListFile() , "r" ) ) { die "Could not open structures list file $! \n" };
		while ( my $structure = <$fh> ) {
    		chomp( $structure );
    		$this->{'providedStructures'}->{$structure} = 0; # hash with provided structure names as keys, will be used later in filterByProvidedLists
    	}
	}
	##
	if ( ( defined $this->geneListFile() ) && ( not defined $this->structureListFile() ) ) {
		$this->{'listOption'} = "GeneListOnly";
		#print "gene list only\n";
	}
	if ( ( not defined $this->geneListFile() ) && ( defined $this->structureListFile() ) ) {
		$this->{'listOption'} = "StructureListOnly";
		warn "HotSpot3D::Cluster::setOptions warning: only the mutations present in the given structures are being considered\n";
		#print "structure list only\n";
	}
	if ( ( defined $this->geneListFile() ) && ( defined $this->structureListFile() ) ) {
		$this->{'listOption'} = "BothLists";
		#print "both lists\n";
	}
	if ( ( not defined $this->geneListFile() ) && ( not defined $this->structureListFile() ) ) {
		$this->{'listOption'} = "None";
		#print "no lists\n";
	}
}

sub launchClustering {
	my $this = shift; # removed my ( $this, $help ) = @_;
	if ( $this->{'clustering'} eq $NETWORK ) {
		# if ( $help ) {
		# 	die $this->help_text();
		# } else {
			TGI::Mutpro::Main::Network->new( $this );
			exit;
		# }
	} elsif ( $this->{'clustering'} eq $DENSITY ) {
		# if ( $help ) {
		# 	die $this->density_help_text();
		# } else{
			TGI::Mutpro::Main::Density->new($this);
			exit;
		# }
	}
	return;
}

sub vertexFilter {
#$this->vertexFilter( $temp_mutations , $temp_distance_matrix , $mutations , $distance_matrix );
	my ( $this , $temp_mutations , $temp_distance_matrix , $mutations , $distance_matrix ) = @_;
	if ( $this->{'vertex_type'} eq $SITE ) {
		print STDOUT "Filtering vertices\n";
#TODO if using a different .maf from search step, then some mutations can be missed
		my $vertexMap = {}; #a hash to map isSameProteinPosition vertices (and others to their selves)-- map=f()
		my @mKeys = sort keys %{$temp_mutations}; #an array to store all the mutation keys

		#filter same vertices and generate a map
		foreach my $mutationKey1 ( @mKeys ) {
			next if not exists $temp_mutations->{$mutationKey1};
			foreach my $mutationKey2 ( @mKeys ) {
				next if not exists $temp_mutations->{$mutationKey2};
				if ( $mutationKey1 eq $mutationKey2 ) {
					$vertexMap->{$mutationKey2} = $mutationKey1;
					# print "ACSW::VertexFilter::Equal $mutationKey2 \=\=\> $mutationKey1\n";
					next;
				}
				elsif ( $this->isSameProteinPosition( $temp_mutations , $mutationKey1 , $mutationKey2 ) ) { #if same site
					$vertexMap->{$mutationKey2} = $mutationKey1;
					print "ACSW::VertexFilter::SameSite $mutationKey2 \=\=\> $mutationKey1\n";
					delete $temp_mutations->{$mutationKey2};
				}
			}
		}
		print "ACSW::VertexFilter::mutations filtering done\n";

		#generate representative annotations
		foreach my $mutationKey ( keys %{$temp_mutations} ) {
			my ( $ra , $pk ) = $this->getARepresentativeAnnotation( $temp_mutations , $mutationKey );
			$mutations->{$mutationKey}->{$ra}->{$pk} = 1;
		}
		print "ACSW::VertexFilter::mutations representative annotation done\n";

		#use mk'=f(mk) vertex map in the distance matrix
		foreach my $structure ( keys %{$temp_distance_matrix} ) {		
			foreach my $mutationKey1 ( keys %{$temp_distance_matrix->{$structure}} ) {
				my $mutationKey1prime = $vertexMap->{$mutationKey1}; # mk1'=f(mk1)
				foreach my $mutationKey2 ( keys %{$temp_distance_matrix->{$structure}->{$mutationKey1}} ) {
					my $mutationKey2prime = $vertexMap->{$mutationKey2}; # mk2'=f(mk2)
					if ( not defined $mutationKey1prime and not defined $mutationKey2prime ) {
						print "ACSW::VertexFilter::Error::undefined prime keys $mutationKey1 \= $mutationKey1prime and $mutationKey2 \= $mutationKey2prime\n";
					}
					if ( exists $distance_matrix->{$structure}->{$mutationKey1prime}->{$mutationKey2prime} ) {
						if ( $distance_matrix->{$structure}->{$mutationKey1prime}->{$mutationKey2prime} ne $temp_distance_matrix->{$structure}->{$mutationKey1}->{$mutationKey2}) {
							print "ACSW::VertexFilter::Error::non-matching distances $mutationKey1 \= $mutationKey1prime and $mutationKey2 \= $mutationKey2prime, d\($mutationKey1,$mutationKey2\) ne d\($mutationKey1prime,$mutationKey2prime\)\n";
						}
					}
					else {
						$distance_matrix->{$structure}->{$mutationKey1prime}->{$mutationKey2prime} = $temp_distance_matrix->{$structure}->{$mutationKey1}->{$mutationKey2};
					}
				}
			}
		}
	} else {
		%{$mutations} = %{$temp_mutations};
		%{$distance_matrix} = %{$temp_distance_matrix};
	}
	$temp_mutations = undef;
	$temp_distance_matrix = undef;
	# print "distance_matrix\n";
	# print Dumper $distance_matrix;
	return;
}

sub getARepresentativeAnnotation {
	my ( $this , $mutations , $mutationKey ) = @_;
	my $ra = ".:.";
	my $pk = ".:p.";
	foreach my $refAlt ( keys %{$mutations->{$mutationKey}} ) {
		$ra = $refAlt;
		foreach my $proteinKey ( keys %{$mutations->{$mutationKey}->{$refAlt}} ) {
			$pk = $proteinKey;
			last;
		}
		last;
	}
	return ( $ra , $pk );
}

sub checkPartners {
	my ( $this , $distance_matrix , $structure , $mutationKey1 , $mutationKey2 , $temp_distance_matrix , $temp_mutations , $mutations ) = @_;
	my $storeIt = 1;
	foreach my $partner ( keys %{$distance_matrix->{$structure}->{$mutationKey1}} ) { #foreach known partner
		if ( $this->isSameProteinPosition( $temp_mutations , $mutationKey2 , $partner ) ) { #is the partner at the same site as our new key
			$storeIt = 0;
		}
	}
	if ( $storeIt ) { #keep this pair
		$this->setRepresentative( $distance_matrix , $structure , $mutationKey1 , $mutationKey2 , $temp_distance_matrix , $temp_mutations , $mutations );
	}
}

sub setRepresentative {
	my ( $this , $distance_matrix , $structure , $mutationKey1 , $mutationKey2 , $temp_distance_matrix , $temp_mutations , $mutations ) = @_;
	$distance_matrix->{$structure}->{$mutationKey1}->{$mutationKey2} = $temp_distance_matrix->{$structure}->{$mutationKey1}->{$mutationKey2};
	foreach my $refAlt ( keys %{$temp_mutations->{$mutationKey1}} ) { #get A refAlt for representation
		foreach my $proteinKey ( sort keys %{$temp_mutations->{$mutationKey1}->{$refAlt}} ) { #get A protein key for representation
			$mutations->{$mutationKey1}->{$refAlt}->{$proteinKey} = $temp_mutations->{$mutationKey1}->{$refAlt}->{$proteinKey};
			print join( "\t" , ( "setRep" , $mutationKey1 , $refAlt , $proteinKey , $mutationKey2 ) )."\n";
			last; #only need one
		}
		last; #only need one
	}
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
	return if ( not defined $this->pairwiseFile() );
	print STDOUT "HotSpot3D::Cluster::readPairwise\n";
	my $fh = new FileHandle;
	unless( $fh->open( $this->pairwiseFile() , "r" ) ) { die "Could not open pairwise file $! \n" };
	my $pi = 0;
	my $headline = $fh->getline(); chomp( $headline );
	my %pairwisecols = map{ ( $_ , $pi++ ) } split( /\t/ , $headline );
	my @required = ( "Gene1" , "Chromosome1" , "Start1" , "Stop1" ,
					 "Mutation1" , "Chain1" , "Position1" , 
					 "Gene2" , "Chromosome2" , "Start2" , "Stop2" ,
					 "Mutation2" , "Chain2" , "Position2" , 
					 "LinearDistance" , "DistanceInfo" );
	unless(	defined( $pairwisecols{"Gene1"} )
		and defined( $pairwisecols{"Chromosome1"} )
		and defined( $pairwisecols{"Start1"} )
		and defined( $pairwisecols{"Stop1"} )
		and defined( $pairwisecols{"Mutation1"} )
		and defined( $pairwisecols{"Chain1"} )
		and defined( $pairwisecols{"Position1"} )
		and defined( $pairwisecols{"Gene2"} )
		and defined( $pairwisecols{"Chromosome2"} )
		and defined( $pairwisecols{"Start2"} )
		and defined( $pairwisecols{"Stop2"} )
		and defined( $pairwisecols{"Mutation2"} )
		and defined( $pairwisecols{"Chain2"} )
		and defined( $pairwisecols{"Position2"} )
		and defined( $pairwisecols{"LinearDistance"} )
		and defined( $pairwisecols{"DistanceInfo"} ) ) {
		die "not a valid pairwise file\n";
	}
	my @wantpairwisecols = (	$pairwisecols{"Gene1"} ,
						$pairwisecols{"Chromosome1"} ,
						$pairwisecols{"Start1"} ,
						$pairwisecols{"Stop1"} ,
						$pairwisecols{"Mutation1"} ,
						$pairwisecols{"Chain1"} ,
						$pairwisecols{"Position1"} ,
						$pairwisecols{"Gene2"} ,
						$pairwisecols{"Chromosome2"} ,
						$pairwisecols{"Start2"} ,
						$pairwisecols{"Stop2"} ,
						$pairwisecols{"Mutation2"} ,
						$pairwisecols{"Chain2"} ,
						$pairwisecols{"Position2"} ,
						$pairwisecols{"LinearDistance"} ,
						$pairwisecols{"DistanceInfo"} );
	my $pdbCount;
	foreach my $line ( <$fh> ) {
		chomp( $line );
		my @line = split /\t/ , $line;
		my ( $gene1 , $chromosome1 , $start1 , $stop1 , $aa_1 , $chain1 , $loc_1 , 
			 $gene2 , $chromosome2 , $start2 , $stop2 , $aa_2 , $chain2 , $loc_2 , 
			 $linearDistance , $infos ) = @line[@wantpairwisecols];

		if ( $this->checkByProvidedLists( $gene1, $gene2, $infos ) ) { # check to filter out mutation pairs according to a given gene list or structures list or both
			if ( defined $this->structureListFile() ) { # if true, get only the pdbIDs which are in this list
				$infos = $this->filterInfosByGivenStructures( $infos );
				#print "new infos= $infos\n";
			}
			if ( $this->checkMeric( $gene1 , $gene2 , $chain1 , $chain2 ) ) {
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
		}
	}
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

sub checkByProvidedLists {
	my ( $this, $gene1, $gene2, $infos ) = @_;
	#print join( "  " , ( $gene1 , $gene2 , $infos , "" ) );
	if ( $this->{'listOption'} eq "None" ) { # no lists given, proceed as usual
		#print "no lists: process as usual\n";
		return 1;

	} elsif ( $this->{'listOption'} eq "GeneListOnly" ) {
		if ( (exists $this->{'providedGenes'}->{$gene1}) && (exists $this->{'providedGenes'}->{$gene2}) ) {
			#print "gene check okay\n";
			return 1;
		} else { return 0; }

	} elsif ( $this->{'listOption'} eq "StructureListOnly" ) {
		my $newInfos = $this->filterInfosByGivenStructures( $infos );
		if ( defined $newInfos ) {
			#print "structure check okay\n";
			return 1;
		}else { return 0; }

	} elsif ( $this->{'listOption'} eq "BothLists" ) {
		my $newInfos = filterInfosByGivenStructures( $this, $infos );
		if ( (defined $newInfos) && (exists $this->{'providedGenes'}->{$gene1}) && (exists $this->{'providedGenes'}->{$gene2}) ) {
			#print "both checks okay\n";
			return 1;
		} else { return 0; }
	}
}

sub filterInfosByGivenStructures { # takes the infos column from pairwise file and filter out the structures that are not present in the structures list
	my ( $this, $infos ) = @_;
	my $newInfos = undef;
	my @infos = split /\|/ , $infos;
	foreach my $info ( @infos ) {
		chomp( $info );
		next if ( $info eq "" ); 
		my ( $distance , $pdbID , $pvalue ) = split / / , $info;
		next if ( not exists $this->{'providedStructures'}->{$pdbID} );
		if ( defined $newInfos ) {
			$newInfos = join ( "\|", $newInfos, $info );
		} else {
			$newInfos = $info;
		}
	}
	return $newInfos;
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
	return if ( not defined $this->mafFile() );
	print STDOUT "HotSpot3D::Cluster::readMAF\n";
	my $fh = new FileHandle;
	die "Could not open .maf file\n" unless( $fh->open( $this->mafFile() , "r" ) );
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
	my $nMutations = 0;
	foreach my $line ( <$fh> ) {
		chomp( $line );
		my @line = split /\t/ , $line;
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
				$nMutations++;
			} #if mutation is missense or in frame
		} #if columns in range
	}
	$fh->close();
	print STDOUT $nMutations." mutations from .maf\n";
	return;
}

sub setMutation {
	my ( $this , $mutations , $mutation , $barID , $weight ) = @_;

	my $mutationKey = $this->makeMutationKey( $mutation , "" );
	my $refAlt = &combine( $mutation->reference() , $mutation->alternate() );
	my $proteinKey = $this->makeProteinKey( $mutation );
	#print "setMutation: ".$refAlt."\n";

	#print $this->printMutations( $mutations , "preset" );
	if ( $refAlt =~ m/$PTM/ ) {
		if ( exists $mutations->{$mutationKey}->{$refAlt}->{$proteinKey} ) {
			#print "PTM AGAIN\n";
			#print join( "\t" , ( $mutationKey , $refAlt , $proteinKey , $barID , $weight ) )."\n";
			return;
		}
	}
	if ( exists $mutations->{$mutationKey}->{$refAlt}->{$proteinKey} ) {
		#print "existing\n";
		if ( $this->{'vertex_type'} ne $WEIGHT ) {
			$mutations->{$mutationKey}->{$refAlt}->{$proteinKey} += 1;
		} else {
			$mutations->{$mutationKey}->{$refAlt}->{$proteinKey} = $weight;
		}
	} else {
		#print "new\n";
		if ( $this->{'vertex_type'} eq $WEIGHT ) {
			$mutations->{$mutationKey}->{$refAlt}->{$proteinKey} = $weight;
		} else {
			$mutations->{$mutationKey}->{$refAlt}->{$proteinKey} += 1;
		}
	} #if mutation exists
	#print $this->printMutations( $mutations , "postset" );
	#my $w = $mutations->{$mutationKey}->{$refAlt}->{$proteinKey};
	#print "ADDED: ".$w."\n";
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

sub writeClustersFile {
	my ( $this , $fh , $results ) = @_;
	#print STDOUT "HotSpot3D::Cluster::writeClustersFile\n";
	#print Dumper( $results );
	foreach my $line ( sort keys %{$results} ) {
		$fh->print( $line."\n" );
	}
	return;
	foreach my $structure ( %{$results} ) {
		if ( exists $results->{$structure} ) {
			my @mutations = uniq( @{$results->{$structure}} );
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
	return;
}

sub generateFilename {
	my $this = shift;
	my @outFilename;
	if ( $this->{'output_prefix'} ) {
		push @outFilename , $this->{'output_prefix'};
	} else {
		if ( $this->mafFile() ) {
			my $maf = basename( $this->mafFile() );
			push @outFilename , $maf;
		}
		if ( $this->pairwiseFile() ne "" ) {
			my $clean = basename( $this->pairwiseFile() );
			if ( $clean ne '' and scalar @outFilename > 1 ) {
				push @outFilename , $clean;
			} elsif ( $clean ne '' ) {
				push @outFilename , $clean;
			}
		}
		if ( $this->sitePairsFile() ne "" ) {
			my $clean = basename( $this->sitePairsFile() );
			if ( $clean ne '' and scalar @outFilename > 1 ) {
				push @outFilename , $clean;
			} elsif ( $clean ne '' ) {
				push @outFilename , $clean;
			}
		}
		if ( $this->musitePairsFile() ne "" ) {
			my $clean = basename( $this->musitePairsFile() );
			if ( $clean ne '' and scalar @outFilename > 1 ) {
				push @outFilename , $clean;
			} elsif ( $clean ne '' ) {
				push @outFilename , $clean;
			}
		}
		if ( $this->drugsCleanFile() ne "" ) {
			my $clean = basename( $this->drugsCleanFile() );
			if ( $clean ne '' and scalar @outFilename > 1 ) {
				push @outFilename , $clean;
			} elsif ( $clean ne '' ) {
				push @outFilename , $clean;
			}
		}
		push @outFilename , $this->{'vertex_type'};
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

sub readSite {
	my ( $this , $mutations ) = @_;
    return if ( not defined $this->siteListFile() ); #if musite pairs included
	print STDOUT "HotSpot3D::Cluster::readSite - ".$this->siteListFile()."\n";
	my $fh = new FileHandle;
	unless( $fh->open( $this->siteListFile() , "r" ) ) {
		die "Could not open site list file $! \n"
	};
	my $si = 0;
	my $headline = $fh->getline(); chomp( $headline );
	my %sitecols = map{ ( $_ , $si++ ) } split( /\t/ , $headline );
	my @required = ( "Hugo_Symbol" , $this->{'transcript_id_header'} , 
					 "Position" , "Feature" );
	unless(	defined( $sitecols{"Hugo_Symbol"} )
		and defined( $sitecols{$this->{'transcript_id_header'}} )
		and defined( $sitecols{"Position"} )
		and defined( $sitecols{"Feature"} ) ) {
		die "not a valid site file\n";
	}
	my @wantsitecols = ( $sitecols{"Hugo_Symbol"} ,
						$sitecols{$this->{'transcript_id_header'}} ,
						$sitecols{"Position"} ,
						$sitecols{"Feature"} );
	my $nSites = 0;
	foreach my $line ( <$fh> ) { 
		chomp( $line );
		my @line = split /\t/ , $line;
		my ( $gene , $transcript , $position , $feature ) = @line[@wantsitecols];
		my $proteinSite = TGI::ProteinVariant->new();
		my $site = TGI::Variant->new(); #for the purpose of making keys, will look like a mutation
		$site->gene( $gene );
		$site->chromosome( $transcript );
		$site->start( $position );
		$site->stop( $NULL );
		$site->reference( $NULL );
		$site->alternate( $PTM );
		$proteinSite->transcript( $transcript );
		$proteinSite->aminoAcidChange( $site );
		$site->addProteinVariant( $proteinSite );
		my $weight = 1;
		$this->setMutation( $mutations , $site , $feature , $weight );
		$nSites++;
	}
	$fh->close();
	print STDOUT $nSites." sites\n";
	return;
}

sub readMutationSitePairs {
	my ( $this , $mutations , $distance_matrix ) = @_;
    return if ( not defined $this->musitePairsFile() ); #if musite pairs included
	print STDOUT "HotSpot3D::Cluster::readMutationSitePairs\n";
	my $fh = new FileHandle;
	unless( $fh->open( $this->musitePairsFile() , "r" ) ) {
		die "Could not open musite pairs data file $! \n"
	};
	my $msi = 0;
	my $headline = $fh->getline(); chomp( $headline );
	my %musitecols = map{ ( $_ , $msi++ ) } split( /\t/ , $headline );
	my @required = ( "Gene1" , "Chromosome1" , "Start1" , "Stop1" , 
					 "Mutation1" , "Chain1" , "Position1" , 
					 "Gene2" , "Transcript2" , "TranscriptPosition2" , "Site2" ,
					 "Chain2" , "Position2" , "Feature2" , 
					 "LinearDistance" , "DistanceInfo" );
	unless(	defined( $musitecols{"Gene1"} )
		and defined( $musitecols{"Chromosome1"} )
		and defined( $musitecols{"Start1"} )
		and defined( $musitecols{"Stop1"} )
		and defined( $musitecols{"Mutation1"} )
		and defined( $musitecols{"Chain1"} )
		and defined( $musitecols{"Position1"} )
		and defined( $musitecols{"Gene2"} )
		and defined( $musitecols{"Transcript2"} )
		and defined( $musitecols{"TranscriptPosition2"} )
		and defined( $musitecols{"Site2"} )
		and defined( $musitecols{"Chain2"} )
		and defined( $musitecols{"Position2"} )
		and defined( $musitecols{"Feature2"} )
		and defined( $musitecols{"LinearDistance"} )
		and defined( $musitecols{"DistanceInfo"} ) ) {
		die "not a valid musite file\n";
	}
	my @wantmscols = (	$musitecols{"Gene1"} ,
						$musitecols{"Chromosome1"} ,
						$musitecols{"Start1"} ,
						$musitecols{"Stop1"} ,
						$musitecols{"Mutation1"} ,
						$musitecols{"Chain1"} ,
						$musitecols{"Position1"} ,
						$musitecols{"Gene2"} ,
						$musitecols{"Transcript2"} ,
						$musitecols{"TranscriptPosition2"} ,
						$musitecols{"Site2"} ,
						$musitecols{"Chain2"} ,
						$musitecols{"Position2"} ,
						$musitecols{"Feature2"} ,
						$musitecols{"LinearDistance"} ,
						$musitecols{"DistanceInfo"} );
	my $pdbCount = {};
	my $nMuSitePairs = 0;
	foreach my $line ( <$fh> ) { 
		chomp( $line );
		my @line = split /\t/ , $line;
		my ( $gene1 , $chromosome1 , $start1 , $stop1 , $aaChange1 , 
			 $chain1 , $position1 , 
			 $gene2 , $transcript2 , $transcriptPosition2 , $site2 , 
			 $chain2 , $position2 , $feature2 , 
			 $linearDistance , $infos ) = @line[@wantmscols];
		my $proteinMutation = TGI::ProteinVariant->new();
		my $mutation = TGI::Variant->new();
		$mutation->gene( $gene1 );
		$mutation->chromosome( $chromosome1 );
		$mutation->start( $start1 );
		$mutation->stop( $stop1 );
		$proteinMutation->aminoAcidChange( $aaChange1 );
		$mutation->addProteinVariant( $proteinMutation );
		my $site = TGI::Variant->new(); #for the purpose of making keys, will look like a mutation
		$site->gene( $gene2 );
		$site->chromosome( $transcript2 );
		$site->start( $transcriptPosition2 );
		$site->stop( $NULL );
		$site->reference( $NULL );
		$site->alternate( $PTM );
		$proteinMutation->transcript( $transcript2 );
		$proteinMutation->aminoAcidChange( $site2 );
		$site->addProteinVariant( $proteinMutation );
		$this->setDistance( $distance_matrix , $site , 
							$mutation , $chain1 , $chain2 , 
							$infos , $pdbCount );
		my $weight2 = 1;
		$this->setMutation( $mutations , $site , $feature2 , $weight2 );
		$nMuSitePairs++;
	}
	$fh->close();
	print STDOUT $nMuSitePairs." mutation-site pairs from .musites\n";
	return;
}

sub readSiteSitePairs {
	my ( $this , $mutations , $distance_matrix ) = @_;
    return if ( not defined $this->sitePairsFile() ); #if sitesite pairs included
	print STDOUT "HotSpot3D::Cluster::readSiteSitePairs\n";
	my $fh = new FileHandle;
	unless( $fh->open( $this->sitePairsFile() , "r" ) ) {
		die "Could not open sitesite pairs data file $! \n"
	};
	my $si = 0;
	my $headline = $fh->getline(); chomp( $headline );
	my %sitesitecols = map{ ( $_ , $si++ ) } split( /\t/ , $headline );
	my @required = ( "Gene1" , "Transcript1" , "TranscriptPosition1" , "Site1" ,
					 "Chain1" , "Position1" , "Feature1" , 
					 "Gene2" , "Transcript2" , "TranscriptPosition2" , "Site2" ,
					 "Chain2" , "Position2" , "Feature2" , 
					 "LinearDistance" , "DistanceInfo" );
	unless(	defined( $sitesitecols{"Gene1"} )
		and defined( $sitesitecols{"Transcript1"} )
		and defined( $sitesitecols{"TranscriptPosition1"} )
		and defined( $sitesitecols{"Site1"} )
		and defined( $sitesitecols{"Chain1"} )
		and defined( $sitesitecols{"Position1"} )
		and defined( $sitesitecols{"Feature1"} )
		and defined( $sitesitecols{"Gene2"} )
		and defined( $sitesitecols{"Transcript2"} )
		and defined( $sitesitecols{"TranscriptPosition2"} )
		and defined( $sitesitecols{"Site2"} )
		and defined( $sitesitecols{"Chain2"} )
		and defined( $sitesitecols{"Position2"} )
		and defined( $sitesitecols{"Feature2"} )
		and defined( $sitesitecols{"LinearDistance"} )
		and defined( $sitesitecols{"DistanceInfo"} ) ) {
		die "not a valid sites file\n";
	}
	my @wantsitecols = (	$sitesitecols{"Gene1"} ,
						$sitesitecols{"Transcript1"} ,
						$sitesitecols{"TranscriptPosition1"} ,
						$sitesitecols{"Site1"} ,
						$sitesitecols{"Chain1"} ,
						$sitesitecols{"Position1"} ,
						$sitesitecols{"Feature1"} ,
						$sitesitecols{"Gene2"} ,
						$sitesitecols{"Transcript2"} ,
						$sitesitecols{"TranscriptPosition2"} ,
						$sitesitecols{"Site2"} ,
						$sitesitecols{"Chain2"} ,
						$sitesitecols{"Position2"} ,
						$sitesitecols{"Feature2"} ,
						$sitesitecols{"LinearDistance"} ,
						$sitesitecols{"DistanceInfo"} );
	my $pdbCount = {};
	my $nSitePairs = 0;
	foreach my $line ( <$fh> ) { 
		chomp( $line );
		my @line = split /\t/ , $line;
		map{ $line =~ s/"//g } @line;
		my ( $gene1 , $transcript1 , $transcriptPosition1 , $site1 , 
			 $chain1 , $position1 , $feature1 , 
			 $gene2 , $transcript2 , $transcriptPosition2 , $site2 , 
			 $chain2 , $position2 , $feature2 , 
			 $linearDistance , $infos ) = @line[@wantsitecols];
		my $proteinSite1 = TGI::ProteinVariant->new();
		my $tsite1 = TGI::Variant->new();
		$tsite1->gene( $gene1);
		$tsite1->chromosome( $transcript1 );
		$tsite1->start( $transcriptPosition1 );
		$tsite1->stop( $NULL );
		$tsite1->reference( $NULL );
		$tsite1->alternate( $PTM );
		$proteinSite1->transcript( $transcript1 );
		$proteinSite1->aminoAcidChange( $site1 );
		$tsite1->addProteinVariant( $proteinSite1 );
		my $proteinSite2 = TGI::ProteinVariant->new();
		my $tsite2 = TGI::Variant->new(); #for the purpose of making keys, will look like a site21
		$tsite2->gene( $gene2 );
		$tsite2->chromosome( $transcript2 );
		$tsite2->start( $transcriptPosition2 );
		$tsite2->stop( $NULL );
		$tsite2->reference( $NULL );
		$tsite2->alternate( $PTM );
		$proteinSite2->aminoAcidChange( $site2 );
		$proteinSite2->transcript( $transcript2 );
		$tsite2->addProteinVariant( $proteinSite2 );
		$this->setDistance( $distance_matrix , $tsite1 , 
							$tsite2 , $chain1 , $chain2 , 
							$infos , $pdbCount );
		my $weight1 = 1;
		$this->setMutation( $mutations , $tsite1 , $feature1 , $weight1 );
		my $weight2 = 1;
		$this->setMutation( $mutations , $tsite2 , $feature2 , $weight2 );
		$nSitePairs++;
	}
	$fh->close();
	print STDOUT $nSitePairs." site-site pairs from .sites\n";
	#print $this->printMutations( $mutations , "sitesite" );
	return;
}

sub getMutationSitePairs {
	#$this->getMutationSitePairs( $distance_matrix );
	my ( $this , $mutations , $distance_matrix ) = @_;
	print STDOUT "HotSpot3D::Cluster::getMutationSitePairs\n";
	$this->readMutationSitePairs( $mutations , $distance_matrix );
	return;
}

sub getSiteSitePairs {
	#$this->getSiteSitePairs( $distance_matrix );
	my ( $this , $mutations , $distance_matrix ) = @_;
	print STDOUT "HotSpot3D::Cluster::getSiteSitePairs\n";
	$this->readSiteSitePairs( $mutations , $distance_matrix );
	return;
}

sub getDrugMutationPairs {
	#$this->getDrugMutationPairs( $distance_matrix );
	my ( $this , $mutations , $distance_matrix ) = @_;
	print STDOUT "HotSpot3D::Cluster::getDrugMutationPairs\n";
	$this->readDrugClean( $mutations , $distance_matrix );
	return;
}

sub readDrugClean {
	#$this->readDrugClean( $distance_matrix );
	my ( $this , $mutations , $distance_matrix ) = @_;
	print STDOUT "HotSpot3D::Cluster::readDrugClean\n";
    return if ( not defined $this->drugsCleanFile() ); #if drug pairs included
	my $fh = new FileHandle;
	unless( $fh->open( $this->drugsCleanFile() , "r" ) ) {
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
	foreach my $line ( <$fh> ) { 
			chomp( $line );
			my @line = split /\t/ , $line;
			map{ $line =~ s/"//g } @line;
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
	}
	$fh->close();
	return;
}

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
				unless ( $refAlt1 =~ m/$PTM/ ) {
					print "...next, no match aaChange1\n";
					next;
				}
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
						unless ( $refAlt2 =~ m/$PTM/ ) {
							print "...next, no match aaChange2\n";
							next;
						}
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

#sub checkProteinPosition {
#	my ( $this , $proteinKey1 , $proteinKey2 ) = @_;
#	#print join( "\t" , ( $aaPosition1 , $aaReference1 , $aaPosition2 , $aaReference2 ) )."\t";
#	return 0;
#}

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
	#print "HotSpot3D::Cluster::setShortestDistance\n";
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
--structure-dependent        Clusters for each structure or across all structures (dependent or independent), default: independent 

--help                       this message

HELP
}

sub help_text {
    my $this = shift;
#--site-file                  The site list file used in proximity search step
        return <<HELP

Usage: hotspot3d cluster [options]

                             REQUIRE AT LEAST ONE
--sites-file                 A .sites file with site-site pairs
--pairwise-file              A .pairwise file with mutation-mutation pairs (provide maf-file)
--drug-clean-file            A .drugs.*target.clean file with mutation-drug pairs (provide maf-file)
--musites-file               A .musites file with mutation-site pairs (provide maf-file and site-file)

                             CONDITIONAL REQUIREMENT
--maf-file                   .maf file used in proximity search step
                                 necessary for pairwise, drug-clean, or musites pair data

                             OPTIONAL
--output-prefix              Output prefix, default: 3D_Proximity
--p-value-cutoff             P_value cutoff (<), default: 0.05 (if 3d-distance-cutoff also not set)
--3d-distance-cutoff         3D distance cutoff (<), default: 100 (if p-value-cutoff also not set)
--linear-cutoff              Linear distance cutoff (> peptides), default: 0
--max-radius                 Maximum cluster radius (max network geodesic from centroid, <= Angstroms), default: 10
--clustering                 Cluster using network or density-based methods (network or density), default: network
--vertex-type                Graph vertex type for network-based clustering (recurrence, unique, site or weight), default: site
                                 recurrence vertices are the genomic mutations for each sample from the given .maf
                                 unique vertices are the specific genomic changes
                                 site vertices are the affected protein positions
                                 weight vertices are the genomic mutations with a numerical weighting
--distance-measure           Pair distance to use (shortest or average), default: average
--structure-dependent        Clusters for each structure or across all structures, default (no flag): independent
--subunit-dependent          Clusters for each subunit or across all subunits, default (no flag): independent
--meric-type                 Clusters for each intra-molecular (both monomers and homomers), monomer, homomer, 
                                 inter-molecular (heteromers), heteromer, multimer (simultaneously homomer & heteromer), or any *mer 
                                 (intra, monomer, homomer, inter, heteromer, multimer, or unspecified), default: unspecified
--transcript-id-header       .maf file column header for transcript id's, default: transcript_name
--amino-acid-header          .maf file column header for amino acid changes, default: amino_acid_change 
--weight-header              .maf file column header for mutation weight, default: weight (used if vertex-type = weight)
--parallel                   Parallelization for structure and subunit dependent runs (none or local), default: none
--max-processes              Set if using parallel type local (CAUTION: make sure you know your max CPU processes)
--gene-list-file             Choose mutations from the genes given in this list
--structure-list-file        Choose mutations from the structures given in this list
--use-JSON                   Use pre-encoded mutations and distance-matrix hashes in json format, default (no flag): do not use json
--mutations-hash-json-file   JSON encoded mutations hash file produced by a previous cluster run
--distance-matrix-json-file  JSON encoded distance-matrix hash file produced by a previous cluster run



--help                       this message

HELP

}

1;
