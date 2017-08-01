package TGI::Mutpro::Main::Visual;
#
#----------------------------------
# $Authors: Adam Scott & Beifang Niu 
# $Date: 2014-01-14 14:34:50 -0500 (Tue Jan 14 14:34:50 CST 2014) $
# $Revision: 0.3 $
# $URL: $
# $Doc: $ clusters visualization
#----------------------------------
#
use strict;
use warnings;

use Carp;
use Cwd;
use Getopt::Long;
use LWP::Simple;
use IO::File;
use FileHandle;
use List::MoreUtils qw( uniq );

my $NETWORK = "network";
my $DENSITY = "density";
my $PTM = "ptm";

sub new {
    my $class = shift;
    my $this = {};
	$this->{_CLUSTERS_FILE} = undef;
    $this->{_PAIRWISE_FILE} = undef;
    $this->{_DRUG_PAIRS_FILE} = undef;
    $this->{_OUTPUT_FILE} = undef;
    $this->{_OUTPUT_DIR} = getcwd;
    $this->{_PDB_DIR} = getcwd;
    $this->{_PYMOL} = '/usr/bin/pymol';
    $this->{_SCRIPT_ONLY} = undef;
    $this->{_MUT_COLOR} = '1,0,0';
    $this->{_SITE_COLOR} = '0,0,1';
    $this->{_COMPOUND_COLOR} = 'util.cbag';
    $this->{_MUT_STYLE} = 'spheres';
    $this->{_SITE_STYLE} = 'spheres';
    $this->{_COMPOUND_STYLE} = 'sticks';
    $this->{_BG_COLOR} = 'white';
    $this->{_NO_LABEL} = undef;
    $this->{_STAT} = undef;
	$this->{_PDB} = undef;
	$this->{mutations} = undef;
	$this->{clusters} = undef;
	$this->{chains} = undef;
	$this->{compounds} = undef;
	$this->{locations} = undef;
	$this->{centroids} = undef;
	$this->{clusters_file_type} = undef;
    bless $this, $class;
    $this->process();
    return $this;
}

sub process {
    my $this = shift;
    my ( $help, $options );
    unless( @ARGV ) { die $this->help_text(); }
    $options = GetOptions (
		'clusters-file=s' => \$this->{_CLUSTERS_FILE},
        'pairwise-file=s' => \$this->{_PAIRWISE_FILE},
        'musites-file=s' => \$this->{_MUSITES_FILE},
        'sites-file=s' => \$this->{_SITES_FILE},
        'drug-pairs-file=s' => \$this->{_DRUG_PAIRS_FILE},
        'output-dir=s' => \$this->{_OUTPUT_DIR},
        'output-file=s' => \$this->{_OUTPUT_FILE},
        'pdb-dir=s' => \$this->{_PDB_DIR},
		'pdb=s' => \$this->{_PDB},
        'pymol=s' => \$this->{_PYMOL},
		'script-only' => \$this->{_SCRIPT_ONLY},
        'mut-color=s' => \$this->{_MUT_COLOR},
        'site-color=s' => \$this->{_SITE_COLOR},
        'compound-color=s' => \$this->{_COMPOUND_COLOR},
        'mut-style=s' => \$this->{_MUT_STYLE},
        'site-style=s' => \$this->{_SITE_STYLE},
        'compound-style=s' => \$this->{_COMPOUND_STYLE},
        'no-label' => \$this->{_NO_LABEL},
        'bg-color=s' => \$this->{_BG_COLOR},
        'clusters-file-type=s' => \$this->{clusters_file_type},
        'help' => \$help
    );
    if ( $help ) { die $this->help_text(); }
    unless( $options ) { die $this->help_text(); }
    if (    ( not $this->{_PAIRWISE_FILE} ) 
		and ( not $this->{_MUSITES_FILE} ) 
		and ( not $this->{_SITES_FILE} ) 
		and ( not $this->{_DRUG_PAIRS_FILE} ) ) { #no input
		warn 'You must provide a pairwise, musites, sites, or drug pairs file! ', "\n";
		die $this->help_text();
	}
    unless( $this->{_CLUSTERS_FILE} ) { warn 'You must provide a cluster file ! ', "\n"; die $this->help_text(); }
    unless( -e $this->{_CLUSTERS_FILE} ) { warn ' cluster file is not exist  ! ', "\n"; die $this->help_text(); }
    unless( $this->{_OUTPUT_FILE} || $this->{_OUTPUT_DIR} ) { warn 'You must provide an output file or directory! ', "\n"; die $this->help_text(); }
	if ( $this->{_PAIRWISE_FILE} ) {
		if ( not -e $this->{_PAIRWISE_FILE} ) {
			warn ' pairwise file does not exist  ! ', "\n";
			die $this->help_text();
		}
	} 
	if ( $this->{_MUSITES_FILE} ) {
		if ( not -e $this->{_MUSITES_FILE} ) {
			warn ' musites file does not exist  ! ', "\n";
			die $this->help_text();
		}
	}
	if ( $this->{_SITES_FILE} ) {
		if ( not -e $this->{_SITES_FILE} ) {
			warn ' sites file does not exist  ! ', "\n";
			die $this->help_text();
		}
	}
	if ( $this->{_DRUG_PAIRS_FILE} ) {
		if ( not -e $this->{_DRUG_PAIRS_FILE} ) {
			warn ' drug pairs file does not exist  ! ', "\n";
			die $this->help_text();
		}
	}
    map{ $this->{_STAT}{$_} = 0; } qw( num_muts pdb pairs );

    if ( not defined $this->{clusters_file_type} ) {
		$this->{clusters_file_type} = $NETWORK;
		warn "HotSpot3D::Visual::setOptions warning: no clusters-file-type option given, setting to default network\n";
	}

	if ( $this->{clusters_file_type} eq $NETWORK ) {
		$this->getMappingLocations(  );

		#$this->checkPDB(  );

		unless( $this->{mutations} ) { warn "No proximal mutation pairs were found on ".$this->{_PDB}."\n"; }
		unless( $this->{sites} ) { warn "No proximal site pairs were found on ".$this->{_PDB}."\n"; }
		$this->getClusters(  );
		unless( $this->{clusters} ) { warn "No cluster to display on ".$this->{_PDB}."\n"; }
		$this->makePyMOLScript(  );

		$this->execute(  );

	    return 1;
	}
	elsif ( $this->{clusters_file_type} eq $DENSITY ) { 
		$this->densityVisual();
	}
	else {
		die "Error: clusters-file-type should be one of the following: network or density (default:network)\n";
	}
}

# check PDBs  
sub addPDB {
    my ( $this ) = @_;
	my $pdbfile = $this->{_PDB_DIR}."/".$this->{_PDB}.".pdb";
	my $download;
	unless( -e $pdbfile ) {
		print STDOUT "Downloading ".$this->{_PDB}.".pdb\n";
        my $pdburl = "http://www.rcsb.org/pdb/files/$this->{_PDB}.pdb";
        $download = get( $pdburl );
		my $fh = new FileHandle;
		$fh->open( $pdbfile , "w" );
		$fh->print( $download );
		$fh->close();
	}
	unless( -e $pdbfile ) { warn "Could not download ".$pdbfile."\n";
	} else {
		#print STDOUT "Compressing ".$pdbfile."\n";
		#system( "gzip ".$pdbfile );
	}
	#TO DO: use residue name to check all possible positions of drugs
    return 1;
}

sub getMappingLocations {
	my ( $this ) = @_;
	$this->readPairwise();
	$this->readMuSites();
	$this->readSites();
	$this->readDrugPairs();
	return 1;
}

sub readPairwise {
	my ( $this ) = @_;
	if ( $this->{_PAIRWISE_FILE} ) {
		my $fh = new FileHandle;
		print STDOUT "Getting mutation-mutation pairs\n";
		unless( $fh->open( $this->{_PAIRWISE_FILE} ) ) { die "Could not open pairwise file\n" };
		my @ccols;
		while ( my $line = <$fh> ) {
			chomp( $line );
			if ( $line =~ /Gene/ ) {
				my $i = 0;
				my %ccols = map{ ( $_ , $i++ ) } split( /\t/ , $line );
				unless(	defined( $ccols{"Gene1"} )
					and defined( $ccols{"Mutation1"} )
					and defined( $ccols{"Chain1"} )
					and defined( $ccols{"Position1"} )
					and defined( $ccols{"Gene2"} )
					and defined( $ccols{"Mutation2"} )
					and defined( $ccols{"Chain2"} )
					and defined( $ccols{"Position2"} ) ) {
					die "Not a valid pairwise file!\n";
				}
				@ccols = (	$ccols{"Gene1"} ,
							$ccols{"Mutation1"} , 
							$ccols{"Chain1"} , 
							$ccols{"Position1"} , 
							$ccols{"Gene2"} ,
							$ccols{"Mutation2"} , 
							$ccols{"Chain2"} , 
							$ccols{"Position2"} );
			} else {
				if ( $line =~ /$this->{_PDB}/ ) {
					my ( $gene1 , $mu1 , $chain1 , $res1 , $gene2 , $mu2 , $chain2 , $res2 ) = ( split /\t/ , $line )[@ccols];
					$chain1 =~ s/\[(\w)\]/$1/g;
					$chain2 =~ s/\[(\w)\]/$1/g;
					$this->{mutations}->{$gene1.":".$mu1}->{$chain1.":".$res1} = 0;
					$this->{mutations}->{$gene2.":".$mu2}->{$chain2.":".$res2} = 0;
					$this->{chains}->{$chain1}->{$gene1} = 0;
					$this->{chains}->{$chain2}->{$gene2} = 0;
					push @{$this->{locations}->{$gene1}->{$chain1}->{$res1}} , $mu1;
					push @{$this->{locations}->{$gene2}->{$chain2}->{$res2}} , $mu2;
				}
			}
		}
		$fh->close();
	}
	return;
}

sub readMuSites {
	my ( $this ) = @_;
	if ( $this->{_MUSITES_FILE} ) {
		my $fh = new FileHandle;
		print STDOUT "Getting mutation-site pairs\n";
		unless( $fh->open( $this->{_MUSITES_FILE} ) ) { die "Could not open musites file\n" };
		my @ccols;
		while ( my $line = <$fh> ) {
			chomp( $line );
			if ( $line =~ /Gene/ ) {
				my $i = 0;
				my %ccols = map{ ( $_ , $i++ ) } split( /\t/ , $line );
				unless(	defined( $ccols{"Gene1"} )
					and defined( $ccols{"Mutation1"} )
					and defined( $ccols{"Chain1"} )
					and defined( $ccols{"Position1"} )
					and defined( $ccols{"Gene2"} )
					and defined( $ccols{"Site2"} )
					and defined( $ccols{"Chain2"} )
					and defined( $ccols{"Position2"} ) ) {
					die "Not a valid musites file!\n";
				}
				@ccols = (	$ccols{"Gene1"} ,
							$ccols{"Mutation1"} , 
							$ccols{"Chain1"} , 
							$ccols{"Position1"} , 
							$ccols{"Gene2"} ,
							$ccols{"Site2"} , 
							$ccols{"Chain2"} , 
							$ccols{"Position2"} );
			} else {
				if ( $line =~ /$this->{_PDB}/ ) {
					my ( $gene1 , $mu1 , $chain1 , $res1 , $gene2 , $site2 , $chain2 , $res2 ) = ( split /\t/ , $line )[@ccols];
					$chain1 =~ s/\[(\w)\]/$1/g;
					$chain2 =~ s/\[(\w)\]/$1/g;
					$this->{mutations}->{$gene1.":".$mu1}->{$chain1.":".$res1} = 0;
					$this->{sites}->{$gene2.":".$site2}->{$chain2.":".$res2} = 0;
					$this->{chains}->{$chain1}->{$gene1} = 0;
					$this->{chains}->{$chain2}->{$gene2} = 0;
					push @{$this->{locations}->{$gene1}->{$chain1}->{$res1}} , $mu1;
					push @{$this->{locations}->{$gene2}->{$chain2}->{$res2}} , $site2;
				}
			}
		}
		$fh->close();
	}
	return;
}

sub readSites {
	my ( $this ) = @_;
	if ( $this->{_SITES_FILE} ) {
		my $fh = new FileHandle;
		print STDOUT "Getting site-site pairs\n";
		unless( $fh->open( $this->{_SITES_FILE} ) ) { die "Could not open sites file\n" };
		my @ccols;
		while ( my $line = <$fh> ) {
			chomp( $line );
			if ( $line =~ /Gene/ ) {
				my $i = 0;
				my %ccols = map{ ( $_ , $i++ ) } split( /\t/ , $line );
				unless(	defined( $ccols{"Gene1"} )
					and defined( $ccols{"Site1"} )
					and defined( $ccols{"Chain1"} )
					and defined( $ccols{"Position1"} )
					and defined( $ccols{"Gene2"} )
					and defined( $ccols{"Site2"} )
					and defined( $ccols{"Chain2"} )
					and defined( $ccols{"Position2"} ) ) {
					die "Not a valid musites file!\n";
				}
				@ccols = (	$ccols{"Gene1"} ,
							$ccols{"Site1"} , 
							$ccols{"Chain1"} , 
							$ccols{"Position1"} , 
							$ccols{"Gene2"} ,
							$ccols{"Site2"} , 
							$ccols{"Chain2"} , 
							$ccols{"Position2"} );
			} else {
				if ( $line =~ /$this->{_PDB}/ ) {
					my ( $gene1 , $site1 , $chain1 , $res1 , $gene2 , $site2 , $chain2 , $res2 ) = ( split /\t/ , $line )[@ccols];
					$chain1 =~ s/\[(\w)\]/$1/g;
					$chain2 =~ s/\[(\w)\]/$1/g;
					$this->{sites}->{$gene1.":".$site1}->{$chain1.":".$res1} = 0;
					$this->{sites}->{$gene2.":".$site2}->{$chain2.":".$res2} = 0;
					$this->{chains}->{$chain1}->{$gene1} = 0;
					$this->{chains}->{$chain2}->{$gene2} = 0;
					push @{$this->{locations}->{$gene1}->{$chain1}->{$res1}} , $site1;
					push @{$this->{locations}->{$gene2}->{$chain2}->{$res2}} , $site2;
				}
			}
		}
		$fh->close();
	}
	return;
}

sub readDrugPairs {
	my ( $this ) = @_;
	if ( $this->{_DRUG_PAIRS_FILE} ) {
		my $dfh = new FileHandle;
		print STDOUT "Getting drug-mutation pairs\n";
		unless( $dfh->open( $this->{_DRUG_PAIRS_FILE} ) ) { die "Could not open drug pairs file\n" };
		my @ccols;
		while ( my $line = <$dfh> ) {
			chomp( $line );
			if ( $line =~ /Gene/ ) {
				my $i = 0;
				my %ccols = map{ ( $_ , $i++ ) } split( /\t/ , $line );
				unless(	defined( $ccols{"Drug"} )
					and defined( $ccols{"Chain1"} )
					and defined( $ccols{"Position1"} )
					and defined( $ccols{"Gene2"} )
					and defined( $ccols{"Mutation2"} )
					and defined( $ccols{"Chain2"} )
					and defined( $ccols{"Position2"} ) ) {
					die "Not a valid drug-clean file!\n";
				}
				@ccols = (	$ccols{"Gene1"} ,
							$ccols{"Site1"} , 
							$ccols{"Chain1"} , 
							$ccols{"Position1"} , 
							$ccols{"Gene2"} ,
							$ccols{"Site2"} , 
							$ccols{"Chain2"} , 
							$ccols{"Position2"} );
			} else {
				if ( $line =~ /$this->{_PDB}/ ) {
					my ( $drug , $chain1 , $res1 , $gene2 , $mu2 , $chain2 , $res2 ) = ( split /\t/ , $line )[@ccols];
					$chain1 =~ s/\[(\w)\]/$1/g;
					$chain2 =~ s/\[(\w)\]/$1/g;
					#$mutations->{$drug.":".$gene2}->{$chain1.":".$res1} = 0;
					$this->{mutations}->{$gene2.":".$mu2}->{$chain2.":".$res2} = 0;
					$this->{compounds}->{$drug}->{$chain1} = $res1;
					$this->{chains}->{$chain1}->{$drug} = 0;
					$this->{chains}->{$chain2}->{$gene2} = 0;
					push @{$this->{locations}->{$drug}->{$chain1}->{$res1}} , $drug;
					push @{$this->{locations}->{$gene2}->{$chain2}->{$res2}} , $mu2;
				}
			}
		}
		$dfh->close();
	}
	return;
}

sub getClusters {
	my ( $this ) = @_;
	my $fh = new FileHandle;
	unless( $fh->open( $this->{_CLUSTERS_FILE} ) ) { die "Could not open clusters file\n" };
	my @ccols;
	print STDOUT "Getting clusters\n";
	while ( my $line = <$fh> ) {
		chomp( $line );
		if ( $line =~ /Cluster/ ) {
			my $i = 0;
			my %ccols = map{ ( $_ , $i++ ) } split( /\t/ , $line );
			unless(	defined( $ccols{"Cluster"} )
				and defined( $ccols{"Gene/Drug"} )
				and defined( $ccols{"Mutation/Gene"} )
				and defined( $ccols{"Closeness_Centrality"} )
				and defined( $ccols{"Geodesic_From_Centroid"} )
				and defined( $ccols{"Alternate"} ) ) {
				die "Not a valid clusters file!\n";
			}
			@ccols = (	$ccols{"Cluster"} ,
						$ccols{"Gene/Drug"} , 
						$ccols{"Mutation/Gene"} , 
						$ccols{"Closeness_Centrality"} , 
						$ccols{"Geodesic_From_Centroid"} ,
						$ccols{"Alternate"} );
		} else {
			my ( $id , $gd , $mp , $cc , $gfc , $alt ) = ( split /\t/ , $line )[@ccols];
			if ( $id =~ m/\d+\.\d+\.\w/ ) {
				next if ( $id !~ m/$this->{_PDB}/ );
			}
			my $variant = $gd.":".$mp;
			if ( $alt eq $PTM ) {
				if ( exists $this->{sites}->{$variant} ) {
					my @positions = keys %{$this->{sites}->{$variant}};
					if ( $gfc == 0 ) {
						push @{$this->{centroids}->{$id}} , $mp;
						foreach my $pos ( @positions ) { $this->{sites}->{$variant}->{$pos} = 1; }
					}
					$this->{clusters}->{$id}->{$gd}->{$mp} = $cc;
				}
			} else {
				if ( exists $this->{mutations}->{$variant} ) {
					my @positions = keys %{$this->{mutations}->{$variant}};
					if ( $gfc == 0 ) {
						push @{$this->{centroids}->{$id}} , $mp;
						foreach my $pos ( @positions ) { $this->{mutations}->{$variant}->{$pos} = 1; }
					}
					$this->{clusters}->{$id}->{$gd}->{$mp} = $cc;
				}
			}
		}
	}
	$fh->close();
	return 1;
}

sub makePyMOLScript {
	my ( $this ) = @_;
	print STDOUT "Generating PyMOL script (.pml)\n";
	my $residue_line = {};
	my $chain_colors = {};
	my $old_gene = "";
	my $centroid_color = "grey";
	my $centroid_site_color = "black";
	my $site_color = "white";
	my $CHAIN_COLORS = ['pink' , 'palegreen' , 'lightblue' , 'lightmagenta' , 'lightorange' , 'lightpink' , 'paleyellow' , 'lightteal' , 'aquamarine' , 'palecyan'];
	my $CLUSTER_COLORS = ['red' , 'blue' , 'green' , 'magenta' , 'forest' , 'brightorange' , 'chocolate' , 'deepblue' , 'deepolive' , 'deeppurple' , 'bluewhite' , 'lime' , 'purple' , 'dash' , 'orange' , 'brown' , 'salmon' , 'oxygen' , 'wheat' , 'violetpurple' , 'limon' , 'sand' , 'raspberry' , 'slate'];
	foreach my $id ( sort keys %{$this->{clusters}} ) {
		print STDOUT "\tWorking on cluster ".$id."\n";
		my $cluster_color = shift @{$CLUSTER_COLORS};
		my $drawing_centroid = 0;
		foreach my $genedrug ( keys %{$this->{clusters}->{$id}} ) {
			#print "\t".$id." has gene ".$genedrug."\n";
			my @cluster_AA_changes = keys %{$this->{clusters}->{$id}->{$genedrug}};
			my $clusters_AA_changes = join( "," , @cluster_AA_changes );
			foreach my $chain ( sort keys %{$this->{chains}} ) {
				#print "\t\t".$id.": ".$genedrug." is on chain ".$chain."\n";
				if ( not exists $chain_colors->{$chain} ) {
					$chain_colors->{$chain} = shift @{$CHAIN_COLORS};
				}
				foreach my $location ( sort {$a <=> $b} keys %{$this->{locations}->{$genedrug}->{$chain}} ) {
					#print "\t\t\t".$id.": ".$genedrug." has a mutation at chain ".$chain.", residue ".$location."\n";
					my @location_AA_changes = uniq( @{$this->{locations}->{$genedrug}->{$chain}->{$location}} );
					my $comment = $genedrug.": ";
					$comment .= join( ", " , @location_AA_changes );
					my $AA_change = shift @location_AA_changes;
					my $residue = $AA_change;
					$residue =~ s/^p\.(.+)\D*$/$1/g;
					my $isSite = 0;
					if ( exists $this->{sites}->{$genedrug.":".$AA_change} ) {
						$isSite = 1;
					}
					if ( ( grep{ $_ eq $residue } @cluster_AA_changes ) || ( grep{ $_ eq $AA_change } @cluster_AA_changes ) ) {
						if ( not exists $residue_line->{$chain}->{$location} ) {
							#print "\t\t\t\tMaking line ".$id.": ".$genedrug.":".join( "," , @location_AA_changes )." - ".$chain.", ".$location."\n";
							my $label = "";
							if ( $isSite ) {
								$label = $cluster_color."_".$genedrug."_ptm".$residue."_".$chain;
							} else {
								$label = $cluster_color."_".$genedrug."_".$residue."_".$chain;
							}
							$residue_line->{$chain}->{$location} = "sele ".$label.", (resi ".$location." and chain ".$chain."); ";
							if ( grep{ $_ eq $AA_change } @{$this->{centroids}->{$id}} ) { #is centroid
								if ( $isSite ) {
									$residue_line->{$chain}->{$location} .= "color ".$centroid_site_color.", ".$label."; ";
								} else {
									$residue_line->{$chain}->{$location} .= "color ".$centroid_color.", ".$label."; ";
								}
							} else { #is not centroid
								if ( $isSite ) {
									$residue_line->{$chain}->{$location} .= "color ".$site_color.", ".$label."; ";
								} else {
									$residue_line->{$chain}->{$location} .= "color ".$cluster_color.", ".$label."; ";
								}
							}
							if ( $isSite ) {
								$residue_line->{$chain}->{$location} .= "show ".$this->{_SITE_STYLE}.", ".$label."; ";
							} else {
								$residue_line->{$chain}->{$location} .= "show ".$this->{_MUT_STYLE}.", ".$label."; ";
							}
							#$residue_line->{$chain}->{$location} .= "set label_font_id, 10; ";
							#$residue_line->{$chain}->{$location} .= "set label_position, (3,2,1); ";
							#$residue_line->{$chain}->{$location} .= "set label_size, 20; ";
							#$residue_line->{$chain}->{$location} .= "label ".$label." and name CA, \"".$residue."\";\n";
							$residue_line->{$chain}->{$location} .= "#".$comment."\n";
						}
					} #if
				} #location
			} #chain
		} #gene/drug
	} #cluster id
	if ( $this->{compounds} ) {
		foreach my $compound ( keys %{$this->{compounds}} ) {
			foreach my $chain ( keys %{$this->{compounds}->{$compound}} ) {
				my $location = $this->{compounds}->{$compound}->{$chain};
				$residue_line->{$chain}->{$location} = "sele ".$compound.", (resi ".$location." and chain ".$chain."); ";
				if ( $this->{_COMPOUND_COLOR} =~ /util/ ) {
					$residue_line->{$chain}->{$location} .= "util.cbag ".$compound.";\n"
				} else {
					$residue_line->{$chain}->{$location} .= "color ".$this->{_COMPOUND_COLOR}.", ".$compound."; ";
					$residue_line->{$chain}->{$location} .= "show ".$this->{_COMPOUND_STYLE}.", ".$compound.";\n";
				}
			}
		}
	}
	my $fh = new FileHandle;
	my $outFilename; 
	if ( defined $this->{_OUTPUT_FILE} ) { 
		$outFilename = $this->{_OUTPUT_FILE};
	} else {
		$this->{_OUTPUT_FILE} = $this->{_CLUSTERS_FILE}.".".$this->{_PDB}.".pml";
		$outFilename = $this->{_OUTPUT_FILE};
	}
	unless( $fh->open( $outFilename , "w" ) ) { die "Could not open output file, $outFilename\n"; }

	#viewer
	$fh->print( "reinitialize everything;\n" );
	#$this->addPDB(  );
	if ( -e $this->{_PDB_DIR}."/".$this->{_PDB} ) {
		$fh->print( "load ".$this->{_PDB_DIR}."/".$this->{_PDB}.".pdb;\n" );
	} else {
		$fh->print( "load http://www.rcsb.org/pdb/files/".$this->{_PDB}.".pdb;\n" );
	}
	$fh->print( "viewport 480,480;\n" );
	$fh->print( "preset.publication(\"".$this->{_PDB}."\");\n" );
	$fh->print( "#show mesh;\n" );
	$fh->print( "\n" );
	$fh->print( "bg_color ".$this->{_BG_COLOR}.";\n" );
	$fh->print( "set ray shadows, off;\n" );
	$fh->print( "set ray_opaque_background, off;\n" );
	$fh->print( "\n" );

	#chains
	foreach my $chain ( sort keys %{$chain_colors} ) {
		my $comment = join( "\," , sort keys %{$this->{chains}->{$chain}} );
		$fh->print( "color ".$chain_colors->{$chain}.", chain ".$chain."; #".$comment."\n" );
	}
	$fh->print( "\n" );

	#residues
	foreach my $chain ( sort keys %{$residue_line} ) {
		$fh->print( "#Residues from chain ".$chain."\n" );
		foreach my $location ( sort {$a <=> $b} keys %{$residue_line->{$chain}} ) {
			$fh->print( $residue_line->{$chain}->{$location} );
		}
		$fh->print( "\n" );
	}
	$fh->print( "\n" );

	#movie
	$fh->print( "\n\n#mset 1\n#mdo 1: turn x,0.2; turn y,0.2; turn z,0.2;\n#mplay" );

	$fh->close();

	return 1;
}

sub execute {
	my ( $this ) = @_;
	if ( not $this->{_SCRIPT_ONLY} ) {
		my $pymol = $this->{_PYMOL};
		my $pml = $this->{_OUTPUT_FILE};
		print STDOUT "HotSpot3D ... Visualize!\n\n";
		system( "$pymol $pml" );
		#pymol <output-file>
		#-x disables external GUI module
	}
}

############################################################################
########################### Density Visual #################################
############################################################################

sub densityVisual {
	my $this = shift;

	my $ClustersFile = $this->{_CLUSTERS_FILE};
	my $PairwiseFile = $this->{_PAIRWISE_FILE};
	my $PDB = $this->{_PDB};

	$this->GetFileName(); # obtain the clusters-file-name only; for file-naming purposes

	########################  Reading from *.clusters  #############################

	my @InputClusters;

	open(IN, "<$ClustersFile") || die "Can't open $ClustersFile: $!";
	while (my $line = <IN>) {
		if ( not $line =~ /Cluster/ ) {
			chomp $line;
			my @tabs3 = split(/\t/,$line);
			push @InputClusters, [$tabs3[0],$tabs3[1],$tabs3[2],$tabs3[7],$tabs3[8],$tabs3[9]]; # ClusterID, Gene, Mutation, Epsilon_prime, Avg_density, Covering_clusters, last_entry_to_mark_processed
		}
	}
	close(IN);

	###############################################################################

	#print Dumper \@InputClusters;

	# my $this = {};
	#my @spectrum = ("blue_yellow","rainbow","green_red","red_yellow","blue_green","cyan_red");
	my $ChainColors = ['pink' , 'palegreen' , 'lightblue' , 'lightmagenta' , 'lightorange' , 'lightpink' , 'paleyellow' , 'lightteal' , 'aquamarine' , 'palecyan'];
	my $SurfaceColors = ['violetpurple','violet','deeppurple', 'purple','lightmagenta', 'blue', 'lightblue', 'bluewhite', 'oxygen' , 'green', 'deepolive' , 'wheat' , 'lime' , 'brown' , 'orange' , 'yellow', 'salmon' , 'red',  'violetpurple' , 'limon' , 'sand' , 'raspberry' , 'slate', 'red','blue', 'green', 'yellow', 'deepolive' , 'deeppurple' , 'bluewhite' , 'lime' , 'purple' , 'dash' , 'orange' , 'brown' , 'salmon' , 'oxygen' , 'wheat' , 'violetpurple' , 'limon' , 'sand' , 'raspberry' , 'slate','red','blue', 'green', 'yellow', 'deepolive' , 'deeppurple' , 'bluewhite' , 'lime' , 'purple' , 'dash' , 'orange' , 'brown' , 'salmon' , 'oxygen' , 'wheat' , 'violetpurple' , 'limon' , 'sand' , 'raspberry' , 'slate'];
	#my $spectrumPalette = shift @spectrum;
	#my $ClusterSurfaceColor = shift $SurfaceColors;
	my $currentColor = shift @{$SurfaceColors};


	for (my $i = 0; $i < scalar @InputClusters; $i++) {
		my @IDs = split(/\./,$InputClusters[$i][0]);
		$this->{ClusterStructure}->{$IDs[0]}->{$IDs[1]} = $IDs[2];
	}
	#print Dumper \%ClusterStructure;

	$this->getMappingLocations( );
	print "Done getting mapping locations.\n";

	#################  Clusters present in the given PDB structure  #####################

	my $OutFile1 = "DensityVisual.$this->{clusters_file_name_only}.$this->{_PDB}.presentClusters";
	open (OUT, ">$OutFile1");

	for (my $i = 0; $i < scalar @InputClusters; $i++) { 
		if ( exists $this->{mutations}->{$InputClusters[$i][1].":".$InputClusters[$i][2]}) {
			print OUT "$InputClusters[$i][0]\t$InputClusters[$i][1]\t$InputClusters[$i][2]\t$InputClusters[$i][3]\t$InputClusters[$i][4]\t$InputClusters[$i][5]\n";
		}
	}
	close (OUT);

	######################################################################################

	print "Preparing the pymol script\n";

	my $fh = new FileHandle;
	my $outFilename = "DensityVisual.$this->{clusters_file_name_only}.$PDB.pml";
	$fh->open( $outFilename , "w"  );
	$fh->print( "reinitialize everything;\n" );
	$fh->print( "load http://www.rcsb.org/pdb/files/".$PDB.".pdb;\n" );

	$fh->print( "viewport 480,480;\n" );
	$fh->print( "preset.publication(\"".$PDB."\");\n" );
	$fh->print( "#show mesh;\n" );
	$fh->print( "\n" );
	$fh->print( "bg_color white;\n" );
	$fh->print( "\n" );
	foreach my $chain ( keys %{$this->{chains}} ) {
		my $chain_color = shift @{$ChainColors};
		$fh->print( "color ".$chain_color.", chain ".$chain.";\n" );
	}

	my ( $CurrentSC, $CurrentLevel, $CurrentSub, $PrevEntry );

	for (my $i = 0; $i < scalar @InputClusters; $i++) { # Find the first variant which exists in the given PDB structure
		if ( exists $this->{mutations}->{$InputClusters[$i][1].":".$InputClusters[$i][2]}) {
			my @IDs = split(/\./,$InputClusters[$i][0]);
			$CurrentSC = $IDs[0];
			$CurrentLevel = $IDs[1];
			$CurrentSub = $IDs[2];
			$PrevEntry = $i;
			last;
		}
	}

	if (not defined $PrevEntry) {
		print "******** None of the variants appear in the given PDB structure *********\n";
	}

	my $ColorIndex = 0; # index to keep track of colors

	for (my $i = $PrevEntry; $i < scalar @InputClusters; $i++) { 
		my @IDs = split(/\./,$InputClusters[$i][0]);
		my $variant = $InputClusters[$i][1].":".$InputClusters[$i][2];

		if ( exists $this->{mutations}->{$InputClusters[$i][1].":".$InputClusters[$i][2]}) {
			if ($IDs[0] == $CurrentSC) {
				if ($IDs[2] == $CurrentSub && $IDs[1] == $CurrentLevel) {
					setColors($variant, $InputClusters[$i], $this, $IDs[1], $fh, $currentColor);
				}
				else {
					if ($CurrentSub != 0) { # previous entry was NOT the last entry of a super cluster 
						if ($InputClusters[$PrevEntry][5] eq 1 || $InputClusters[$PrevEntry][5] eq 2) { # No clusters covered below by this sub
							$fh->print( "create S.".$currentColor."_".$InputClusters[$PrevEntry][0].", ".$currentColor."_".$InputClusters[$PrevEntry][0]."\n" ); 
							$fh->print( "set surface_color, ".$currentColor.", S.".$currentColor."_".$InputClusters[$PrevEntry][0]."\n" );
							$fh->print( "show surface, S.".$currentColor."_".$InputClusters[$PrevEntry][0]."\n" ); 
						}
						else {
							foreach my $sub (split(":",$InputClusters[$PrevEntry][5])) {
								$fh->print( "sele ".$currentColor."_".$InputClusters[$PrevEntry][0].", (".$currentColor."_".$InputClusters[$PrevEntry][0].", ".$this->{ClusterColors}->{$sub}."_".$sub.");\n" ); 
							}
							$fh->print( "create S.".$currentColor."_".$InputClusters[$PrevEntry][0].", ".$currentColor."_".$InputClusters[$PrevEntry][0]."\n" ); 
							$fh->print( "set surface_color, ".$currentColor.", S.".$currentColor."_".$InputClusters[$PrevEntry][0]."\n" );
							$fh->print( "show surface, S.".$currentColor."_".$InputClusters[$PrevEntry][0]."\n" );							
						}
						$CurrentLevel = $IDs[1];
						$CurrentSub = $IDs[2];
					}
					else { # previous entry WAS the last entry of a super cluster 
						$CurrentLevel = $IDs[1];
						$CurrentSub = $IDs[2];
					}
					# Done with taking care of the end of previous sub cluster. Now change the color and record the new sub:
					$ColorIndex++;
					$currentColor = $SurfaceColors->[$ColorIndex]; 
					setFirstColors($variant, $InputClusters[$i], $this, $IDs[1], $fh, $currentColor); # set color for the current entry
				}
				
			}
			else {
				$ColorIndex = 0; # Start of a new super cluster. Start from the begining
				$currentColor = $SurfaceColors->[$ColorIndex]; 
				setFirstColors($variant, $InputClusters[$i], $this, $IDs[1], $fh, $currentColor); 
				$CurrentSC = $IDs[0];
				$CurrentLevel = $IDs[1];
				$CurrentSub = $IDs[2];
			}
		$PrevEntry = $i;
		}	
	}
	$fh->close();
	print "pymol script written\n";
	print "Done.\n";

	#print Dumper \@InputClusters;
}

#########
### Density Visual functions
#########

sub test3 {
	my ($fh, $ref) = @_;
	$fh->print($ref->[0].";\n");
}

sub setColors {
	my ( $variant, $ClusterArrayRef, $this, $SecondDigit, $fh, $currentColor ) = @_; # $ClusterArrayRef = $InputClusters[$i]
	$variant =~ /(\w+)\:\D\.(\D+\d+)\D/g;
	my $GeneName = $1;
	my $MutName = $2;

	foreach my $key (keys %{$this->{mutations}->{$variant}}) {
		my @ChainRes = split(":", $key);
		my $chain = shift @ChainRes;
		my $res = shift @ChainRes;

		if ($SecondDigit == 0) { # super cluster
			# #$fh->print( "sele ".$currentColor."_".$GeneName."_".$MutName."_".$chain.", (resi ".$res." and chain ".$chain.");\n" ); #CHECK whether needed
			# $fh->print( "color ".$currentColor.", (resi ".$res." and chain ".$chain.");\n" );
			# #$fh->print( "show spheres, (resi ".$res." and chain ".$chain.");\n" );

			### ***** The following three lines show spheres with sand color for all the variants in the super cluster. Some of these variants get assigned another color later.
			### TO DO: select only outer variants in the super cluster; Add a surface to the last sub cluster
			$fh->print( "sele sand_".$GeneName."_".$MutName."_".$chain.", (resi ".$res." and chain ".$chain.");\n" ); #CHECK whether needed
			$fh->print( "color sand, (resi ".$res." and chain ".$chain.");\n" );
			$fh->print( "show spheres, (resi ".$res." and chain ".$chain.");\n" );
		}
		elsif ($SecondDigit == 1) { # first level
			$fh->print( "color ".$currentColor.", (resi ".$res." and chain ".$chain.");\n" );
			$fh->print( "show spheres, (resi ".$res." and chain ".$chain.");\n" );
			$fh->print( "sele ".$currentColor."_".$GeneName."_".$MutName."_".$chain.", (resi ".$res." and chain ".$chain.");\n" );
			$fh->print( "sele ".$currentColor."_".$ClusterArrayRef->[0].", (".$currentColor."_".$ClusterArrayRef->[0].", ".$currentColor."_".$GeneName."_".$MutName."_".$chain.");\n" ); # add to the object named by the cluster ID
			$this->{mutations}->{$variant}->{$chain.":".$res} = 1; #$ClusterArrayRef->[6] = 1; # set as processed
			$this->{ClusterColors}->{$ClusterArrayRef->[0]} = $currentColor; # To keep track of the color of a given cluster ID. Usefull in setting surfaces
		}
		elsif ($SecondDigit > 1) { # higher levels   
			if ( $ClusterArrayRef->[5] eq 2 ) { # this subcluster doesn't cover anything below it (no need to check PROCESSED)
				$fh->print( "color ".$currentColor.", (resi ".$res." and chain ".$chain.");\n" );
				$fh->print( "show spheres, (resi ".$res." and chain ".$chain.");\n" );
				$fh->print( "sele ".$currentColor."_".$GeneName."_".$MutName."_".$chain.", (resi ".$res." and chain ".$chain.");\n" );
				$fh->print( "sele ".$currentColor."_".$ClusterArrayRef->[0].", (".$currentColor."_".$ClusterArrayRef->[0].", ".$currentColor."_".$GeneName."_".$MutName."_".$chain.");\n" ); # add to the object named by the cluster ID
				$this->{mutations}->{$variant}->{$chain.":".$res} = 1; #$ClusterArrayRef->[6] = 1; # set as processed
				$this->{ClusterColors}->{$ClusterArrayRef->[0]} = $currentColor; # To keep track of the color of a given cluster ID. Usefull in setting surfaces
			}
			elsif ( $ClusterArrayRef->[5] ne 2 ) {
				if ($this->{mutations}->{$variant}->{$chain.":".$res} == 0) {
					$fh->print( "color ".$currentColor.", (resi ".$res." and chain ".$chain.");\n" );
					$fh->print( "show spheres, (resi ".$res." and chain ".$chain.");\n" );
					$fh->print( "sele ".$currentColor."_".$GeneName."_".$MutName."_".$chain.", (resi ".$res." and chain ".$chain.");\n" );
					$fh->print( "sele ".$currentColor."_".$ClusterArrayRef->[0].", (".$currentColor."_".$ClusterArrayRef->[0].", ".$currentColor."_".$GeneName."_".$MutName."_".$chain.");\n" ); # add to the object named by the cluster ID
					$this->{mutations}->{$variant}->{$chain.":".$res} = 1; #$ClusterArrayRef->[6] = 1; # set as processed
					$this->{ClusterColors}->{$ClusterArrayRef->[0]} = $currentColor; # To keep track of the color of a given cluster ID. Usefull in setting surfaces
				}
				else {$fh->print( "sele ".$currentColor."_".$ClusterArrayRef->[0].", (".$currentColor."_".$ClusterArrayRef->[0].", resi ".$res." and chain ".$chain.");\n" );} # add to the object named by the cluster ID	
				$this->{ClusterColors}->{$ClusterArrayRef->[0]} = $currentColor; # To keep track of the color of a given cluster ID. Usefull in setting surfaces			
			}
		}		
	}	
}

sub setFirstColors {
	my ( $variant, $ClusterArrayRef, $this, $SecondDigit, $fh, $currentColor ) = @_; # $ClusterArrayRef = $InputClusters[$i]
	$variant =~ /(\w+)\:\D\.(\D+\d+)\D/g;
	my $GeneName = $1;
	my $MutName = $2;

	foreach my $key (keys %{$this->{mutations}->{$variant}}) {
		my @ChainRes = split(":", $key);
		my $chain = shift @ChainRes;
		my $res = shift @ChainRes;

		if ($SecondDigit == 0) { # super cluster
			$fh->print( "\n" );
			$fh->print( "color ".$currentColor.", (resi ".$res." and chain ".$chain.");\n" );
		}
		elsif ($SecondDigit == 1) { # first level
			$fh->print( "\n" );
			$fh->print( "color ".$currentColor.", (resi ".$res." and chain ".$chain.");\n" );
			$fh->print( "show spheres, (resi ".$res." and chain ".$chain.");\n" );
			$fh->print( "sele ".$currentColor."_".$GeneName."_".$MutName."_".$chain.", (resi ".$res." and chain ".$chain.");\n" );
			$fh->print( "sele ".$currentColor."_".$ClusterArrayRef->[0].", ".$currentColor."_".$GeneName."_".$MutName."_".$chain.";\n" ); # add to the object named by the cluster ID
			$this->{mutations}->{$variant}->{$chain.":".$res} = 1; #$ClusterArrayRef->[6] = 1; # set as processed
			$this->{ClusterColors}->{$ClusterArrayRef->[0]} = $currentColor; # To keep track of the color of a given cluster ID. Usefull in setting surfaces
		}
		elsif ($SecondDigit > 1) { # higher levels 
			$fh->print( "\n" );  
			if ( $ClusterArrayRef->[5] eq 2 ) { # this subcluster doesn't cover anything below it (no need to check PROCESSED)
				$fh->print( "color ".$currentColor.", (resi ".$res." and chain ".$chain.");\n" );
				$fh->print( "show spheres, (resi ".$res." and chain ".$chain.");\n" );
				$fh->print( "sele ".$currentColor."_".$GeneName."_".$MutName."_".$chain.", (resi ".$res." and chain ".$chain.");\n" );
				$fh->print( "sele ".$currentColor."_".$ClusterArrayRef->[0].", ".$currentColor."_".$GeneName."_".$MutName."_".$chain.";" );  # add to the object named by the cluster ID
				$this->{mutations}->{$variant}->{$chain.":".$res} = 1; #$ClusterArrayRef->[6] = 1; # set as processed
				$this->{ClusterColors}->{$ClusterArrayRef->[0]} = $currentColor; # To keep track of the color of a given cluster ID. Usefull in setting surfaces
			}
			elsif ( $ClusterArrayRef->[5] ne 2 ) {
				if ($this->{mutations}->{$variant}->{$chain.":".$res} == 0) {
					$fh->print( "color ".$currentColor.", (resi ".$res." and chain ".$chain.");\n" );
					$fh->print( "show spheres, (resi ".$res." and chain ".$chain.");\n" );
					$fh->print( "sele ".$currentColor."_".$GeneName."_".$MutName."_".$chain.", (resi ".$res." and chain ".$chain.");\n" );
					$fh->print( "sele ".$currentColor."_".$ClusterArrayRef->[0].", ".$currentColor."_".$GeneName."_".$MutName."_".$chain."\n" ); # add to the object named by the cluster ID
					$this->{mutations}->{$variant}->{$chain.":".$res} = 1; #$ClusterArrayRef->[6] = 1; # set as processed
					$this->{ClusterColors}->{$ClusterArrayRef->[0]} = $currentColor; # To keep track of the color of a given cluster ID. Usefull in setting surfaces
				}
				else {$fh->print( "sele ".$currentColor."_".$ClusterArrayRef->[0].", ( resi ".$res." and chain ".$chain.");\n" );} # add to the object named by the cluster ID
				$this->{ClusterColors}->{$ClusterArrayRef->[0]} = $currentColor; # To keep track of the color of a given cluster ID. Usefull in setting surfaces
				
			}
		}		
	}	
}

sub GetFileName {
    my $this = shift;

    my @tempArray = split( "/",$this->{_CLUSTERS_FILE}) ;
    $this->{clusters_file_name_only} = pop @tempArray;

    return $this;
}

sub help_text{
    my $this = shift;
        return <<HELP

Usage: hotspot3d visual [options]

                             REQUIRED
--clusters-file              Clusters file
--pdb                        PDB ID on which to view clusters

                             AT LEAST ONE
--pairwise-file              A .pairwise file
--musites-file               A .musites file
--sites-file                 A .sites file
--drug-pairs-file            A .drug*clean file (either target or nontarget)

                             OPTIONAL
--output-file                Output filename for single PyMol script, default: hotspot3d.visual.pml
--pymol                      PyMoL program location, default: /usr/bin/pymol
--output-dir                 Output directory for multiple PyMol scripts, current working directory
--pdb-dir                    PDB file directory, default: current working directory
--bg-color                   background color, default: white
--mut-color                  mutation color, default: red
--mut-style                  mutation style, default: spheres
--site-color                 site color, default: blue
--site-style                 site style, default: sticks
--compound-color             compound color, default: util.cbag
--compound-style             compound style, default: sticks if compound-color, util.cbag otherwise
--script-only                If included (on), pymol is not run with new <output-file> when finished, default: off
--clusters-file-type         which clustering module created your clusters-file? network or density, default: network

--help                       this message

Tip: To run an already created .pml file, run pymol <your output-file>

HELP

}

1;

