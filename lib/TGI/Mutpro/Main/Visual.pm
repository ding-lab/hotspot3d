package TGI::Mutpro::Main::Visual;
#
#----------------------------------
# $Authors: Adam Scott & Beifang Niu 
# $Date: 2014-01-14 14:34:50 -0500 (Tue Jan 14 14:34:50 CST 2014) $
# $Revision:  $
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
    $this->{_COMPOUND_COLOR} = 'util.cbag';
    $this->{_MUT_STYLE} = 'spheres';
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
        'drug-pairs-file=s' => \$this->{_DRUG_PAIRS_FILE},
        'output-dir=s' => \$this->{_OUTPUT_DIR},
        'output-file=s' => \$this->{_OUTPUT_FILE},
        'pdb-dir=s' => \$this->{_PDB_DIR},
		'pdb=s' => \$this->{_PDB},
        'pymol=s' => \$this->{_PYMOL},
		'script-only' => \$this->{_SCRIPT_ONLY},
        'mut-color=s' => \$this->{_MUT_COLOR},
        'compound-color=s' => \$this->{_COMPOUND_COLOR},
        'mut-style=s' => \$this->{_MUT_STYLE},
        'compound-style=s' => \$this->{_COMPOUND_STYLE},
        'no-label' => \$this->{_NO_LABEL},
        'bg-color=s' => \$this->{_BG_COLOR},
        'help' => \$help
    );
    if ( $help ) { warn help_text(); exit 0; }
    unless( $options ) { die $this->help_text(); }
    if ( ( not $this->{_PAIRWISE_FILE} ) && ( not $this->{_DRUG_PAIRS_FILE} ) ) { #no input
		warn 'You must provide a pairwise file or a drug pairs file! ', "\n";
		die $this->help_text();
	}
	if ( not $this->{_DRUG_PAIRS_FILE} ) { #have drug pairs input
		if ( not -e $this->{_PAIRWISE_FILE} ) {
			warn ' pairwise file does not exist  ! ', "\n";
			die $this->help_text();
		}
	} elsif ( not $this->{_PAIRWISE_FILE} ) {
		if ( not -e $this->{_DRUG_PAIRS_FILE} ) {
			warn ' drug pairs file does not exist  ! ', "\n";
			die $this->help_text();
		}
	}
    unless( $this->{_CLUSTERS_FILE} ) { warn 'You must provide a cluster file ! ', "\n"; die $this->help_text(); }
    unless( -e $this->{_CLUSTERS_FILE} ) { warn ' cluster file is not exist  ! ', "\n"; die $this->help_text(); }
    unless( $this->{_OUTPUT_FILE} || $this->{_OUTPUT_DIR} ) { warn 'You must provide an output file or directory! ', "\n"; die $this->help_text(); }
    map{ $this->{_STAT}{$_} = 0; } qw( num_muts pdb pairs );
	
	$this->getMappingLocations(  );

	#$this->checkPDB(  );

	unless( $this->{mutations} ) { warn "No proximal pairs were found on ".$this->{_PDB}."\n"; }
	$this->getClusters(  );
	unless( $this->{clusters} ) { warn "No cluster to display on ".$this->{_PDB}."\n"; }
	$this->makePyMOLScript(  );

	$this->execute(  );

    return 1;
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
	if ( $this->{_PAIRWISE_FILE} ) {
		my $fh = new FileHandle;
		print STDOUT "Getting mutation-mutation pairs\n";
		unless( $fh->open( $this->{_PAIRWISE_FILE} ) ) { die "Could not open pairwise file\n" };
		while ( my $line = <$fh> ) {
			chomp( $line );
			if ( $line =~ /$this->{_PDB}/ ) {
				my ( $gene1 , $mu1 , $chain1 , $res1 , $gene2 , $mu2 , $chain2 , $res2 ) = ( split /\t/ , $line )[0,4,5,6,9,13,14,15];
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
		$fh->close();
	}
	if ( $this->{_DRUG_PAIRS_FILE} ) {
		my $dfh = new FileHandle;
		print STDOUT "Getting drug-mutation pairs\n";
		unless( $dfh->open( $this->{_DRUG_PAIRS_FILE} ) ) { die "Could not open drug pairs file\n" };
		while ( my $line = <$dfh> ) {
			chomp( $line );
			if ( $line =~ /$this->{_PDB}/ ) {
				my ( $drug , $chain1 , $res1 , $gene2 , $mu2 , $chain2 , $res2 ) = ( split /\t/ , $line )[0,3,4,5,9,10,11];
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
		$dfh->close();
	}
	return 1;
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
				and defined( $ccols{"Geodesic_From_Centroid"} ) ) {
				die "Not a valid clusters file!\n";
			}
			@ccols = (	$ccols{"Cluster"} ,
						$ccols{"Gene/Drug"} , 
						$ccols{"Mutation/Gene"} , 
						$ccols{"Closeness_Centrality"} , 
						$ccols{"Geodesic_From_Centroid"} );
		} else {
			my ( $id , $gd , $mp , $cc , $gfc ) = ( split /\t/ , $line )[@ccols];
			my $variant = $gd.":".$mp;
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
	return 1;
}

sub makePyMOLScript {
	my ( $this ) = @_;
	print STDOUT "Generating PyMOL script (.pml)\n";
	my $residue_line = {};
	my $chain_colors = {};
	my $old_gene = "";
	my $centroid_color = "grey";
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
					$residue =~ s/^p\.(.+)\D+$/$1/g;
					if ( ( grep{ $_ eq $residue } @cluster_AA_changes ) || ( grep{ $_ eq $AA_change } @cluster_AA_changes ) ) {
						if ( not exists $residue_line->{$chain}->{$location} ) {
							#print "\t\t\t\tMaking line ".$id.": ".$genedrug.":".join( "," , @location_AA_changes )." - ".$chain.", ".$location."\n";
							my $label = $cluster_color."_".$genedrug."_".$residue."_".$chain;
							$residue_line->{$chain}->{$location} = "sele ".$label.", (resi ".$location." and chain ".$chain."); ";
							if ( grep{ $_ eq $AA_change } @{$this->{centroids}->{$id}} ) { #is centroid
								$residue_line->{$chain}->{$location} .= "color ".$centroid_color.", ".$label."; ";
							} else { #is not centroid
								$residue_line->{$chain}->{$location} .= "color ".$cluster_color.", ".$label."; ";
							}
							$residue_line->{$chain}->{$location} .= "show ".$this->{_MUT_STYLE}.", ".$label."; ";
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
	my $outFilename = $this->{_CLUSTERS_FILE}.".".$this->{_PDB}.".pml";
	unless( $fh->open( $outFilename , "w" ) ) { die "Could not open output file, $outFilename\n"; }

	#viewer
	$fh->print( "reinitialize everything;\n" );
	$this->addPDB(  );
	#if ( -e $this->{_PDB_DIR}.$this->{_PDB} ) {
		$fh->print( "load ".$this->{_PDB}.".pdb;\n" );
	#	$fh->print( "#load http://www.rcsb.org/pdb/files/".$this->{_PDB}.".pdb.gz;\n" );
	#} else {
	#	$fh->print( "#load ".$this->{_PDB}.".pdb;\n" );
	#	$fh->print( "load http://www.rcsb.org/pdb/files/".$this->{_PDB}.".pdb.gz;\n" );
	#}
	$fh->print( "viewport 480,480;\n" );
	$fh->print( "preset.publication(\"".$this->{_PDB}."\");\n" );
	$fh->print( "#show mesh;\n" );
	$fh->print( "\n" );
	$fh->print( "bg_color ".$this->{_BG_COLOR}.";\n" );
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
		print "HotSpot3D ... Visualize!\n\n";
		system( "$pymol $pml" );
		#pymol <output-file>
		#-x disables external GUI module
	}
}

sub help_text{
    my $this = shift;
        return <<HELP

Usage: hotspot3d visual [options]

                             REQUIRED
--clusters-file              Clusters file
--pdb                        PDB ID on which to view clusters

                             AT LEAST ONE
--pairwise-file              Pairwise file
--drug-pairs-file            Drug pairs file (target/nontarget/hs3dd)

                             OPTIONAL
--output-file                Output filename for single PyMol script, default: hotspot3d.visual.pml
--pymol                      PyMoL program location, default: /usr/bin/pymol
--output-dir                 Output directory for multiple PyMol scripts, current working directory
--pdb-dir                    PDB file directory, default: current working directory
--bg-color                   background color, default: white
--mut-color                  mutation color, default: red
--mut-style                  mutation style, default: spheres
--compound-color             compound color, default: util.cbag
--compound-style             compound style, default: sticks if compound-color, util.cbag otherwise
--script-only                If included (on), pymol is not run with new <output-file> when finished, default: off

--help                       this message

Tip: To run an already created .pml file, run pymol <your output-file>

HELP

}

1;

