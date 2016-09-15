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

use Data::Dumper;

my $WEIGHT = "weight";
my $RECURRENCE = "recurrence";
my $UNIQUE = "unique";

sub new {
    my $class = shift;
    my $this = {};
    $this->{'pairwise_file'} = '3D_Proximity.pairwise';
    $this->{'collapsed_file'} = '3D_Proximity.pairwise.singleprotein.collapsed';
    $this->{'drug_clean_file'} = undef;
    $this->{'output_prefix'} = undef;
    $this->{'p_value_cutoff'} = 0.05;
    $this->{'linear_cutoff'} = 0;
	$this->{'max_radius'} = 10;
	$this->{'vertex_type'} = $RECURRENCE;
    $this->{'maf_file'} = undef;
    $this->{'amino_acid_header'} = "amino_acid_change";
    $this->{'transcript_id_header'} = "transcript_name";
    $this->{'weight_header'} = $WEIGHT;
    bless $this, $class;
    $this->process();
    return $this;
}

sub process {
    my $this = shift;
    my ( $help, $options );
    unless( @ARGV ) { die $this->help_text(); }
    $options = GetOptions (
        'output-prefix=s' => \$this->{'output_prefix'},
        'pairwise-file=s' => \$this->{'pairwise_file'},
        'collapsed-file=s' => \$this->{'collapsed_file'},
        'drug-clean-file=s' => \$this->{'drug_clean_file'},
        'p-value-cutoff=f' => \$this->{'p_value_cutoff'},
        'linear-cutoff=f' => \$this->{'linear_cutoff'},
        'max-radius=f' => \$this->{'max_radius'},
        'vertex-type=s' => \$this->{'vertex_type'},
        'maf-file=s' => \$this->{'maf_file'},
        'amino-acid-header=s' => \$this->{'amino_acid_header'},
        'transcript-id-header=s' => \$this->{'transcript_id_header'},		
        'weight-header=s' => \$this->{'weight_header'},		
        'help' => \$help,
    );
    if ( $help ) { print STDERR help_text(); exit 0; }
    unless( $options ) { die $this->help_text(); }
    if ( ( not defined $this->{'collapsed_file'} ) and ( not defined $this->{'drug_clean_file'} ) ) {
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
    unless( $this->{'pairwise_file'} ) { warn 'You must provide a pairwise file! ', "\n"; die $this->help_text(); }
    unless( -e $this->{'pairwise_file'} ) { warn "The input pairwise file (".$this->{'pairwise_file'}.") does not exist! ", "\n"; die $this->help_text(); }
	if ( $this->{'vertex_type'} ne $RECURRENCE and $this->{'vertex_type'} ne $UNIQUE and $this->{'vertex_type'} ne $WEIGHT ) {
		warn "vertex_type option not recognized as \'recurrence\', \'unique\', or \'weight\'\n";
		warn "Using default vertex_type = \'recurrence\'\n";
		$this->{'vertex_type'} = $RECURRENCE;
	}
	if ( $this->{'vertex_type'} ne $UNIQUE ) {
		unless( $this->{'maf_file'} ) { warn 'You must provide a .maf file if not using unique vertex type! ', "\n"; die $this->help_text(); }
		unless( -e $this->{'maf_file'} ) { warn "The input .maf file )".$this->{'maf_file'}.") does not exist! ", "\n"; die $this->help_text(); }
	}
    ## processing procedure
	my ( %clusterings , %distance_matrix , %pdb_loc, %aa_map, %master, %mut_chrpos , %locations , %Variants , $WEIGHT );
	$WEIGHT = "weight";
#####
#	drug-mutation pairs
#####
    if ( $this->{'drug_clean_file'} ) { #if drug pairs included
		print STDOUT "\nWorking on HotSpot3D drug pairs ..... \n";
		my $fh = new FileHandle;
		unless( $fh->open( $this->{'drug_clean_file'} , "r" ) ) { die "Could not open drug pairs data file $! \n" };
		my $rxi = 0;
		my $headline = $fh->getline(); chomp( $headline );
		my %drugcols = map{ ( $_ , $rxi++ ) } split( /\t/ , $headline );
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
				$mut_chrpos{$gene2.":".$m2}{$chr."_".$start."_".$stop} = 1;
				$pval =~ s/"//g;
				my $first = $drug.":".$gene2;
				my $second = $gene2.":".$m2;
				my $info = $dist.":".$pval;
				$master{$first}{$second} = $info; #store necessary pair info
				$this->redundant(\%pdb_loc, \%aa_map, \%locations , $gene2, $m2, $loc); #filter transcripts
		} $fh->getlines; 
		$fh->close();
		#Pick longest transcript representation of unique
		foreach my $gene ( keys %aa_map ) {
			foreach my $aa ( @{$aa_map{$gene}} ) {
				if ( $aa ne 'NA' ) {
					if ( exists $locations{$gene}{$aa} ) {
						my @locs = @{$locations{$gene}{$aa}};
						my $orig_loc = $aa;
						my $orig_letters = $aa;
						$orig_letters =~ s/[^A-Z]//g;
						$orig_loc =~ s/\D//g;
						my @pdb_locations=@{$pdb_loc{$gene}};
						foreach my $current_loc ( @locs ) {
							if ( grep{ $_ eq $current_loc } @pdb_locations ) {
								my @idx = grep{ $pdb_locations[$_] eq $current_loc } 0..$#pdb_locations;
								foreach my $current_idx ( @idx ) {
									my $current = $aa_map{$gene}[$current_idx];
									my $current_int = $current;
									my $current_letters = $current;
									$current_int =~ s/\D//g;
									$current_letters =~ s/[^A-Z]//g;
									if ( $orig_loc > $current_int && $orig_letters eq $current_letters ) {
										$aa_map{$gene}[$current_idx] = $aa;
									} #foreach good AA, keep track
								}
							}
						} #foreach location
					} #if location exists
				} #if not na
			} #foreach aa
		} #foreach gene
		#cluster drug-mutation pairs and build distance matrix
		foreach my $first ( keys %master ) {
			foreach my $second ( keys %{$master{$first}} ) {
				my ( $gene2 , $m2 ) = split ":" , $second;
				if ( grep{ $_ eq $m2 } @{$aa_map{$gene2}} ) { 
					my @mutations = ();
					push @mutations , $first;
					if ( $this->{'vertex_type'} eq $UNIQUE ) { $m2 =~ s/\D+(\d+)\D+/$1/g; }
					$second = $gene2.":".$m2;
					push @mutations , $second; #@mus2;
					my ( $dist , $pval ) = split ":" , $master{$first}{$second};
					$this->AHC( $pval , $this->{'p_value_cutoff'} , \%clusterings , \@mutations );
					if ( $pval < $this->{'p_value_cutoff'} ) {
						$distance_matrix{$first}{$second} = $dist;
						$distance_matrix{$second}{$first} = $dist;
					}
				}	
			}
		} 
	} #if drug pairs included
#####
#	pairwise data
#####
	print STDOUT "\nReading in pairwise data ... \n";
	my $fh = new FileHandle;
	unless( $fh->open( $this->{'pairwise_file'} , "r" ) ) { die "Could not open pairwise file $! \n" };
	map {
		my ( $gene1 , $chr1 , $start1 , $stop1 , $aa_1 , $loc_1 , $gene2 , $chr2 , $start2 , $stop2 , $aa_2 , $loc_2 ) = (split /\t/)[0..4,6,9..13,15];
		$mut_chrpos{$gene1.":".$aa_1}{$chr1."_".$start1."_".$stop1} = 1;
		$mut_chrpos{$gene2.":".$aa_2}{$chr2."_".$start2."_".$stop2} = 1;
		$this->redundant(\%pdb_loc, \%aa_map, \%locations , $gene1, $aa_1, $loc_1);
		$this->redundant(\%pdb_loc, \%aa_map, \%locations , $gene2, $aa_2, $loc_2);		
	} $fh->getlines;
	$fh->close();
	#Pick longest transcript representation of unique
	foreach my $gene ( keys %aa_map ) {
		foreach my $aa ( @{$aa_map{$gene}} ) {
			if ( $aa ne 'NA' ) {
				if ( exists $locations{$gene}{$aa} ) {
					my @locs = @{$locations{$gene}{$aa}};
					my $orig_loc = $aa;
					my $orig_letters = $aa;
					$orig_letters =~ s/[^A-Z]//g;
					$orig_loc =~ s/\D//g;
					my @pdb_locations=@{$pdb_loc{$gene}};
					foreach my $current_loc ( @locs ) {
						if ( grep{ $_ eq $current_loc } @pdb_locations ) {
							my @idx = grep{ $pdb_locations[$_] eq $current_loc } 0..$#pdb_locations;
							foreach my $current_idx ( @idx ) {
								my $current = $aa_map{$gene}[$current_idx];
								my $current_int = $current;
								my $current_letters = $current;
								$current_int =~ s/\D//g;
								$current_letters =~ s/[^A-Z]//g;
								if ( $orig_loc > $current_int && $orig_letters eq $current_letters ) {
									$aa_map{$gene}[$current_idx] = $aa;
								} #foreach good AA, keep track
							}
						}
					} #foreach location
				}
			} #if not na
		} #foreach aa
	} #foreach gene
#####
#	collapsed pairwise data
#####
    print STDOUT "\nWorking on collapsed data ... \n";
	if ( $this->{'collapsed_file'} ) {
		unless( $fh->open( $this->{'collapsed_file'} , "r" ) ) { die "Could not open collapsed file $! \n" };
		map {
			chomp;
			my ( $gene1, $m1, $gene2, $m2, $lindists , $dist, $pdb, $pval ) = ( split /\t/ )[0,1,5,6,11..14];
			my @mus1 = split /,/, $m1;
			my @mus2 = split /,/, $m2;
			my @gm1 = (); my @gm2 = ();
			for ( my $i = 0 ; $i < @mus1 ; $i++ ) { 
				my $mut1 = $mus1[$i];
				my $mut2 = $mus2[$i];
				if ( ( grep{ $_ eq $mut1 } @{$aa_map{$gene1}} ) && ( grep{ $_ eq $mut2 } @{$aa_map{$gene2}} ) ) {
					my $first = $gene1.":".$mut1;
					my $second = $gene2.":".$mut2;
					if ( $pval < $this->{'p_value_cutoff'} ) {
						my @lindists = split /\,/ , $lindists;
						my $lin_dist = $lindists[0];
						if ( $lin_dist eq "N/A" ) {
							push @gm1 , ( $first ); #unique mutation
							push @gm2 , ( $second ); #unique mutation
							$distance_matrix{$first}{$second} = $dist;
							$distance_matrix{$second}{$first} = $dist;
						} elsif ( $lin_dist > $this->{'linear_cutoff'} ) {
							push @gm1 , ( $first ); #unique mutation
							push @gm2 , ( $second ); #unique mutation
							$distance_matrix{$first}{$second} = $dist;
							$distance_matrix{$second}{$first} = $dist;
						} #if linear distance okay
					} #if spatial significance okay
				} #if okay transcript representations
			} #foreach transcript representation of mutations
			my @mutations = @gm1;
			push @mutations , @gm2;
			$this->AHC( $pval , $this->{'p_value_cutoff'} , \%clusterings , \@mutations );
		} $fh->getlines;
		$fh->close();
	} #if using collapsed pairs file
#####
#	clean up clusters
#####
    my $i = 0;
    ## REASSIGN CLUSTER IDS BY +1 INCREMENTS (REMOVE GAPS IN ID LIST)
	foreach ( keys %clusterings ) { if ( scalar @{$clusterings{$_}} == 0 ) { delete $clusterings{$_}; } } #assure empty clusters (if any) are removed
    my @keys = sort { $a <=> $b } keys %clusterings;
    while ( $i < scalar keys %clusterings ) {
        if ( $keys[$i] != $i ) {
            $clusterings{$i} = $clusterings{$keys[$i]};
            delete $clusterings{$keys[$i]};
        }
        $i++;
    }
    @keys = sort { $a <=> $b } keys %clusterings;
#####
#	finalize cluster data
#####
    ## SHAVE CLUSTERS TO CORE
    #use distance_matrix matrix
    my %degree_connectivity = ();
    foreach ( keys %distance_matrix ) {
        $degree_connectivity{$_} = scalar keys %{$distance_matrix{$_}};
    }
	if ( $this->{'vertex_type'} ne $UNIQUE ) {
		print STDOUT "\nPreparing to get recurrence or weight ...\n";
		my %variants_from_pairs;
		foreach my $id ( keys %clusterings ) {
			foreach my $gene_mut ( @{$clusterings{$id}} ) {
				my ( $gene , $mutation ) = split /\:/ , $gene_mut;
				if ( exists $mut_chrpos{$gene_mut} ) {
					foreach my $chr_start_stop ( keys %{$mut_chrpos{$gene_mut}} ) {
						my $variant = join( "_" , ( $gene , $mutation , $chr_start_stop ) );
						$variants_from_pairs{$variant} = 1;
					}
				}
			}
		}
		##Mutation recurrence or weight from MAF
		my %mutations;
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
			unless( defined( $mafcols{$this->{"weight_header"}} ) ) { die "\n"; };
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
				if ( exists $variants_from_pairs{$variant} ) {
					my $gene_aachange = $gene.":".$aachange;
					if ( exists $Variants{$gene_aachange} ) {
						if ( $this->{'vertex_type'} ne $WEIGHT ) {
							if ( not exists $mutations{$variant}{$barID} ) {
								$Variants{$gene_aachange} += $weight;
								$mutations{$variant}{$barID} += $weight;
							}
						}
					} else {
						$Variants{$gene_aachange} = $weight;
						$mutations{$variant}{$barID} = $weight;
					}
				}
			}
		} $fh->getlines;
		$fh->close();
	} #if vertex_type
#####
#	write cluster output
#####
	my @outFilename;
	if ( $this->{'output_prefix'} ) {
		push @outFilename , $this->{'output_prefix'};
	} else {
		if ( $this->{'maf_file'} ) {
			push @outFilename , $this->{'maf_file'};
		}
		if ( $this->{'collapsed_file'} ) {
			push @outFilename , $this->{'collapsed_file'};
		}
		if ( defined $this->{'drug_clean_file'} ) {
			if ( $this->{'drug_clean_file'} ne '' and scalar @outFilename > 1 ) {
				push @outFilename , $this->{'drug_clean_file'};
			} elsif ( $this->{'drug_clean_file'} ne '' ) {
				push @outFilename , $this->{'drug_clean_file'};
			}
		}
		push @outFilename , $this->{'linear_cutoff'};
		push @outFilename , $this->{'p_value_cutoff'};
		push @outFilename , $this->{'max_radius'};
	}
	push @outFilename , "clusters";
	my $outFilename = join( "." , @outFilename );
    die "Could not create clustering output file\n" unless( $fh->open( $outFilename , "w" ) );
    $fh->print( "Cluster\tGene/Drug\tMutation/Gene\tDegree_Connectivity\tCloseness_Centrality\tGeodesic_From_Centroid\tRecurrence\n" );
	print STDOUT "Getting Cluster ID's & Centroids\n";
    foreach my $clus_num ( keys %clusterings ) {
		my @clus_mut = @{$clusterings{$clus_num}};
		$this->centroid(\%Variants,\%distance_matrix,\%degree_connectivity,$clus_num,\@clus_mut,$fh,0, 1);
    } #foreach cluster
    my $numclusters = scalar keys %clusterings;
    print STDOUT "Found $numclusters clusters\n";

    return 1;
}
#####
#	sub functions
#####
sub centroid{
		my ($this, $Variants,$distance_matrix,$degree_connectivity,$clus_num,$clus_mut,$fh,$recluster,$counter)=@_;
		my %dist = ();
		foreach my $mut1 ( @{$clus_mut} ) { #initialize geodesics
			my @mu1 = split( ":" , $mut1 );
			foreach my $mut2 ( @{$clus_mut} ) {
				my @mu2 = split( ":" , $mut2 );
				if ( $mu1[1] =~ /p\./ ) { $mu1[1] =~ s/\D//g; }
				if ( $mu2[1] =~ /p\./ ) { $mu2[1] =~ s/\D//g; }
				if ( exists $distance_matrix->{$mut1}->{$mut2} ) {
					$dist{$mut1}{$mut2} = $distance_matrix->{$mut1}->{$mut2};
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
						#print "$current\t$other\t$dist{$current}{$other}\n";	
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
		#print"$count\t$cluster_size\t$centroid\n";
		$count=0;
		#print STDOUT "$clus_num\t$centroid\n";
		if ( exists $dist{$centroid} ) {
			foreach $other ( keys %{$dist{$centroid}} ) {
				my $geodesic = $dist{$centroid}{$other};
				my $degrees = $degree_connectivity->{$other};
				my $closenesscentrality = $centrality{$clus_num}{$other};
				my ( $gene , $mutation ) = split /\:/ , $other;
				$weight = 1;

				if ( $geodesic <= $this->{'max_radius'} ) {
					if ( exists $Variants->{$other} ) {
						$weight = $Variants->{$other};
					}
					$count+=1;
					if ($recluster==1){
						$fh->print( join( "\t" , ( "$clus_num.$counter" , $gene , $mutation , $degrees , $closenesscentrality , $geodesic , $weight ) )."\n" );
						my $index=0;
						$index++ until $clus_mut->[$index] eq $other;
						splice(@{$clus_mut}, $index, 1);
					}
					else {
						$fh->print( join( "\t" , ( $clus_num , $gene , $mutation , $degrees , $closenesscentrality , $geodesic , $weight ) )."\n" );
					}
				}
			} #foreach other vertex in network
		} #if dist for centroid

		if ($count<2)
		{
			return;
		}

		if ($recluster==1) {
			$counter+=1;
			#print STDOUT "$counter\n";
			$this->centroid($Variants,$distance_matrix,$degree_connectivity,$clus_num,$clus_mut,$fh,$recluster, $counter);
		}
}


## CLUSTERING FUNCTION - AGGLOMERATIVE HIERARCHICAL CLUSTERING (AHC)
sub AHC {
    my ( $this, $pval , $pthreshold , $clusterings , $mutations ) = @_;
    if ( $pval < $pthreshold ) { #meets desired significance
        my ( @temp, @found, @combine ); 
        my ( @uniq, $c );
        foreach $c ( keys %{$clusterings} ) { #each cluster
            my @mus_in_cluster = @{$clusterings->{$c}};
            foreach my $mu ( @{$mutations} ) { foreach ( @mus_in_cluster ) { if ( $mu eq $_ ) { push @combine , $c; } } }
        }
        my @uniqcombo = uniq @combine; #cluster types
        if ( scalar @uniqcombo > 0 ) { #collapse clusters into one
            my $collapse_to = min @uniqcombo; #cluster type
            my $j = 0; #iterator type
            while ( $j < scalar @uniqcombo ) {
                push @{$clusterings->{$collapse_to}} , @{$clusterings->{$uniqcombo[$j]}}; #mutation types
                push @{$clusterings->{$collapse_to}} , @{$mutations}; #
                if ( $collapse_to ne $uniqcombo[$j] ) { delete $clusterings->{$uniqcombo[$j]}; }
                $j++;
            }
            @{$clusterings->{$collapse_to}} = uniq @{$clusterings->{$collapse_to}};
        } else { #new cluster
            if ( scalar keys %{$clusterings} > 0 ) {
                $c = ( max keys %{$clusterings} ) + 1;
            } else { $c = 0; }
            push @temp , @{$mutations};
            @uniq = uniq @temp; #mutation types
            $clusterings->{$c} = \@uniq;
        }
        $c = scalar keys %{$clusterings};
        #print STDOUT "New cluster $c\n";
    } #if pval significant

    return 1; 
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
--p-value-cutoff             P_value cutoff (<), default: 0.05
--linear-cutoff              Linear distance cutoff (> peptides), default: 20
--max-radius                 Maximum cluster radius (max network geodesic from centroid, <= Angstroms), default: 10
--vertex-type                Graph vertex type (recurrence, unique, or weight), default: recurrence
--maf-file                   .maf file used in proximity search step (used if vertex-type = recurrence)
--transcript-id-header       .maf file column header for transcript id's, default: transcript_name
--amino-acid-header          .maf file column header for amino acid changes, default: amino_acid_change 
--weight-header              .maf file column header for mutation weight, default: weight (used if vertex-type = weight)

--help                       this message

HELP

}

1;
