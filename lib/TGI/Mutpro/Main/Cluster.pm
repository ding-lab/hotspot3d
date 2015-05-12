package TGI::Mutpro::Main::Cluster;
#
#----------------------------------
# $Authors: Beifang Niu & Adam Scott
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

sub new {
    my $class = shift;
    my $this = {};
    $this->{'inter_intra_proximity_file'} = undef;
    $this->{'target_nontarget_file'} = undef;
    $this->{'p_value_cutoff'} = 0.05;
    $this->{'linear_cutoff'} = 20;
    $this->{'pairwise_file'} = undef;
    $this->{'maf_file'} = undef;
	$this->{'vertex_type'} = 'recurrence';
	$this->{'max_radius'} = 10;
    $this->{'output_file'} = 'HotSpot3D_results.clusters';
    $this->{'stat'} = undef;
    bless $this, $class;
    $this->process();
    return $this;
}

sub process {
    my $this = shift;
    my ( $help, $options );
    unless( @ARGV ) { die $this->help_text(); }
    $options = GetOptions (
        'inter-intra-proximity-file=s' => \$this->{'inter_intra_proximity_file'},
        'target-nontarget-file=s' => \$this->{'target_nontarget_file'},
        'pairwise-file=s' => \$this->{'pairwise_file'},
        'maf-file=s' => \$this->{'maf_file'},
        'p-value-cutoff=f' => \$this->{'p_value_cutoff'},
        'linear-cutoff=f' => \$this->{'linear_cutoff'},
        'vertex-type=s' => \$this->{'vertex_type'},
        'max-radius=f' => \$this->{'max_radius'},
        'output-file=s' =>\$this->{'output_file'},
        'help' => \$help,
    );
    if ( $help ) { print STDERR help_text(); exit 0; }
    unless( $options ) { die $this->help_text(); }
    unless( $this->{'inter_intra_proximity_file'} ) { warn 'You must provide a collapsed pairwise file ! ', "\n"; die $this->help_text(); }
    unless( -e $this->{'inter_intra_proximity_file'} ) { warn ' collapsed pairwise file is not exist  ! ', "\n"; die $this->help_text(); }
    unless( $this->{'pairwise_file'} ) { warn 'You must provide location data file ! ', "\n"; die $this->help_text(); }
    unless( -e $this->{'pairwise_file'} ) { warn ' location data file does not exist  ! ', "\n"; die $this->help_text(); }
    unless( $this->{'maf_file'} and $this->{'vertex_type'} ne 'residue' ) { warn 'You must provide MAF file if not using residue vertex type ! ', "\n"; die $this->help_text(); }
    unless( -e $this->{'maf_file'} and $this->{'vertex_type'} ne 'residue' ) { warn ' MAF file does not exist  ! ', "\n"; die $this->help_text(); }
    ## processing procedure
	my ( %clusterings , %distance_matrix , %pdb_loc, %aa_map, @master, %mut_chrpos , %locations , %Variants );
	my $type = 1; #1 = recurrence, 0 = residue
#####
#	drug-mutation pairs
#####
    if ( $this->{'target_nontarget_file'} and -e $this->{'target_nontarget_file'} ) { #if drug pairs included
        print STDOUT "\nWorking on HotSpot3D drug pairs ..... \n";
        my $fh = new FileHandle;
        unless( $fh->open( $this->{'target_nontarget_file'} , "r" ) ) { die "Could not open drug pairs data file $! \n" };
        map { 
            unless( /^"Drugport"/ ) {
                chomp;
                my @t = split /\t/; 
                map{ $_ =~ s/"//g } @t;
                my ( $drug, $pdb, $gene2, $chr , $start , $stop , $m2, $loc, $infor_3d ) = @t[0,2,5,6,7,8,9,11,15];
                my ( $dist, $pdb2, $pval ) = split / /, $infor_3d;
				$mut_chrpos{$gene2.":".$m2}{$chr."_".$start."_".$stop} = 1;
                $pval =~ s/"//g;
                push @master, ( $drug , $pdb , $gene2 , $m2 , $pdb2 , $dist , $pval ); #store necessary pari info
                $this->redundant(\%pdb_loc, \%aa_map, $gene2, $m2, $loc); #filter transcripts
            }
        } $fh->getlines; 
        $fh->close();
		#Pick longest transcript representation of residue
		foreach my $gene ( keys %aa_map ) {
			foreach my $aa ( @{$aa_map{$gene}} ) {
				if ( $aa ne 'NA' ) {
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
				} #if not na
			} #foreach aa
		} #foreach gene
		#cluster drug-mutation pairs and build distance matrix
		map {
			my ( $drug , $pdb , $gene2 , $m2 , $pdb2 , $dist , $pval ) = @{$_};
			if ( grep{ $_ eq $m2 } @{$aa_map{$gene2}} ) { 
				my @mutations = ();
				my $first = $drug.":".$pdb;
				push @mutations , $first;
				my $second = $gene2.":".$m2;
				push @mutations , $second; #@mus2;
				$this->AHC( $pval , $this->{'p_value_cutoff'} , \%clusterings , \@mutations );
				if ( $pval < $this->{'p_value_cutoff'} ) {
					$distance_matrix{$first}{$second} = $dist;
					$distance_matrix{$second}{$first} = $dist;
				}
			}	
		} @master;
    } #if drug pairs included
#####
#	pairwise data
#####
	print STDOUT "\nReading in pairwise data ... \n";
    my $fh = new FileHandle;
    unless( $fh->open( $this->{'pairwise_file'} , "r" ) ) { die "Could not open pairwise file $! \n" };
    map {
        my ( $gene1 , $chr1 , $start1 , $stop1 , $aa_1 , $loc_1 , $gene2 , $chr2 , $start2 , $stop2 , $aa_2 , $loc_2 ) = (split /\t/)[0,1,2,3,4,6,9,10,11,12,13,15];
		$mut_chrpos{$gene1.":".$aa_1}{$chr1."_".$start1."_".$stop1} = 1;
		$mut_chrpos{$gene2.":".$aa_2}{$chr2."_".$start2."_".$stop2} = 1;
        $this->redundant(\%pdb_loc, \%aa_map, $gene1, $aa_1, $loc_1);
        $this->redundant(\%pdb_loc, \%aa_map, $gene2, $aa_2, $loc_2);		
    } $fh->getlines;
    $fh->close();
	#Pick longest transcript representation of residue
	foreach my $gene ( keys %aa_map ) {
		foreach my $aa ( @{$aa_map{$gene}} ) {
			if ( $aa ne 'NA' ) {
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
			} #if not na
		} #foreach aa
	} #foreach gene
#####
#	collapsed pairwise data
#####
    print STDOUT "\nWorking on collapsed data ... \n";
    unless( $fh->open( $this->{'inter_intra_proximity_file'} , "r" ) ) { die "Could not open collapsed file $! \n" };
    map {
        chomp;
        my ( $gene1, $m1, $gene2, $m2, $lindists , $dist, $pdb, $pval ) = ( split /\t/ )[0,1,5,6,11,12,13,14];
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
					} elsif ( $lin_dist < $this->{'linear_cutoff'} ) {
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
	if ( $this->{'vertex_type'} ne 'recurrence' and $this->{'vertex_type'} ne 'residue' ) {
		print STDERR "vertex_type option not recognized as \'recurrence\' or \'residue\'\n";
		print STDERR "Using default vertex_type = \'recurrence\'\n";
	} elsif ( $this->{'vertex_type'} eq 'residue' ) {
		$type = 0;
	}
	if ( $type == 1 ) {
		my %variants_from_pairs;
		foreach my $id ( keys %clusterings ) {
			foreach my $gene_mut ( @{$clusterings{$id}} ) {
				my ( $gene , $mutation ) = split /\:/ , $gene_mut;
				if ( exists $mut_chrpos{$gene_mut} ) {
					foreach my $css ( keys %{$mut_chrpos{$gene_mut}} ) {
						my @vinfo = split( "_" , $css ); #chr_start_stop
						my $variant = join( "_" , ( $gene , $mutation , @vinfo ) );
						$variants_from_pairs{$variant} = 1;
					}
				}
			}
		}
		##Mutation recurrence from MAF
		my %mutations;
		die "Could not open MAF file\n" unless( $fh->open( "<$this->{'maf_file'}" , "r" ) );
		map {
			chomp;
			if ( /missense/ or /in_frame/ ) {
				my ( $gene , $chr , $start , $stop , $barID , $aachange ) = @{$_}[0,4,5,6,15,47];
				my $variant = join( "_" , ( $gene , $aachange , $chr , $start , $stop ) );
				if ( exists $variants_from_pairs{$variant} ) {
					my $gene_aachange = $gene.":".$aachange;
					if ( exists $Variants{$gene_aachange} ) {
						if ( not exists $mutations{$variant}{$barID} ) {
							$Variants{$gene_aachange}++;
							$mutations{$variant}{$barID}++;
						}
					} else {
						$Variants{$gene_aachange} = 1;
						$mutations{$variant}{$barID} = 1;
					}
				}
			}
		} $fh->getline;
		$fh->close();
	} #if vertex_type
#####
#	write cluster output
#####
    die "Could not create clustering output file\n" unless( $fh->open( ">$this->{'output_file'}" , "w" ) );
    $fh->print( "Cluster\tGene/Drug\tMutation/Gene\tDegree_Connectivity\tCloseness_Centrality\tGeodesic_From_Centroid\tFrequency\n" );
	my %centrality = ();
	print STDOUT "Cluster ID & Centroid\n";
    foreach my $clus_num ( keys %clusterings ) {
		my $max = 0;
		my $centroid;
		my %dist = ();
		my @clus_mut = @{$clusterings{$clus_num}};
		#initialize geodesics
		foreach my $mut1 ( @clus_mut ) {
			my @mu1 = spli( ":" , $mut1 );
			foreach my $mut2 ( @clus_mut ) {
				my @mu2 = split( ":" , $mut2 );
				if ( $mu1[1] =~ /p\./ ) { $mu1[1] =~ s/\D//g; }
				if ( $mu2[1] =~ /p\./ ) { $mu2[1] =~ s/\D//g; }
				if ( exists $distance_matrix{$mut1}{$mut2} ) {
					$dist{$mut1}{$mut2} = $distance_matrix{$mut1}{$mut2};
				} elsif ( ( $mu1[0] eq $mu2[0] ) && ( $mu1[1] eq $mu2[1] ) ) {
					$dist{$mut1}{$mut2} = 0;
				} else {
					$dist{$mut1}{$mut2} = 1000000;
				}
			}
		}
		&floydwarshall( \%dist , \@clus_mut ); #get geodesics
		#calculate closeness centralities
		foreach my $current ( keys %dist ) {
			my $C = 0;
			foreach my $other ( keys %{$dist{$current}} ) {
				my $weight = 1; #stays as 1 if vertex_type eq 'residue'
				if ( $type ==1 ) {
					if ( exists $Variants{$other} ) {
						$weight = $Variants{$other};
					}
				}
				if ( $current ne $other ) {
					my $geodesic = $dist{$current}{$other};
					if ( $geodesic <= $this->{'linear_cutoff'} ) {
						$C += $weight/( 2**$geodesic );
					}
				} else {
					$C += $weight -1;
				}
			}
			$centrality{$clus_num}{$current} = $C;
			if ( $C > $max ) {
				$max = $C;
				$centroid = $current;
			}
		}
		print STDOUT "$clus_num\t$centroid\n";
		if ( exists $dist{$centroid} ) {
			foreach my $other ( keys %{$dist{$centroid}} ) {
				my $geodesic = $dist{$centroid}{$other};
				my $degrees = $degree_connectivity{$other};
				my $closenesscentrality = $centrality{$clus_num}{$other};
				if ( $geodesic <= $this->{'max_radius'} ) {
					my ( $gene , $mutation ) = split /\:/ , $other;
					my $weight = 1;
					if ( $this->{'vertex_type'} == 1 ) {
						if ( exists $Variants{$other} ) {
							$weight = $Variants{$other};
						}
					}
					$fh->print( join( "\t" , ( $clus_num , $gene , $mutation , $degrees , $closenesscentrality , $geodesic , $weight ) )."\n" );
				}
			} #foreach other vertex in network
        } #if dist for centroid
    } #foreach cluster
    my $numclusters = scalar keys %clusterings;
    print STDOUT "Found $numclusters clusters\n";

    return 1;
}
#####
#	sub functions
#####
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
        print STDOUT "New cluster $c\n";
    } #if pval significant

    return 1; 
}

sub redundant {
    my ( $this, $pdb_loc, $aa_map, $gene, $aa, $loc ) = @_;
    my $aa_orig = $aa;
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
            } elsif ( abs($prev_aa-$loc)==abs($aa-$loc) && !(grep(/^$aa_orig/,@{$aa_map->{$gene}})) ) {
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

Usage: hotspot3d cluter [options]

--inter-intra-proximity-file     Both inter & intra molecular 3D proximity interactions
--pairwise-file                  Location data file
--output-file                    Output file
--target-nontarget-file          Both target & nontarget drug data file (optional)
--p-value-cutoff                 P_value cutoff, default <0.05
--linear-cutoff                  Linear distance cutoff, default >20 residues
--max-radius                     Maximum cluster radius (max network geodesic from centroid), default <=10 Angstroms
--vertex-type                    Graph vertex type (recurrence or residue), default recurrence
--maf-file                       MAF file used in proximity search step (used if vertex-type = recurrence)

--help                           this message

NOTE: At least one of two pair files are needed from inter-intra-proximity-file/target-nontarget-file.

HELP

}

1;
