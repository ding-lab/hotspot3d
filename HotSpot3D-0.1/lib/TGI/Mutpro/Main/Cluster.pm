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
    $this->{'data_location_file'} = undef;
    $this->{'output_file'} = 'HotSpot3D_results_cluter';
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
        'data-location-file=s' => \$this->{'data_location_file'},
        'p-value-cutoff=f' => \$this->{'p_value_cutoff'},
        'output-file=s' =>\$this->{'output_file'},
        'help' => \$help,
    );
    if ( $help ) { print STDERR help_text(); exit 0; }
    unless( $options ) { die $this->help_text(); }
    unless( $this->{'inter_intra_proximity_file'} ) { warn 'You must provide a conbined proximity file ! ', "\n"; die $this->help_text(); }
    unless( -e $this->{'inter_intra_proximity_file'} ) { warn ' combined proximity file is not exist  ! ', "\n"; die $this->help_text(); }
    unless( $this->{'data_location_file'} ) { warn 'You must provide location data file ! ', "\n"; die $this->help_text(); }
    unless( -e $this->{'data_location_file'} ) { warn ' location data file is not exist  ! ', "\n"; die $this->help_text(); }
    ## processing procedure
    my %clusterings = (); my %distance_matrix = ();
    if ( $this->{'target_nontarget_file'} and -e $this->{'target_nontarget_file'} ) {
        print STDOUT "\nWorking on HotSpot3D druggable ..... \n";
        my $fh = new FileHandle;
        unless( $fh->open( $this->{'target_nontarget_file'} ) ) { die "Could not open drug data file $! \n" };
        my ( %pdb_loc, %aa_map, @master, ); @master = ();
        map { 
            unless( /^"Drugport"/ ) {
                chomp;
                my @t = split /\t/; 
                map{ $_ =~ s/"//g } @t;
                my ( $drug, $pdb, $gene2, $m2, $loc, $infor_3d, ) = @t[0,2,5,9,11,15];
                my ( $dist, $pdb2, $pval, ) = split / /, $infor_3d;
                $pval =~ s/"//g;
                my @info = ($drug, $pdb, $gene2, $m2, $pdb2, $dist, $pval); #necessary info
                push @master, \@info; #store info
                $this->redundant(\%pdb_loc, \%aa_map, $gene2, $m2, $loc); #filter transcripts
                map {
                    my ( $drug , $pdb , $gene2 , $m2 , $pdb2 , $dist , $pval ) = @{$_};
                    if ( grep{$_ eq $m2} @{$aa_map{$gene2}} ) { 
                        my @mutations = ();
                        my $first = $drug.":".$pdb;
                        push @mutations , $first;
                        my $second = $gene2.":".$m2;
                        push @mutations , $second; #@mus2;
                        $this->AHC( $pval , $this->{'p_value_cutoff'} , \%clusterings , \@mutations );
                        $this->build_distance_matrix( $first , $second , $dist , \%distance_matrix );
                        $this->build_distance_matrix( $second , $first , $dist , \%distance_matrix );
                    }	
                } @master;
            }
        } $fh->getlines; 
        $fh->close();
    }
    my ( %pdb_loc, %aa_map, );
    my $fh = new FileHandle;
    unless( $fh->open( $this->{'data_location_file'} ) ) { die "Could not open data location file $! \n" };
    %pdb_loc = (); %aa_map = ();
    map {
        my ( $gene1, $aa_1, $loc_1, $gene2, $aa_2, $loc_2, ) = (split /\t/)[0,4,6,9,13,15];
        $this->redundant(\%pdb_loc, \%aa_map, $gene1, $aa_1, $loc_1);
        $this->redundant(\%pdb_loc, \%aa_map, $gene2, $aa_2, $loc_2);		
    } $fh->getlines;
    $fh->close();
    print STDOUT "\nWorking on HotSpot3D collapsed ... \n";
    unless( $fh->open( $this->{'inter_intra_proximity_file'} ) ) { die "Could not open proximity file $! \n" };
    map {
        chomp;
        my ( $gene1, $m1, $gene2, $m2, $dist, $pdb, $pval, ) = ( split /\t/ )[0,1,5,6,12,13,14];
        my @mus1 = split /,/, $m1;
        my @mus2 = split /,/, $m2;
        my @gm1 = (); my @gm2 = ();
        for ( my $i = 0 ; $i < @mus1 ; $i++ ) { 
            my $mut1 = $mus1[$i];
            my $mut2 = $mus2[$i];
            if ( (grep{$_ eq $mut1} @{$aa_map{$gene1}}) && (grep{$_ eq $mut2} @{$aa_map{$gene2}}) ) {
                my $first = $gene1.":".$mut1;
                my $second = $gene2.":".$mut2;
                push @gm1 , ( $first ); #unique mutation
                push @gm2 , ( $second ); #unique mutation
                $this->build_distance_matrix( $first , $second , $dist , \%distance_matrix );
                $this->build_distance_matrix( $second , $first , $dist , \%distance_matrix );
            }
        }
        my @mutations = @gm1;
        push @mutations , @gm2;
        $this->AHC( $pval , $this->{'p_value_cutoff'} , \%clusterings , \@mutations );
    } $fh->getlines;
    $fh->close();
    
    my $i = 0;
    my @keys = sort { $a <=> $b } keys %clusterings;
    ## REASSIGN CLUSTER IDS BY +1 INCREMENTS (REMOVE GAPS IN ID LIST)
    while ( $i < scalar keys %clusterings ) {
        if ( $keys[$i] != $i ) {
            $clusterings{$i} = $clusterings{$keys[$i]};
            delete $clusterings{$keys[$i]};
        }
        $i++;
    }
    
    @keys = sort { $a <=> $b } keys %clusterings;
    ## SHAVE CLUSTERS TO CORE
    #use distance_matrix matrix
    my %degree_connectivity = ();
    foreach ( keys %distance_matrix ) {
        $degree_connectivity{$_} = scalar keys %{$distance_matrix{$_}};
    }
    ## PREP THE OUTPUT FILENAME
    ## WRITE THE CLUSTER OUTPUT
    die "Could not create clustering output file\n" unless( $fh->open(">$this->{'output_file'}") );
    $fh->print( "Cluster\tGene/Drug\tMutation/PDB\tDegree_Connectivity\n" );
    foreach my $key ( keys %clusterings ) {
        foreach ( @{$clusterings{$key}} ) {
            my ( $gene , $mutation ) = split ( /\:/ , $_ );
            if ( defined $mutation ) {
                $fh->print( $key."\t".$gene."\t".$mutation."\t".$degree_connectivity{$_}."\n" );
            } else { $fh->print( "$key\t$gene\n" ); }
        }
    }
    my $numclusters = scalar keys %clusterings;
    print STDOUT "Found $numclusters clusters\n";

    return 1;
}

## CLUSTERING FUNCTION - AGGLOMERATIVE HIERARCHICAL CLUSTERING (AHC)
sub AHC {
    my ( $this, $pval , $pthreshold , $clusterings , $mutations ) = @_;
    if ( $pval < $pthreshold ) { #meets desired significance
        my ( @temp, @found, @combine, ); 
        foreach my $c ( keys %{$clusterings} ) { #each cluster
            my @mus_in_cluster = @{$clusterings->{$c}};
            foreach my $mu ( @{$mutations} ) { foreach ( @mus_in_cluster ) { if ( $mu eq $_ ) { push @combine , $c; } } }
        }
        my @uniqcombo = uniq @combine; #cluster types
        my ( @uniq, $c, );
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
    my ( $this, $pdb_loc, $aa_map, $gene, $aa, $loc, ) = @_;
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

## BUILD DISTANCE MATRIX - CAN USE ENTRY STATUS AS ADJACENCY MATRIX
sub build_distance_matrix {
    my ( $this, $m1 , $m2 , $dist , $distance_matrix ) = @_;
    $distance_matrix->{$m1}->{$m2} = $dist;

    return 1;
}

sub help_text{
    my $this = shift;
        return <<HELP

Usage: hotspot3d cluter [options]

--inter-intra-proximity-file     Both inter & intra molecular 3D proximity interactions
--data-location-file             Location data file
--output-file                    Output file
--target-nontarget-file          Both target & nontarget drug data (optional)
--p-value-cutoff		 P_value cutoff, default 0.05

--help			this message

HELP

}

1;
