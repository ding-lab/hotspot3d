package TGI::Mutpro::Main::Density;
#
#----------------------------------
# $Authors: Amila Weerasinghe
# $Date: 2016-11-28 14:34:50 -0500 (Mon Nov 28 14:34:50 CST 2016) $
# $Revision: 2017-09-28 by A.W. (adding the recurrence and weight based clustering) 
# $URL: $
# $Doc: $ Determine density-based clusters from Hotspot3D pairwise data.
# 
#----------------------------------
#
use strict;
use warnings;

use LWP::Simple;
use FileHandle;
use Data::Dumper qw(Dumper);
use List::Util qw[min max shuffle];
use Getopt::Long;
use JSON;

our @ISA = qw( TGI::Mutpro::Main::Cluster );

my $EPSILONDEFAULT = 10;
my $MINPTSDEFAULT = 4;
my $NUMRUNSDEAFAULT = 10;
my $PROBCUTOFFDEFAULT = 100;

my $SHORTESTDISTANCE = "shortest";
my $AVERAGEDISTANCE = "average";
my $INDEPENDENT = "independent";
my $DEPENDENT = "dependent";
my $LOCAL = "local";
my $NONE = "none";

my $WEIGHT = "weight";
my $RECURRENCE = "recurrence";
my $UNIQUE = "unique";
my $SITE = "site";

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
    my $WEIGHT = "weight";

    if ( not defined $this->{'Epsilon'} ) {
    	warn "HotSpot3D::Cluster warning: no Epsilon value given, setting to default Epsilon value = 10\n";
			$this->{'Epsilon'} = $EPSILONDEFAULT;
    }
    if ( not defined $this->{'MinPts'} ) {
    	warn "HotSpot3D::Cluster warning: no MinPts value given, setting to default MinPts value = 4\n";
			$this->{'MinPts'} = $MINPTSDEFAULT;
    }
    if ( not defined $this->{'number_of_runs'} ) {
    	warn "HotSpot3D::Cluster warning: no number-of-runs given, setting to default value = 10\n";
			$this->{'number_of_runs'} = $NUMRUNSDEAFAULT;
    }
    if ( not defined $this->{'probability_cut_off'} ) {
    	warn "HotSpot3D::Cluster warning: no probability-cut-off given, setting to default value = 100\n";
			$this->{'probability_cut_off'} = $PROBCUTOFFDEFAULT;
    }
    print STDOUT "Preparing Density-based Clusters...\n";

    #####################################################
    
    $this->GetFileName(); # just to get only the name from the pairwise file destination
    $this->setRWD(); # to retrieve the path to TGI::Mutpro::Main. Used to call RScripts using perl

    #####################################################

	if ( $this->{'JSON_status'} ) {
		$mutations = readHash( $this->{'mutations_json_hash'} );
		$distance_matrix = readHash( $this->{'distance_matrix_json_hash'} );
        $siteVertexMap = readHash( $this->{'siteVertexMap_json_hash'} );	
	}
	else {
		unless( $this->{'pairwise_file'} ) { warn 'You must provide a pairwise file! ', "\n"; die $this->help_text(); }
    	unless( -e $this->{'pairwise_file'} ) { warn "The input pairwise file (".$this->{'pairwise_file'}.") does not exist! ", "\n"; die $this->help_text(); }

		$this->readMAF( $temp_mutations );
        # print "temp_mutations\n";
        # print Dumper $temp_mutations;
		$this->getDrugMutationPairs( $temp_distance_matrix );
		$this->getMutationMutationPairs( $temp_distance_matrix );
		$this->vertexFilter( $temp_mutations , $temp_distance_matrix , $mutations , $distance_matrix, $siteVertexMap );
        # $this->setSameSiteDistancesToZero( $distance_matrix, $mutations );
  
		printHash( $distance_matrix, "distance_matrix_hash" );
		printHash( $mutations, "mutations_hash" );
        printHash( $siteVertexMap, "siteVertexMap_hash" );

        # print "mutations\n";
        # print Dumper $mutations;
        # print "distance_matrix\n";
        # print Dumper $distance_matrix;
        # print "siteVertexMap\n";
        # print Dumper $siteVertexMap;
	}

    $this->{"siteVertexMap"} = $siteVertexMap; # store the reference to siteVertexMap
    $this->{"mutationsHashRef"} = $mutations; # store the reference to mutations 

    #print "distance matrix\n";
    #print Dumper $distance_matrix;
    #print "mutations\n";
    #print Dumper $mutations;
    #printHash( $distance_matrix, "distance_matrix_hash" );
    #printHash( $mutations, "mutations_hash" );
    #print "read from JSON\n";
    #readHash( $distance_matrix, "distance_matrix_hash.json" );
    #readHash( $mutations, "mutations_hash.json" );
    # print "distance matrix from JSON\n";
    # print Dumper $distance_matrix;
    # print "mutations from JSON\n";
    # print Dumper $mutations;

    $this->densityClustering( $mutations, $distance_matrix );

    my $numStructure = scalar keys %{ $distance_matrix };
    print "Density-based clusters are being calculated for $numStructure structures.\n\n";
   
}

#####
#   Functions
#####

sub densityClustering {
    my ( $this, $mutations, $distance_matrix ) = @_;

    print STDOUT "HotSpot3D::Cluster::Density::densityClustering\n";

    if ( $this->{'parallel'} eq $LOCAL ) {
        $this->localParallelDensityClustering( $mutations , $distance_matrix );
    #} elsif ( $this->{'parallel'} eq $BSUB ) {
    #   $this->bsubParallelNetworkClustering( $clusterings , $mutations , $distance_matrix );
    } else {
        $this->noParallelDensityClustering( $mutations , $distance_matrix );
    }
    return;
}

sub localParallelDensityClustering {
    my ( $this , $mutations , $distance_matrix ) = @_;

    my $pairwiseFN = "$this->{'distance_measure'}.$this->{pairwise_file_name_only}"; 
    print STDOUT "\tParallel clustering over structures with up to ".$this->{'max_processes'}." processes\n";
    my $pm = Parallel::ForkManager->new( $this->{'max_processes'} );

    DATA_LOOP:
    foreach my $structure ( keys %{$distance_matrix} ) { # run the density calculation for each available structure
        #print "Structure= $structure\n";
        # name output files as *.$pdbID.structure.*
        $this->{pairwise_file_name_only} = "$structure.Structure.$pairwiseFN";

        # call the SetOfNodes hash for each structure
        $this->{"CurrentSetOfNodes"} = $distance_matrix->{$structure};

        ###### Reference run: start
        $this->MainOPTICS( $structure , $distance_matrix, $mutations ); # perform OPTICS for the first time
        $this->RunSuperClustersID(); # perform Clustering for the reference run

        print "Reference run: Done.\nStart probability calculation for $structure\n";

        $this->getClusterProbabilities( $structure , $distance_matrix, $mutations ); # perform cluster-membership probability calculation

        print "\nProbability Calculation is Done for $structure.\n\n";
    }
    $pm->wait_all_children;

    return;
}

sub noParallelDensityClustering {
    my ( $this , $mutations , $distance_matrix ) = @_;
    print "Serially clustering over structures\n";

    my $pairwiseFN = "$this->{'distance_measure'}.$this->{pairwise_file_name_only}"; 
    foreach my $structure ( keys %{$distance_matrix} ) { # run the density calculation for each available structure
        #print "Structure= $structure\n";
        # name output files as *.$pdbID.structure.*
        $this->{pairwise_file_name_only} = "$structure.Structure.$pairwiseFN";

        # call the SetOfNodes hash for each structure
        $this->{"CurrentSetOfNodes"} = $distance_matrix->{$structure};

        ###### Reference run: start
        $this->MainOPTICS( $structure , $distance_matrix, $mutations ); # perform OPTICS for the first time
        $this->RunSuperClustersID(); # perform Clustering for the reference run

        print "Reference run: Done.\nStart probability calculation for $structure\n";

        $this->getClusterProbabilities( $structure , $distance_matrix, $mutations ); # perform cluster-membership probability calculation

        print "\nProbability Calculation is Done for $structure.\n\n";
    }
}

sub MainOPTICS {
    my $this = shift;
	my $structure = shift;
    my $distance_matrix = shift;
    my $mutations = shift;

    my $Epsilon = $this->{Epsilon};
    my $MinPts = $this->{MinPts};

    my $SetOfNodes = $this->{CurrentSetOfNodes};
	$this->resetProcessed( $structure , $SetOfNodes );# important in prob. runs; use the same set of nodes but start from the begining by setting processinfo to false

    #my $SetOfNodesRef = \%SetOfNodes; # just because Mac and Linux compatibility

    my @SetOfCores;
    my @SetOfEdges;
    my $CoreHash = {};

    foreach my $key ( keys %{$SetOfNodes} ) {
        if ( isThisAcore( $key, $SetOfNodes, $mutations, $Epsilon, $MinPts ) ) { # 1-core, 0-edge
            push @SetOfCores, $key;
            $CoreHash->{$key} = 0;
        }
        else {
            push @SetOfEdges, $key;
        }
    }
    @SetOfCores = shuffle @SetOfCores;
    @SetOfEdges = shuffle @SetOfEdges;
    my @SetOfCoresThenEdges = ( @SetOfCores, @SetOfEdges );

    # print "Set of Cores\n";
    # print Dumper \@SetOfCores;

    # my $neighbors1 = GetNeighbors( "SMAD4:18:48604788:48604788", $Epsilon, $SetOfNodes );
    # print "neighbors\n";
    # print Dumper $neighbors1;
    # my $coredist1 = GetCoreDistance( "SMAD4:18:48604788:48604788", $neighbors1, $MinPts, $mutations );
    # print "coredist = $coredist1\n";
    # if ( isThisAcore("SMAD4:18:48604788:48604788", $SetOfNodes, $mutations, $Epsilon, $MinPts )) {
    #     print "is a core\n";
    # }
    # else { print "not a core\n"};
    # print "distance_matrix\n";
    # print Dumper $distance_matrix;
    # print "mutations\n";
    # print Dumper $mutations;

    ###########################################################

    my @OrderedNodes;

    ################# Main OPTICS function ####################

    foreach my $p ( @SetOfCoresThenEdges ) { # p - a node(object)
        #print "first p=$p\n";
        if ( not $this->hasBeenProcessed($structure, $p) ) { # not processed yet
            ########## Expand Cluster Order ###########
            my %neighbors; # is a hash with keys neigbor indices whose values are mutual separations
            my %OrderSeeds; # is a hash to add seeds
            %neighbors = %{GetNeighbors($p, $Epsilon, $SetOfNodes)};
            $this->setProcessStatus($structure, $p, 1); # set as processed
            my $RD = undef; # reachability distance
            my $CD; # core distance
            $CD = GetCoreDistance($p, \%neighbors, $MinPts, $mutations);
            # print "p=$p and ";
            # print "CD=$CD\n";
            # print "neighbors\n";
            # hashQC($this, $mutations, \%neighbors);
            pushToOrderedNodesArray($p, $RD, $CD, $mutations, \@OrderedNodes, $this); # write to the file 

            if ( defined $CD ) {
                OrderSeedsUpdate($this, $structure, \%neighbors, $CD, \%OrderSeeds);
                # print "For p=$p, OrderSeeds= \n";
                # hashQC($this, $mutations, \%OrderSeeds);

                my $PrevObj = $p; # used to get the current obj. (To check whether variants are at the same location)
                while (scalar keys %OrderSeeds != 0) {
                    my @SeedKeys = sort { $OrderSeeds{$a} <=> $OrderSeeds{$b} } keys %OrderSeeds;
                    my @SeedValues = @OrderSeeds{@SeedKeys};

                    my $CurrentObject = GetCurrentObject(\@SeedValues, \@SeedKeys, $PrevObj, $CoreHash, $this, $mutations);
                    $PrevObj = $CurrentObject;
                    #print "\n\n current object= $CurrentObject\t neighbors=\n";
                    %neighbors = %{GetNeighbors($CurrentObject, $Epsilon, $SetOfNodes)};
                    #hashQC($this, $mutations, \%neighbors);

                    $this->setProcessStatus($structure, $CurrentObject, 1); # set as processed
                    $RD = $SeedValues[0];
                    $CD = GetCoreDistance($CurrentObject, \%neighbors, $MinPts, $mutations);
                    pushToOrderedNodesArray($CurrentObject, $RD, $CD, $mutations, \@OrderedNodes, $this); # write to the file 
                    delete $OrderSeeds{$CurrentObject};
                    if (defined $CD) {
                        #print "\tCurrent object is a core.(CD=$CD)\n Updated Order seeds list\n";
                        OrderSeedsUpdate($this, $structure, \%neighbors, $CD, \%OrderSeeds);
                        #hashQC($this, $mutations, \%OrderSeeds);
                    }
                }
            }
            # print "p=$p,(undefined CD), go to next\n";
        }
    }

    ### Replacing undefined RD by 10
    for (my $i = 0; $i < scalar @OrderedNodes; $i++) {
        if (not defined $OrderedNodes[$i][1]) {
            $OrderedNodes[$i][1] = 10;
        }
    }
    $this->{"CurrentRDarray"} = \@OrderedNodes;

    # my $OrderedFile = "./RD.$Epsilon.$MinPts.$this->{pairwise_file_name_only}";
    # open (OUT, ">$OrderedFile");
    # foreach my $x (1...scalar keys %SetOfNodes) {
    #     if (defined $OrderedNodes[$x-1][1]) {
    #         print OUT "$OrderedNodes[$x-1][0]\t $OrderedNodes[$x-1][1]\n";
    #     }
    #     else {
    #         print OUT "$OrderedNodes[$x-1][0]\t 10\n";
    #     }
    # }
    # close (OUT);


    # print "Current RD array\n";
    # print Dumper $this->{CurrentRDarray};

    return $this;
}

#####
##   Optics sub-methods
#####

sub isThisAcore {
    my ( $Obj, $SetOfNodes, $mutations, $Epsilon, $MinPts ) = @_;

    my $neighbors = GetNeighbors( $Obj, $Epsilon, $SetOfNodes );
    my $CoreDist = GetCoreDistance ( $Obj, $neighbors, $MinPts, $mutations );
    if ( defined $CoreDist ) {
        return 1;
    }
    else {
        return 0;
    }

}

sub GetNeighbors { # retrieves the epsilon-neighborhood
    my ( $Obj, $Epsilon, $SetOfNodes ) = @_;
    my $neighbors = {};
    foreach my $mutation_key ( keys %{$SetOfNodes->{$Obj}} ) {
        if ( $SetOfNodes->{$Obj}->{$mutation_key} <= $Epsilon ) {
            $neighbors->{$mutation_key} = $SetOfNodes->{$Obj}->{$mutation_key};
        }
    }
    return $neighbors;
}

sub GetCoreDistance {
    my ( $Obj, $neighbors, $MinPts, $mutations ) = @_;

    $neighbors->{$Obj} = 0; # adding the object so that it is counted towards the MinPts
    my $CoreDist = undef;
    my $neighbors_count = 0;
    my @keys = sort { $neighbors->{$a} <=> $neighbors->{$b} } keys %{$neighbors}; # sort keys according to distances
    # print "keys=\n";
    # print Dumper \@keys;

    for (my $i = 0; $i < scalar @keys; $i++) {
        my $mk = $keys[$i]; # neighboring mutation key
        foreach my $refAlt ( keys %{$mutations->{$mk}} ) {
            # print "$keys[$i]\n";
            # print "\tneighbors_count prev = $neighbors_count\n";
            my $pk = ( keys %{$mutations->{$mk}->{$refAlt}} )[0]; # get the first protein key corresponding to the ref:alt (there will be just one)
            my $weight = $mutations->{$mk}->{$refAlt}->{$pk}; # get the weight (for SITE=1, RECURRENCE=total frequency, WEIGHT=total weight)
            $neighbors_count += $weight; # add the number of "mutations" to the count (the number depends on what you're clustering on; as explained above)
        }
        #print "\tneighbors_count after = $neighbors_count\n";
        if ( $neighbors_count >= $MinPts ) {
            $CoreDist = $neighbors->{$keys[$i]};
            #print "\tCore dist= $CoreDist\n";
            last;
        }
    }
    return $CoreDist;
}

sub OrderSeedsUpdate {
    my ( $this, $structure, $neighbors, $CD, $OrderSeeds ) = @_;

    foreach my $mutation_key ( keys %{$neighbors} ) {
        if ( not $this->hasBeenProcessed( $structure, $mutation_key ) ) {
            my $new_r_dist = max ( $CD, $neighbors->{$mutation_key} );
            if ( exists $OrderSeeds->{$mutation_key} ) {
                if ( $new_r_dist < $OrderSeeds->{$mutation_key} ) {
                    $OrderSeeds->{$mutation_key} = $new_r_dist;
                }
            }
            else {
                $OrderSeeds->{$mutation_key} = $new_r_dist;
            }
        }
    }
}

sub CombineWords {
    my ($word1,$word2)=@_;
    return $word1.":".$word2;
}

sub GetCurrentObject { # To check whether variants are at the same location
    my ( $ValueSet, $KeySet, $PrevObj, $CoreHash, $this, $mutations ) = @_;
    my @SmallestKeys;
    for ( my $i = 0; $i < scalar @$ValueSet; $i++ ) {
        if ( $ValueSet->[$i] == $ValueSet->[0] ) {
            push @SmallestKeys, $KeySet->[$i];
        }
    }

    if ( scalar @SmallestKeys > 1 ) { # more than one variant has the smallest RD
        # if exists a core, get it first
        for ( my $i = 0; $i < scalar @SmallestKeys; $i++ ) {
            if ( exists $CoreHash->{$SmallestKeys[$i]} ) {
                my $movingKey = $SmallestKeys[$i];
                splice @SmallestKeys, $i, 1;
                unshift @SmallestKeys, $movingKey;
                last;
            }
        }

        # if there's a variant at the same position get it even before that
        foreach my $mutation_key ( @SmallestKeys ) {
            if ( $this->isSameProteinPosition( $mutations , $PrevObj , $mutation_key ) == 1 ) {
                unshift @SmallestKeys, $mutation_key;
                last;
            }
        }
    }
    return shift @SmallestKeys;
}

sub pushToOrderedNodesArray {
    my ( $Obj, $RD, $CD, $mutations, $OrderedNodes, $this ) = @_;

    my ( $gene, $chromosome, $start, $stop ) = @{$this->splitMutationKey($Obj)};

    foreach my $RefAlt ( sort keys %{$mutations->{$Obj}} ) {
        my @proteinKeys = sort keys %{$mutations->{$Obj}->{$RefAlt}};
        my $proteinKey = shift @proteinKeys; # take the first key, should be the same except transcripts
        my $weight = $mutations->{$Obj}->{$RefAlt}->{$proteinKey};
        my $altTranscriptColumn = join( "|", @proteinKeys);
        #print "To File: $proteinKey, $RD, $CD\n";
        my ( $ref, $alt ) = @{$this->splitRefAltKey($RefAlt)};
        my ( $transcript, $aa ) = @{$this->splitProteinKey($proteinKey)};

        my $gene_aa = join ( ":", $gene, $aa );
        my $genomicData = join ( "\t", $chromosome, $start, $stop, $ref, $alt, $transcript );

        my $genomicData_altTransColumn = join ( "\t", $genomicData, $altTranscriptColumn );

        push @{$OrderedNodes}, [$gene_aa, $RD, $CD, $genomicData_altTransColumn, $weight];
    }
}

sub hashQC { # used for QC'ing OPTICS. Works for both OrderedSeeds and neigbors hashes; prints out: mutation_keys,distances/RDval,protein_keys    
    my ( $this, $mutations, $OrderSeeds ) = @_;

    my @SeedKeys = sort { $OrderSeeds->{$a} <=> $OrderSeeds->{$b} } keys %{$OrderSeeds};
    
    foreach my $key (@SeedKeys) {
        my @proteinKeys;
        foreach my $refAlt ( sort keys %{$mutations->{$key}} ) {
            push @proteinKeys, (sort keys %{$mutations->{$key}->{$refAlt}})[0];
        }
        my $proteinKey = join ( "|", @proteinKeys);
        print "$key\t$OrderSeeds->{$key}\t$proteinKey\t$this->{processed}->{$key}\n";
    }
    print "\n";
}

sub GetFileName {
    my $this = shift @_;

    my @tempArray = split( "/",$this->{'pairwise_file'}) ;

    $this->{'pairwise_file_name_only'} = pop @tempArray;

    if ( not defined $this->{'pairwise_file_name_only'} or $this->{'pairwise_file_name_only'} eq '' ) {
        $this->{'pairwise_file_name_only'} = 'HS3D.DensityClustering';
    }

    return $this;
}

####################################################################################
############################### Clustering Method ##################################
####################################################################################

sub RunSuperClustersID {
    my $this = shift @_;

    my $MinPts = $this->{MinPts};
    my @InitialSet = @{$this->{CurrentRDarray}};
    my %Clusters;
    my @ClusterArray;


    # Identify super clusters:

    my $scs = 0; # super cluster start
    for (my $i = 1; $i < scalar @InitialSet; $i++) {
        #print "i=$i\n";
        if ( $InitialSet[$i][1] == 10 ) {
            #print "RD($i)=10\n";

            # if ( $InitialSet[$i-1][1] == 10 ) {
            #     $scs = $i;
            # }
            # else {
                $Clusters{SuperClusters}{$scs} = $i;
                $scs = $i;
            # }       
        }
    }
    #print Dumper \%Clusters;

    # Go inside each super cluster and start scanning from the smallest r(x):
    my $idSuperCluster = -1;

    foreach my $superC ( sort { $Clusters{SuperClusters}{$a} <=> $Clusters{SuperClusters}{$b} } keys %{$Clusters{SuperClusters}} ) {
        #print "Start of SC = $superC, $InitialSet[$superC][0]\n";
        my $MinRDat = $superC; # Lowest r(x) is at this x value.
        for (my $i = $superC; $i <= $Clusters{SuperClusters}{$superC} ; $i++) {
            if ($InitialSet[$i][1] < $InitialSet[$MinRDat][1]) {
                $MinRDat = $i;
            }
        }
        #print "\tMinRD=$InitialSet[$MinRDat][0],$InitialSet[$MinRDat][1]\n";
        my $MaxRD = 10;
        for (my $i = $superC+2; $i < $Clusters{SuperClusters}{$superC} ; $i++) {
            if ($InitialSet[$i-1][1] < $InitialSet[$i+1][1]) {
                $MaxRD = $InitialSet[$i][1];
            }
        }

        # Recording the Super Cluster
        my $superClusterCount = 0; # a counter to record the number of weights (recurrence/pathogenicity) inside a super cluster
        for (my $i = $superC; $i <= $Clusters{SuperClusters}{$superC}-1 ; $i++) {
            $superClusterCount += $InitialSet[$i][4]; # add the weight (for SITE=1, RECURRENCE=frequency, WEIGHT=weight)
        }
        # check if the current super cluster satisfies the MinPts condition
        if ( $superClusterCount >= $MinPts ) {
            $idSuperCluster++;
            my $ClusTot = 0;
            for (my $t = $superC+1; $t <= $Clusters{SuperClusters}{$superC}-1; $t++) { # start+1 b/c I want to get rid of the first tall peak, (-1 b/c end has included the start of the next)
                $ClusTot += $InitialSet[$t][1];
            }
            if ( $ClusTot == 0 ) { # according to the above calculation singletons will have a defined $ClusTot value equals to zero
                $ClusTot = 0.1 * $superClusterCount; # this way, singleton super clusters will have avg cluster density = 0.1
            }
            my $ClusAvg = $ClusTot/$superClusterCount;
            push @ClusterArray, [$superC,$Clusters{SuperClusters}{$superC}-1,9.9,$ClusAvg,$idSuperCluster."."."0"."."."0"]; # Artificially recording the super cluster. Randomly picked 9.9. (-1 b/c end has included the start of the next)
        }

        # Start scanning from epsilon prime= minimum r(x) to maximum r(x):

        my $nsubc = 0; # number of sub-clusters
        my $nsubcPre = 0; # previous number of sub-clusters
        my $TempSubClusterRef;
        my $idLevel = 1;

        for (my $Cutoff = $InitialSet[$MinRDat][1] ; $Cutoff <=  $MaxRD + 0.5 ; $Cutoff += 0.1) {
            my $idSubCluster = 1;
            my $n = $superC;
            my $counter = 0;
            my $s = $superC;
            while ($n < $Clusters{SuperClusters}{$superC}) { # n < (end of the super cluster) bc we don't want the last bar in the super cluster in a sub cluster
                if ($InitialSet[$n][1] < $Cutoff) {
                    $counter += $InitialSet[$n][4]; # add the weight (for SITE=1, RECURRENCE=frequency, WEIGHT=weight)
                    # $counter++;
                    if ($counter >= $MinPts) {
                        $Clusters{SubClusters}{$superC}{$s} = $n;
                    }
                    $n++;
                }
                else {
                    if ( $Cutoff == $InitialSet[$MinRDat][1] and $InitialSet[$n][4] >= $MinPts ) { # allowing singleton sub clusters to form at the lowest level
                        $Clusters{SubClusters}{$superC}{$n} = $n; 
                    }
                    # regular sub cluster stuff
                    $s = $n;
                    $counter = $InitialSet[$n][4]; # reset $counter with the next bar having RD>Cutoff
                    $n++;
                }
            }
            
            if (exists $Clusters{SubClusters}{$superC}) {  # if at leat one sub-cluster is found in the corresponding super cluster:
                $nsubc = scalar keys %{$Clusters{SubClusters}{$superC}};
                #print "At cutoff=$Cutoff\tnsubc = $nsubc\t nsubc_previous=$nsubcPre\n";
                #print "\tnew cluster\n\t";
                #print Dumper $Clusters{SubClusters}{$superC};
                
                            
                if ($nsubc != $nsubcPre) {
                    my @OrderedS = sort { $Clusters{SubClusters}{$superC}{$a} <=> $Clusters{SubClusters}{$superC}{$b} } keys %{$Clusters{SubClusters}{$superC}};
                    my @OrderedE = @{$Clusters{SubClusters}{$superC}}{@OrderedS};
                    foreach my $start (@OrderedS) {
                        my $ClusTot = 0;
                        my $tempCounter = 0;
                        if ( $start == $Clusters{SubClusters}{$superC}{$start} ) { # this sub cluster is a singleton
                            $ClusTot = 0.1; # singletons are given cluster avg = 0.1 #$InitialSet[$start][1];
                            $tempCounter = 1;
                        }
                        else { # this sub cluster is NOT a singleton
                            for (my $t = $start+1; $t <= $Clusters{SubClusters}{$superC}{$start}; $t++) { # start+1 b/c I want to get rid of the first tall peak
                                $ClusTot = $ClusTot + ($InitialSet[$t][1] * $InitialSet[$t][4]);
                                $tempCounter += $InitialSet[$t][4];
                            }
                        }
                        my $ClusAvg = $ClusTot/$tempCounter;
                        push @ClusterArray, [$start,$Clusters{SubClusters}{$superC}{$start},$Cutoff,$ClusAvg,$idSuperCluster.".".$idLevel.".".$idSubCluster];
                        $idSubCluster++;
                    }
                    $nsubcPre = $nsubc;
                    $idLevel++;
                    $TempSubClusterRef = $Clusters{SubClusters}{$superC}; # we store the current subclusters hash to detect overlapping and covering
                }
                elsif ($nsubc == $nsubcPre && defined $TempSubClusterRef) {
                    # print "\tinside the second loop\n";
                    # print "\tTemporary cluster\n\t";
                    # print Dumper $TempSubClusterRef;
                    my $OverlapCheck = 1;
                    foreach my $subC (keys %{$Clusters{SubClusters}{$superC}}) {
                        foreach my $subCpre (keys %{$TempSubClusterRef}) {
                            my ($ss, $se, $es, $ee); # s-start, e-end ; new cluster comes first
                            $ss = $subC - $subCpre;
                            $se = $subC - $TempSubClusterRef->{$subCpre};
                            $es = $Clusters{SubClusters}{$superC}{$subC} - $subCpre;
                            $ee = $Clusters{SubClusters}{$superC}{$subC} - $TempSubClusterRef->{$subCpre};
                            #print "\t\tFor subC=$subC, subCpre=$subCpre: ($ss, $es, $se, $ee)";
                            if ( $Clusters{SubClusters}{$superC}{$subC} - $subC != 0 and $TempSubClusterRef->{$subCpre} - $subCpre != 0 ) { # singletons are not considered for the overlap check
                                if (($ss*$es)> 0 && ($se*$ee) > 0) {
                                    $OverlapCheck = 0; # one of new sub clusters does not overlap with the previous sub cluster in consideration, check the next previous sub cluster
                                    #print "\tOverlapCheck = $OverlapCheck  continue checking\n";
                                }
                                else { # one of new sub clusters overlaps with one of the previous sub clusters (if check is necessary to avoid dividing by zeros)
                                    if ( (($Clusters{SubClusters}{$superC}{$subC} - $subC) / ($TempSubClusterRef->{$subCpre} - $subCpre)) >= 2 ) { # the expansion of the cluster is more than 100% --> record
                                        $OverlapCheck = 0;
                                    }
                                    else {
                                        $OverlapCheck = 1;
                                        #print "\tOverlapCheck = $OverlapCheck  stop checking, found one overlap\n";
                                    }
                                    last; # one of new sub clusters overlaps with one of the previous sub clusters (or its size increased by more than 100%), stop checking
                                }
                            }
                        }
                        if ($OverlapCheck == 0) {
                            #print "\tOverlapCheck=0 passed subC=$subC\n";
                            last; # one of new sub clusters does not overlap with any of the previous sub clusters (or there's an overlap, but its size increased by more than 100%), stop checking for other non-overlapping new clusters. One is enough!
                        }
                    }
                    if ($OverlapCheck == 0) {
                        my @OrderedS = sort { $Clusters{SubClusters}{$superC}{$a} <=> $Clusters{SubClusters}{$superC}{$b} } keys %{$Clusters{SubClusters}{$superC}};
                        my @OrderedE = @{$Clusters{SubClusters}{$superC}}{@OrderedS};
                        foreach my $start (@OrderedS) {
                            my $ClusTot = 0;
                            my $tempCounter = 0;
                            if ( $start == $Clusters{SubClusters}{$superC}{$start} ) { # this sub cluster is a singleton
                                $ClusTot = 0.1; # singletons are given cluster avg = 0.1 #$InitialSet[$start][1];
                                $tempCounter = 1;
                            }
                            else { # this sub cluster is NOT a singleton
                                for (my $t = $start+1; $t <= $Clusters{SubClusters}{$superC}{$start}; $t++) { # start+1 b/c I want to get rid of the first tall peak
                                    $ClusTot = $ClusTot + ($InitialSet[$t][1] * $InitialSet[$t][4]);
                                    $tempCounter += $InitialSet[$t][4];
                                }
                            }
                            my $ClusAvg = $ClusTot/$tempCounter;
                            push @ClusterArray, [$start,$Clusters{SubClusters}{$superC}{$start},$Cutoff,$ClusAvg,$idSuperCluster.".".$idLevel.".".$idSubCluster];
                            $idSubCluster++;
                        }
                        $idLevel++;
                        $TempSubClusterRef = $Clusters{SubClusters}{$superC}; # we store the current subclusters hash to detect overlapping and covering
                    }
                }
            }
            
            delete $Clusters{SubClusters}{$superC}; # since we track the number of subclusters by nsubc and nsubcPre we delete the hash--> important to correctly track the number
            #print "\n\n";
        }
        
        #print Dumper \%Clusters;   
    }
    #print Dumper \@ClusterArray;

    #################### Check for covering clusters #########################

    my @OrderedClusters;
    my @IDsArray; # 3D array with SuperCluster, level and SubCluster indices (Note that level and subcluster index starts from 1; @IDsArray[i][0][0] will correspond to ith SuperCluster)
    for (my $i = 0; $i < scalar @ClusterArray; $i++) {
        my @IDs = split(/\./,$ClusterArray[$i][4]);
        #print "$IDs[0],$IDs[1],$IDs[2]\t$ClusterArray[$i][4]\n";
        if ($IDs[2] != 0) {
            $IDsArray[$IDs[0]][$IDs[1]][$IDs[2]-1] = $i;
        }
        elsif ($IDs[2] == 0 && $IDs[1] == 0) {
            $IDsArray[$IDs[0]][$IDs[1]][$IDs[2]] = $i;
        }   
    }

    #print Dumper \@IDsArray;
    # my $pr1 = scalar @{$IDsArray[0]};
    # print "number of levels = $pr1\n";

    for (my $i = 0; $i < scalar @IDsArray; $i++) {
        for (my $j = 0; $j < scalar @{$IDsArray[$i]}; $j++) {
            for (my $k = 0; $k < scalar @{$IDsArray[$i][$j]}; $k++) {
                #print "$i, $j, $k\t$IDsArray[$i][$j][$k]\n";
                if ($j == 0) { # zeroth level is the supercluster
                    $ClusterArray[$IDsArray[$i][$j][$k]][5] = 0;
                }
                elsif ($j == 1) { # first level is not covering any
                    $ClusterArray[$IDsArray[$i][$j][$k]][5] = 1;
                }
                elsif ($j > 1) { # check for covering
                    my $CurrentS = $ClusterArray[$IDsArray[$i][$j][$k]][0];
                    my $CurrentE = $ClusterArray[$IDsArray[$i][$j][$k]][1];
                    my $FifthCol;

                    for (my $l = 0; $l < scalar @{$IDsArray[$i][$j-1]}; $l++) { # check whether the current sub covers any of the clusters in the previous level 
                        my $PrevS = $ClusterArray[$IDsArray[$i][$j-1][$l]][0];
                        my $PrevE = $ClusterArray[$IDsArray[$i][$j-1][$l]][1];

                        if ($PrevS >= $CurrentS && $PrevE <= $CurrentE) {
                            if (defined $FifthCol) {
                                $FifthCol = $FifthCol.":".$ClusterArray[$IDsArray[$i][$j-1][$l]][4];
                            }
                            else {
                                $FifthCol = $ClusterArray[$IDsArray[$i][$j-1][$l]][4];
                            }
                        }
                    }
                    if (defined $FifthCol) {
                        $ClusterArray[$IDsArray[$i][$j][$k]][5] = $FifthCol;
                    }
                    else {
                        $ClusterArray[$IDsArray[$i][$j][$k]][5] = 2;
                    }

                }
            }
        }
    }

    #print Dumper \@ClusterArray;

    ########################  Writing to files  ##############################

    my $OrderedFile = "./RD.$this->{Epsilon}.$this->{MinPts}.$this->{pairwise_file_name_only}.out";
    open (OUT, ">$OrderedFile");
    for (my $i = 0; $i < scalar @{$this->{CurrentRDarray}}; $i++) {
        my $ele1 = ${$this->{CurrentRDarray}}[$i][0]; # Gene:AAchange
        my $ele2 = ${$this->{CurrentRDarray}}[$i][1]; # RD
        my $genomicAnnotation = ${$this->{CurrentRDarray}}[$i][3];
        my $weight = ${$this->{CurrentRDarray}}[$i][4]; # weight
        print OUT "$ele1\t$ele2\t$genomicAnnotation\t$weight\n";
    }
    close (OUT);

    ########################################

    my $OutFile1 = "RD.$this->{Epsilon}.$this->{MinPts}.$this->{pairwise_file_name_only}.clusters.plot";
    open (OUT, ">$OutFile1");

    for (my $i = 0; $i < scalar @ClusterArray; $i++) {
        print OUT "$ClusterArray[$i][0]\t$ClusterArray[$i][1]\t$ClusterArray[$i][2]\t$ClusterArray[$i][3]\t$ClusterArray[$i][4]\n";
        # start--stop--cutoff--clusterAvg--clusterID
    }
    close (OUT);

    ########################################

    my $OutFile2 = "RD.$this->{Epsilon}.$this->{MinPts}.$this->{pairwise_file_name_only}.clusters";

    open (OUT, ">$OutFile2");

    print OUT "Cluster\tGene/Drug\tMutation/Gene\tDegree_Connectivity\tCloseness_Centrality\tGeodesic_From_Centroid\tRecurrence\tEpsilon_prime\tAvg_density\tCovering_clusters\tChromosome\tStart\tStop\tReference\tAlternate\tTranscript\tAlternative_Transcripts\n";

    for (my $i = 0; $i < scalar @ClusterArray; $i++) {
        for (my $j = $ClusterArray[$i][0]; $j <= $ClusterArray[$i][1]; $j++) {
            my @characters = split(":",$InitialSet[$j][0]);
            my $Gene = $characters[0];
            my $Mutation = $characters[1];
            my $genomicAnnotation = $InitialSet[$j][3];
            print OUT "$ClusterArray[$i][4]\t$Gene\t$Mutation\t0\t0\t1\t0\t$ClusterArray[$i][2]\t$ClusterArray[$i][3]\t$ClusterArray[$i][5]\t$genomicAnnotation\n";
        }
    }

    close (OUT);

    ########################################
    # Shiny File

    my $OutFile3 = new FileHandle;
    my $outFilename3 = "RD.$this->{Epsilon}.$this->{MinPts}.$this->{pairwise_file_name_only}.clusters.shiny.R";
    $OutFile3->open( $outFilename3 , "w"  );

    $OutFile3->print( "library(shiny)\n" );
    $OutFile3->print( "ui <- basicPage(
      plotOutput(\"plot1\",click = \"plot_click\", hover = \"plot_hover\", brush = \"plot_brush\" ,  height = 900),
      verbatimTextOutput(\"info\")
    )\n" );
    $OutFile3->print( "y = read.table(\"./RD.$this->{Epsilon}.$this->{MinPts}.$this->{pairwise_file_name_only}.out\")\n
    z = read.table(\"./RD.$this->{Epsilon}.$this->{MinPts}.$this->{pairwise_file_name_only}.clusters.plot\")\n" );
    $OutFile3->print( "RD<-y[[2]]\nID<-y[[1]]\nx0<-z[[1]]\ny0<-z[[3]]\nx1<-z[[2]]+1\ny1<-z[[3]]\nCluster<-z[[5]]\n" );
    $OutFile3->print( "server <- function(input, output) {\n  output\$plot1 <- renderPlot({\n   barplot(RD,names.arg=ID,main=\"Reachability Plot:Epsilon=8 MinPts=4\",col=\"Red\", cex.names=0.6,border=NA, space=0, las=2, ylab=\"Reachabilty Distance (A)\")\n    segments (x0,y0,x1,y1)\n
        text(x1+2,y0,Cluster, cex=1)\n
      })\n" );

    $OutFile3->print( "  output\$info <- renderText({\n
        RDID_str <- function(e) {\n
          if (is.null(e\$x)) return(\"\")\n
          name <- ID[round(e\$x+1)]\n
          RDval <- RD[round(e\$x+1)]\n
          paste0(\"ID=\", name, \"  RD=\", RDval)\n
        }\n
        RDID_range_str <- function(e) {\n
          if(is.null(e)) return(\"NULL\\n\")\n
          selrange <- ID[c(round(e\$xmin+1):round(e\$xmax+1))]\n
          paste0(selrange)\n
        }\n
        paste0(\n
          \"   hover:\", RDID_str(input\$plot_hover),\n
          \"\\n\", RDID_range_str(input\$plot_brush)\n
        )\n
      })\n
    }\n" );

    $OutFile3->print( "shinyApp(ui, server)" );

    $OutFile3->close();

    ############################  Plotting  #################################

    #system ("Rscript ClustersLines.R $inputFile ./Results/$ARGV[0].SuperClustersID.plot SuperClustersID.$ARGV[0].pdf $Epsilon $MinPts");
    system ("Rscript $this->{'R_path'}/HorizClustersLines.R RD.$this->{Epsilon}.$this->{MinPts}.$this->{pairwise_file_name_only}.out RD.$this->{Epsilon}.$this->{MinPts}.$this->{pairwise_file_name_only}.clusters.plot RD.$this->{Epsilon}.$this->{MinPts}.$this->{pairwise_file_name_only}.Horiz.pdf $this->{Epsilon} $this->{MinPts}");
    system ("Rscript $this->{'R_path'}/ColorScore.R RD.$this->{Epsilon}.$this->{MinPts}.$this->{pairwise_file_name_only}.out RD.$this->{Epsilon}.$this->{MinPts}.$this->{pairwise_file_name_only}.clusters.plot RD.$this->{Epsilon}.$this->{MinPts}.$this->{pairwise_file_name_only}.Horiz.weights.pdf $this->{Epsilon} $this->{MinPts}");
    #system ("Rscript EasyClustersLines.R $inputFile ./Results/$ARGV[0].SuperClustersID.plot SuperClustersID.$ARGV[0].Easy.pdf $Epsilon $MinPts");
    #print "Done.\n";
}


####################################################################################
############################ Clustering Probabilities ##############################
####################################################################################

sub getClusterProbabilities{
    my $this = shift;
	my $structure = shift;
    my $distance_matrix = shift;
    my $mutations = shift;

    my $Epsilon = $this->{'Epsilon'};
    my $MinPts = $this->{'MinPts'};
    my $NumRuns = $this->{'number_of_runs'};
    my $ProbCutOff = ($this->{'probability_cut_off'}/100)*$NumRuns;
    my @InitialCuts;
    my @InitialRD = @{$this->{CurrentRDarray}};
    delete $this->{CurrentRDarray}; # clean-up because this key will be used in other OPTICS runs below
    my $clustersFile = "./RD.$this->{Epsilon}.$this->{MinPts}.$this->{pairwise_file_name_only}.clusters";

    open(IN, "<$clustersFile") || die "Can't open $clustersFile: $!";
    while (my $line = <IN>) {
        if ( not $line =~ /Cluster/ ) {
            chomp $line;
            my @tabs2 = split(/\t/,$line);
            my @annotations = splice @tabs2, 10, 7;
            my $genomicAnnotation = join ( "\t", @annotations);
            push @InitialCuts, [$tabs2[0], $tabs2[1], $tabs2[2], $tabs2[7], $tabs2[8], $tabs2[9], $genomicAnnotation];
            # Cluster   Gene/Drug   Mutation/Gene   Epsilon_prime   Avg_density  tCovering_clusters   annotations(several columns)
        }
    }

    ###############################################################################
    # print "InitialCuts\n";
    # print Dumper \@InitialCuts;

    for (my $i = 0; $i < scalar @InitialCuts; $i++) {
        $InitialCuts[$i][0] =~ /(\d+)\.(\d+)\.(\d+)/g;
        if ($2 != 0) { # sub clusters
            $this->{"InitialCuts"}->{$1}->{$2} = $InitialCuts[$i][3];
        }
        else { # super clusters, put epsilon'=10
            $this->{"InitialCuts"}->{$1}->{$2} = 10;
        }
        $this->{"Variants"}->{$InitialCuts[$i][1].":".$InitialCuts[$i][2]}->{"run0"}->{$1}->{$2}->{$3}->{"AverageDensity"} = $InitialCuts[$i][4];
        $this->{"Variants"}->{$InitialCuts[$i][1].":".$InitialCuts[$i][2]}->{"run0"}->{$1}->{$2}->{$3}->{"CoveringClusters"} = $InitialCuts[$i][5];
        $this->{"Variants"}->{$InitialCuts[$i][1].":".$InitialCuts[$i][2]}->{"run0"}->{$1}->{$2}->{$3}->{"Annotations"} = $InitialCuts[$i][6];
        $this->{"Memberships"}->{$1}->{$2}->{$3}->{$InitialCuts[$i][1].":".$InitialCuts[$i][2]} = 1;
    }

    # print Dumper $this;

    ################################################################################

    #print "\nOrderedNodes array=\n";

    #print Dumper $this->{CurrentRDarray};

    for (my $run = 1; $run < $this->{'number_of_runs'}; $run++) {

        $this->MainOPTICS( $structure , $distance_matrix, $mutations ); # Generate a random ordred RD set (random OPTICS)

        $this->GetSuperClusters($this->{CurrentRDarray}, $MinPts); # Identify super clusters
        # print "CurrentRDarray\n";
        # print Dumper $this->{CurrentRDarray};
        # print "CurrentSuperClusters\n";
        # print Dumper $this->{"CurrentSuperClusters"};
        $this->GetSuperClusterMapping(); # Map new super clusters to the ones in run0

        $this->GetSubClusters($MinPts, $run); # Perform clustering at initial epsilon cuts
        $this->GetSubClusterMapping(); # Map new subclusters to the ones in run0

        $this->RecordMemberships(); # if exists increase the number, else add to the list

        #writing to file and plotting
        my $OrderedFile1 = "./$run.RD.out";
        open (OUT, ">$OrderedFile1");
        for (my $i = 0; $i < scalar @{$this->{CurrentRDarray}}; $i++) {
            my $ele1 = ${$this->{CurrentRDarray}}[$i][0];
            my $ele2 = ${$this->{CurrentRDarray}}[$i][1];
            print OUT "$ele1\t$ele2\n";
        }
        close (OUT);

        # print "SubClusters=\n";
        # print Dumper $this->{SubClusters}->{0};
        # print "SubCluster Matching=\n";
        # print Dumper $this->{SubClusterMatching}->{0};
        # print "SubCluster Mapping=\n";
        # print Dumper $this->{SubClusterMap}->{0};
        # print "\n-----------------------------------------------------------------\n";

        # Clean up hashes specific to the current run
        delete $this->{CurrentRDarray};
        delete $this->{CurrentSuperClusters};
        delete $this->{SuperClusterMatching};
        delete $this->{SuperClusterMap};
        delete $this->{SubClusters};
        delete $this->{SubClusterMap};
        delete $this->{SubClusterMatching};

        system ("Rscript $this->{'R_path'}/BruteForceClustersLines.R ./$run.RD.out ./$run.clusters ./$run.pdf $Epsilon $MinPts");


    } # end of run

    # print "memberships\n";
    # print Dumper $this->{"Memberships"};

    # Data to generate the Cluster Membership Probability Plot
    my $FinalDataFile1 = "./ProbabilityData.$this->{'pairwise_file_name_only'}";
    open (OUT, ">$FinalDataFile1");
        print OUT "Variant\tProbability\tClusterID\tSuperClusterID\n";

        for (my $i = 0; $i < scalar @InitialRD; $i++) {
            my $variant1 = $InitialRD[$i][0];

            foreach my $SCID (keys %{$this->{Memberships}}) {
                foreach my $levelID (keys %{$this->{Memberships}->{$SCID}}) {  
                    foreach my $SubID (keys %{$this->{Memberships}->{$SCID}->{$levelID}}) {
                        if (exists $this->{Memberships}->{$SCID}->{$levelID}->{$SubID}->{$variant1}) {
                            my $Occurance1 = $this->{Memberships}->{$SCID}->{$levelID}->{$SubID}->{$variant1};
                            $Occurance1 = $Occurance1/$NumRuns;
                            print OUT "$variant1\t$Occurance1\t$SCID.$levelID.$SubID\t$SCID\n";
                        }
                    }
                }
            }
        }
    close (OUT);

    # Clusters output file (will be used in the visual)
    my $FinalDataFile2 = "./$Epsilon.$MinPts.$this->{'pairwise_file_name_only'}.Prob.$this->{'probability_cut_off'}.clusters";
    open (OUT, ">$FinalDataFile2");
        print OUT "Cluster\tGene/Drug\tMutation/Gene\tDegree_Connectivity\tCloseness_Centrality\tGeodesic_From_Centroid\tRecurrence/Weight\tEpsilon_prime\tAvg_density\tCovering_clusters\tChromosome\tStart\tStop\tReference\tAlternate\tTranscript\tAlternative_Transcripts\n";

        # print "this->variants\n";
        # print Dumper $this->{Variants};
        for (my $SCID = 0; $SCID < scalar keys %{$this->{Memberships}}; $SCID++) {
            for (my $levelID = 0; $levelID < scalar keys %{$this->{Memberships}->{$SCID}}; $levelID++) {
                for (my $SubID = 0; $SubID < scalar keys %{$this->{Memberships}->{$SCID}->{$levelID}}; $SubID++) {
                    foreach my $variant (keys %{$this->{Memberships}->{$SCID}->{$levelID}->{$SubID}}) {
                        if ($this->{Memberships}->{$SCID}->{$levelID}->{$SubID}->{$variant} >= $ProbCutOff) {
                            my $CurrentEpsilon = $this->{InitialCuts}->{$SCID}->{$levelID};
                            my $CurrentAvgDensity = $this->{Variants}->{$variant}->{run0}->{$SCID}->{$levelID}->{$SubID}->{AverageDensity};
                            my $CoveringClusters =  $this->{Variants}->{$variant}->{run0}->{$SCID}->{$levelID}->{$SubID}->{CoveringClusters};
                            my $genomicAnnotation = $this->{Variants}->{$variant}->{run0}->{$SCID}->{$levelID}->{$SubID}->{Annotations};
                            $variant =~ /(\w+)\:(\D\.\D+\d+\D+)/g;
                            # print OUT "$SCID.$levelID.$SubID\t$1\t$2\t0\t0\t0\t0\t$CurrentEpsilon\t$CurrentAvgDensity\t$CoveringClusters\t$genomicAnnotation\n";

                            if ( scalar keys %{$this->{Memberships}->{$SCID}->{$levelID}->{$SubID}} == 1 and $CurrentAvgDensity == 0.1 ) { # artificially reset epsilon-prime to 0.1 for singleton clusters
                                $CurrentEpsilon = 0.1;
                            }

                            # print clusters with other mutations (not only representative vertices, but all)
                            
                            my $gene = $1;
                            my ( $chromosome, $start, $stop, $ref, $alt, $transcript ) = ( split /\t/, $genomicAnnotation )[0,1,2,3,4,5];
                            my $mKey1 = join( ":", $gene, $chromosome, $start, $stop );
                            my $refAlt1 = join( ":", $ref, $alt );
                            my $pKey1 = join( ":", $transcript, $2 );
                            my $recurrence1 = $this->{"siteVertexMap"}->{$mKey1}->{$mKey1}->{$refAlt1}->{$pKey1};
                            print OUT "$SCID.$levelID.$SubID\t$1\t$2\t0\t0\t0\t$recurrence1\t$CurrentEpsilon\t$CurrentAvgDensity\t$CoveringClusters\t$genomicAnnotation\n"; # representative vertex

                            foreach my $mKey2 ( keys %{$this->{"siteVertexMap"}->{$mKey1}} ) {
                                my ( $chromosome2, $start2, $stop2 ) = (split /:/, $mKey2)[1,2,3];
                                
                                foreach my $refAlt2 ( keys %{$this->{"siteVertexMap"}->{$mKey1}->{$mKey2}} ) {
                                    my ( $ref2, $alt2 ) = split /:/, $refAlt2;
                                    
                                    foreach my $pKey2 ( keys %{$this->{"siteVertexMap"}->{$mKey1}->{$mKey2}->{$refAlt2}} ) {
                                        
                                        my ( $transcript2, $aa2 ) =  split /:/, $pKey2;
                                        my $recurrence2 = $this->{"siteVertexMap"}->{$mKey1}->{$mKey2}->{$refAlt2}->{$pKey2};
                                        my $genomicAnnotation2 = join( "\t", $chromosome2, $start2, $stop2, $ref2, $alt2, $transcript2, "1" ); # included "1" in the alternateTrans field
                                        
                                        if ( $mKey2 eq $mKey1 and $refAlt1 eq $refAlt2 ) { # if two mutation keys and ref:alt's are equal, output all the other pKeys
                                            if ( $pKey2 ne $pKey1 ) {
                                                $recurrence1 = $this->{"siteVertexMap"}->{$mKey1}->{$mKey2}->{$refAlt2}->{$pKey2};
                                                print OUT "$SCID.$levelID.$SubID\t$1\t$2\t0\t0\t0\t$recurrence1\t$CurrentEpsilon\t$CurrentAvgDensity\t$CoveringClusters\t$genomicAnnotation2\n";
                                            }
                                        }
                                        else { # other mutations which are at the same site (alternate ones)
                                            print OUT "$SCID.$levelID.$SubID\t$gene\t$aa2\t0\t0\t0\t$recurrence2\t$CurrentEpsilon\t$CurrentAvgDensity\t$CoveringClusters\t$genomicAnnotation2\n";
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
    close (OUT);

    # # print "SC matching=\n";
    # # print Dumper $this->{SuperClusterMatching};
    # print "SC map=\n";
    # print Dumper $this->{SuperClusterMap};
    # print "SubClusters=\n";
    # print Dumper $this->{SubClusters};
    # # print "Initial cut=\n";
    # # print Dumper $this->{InitialCuts};
    # print "SubCluster Matching=\n";
    # print Dumper $this->{SubClusterMatching};
    # print "SubCluster Mapping=\n";
    # print Dumper $this->{SubClusterMap};

    # print "Memberships=\n";
    # print Dumper $this->{Memberships};

    system ("Rscript $this->{'R_path'}/MembershipProbability.R $FinalDataFile1 $NumRuns $this->{'pairwise_file_name_only'}");

    # clean-up : This is important because in structure-wise clustering these hash names will be re-used.
    delete $this->{Memberships};
    delete $this->{InitialCuts};
    delete $this->{Variants};
    #print "Done.\n";
}

#####
##   Probability sub-methods
#####

sub GetSuperClusters { # Identify super clusters: (singleton super clusters are possible to form)
    my ($this, $arrayRef, $MinPts) = @_;

    my $scs = 0; # super cluster start
    my $superClusterCount = ${$arrayRef}[$scs][4]; # count the weight (for SITE=1, RECURRENCE=frequency, WEIGHT=weight)
    for (my $i = 1; $i < scalar @{$arrayRef}; $i++) {
        #print "i=$i\n";
        if ( ${$arrayRef}[$i][1] == 10 ) {
            #print "RD($i)=10\n";
            if ( $superClusterCount >= $MinPts ) {
                $this->{"CurrentSuperClusters"}->{$scs} = $i-1;
            }
            $scs = $i; # reset the super cluster starting point
            $superClusterCount = ${$arrayRef}[$scs][4]; # reset the count  
        }
        else {
            $superClusterCount += ${$arrayRef}[$i][4]; # add the weight (for SITE=1, RECURRENCE=frequency, WEIGHT=weight)
        }   
    }
    return 1;
}

sub GetSuperClusterMapping {
    my $this = shift @_;
    foreach my $CurrentSCstart (keys %{$this->{CurrentSuperClusters}}) { # finding matching super clusters
        #if ($this->{CurrentSuperClusters}->{$CurrentSCstart} - $CurrentSCstart >= $MinPts) {
            for (my $i = $CurrentSCstart; $i <= $this->{CurrentSuperClusters}->{$CurrentSCstart}; $i++) {
                # print "current variant=${$this->{CurrentRDarray}}[$i][0]\n";
                # print "In Old=\n";
                # print Dumper keys %{$this->{Variants}->{${$this->{CurrentRDarray}}[$i][0]}->{run0}};
                my @SCArray = keys %{$this->{Variants}->{${$this->{CurrentRDarray}}[$i][0]}->{run0}}; #SCArray will only contain one element
                if (defined $SCArray[0]) {
                    my $TotHits;
                    if (exists $this->{SuperClusterMatching}->{$CurrentSCstart}->{$SCArray[0]}) {
                        $TotHits = $this->{SuperClusterMatching}->{$CurrentSCstart}->{$SCArray[0]};
                    }
                    else {
                        $TotHits = 0;
                    }
                    $TotHits++;
                    $this->{"SuperClusterMatching"}->{$CurrentSCstart}->{$SCArray[0]} = $TotHits;   
                }
            }           
        #}
        
        my @SCmatchArray = sort { $this->{SuperClusterMatching}->{$CurrentSCstart}->{$a} <=> $this->{SuperClusterMatching}->{$CurrentSCstart}->{$b} } keys %{$this->{SuperClusterMatching}->{$CurrentSCstart}};
        my $SCmatch = pop @SCmatchArray; # this is necessary because sometimes super cluster membership changes at low densities
            #print "SC map= $CurrentSCstart\t$SCmatch\n";
        if (defined $SCmatch && $SCmatch ne '') {
            $this->{"SuperClusterMap"}->{$SCmatch}->{$CurrentSCstart} = $this->{"CurrentSuperClusters"}->{$CurrentSCstart};
        }
    }

}

sub GetSubClusters {
    my ($this, $MinPts, $run) = @_;
    my $runClusterFile = "./$run.clusters";
    open (OUT, ">$runClusterFile");
    foreach my $SCID (keys %{$this->{InitialCuts}}) {
        foreach my $levelID (keys %{$this->{InitialCuts}->{$SCID}}) {
            my $SubID = 1;
            my $EpsilonCut = $this->{InitialCuts}->{$SCID}->{$levelID};
            my @nStartArray = keys %{$this->{SuperClusterMap}->{$SCID}};
            my $nStart = shift @nStartArray;
            my $n = $nStart;
            my $counter = 0;
            my $s = $nStart;
            while ($n <= $this->{SuperClusterMap}->{$SCID}->{$nStart}) {
                if (${$this->{CurrentRDarray}}[$n][1] < $EpsilonCut) {
                    $counter += ${$this->{CurrentRDarray}}[$n][4]; # add the weight (for SITE=1, RECURRENCE=frequency, WEIGHT=weight)
                    # $counter++;
                    if ($counter >= $MinPts) {
                        $this->{SubClusters}->{$SCID}->{$levelID}->{$s} = $n;
                    }
                    $n++;
                }
                else {
                    if ( $levelID == 1 and ${$this->{CurrentRDarray}}[$n][4] >= $MinPts ) { # allowing singleton sub clusters to form at the lowest level
                        $this->{SubClusters}->{$SCID}->{$levelID}->{$n} = $n; 
                    }
                    # regular sub cluster stuff
                    $s = $n;
                    $counter = ${$this->{CurrentRDarray}}[$s][4]; # reset the count
                    $n++;
                }
            }
            foreach my $start1 (keys %{$this->{SubClusters}->{$SCID}->{$levelID}}) {
                my $stop1 = $this->{SubClusters}->{$SCID}->{$levelID}->{$start1};
                print OUT "$start1\t$stop1\t$EpsilonCut\n";
            }
        }
    }
    close (OUT);
}

sub GetSubClusterMapping {
    my $this = shift @_;
    
    foreach my $SCID (keys %{$this->{SubClusters}}) {
        foreach my $levelID (keys %{$this->{SubClusters}->{$SCID}}) {
            my @DummyArray = keys %{$this->{SubClusters}->{$SCID}->{$levelID}}; # contains the start points of sub clusters at levelID
            foreach my $nStart (@DummyArray) {
                #my $nStart = shift @DummyArray;
                my $nStop = $this->{SubClusters}->{$SCID}->{$levelID}->{$nStart};
                for (my $i = $nStart; $i <= $nStop; $i++) {
                    my @DummyArray2 = keys %{$this->{Variants}->{${$this->{CurrentRDarray}}[$i][0]}->{run0}->{$SCID}->{$levelID}}; # contains only one element (subID)
                    if (defined $DummyArray2[0]) {
                        my $TotHits;
                        if (exists $this->{SubClusterMatching}->{$SCID}->{$levelID}->{$nStart}->{$DummyArray2[0]}) {
                            $TotHits = $this->{SubClusterMatching}->{$SCID}->{$levelID}->{$nStart}->{$DummyArray2[0]};
                        }
                        else {
                            $TotHits = 0;
                        }
                        $TotHits++;
                        $this->{SubClusterMatching}->{$SCID}->{$levelID}->{$nStart}->{$DummyArray2[0]} = $TotHits;
                    }
                }

                my @SubMatchArray = sort { $this->{SubClusterMatching}->{$SCID}->{$levelID}->{$nStart}->{$a} <=> $this->{SubClusterMatching}->{$SCID}->{$levelID}->{$nStart}->{$b} } keys %{$this->{SubClusterMatching}->{$SCID}->{$levelID}->{$nStart}};
                my $SubMatch = pop @SubMatchArray; # this is necessary because sometimes sub cluster membership changes at low densities
                    #print "SC map= $CurrentSCstart\t$SCmatch\n";
                if (defined $SubMatch && $SubMatch ne '') {
                    $this->{"SubClusterMap"}->{$SCID}->{$levelID}->{$SubMatch}->{$nStart} = $this->{"SubClusters"}->{$SCID}->{$levelID}->{$nStart};
                }
            }   
        }
    }
}

sub RecordMemberships {
    my $this = shift @_;

    foreach my $SCID (keys %{$this->{SubClusterMap}}) {
        foreach my $levelID (keys %{$this->{SubClusterMap}->{$SCID}}) {    
            foreach my $SubID (keys %{$this->{SubClusterMap}->{$SCID}->{$levelID}}) {
                my @DummyArray = keys %{$this->{SubClusterMap}->{$SCID}->{$levelID}->{$SubID}};
                my $nStart = shift @DummyArray;
                my $nStop = $this->{SubClusterMap}->{$SCID}->{$levelID}->{$SubID}->{$nStart};
                #print "$SCID\t$levelID\t$SubID\$nStart\n";
                for (my $i = $nStart; $i <= $nStop; $i++) {
                    my $Occurance;
                    if (exists $this->{Memberships}->{$SCID}->{$levelID}->{$SubID}->{${$this->{CurrentRDarray}}[$i][0]}) {
                        $Occurance = $this->{Memberships}->{$SCID}->{$levelID}->{$SubID}->{${$this->{CurrentRDarray}}[$i][0]};
                        $Occurance++;
                        $this->{Memberships}->{$SCID}->{$levelID}->{$SubID}->{${$this->{CurrentRDarray}}[$i][0]} = $Occurance;
                    }
                    else {
                        $Occurance = 1;
                        $this->{Memberships}->{$SCID}->{$levelID}->{$SubID}->{${$this->{CurrentRDarray}}[$i][0]} = $Occurance;
                    }
                }
            }
        }
    }
}

sub setRWD {
    my $this = shift;

    my $pathDensity = $INC{"TGI/Mutpro/Main/Density.pm"};
    my @pathArray = split("/",$pathDensity);

    my $pathWD = $pathArray[0];
    for (my $i = 1; $i < scalar @pathArray -1; $i++) {
        $pathWD = $pathWD."/".$pathArray[$i];
    }
    $this->{'R_path'} = $pathWD;
    #print "$this->{'R_path'}\n";


}

sub printHash {
    my ( $hashRef, $fileName ) = @_;
    my $outfile = join( "", $fileName, ".json" ) ;
    my $OUT = FileHandle->new( "$outfile" , "w" );

    my $json = encode_json $hashRef;
    $OUT->print( $json );
    return;
}

sub readHash {
    my $fileName = shift @_;
    my $hashRef = {};
    print "read JSON hash: $fileName\n";
    my $fh = new FileHandle;
    die "Could not open $fileName\n" unless( $fh->open( $fileName , "r" ) );

    my $json = <$fh>;
    $fh->close();

    $hashRef = decode_json( $json );

    return $hashRef;
}

sub help_text{
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

1;
