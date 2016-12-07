use strict;
use warnings;

use LWP::Simple;
use FileHandle;
use Data::Dumper qw(Dumper);
use Scalar::Util qw(looks_like_number);
use List::Util qw[min max shuffle];

######################## Command Line Arguments ################################

# quit unless we have the correct number of command-line args
my $num_args = $#ARGV + 1;
if ($num_args != 6) {
    print "\nUsage: HardClusters.pl RD.*.clusters_in_./Results Epsilon MinPts Pairwisefile_in_./Test Number_of_runs_needed Probability_cut_off_percentage\n\n";
    exit;
}

my $inputFile1 = "./Results/$ARGV[0]";
my $Epsilon = $ARGV[1];
my $MinPts = $ARGV[2];
my $inputRD = "./Results/RD.$Epsilon.$MinPts.$ARGV[3]";
my $NumRuns = $ARGV[4];
my $ProbCutOff = ($ARGV[5]/100)*$NumRuns;

##########################  Reading from RD.out  #############################

my $this = {};
my @InitialCuts;
my @InitialRD;

open(IN, "<$inputFile1") || die "Can't open $inputFile1: $!";
while (my $line = <IN>) {
	if ( not $line =~ /Cluster/ ) {
		chomp $line;
		my @tabs2 = split(/\t/,$line);
		push @InitialCuts, [$tabs2[0],$tabs2[1],$tabs2[2],$tabs2[7],$tabs2[8],$tabs2[9]];
		# Cluster	Gene/Drug	Mutation/Gene	Epsilon_prime	Avg_density	 Covering_clusters
	}
}

open(IN, "<$inputRD") || die "Can't open $inputRD: $!";
while (my $line = <IN>) {
	
	chomp $line;
	my @tabs3 = split(/\t/,$line);
	push @InitialRD, [$tabs3[0],$tabs3[1]];
	# variant  RD
	
}

###############################################################################

#print Dumper \@InitialCuts;

for (my $i = 0; $i < scalar @InitialCuts; $i++) {
	$InitialCuts[$i][0] =~ /(\d+)\.(\d+)\.(\d+)/g;
	if ($2 != 0) {
		$this->{"InitialCuts"}->{$1}->{$2} = $InitialCuts[$i][3];
	}
	else {
		$this->{"InitialCuts"}->{$1}->{$2} = 10;
	}
	$this->{"Variants"}->{$InitialCuts[$i][1].":".$InitialCuts[$i][2]}->{"run0"}->{$1}->{$2}->{$3}->{"AverageDensity"} = $InitialCuts[$i][4];
	$this->{"Variants"}->{$InitialCuts[$i][1].":".$InitialCuts[$i][2]}->{"run0"}->{$1}->{$2}->{$3}->{"CoveringClusters"} = $InitialCuts[$i][5];
	$this->{"Memberships"}->{$1}->{$2}->{$3}->{$InitialCuts[$i][1].":".$InitialCuts[$i][2]} = 1;
}

#print Dumper $this;

################################################################################

#print "\nOrderedNodes array=\n";

#print Dumper $this->{CurrentRDarray};



for (my $run = 1; $run < $NumRuns; $run++) {

	MainOPTICS($this, "$ARGV[3]"); # Generate a random ordred RD set (random OPTICS)

	GetSuperClusters($this, $this->{CurrentRDarray}); # Identify super clusters
	GetSuperClusterMapping($this); # Map new super clusters to the ones in run0

	GetSubClusters($this, $MinPts, $run); # Perform clustering at initial epsilon cuts
	GetSubClesterMapping($this); # Map new subclusters to the ones in run0

	RecordMemberships($this); # if exists increase the number, else add to the list

	#writing to file and plotting
	my $OrderedFile1 = "./Results/runs/$run.RD.out";
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

	system ("Rscript BruteForceClustersLines.R ./Results/runs/$run.RD.out ./Results/runs/$run.clusters ./Results/runs/$run.pdf $Epsilon $MinPts");


} # end of run

# Data to generate the Cluster Membership Probability Plot
my $FinalDataFile1 = "./Results/ProbabilityData.$ARGV[3]";
open (OUT, ">$FinalDataFile1");
	print OUT "Variant\tProbability\tClusterID\n";

	for (my $i = 0; $i < scalar @InitialRD; $i++) {
		my $variant1 = $InitialRD[$i][0];

		foreach my $SCID (keys $this->{Memberships}) {
			foreach my $levelID (keys $this->{Memberships}->{$SCID}) {	
				foreach my $SubID (keys $this->{Memberships}->{$SCID}->{$levelID}) {
					if (exists $this->{Memberships}->{$SCID}->{$levelID}->{$SubID}->{$variant1}) {
						my $Occurance1 = $this->{Memberships}->{$SCID}->{$levelID}->{$SubID}->{$variant1};
						$Occurance1 = $Occurance1/$NumRuns;
						print OUT "$variant1\t$Occurance1\t$SCID.$levelID.$SubID\n";
					}
				}
			}
		}
	}
close (OUT);

# Clusters output file (will be used in the visual)
my $FinalDataFile2 = "./Results/$Epsilon.$MinPts.$ARGV[3].Prob.$ARGV[5].clusters";
open (OUT, ">$FinalDataFile2");
	print OUT "Cluster\tGene/Drug\tMutation/Gene\tDegree_Connectivity\tCloseness_Centrality\tGeodesic_From_Centroid\tRecurrence\tEpsilon_prime\tAvg_density\tCovering_clusters\n";

	for (my $SCID = 0; $SCID < scalar keys $this->{Memberships}; $SCID++) {
		for (my $levelID = 0; $levelID < scalar keys $this->{Memberships}->{$SCID}; $levelID++) {
			for (my $SubID = 0; $SubID < scalar keys $this->{Memberships}->{$SCID}->{$levelID}; $SubID++) {
				foreach my $variant (keys %{$this->{Memberships}->{$SCID}->{$levelID}->{$SubID}}) {
					if ($this->{Memberships}->{$SCID}->{$levelID}->{$SubID}->{$variant} >= $ProbCutOff) {
						my $CurrentEpsilon = $this->{InitialCuts}->{$SCID}->{$levelID};
						my $CurrentAvgDensity = $this->{Variants}->{$variant}->{run0}->{$SCID}->{$levelID}->{$SubID}->{AverageDensity};
						my $CoveringClusters =  $this->{Variants}->{$variant}->{run0}->{$SCID}->{$levelID}->{$SubID}->{CoveringClusters};
						$variant =~ /(\w+)\:(\D\.\D+\d+\D)/g;
						print OUT "$SCID.$levelID.$SubID\t$1\t$2\t0\t0\t0\t0\t$CurrentEpsilon\t$CurrentAvgDensity\t$CoveringClusters\n";
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

print "Memberships=\n";
print Dumper $this->{Memberships};

system ("Rscript MembershipProbability.R $FinalDataFile1 $NumRuns $ARGV[3]");

print "Done.\n";

####################################################################
##########################  Functions  #############################
####################################################################

sub RecordMemberships {
	my $this = shift @_;

	foreach my $SCID (keys $this->{SubClusterMap}) {
		foreach my $levelID (keys $this->{SubClusterMap}->{$SCID}) {	
			foreach my $SubID (keys $this->{SubClusterMap}->{$SCID}->{$levelID}) {
				my @DummyArray = keys $this->{SubClusterMap}->{$SCID}->{$levelID}->{$SubID};
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

sub GetSubClesterMapping {
	my $this = shift @_;
	
	foreach my $SCID (keys $this->{SubClusters}) {
		foreach my $levelID (keys $this->{SubClusters}->{$SCID}) {
			my @DummyArray = keys %{$this->{SubClusters}->{$SCID}->{$levelID}}; # contains the start points of sub clusters at levelID
			foreach my $nStart (@DummyArray) {
				#my $nStart = shift @DummyArray;
				my $nStop = $this->{SubClusters}->{$SCID}->{$levelID}->{$nStart};
				for (my $i = $nStart; $i <= $nStop; $i++) {
					my @DummyArray2 = keys %{$this->{Variants}->{${$this->{CurrentRDarray}}[$i][0]}->{run0}->{$SCID}->{$levelID}}; # contains only one element (subID)
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

				my @SubMatchArray = sort { $this->{SubClusterMatching}->{$SCID}->{$levelID}->{$nStart}->{$a} <=> $this->{SubClusterMatching}->{$SCID}->{$levelID}->{$nStart}->{$b} } keys %{$this->{SubClusterMatching}->{$SCID}->{$levelID}->{$nStart}};
				my $SubMatch = pop @SubMatchArray; # this is necessary because sometimes sub cluster membership changes at low densities
					#print "SC map= $CurrentSCstart\t$SCmatch\n";
				if ($SubMatch ne '') {
					$this->{"SubClusterMap"}->{$SCID}->{$levelID}->{$SubMatch}->{$nStart} = $this->{"SubClusters"}->{$SCID}->{$levelID}->{$nStart};
				}
			}	
		}
	}
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
		#}
		
		my @SCmatchArray = sort { $this->{SuperClusterMatching}->{$CurrentSCstart}->{$a} <=> $this->{SuperClusterMatching}->{$CurrentSCstart}->{$b} } keys %{$this->{SuperClusterMatching}->{$CurrentSCstart}};
		my $SCmatch = pop @SCmatchArray; # this is necessary because sometimes super cluster membership changes at low densities
			#print "SC map= $CurrentSCstart\t$SCmatch\n";
		if ($SCmatch ne '') {
			$this->{"SuperClusterMap"}->{$SCmatch}->{$CurrentSCstart} = $this->{"CurrentSuperClusters"}->{$CurrentSCstart};
		}
	}

}

sub GetSubClusters {
	my ($this, $MinPts, $run) = @_;
	my $OrderedFile2 = "./Results/runs/$run.clusters";
	open (OUT, ">$OrderedFile2");
	foreach my $SCID (keys $this->{InitialCuts}) {
		foreach my $levelID (keys $this->{InitialCuts}->{$SCID}) {
			my $SubID = 1;
			my $EpsilonCut = $this->{InitialCuts}->{$SCID}->{$levelID};
			my @nStartArray = keys $this->{SuperClusterMap}->{$SCID};
			my $nStart = shift @nStartArray;
			my $n = $nStart;
			my $counter = 0;
			my $s = $nStart;
			while ($n <= $this->{SuperClusterMap}->{$SCID}->{$nStart}) {
				if (${$this->{CurrentRDarray}}[$n][1] < $EpsilonCut) {
					$counter++;
					if ($counter >= $MinPts) {
						$this->{SubClusters}->{$SCID}->{$levelID}->{$s} = $n;
					}
					$n++;
				}
				else {
					$s = $n;
					$counter = 0;
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

sub GetSuperClusters { # Identify super clusters:
	my ($this, $arrayRef) = @_;

	my $scs = 0; # super cluster start
	for (my $i = 1; $i < scalar @{$arrayRef}; $i++) {
		#print "i=$i\n";
		if ( ${$arrayRef}[$i][1] == 10 ) {
			#print "RD($i)=10\n";
			$scs = $i;	
		}
		else {
			$this->{"CurrentSuperClusters"}->{$scs} = $i;
		}	
	}
	return 1;
}

sub MainOPTICS {
	my ($this, $Pairwisefile)=@_;

	my %SetOfNodes;

	my $file = "./Test/$Pairwisefile";
	open(IN, "<$file") || die "Can't open $file: $!";

	while (my $line = <IN>) {
		chomp $line;
		my @tabs = split(/\t/,$line);
		my @char19 = split("",$tabs[19]);
		my $dis = $char19[0].$char19[1].$char19[2].$char19[3].$char19[4];
		my $key1 = CombineWords($tabs[0],$tabs[4]);
		my $value1 = CombineWords($tabs[9],$tabs[13]);

		$SetOfNodes{$key1}{distances}{$value1} = $dis;
		$SetOfNodes{$value1}{distances}{$key1} = $dis;
	}
	###### For variants in the same residue and chain 

	foreach my $key ( keys %SetOfNodes ) {
		#print "key= $key\n";
		$key =~ /(\w+)\:\D\.(\D+\d+)\D/g;
		my $keyGene = $1;
		my $keyRes = $2;
		my @hits = grep(/$keyGene\:\D\.$keyRes\D/g, keys %SetOfNodes);
		#print Dumper \@hits;
		foreach my $hit (@hits) {
			if ( $hit ne $key ) {
				$SetOfNodes{$key}{distances}{$hit} = "0";
				$SetOfNodes{$hit}{distances}{$key} = "0";
			}
		}
	}

	foreach my $i (keys %SetOfNodes) {
		$SetOfNodes{$i}{processInfo} = "False";
	}

	#print Dumper \%SetOfNodes;
	print "Number of Objects = ";
	print scalar keys %SetOfNodes;
	print "\n";

	my @SetOfCores;
	my @SetOfEdges;
	foreach my $key ( keys %SetOfNodes ) {
		if ( scalar keys $SetOfNodes{$key}{distances} >= $MinPts ) {
			push @SetOfCores, $key;
		}
		else {
			push @SetOfEdges, $key;
		}
	}
	@SetOfCores = shuffle @SetOfCores;
	@SetOfEdges = shuffle @SetOfEdges;
	my @SetOfCoresThenEdges = ( @SetOfCores, @SetOfEdges );

	###########################################################

	my @OrderedNodes;

	################# Main OPTICS function ####################

	foreach my $p ( @SetOfCoresThenEdges ) {
		#print "first p=$p\n";
		if ($SetOfNodes{$p}{processInfo} =~ "False") {
			########## Expand Cluster Order ###########
			my %neighbors; # is a hash with keys neigbor indices whose values are mutual separations
			my %OrderSeeds; # is a hash to add seeds
			%neighbors = %{GetNeighbors($p,$Epsilon,\%SetOfNodes)};
			$SetOfNodes{$p}{processInfo} = "True"; # set as processed
			my $RD = undef;
			my $CD;
			$CD = GetCoreDistance(\%neighbors,$MinPts);
			# print "p=$p and ";
			# print "CD=$CD\n";
			push @OrderedNodes, [$p,$RD,$CD]; # write to the file 
			if (defined $CD) {
				OrderSeedsUpdate(\%neighbors,$p,$CD, \%OrderSeeds, \%SetOfNodes);
				# print "For p=$p, OrderSeeds= \n";
				# print Dumper \%OrderSeeds;
				my $PrevObj = $p; # used to get the current obj. (To check whether variants are at the same location)
				while (scalar keys %OrderSeeds != 0) {
					my @SeedKeys = sort { $OrderSeeds{$a} <=> $OrderSeeds{$b} } keys %OrderSeeds;
					my @SeedValues = @OrderSeeds{@SeedKeys};
					#my $CurrentObject =  $SeedKeys[0]; # CurrentObject is the object having the least RD in OrderSeeds
					my $CurrentObject = GetCurrentObject(\@SeedValues, \@SeedKeys, $PrevObj);
					$PrevObj = $CurrentObject;
					#print "\n\n current object= $CurrentObject\t neighbors=";
					%neighbors = %{GetNeighbors($CurrentObject,$Epsilon,\%SetOfNodes)};
					#print Dumper \%neighbors;
					#print Dumper $SetOfNodes{$CurrentObject}{distances};
					$SetOfNodes{$CurrentObject}{processInfo} = "True"; # set as processed
					$RD = $SeedValues[0];
					$CD = GetCoreDistance(\%neighbors,$MinPts);
					push @OrderedNodes, [$CurrentObject,$RD,$CD]; # write to the file 
					delete $OrderSeeds{$CurrentObject};
					if (defined $CD) {
						#print "\tCurrent object is a core.(CD=$CD)\n Updated Order seeds list\n\t";
						OrderSeedsUpdate(\%neighbors,$CurrentObject,$CD, \%OrderSeeds, \%SetOfNodes);
						#print Dumper \%OrderSeeds;
					}
				}
			}
			# print "p=$p,(undefined CD) OrderedNodes= \n";
			# print Dumper \@OrderedNodes;
		}
	}

	### Replacing undefined RD by 10
	for (my $i = 0; $i < scalar @OrderedNodes; $i++) {
		if (not defined $OrderedNodes[$i][1]) {
			$OrderedNodes[$i][1] = 10;
		}
	}
	$this->{"CurrentRDarray"} = \@OrderedNodes;
	return $this;
}

sub GetNeighbors {
	my ($Obj, $Epsilon, $Set_ref)=@_;
	my %neighborHash;
	foreach my $i (keys %{$Set_ref->{$Obj}->{distances}}) {
			$neighborHash{$i} = "$Set_ref->{$Obj}->{distances}->{$i}";
	}
	return \%neighborHash;
}

sub GetCoreDistance {
	my ($neighbors_ref, $MinPts)=@_;
	my @keys = sort { $neighbors_ref->{$a} <=> $neighbors_ref->{$b} } keys %{$neighbors_ref}; # sort keys according to distances
	my @vals = @{$neighbors_ref}{@keys};
	my $CoreDist;
	if (scalar keys %{$neighbors_ref} >= $MinPts){
			$CoreDist = $vals[$MinPts-1]; # MinPt^th-distance
		}
	else {
		$CoreDist = undef;
	}
	return $CoreDist;
}

sub OrderSeedsUpdate {
	my ($neighbors_ref, $CenterObject, $CD, $OrderSeeds_ref, $Set_ref) = @_;
	my $c_dist = $CD; 
	my %neighborsHash = % { $neighbors_ref };
	my %OrderSeedsHash = % { $OrderSeeds_ref};
	foreach my $q (keys %{$neighbors_ref}) {
		if (${$Set_ref}{$q}{processInfo} =~ "False") {
			my $new_r_dist = max ($c_dist,${$neighbors_ref}{$q});
			if (exists ${$OrderSeeds_ref}{$q}) {
				if ($new_r_dist < ${$OrderSeeds_ref}{$q}) {
					${$OrderSeeds_ref}{$q}="$new_r_dist";
				}
			}
			else {
					${$OrderSeeds_ref}{$q}="$new_r_dist";
				}
		}
	}
}

sub CombineWords {
	my ($word1,$word2)=@_;
	return $word1.":".$word2;
}

sub GetCurrentObject { # To check whether variants are at the same location
	my ($ValueSet, $KeySet, $PrevObj)=@_;
	my @SmallestKeys;
	for (my $i = 0; $i < scalar @$ValueSet; $i++) {
		if ($ValueSet->[$i] == $ValueSet->[0]) {
			push @SmallestKeys, $KeySet->[$i];
		}
	}
	if (scalar @SmallestKeys > 1) { # more than one variant has the smallest RD
		$PrevObj =~ /(\w+)\:\D\.(\D+\d+)\D/g;
		my $keyGene = $1;
		my $keyRes = $2;
		my @hits = grep(/$keyGene\:\D\.$keyRes\D/g, @SmallestKeys);

		if (scalar @hits > 0) {
			unshift @SmallestKeys, $hits[0];
		}
	}
	return shift @SmallestKeys;
}
