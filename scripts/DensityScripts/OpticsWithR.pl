use strict;
use warnings;

use LWP::Simple;
use FileHandle;
use Data::Dumper qw(Dumper);
use Scalar::Util qw(looks_like_number);
use List::Util qw[min max];

######################## Command Line Arguments ################################

# quit unless we have the correct number of command-line args
my $num_args = $#ARGV + 1;
if ($num_args != 3) {
	print "\nUsage: OpticsWithR.pl epsilon MinPts InputFileName_in_./Test/\n\n";
    exit;
}

my $Epsilon = $ARGV[0];
my $MinPts = $ARGV[1];
print "Preparing Reachability Data\n";

############# Extracting the set of objects ##############

my %SetOfNodes;

my $file = "./Test/$ARGV[2]";
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
my @SetOfCoresThenEdges = ( @SetOfCores, @SetOfEdges );

# print "SetofCores=";
# print Dumper \@SetOfCores;
# print "SetofEdges=";
# print Dumper \@SetOfEdges;
# #push @SetofCores, @SetofEdges;
# print "both=";
# print Dumper \@SetOfCoresThenEdges;
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
		#print "p=$p and ";
		#print "CD=$CD\n";
		push @OrderedNodes, [$p,$RD,$CD]; # write to the file 
		if (defined $CD) {
			OrderSeedsUpdate(\%neighbors,$p,$CD, \%OrderSeeds, \%SetOfNodes);
			#print "For p=$p, OrderSeeds= \n";
			#print Dumper \%OrderSeeds;
			while (scalar keys %OrderSeeds != 0) {
				my @SeedKeys = sort { $OrderSeeds{$a} <=> $OrderSeeds{$b} } keys %OrderSeeds;
				my @SeedValues = @OrderSeeds{@SeedKeys};
				my $CurrentObject =  $SeedKeys[0]; # CurrentObject is the object having the least RD in OrderSeeds
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
					#print "\tCurrent object is a core.\n Updated Order seeds list\n\t";
					OrderSeedsUpdate(\%neighbors,$CurrentObject,$CD, \%OrderSeeds, \%SetOfNodes);
					#print Dumper \%OrderSeeds;
				}
			}
		}
		#print "p=$p, OrderedNodes= \n";
		#print Dumper \@OrderedNodes;
	}
}

######################  Reachability Plot   ########################

my @data;
my @dataX;
my @dataY;
foreach my $x (1...scalar keys %SetOfNodes) {
	push @dataX, "$OrderedNodes[$x-1][0]";
	if (defined $OrderedNodes[$x-1][1]) {
		push @dataY, "$OrderedNodes[$x-1][1]";
	}
	else {
		push @dataY, "10"; # just a large number to denote infinity
	}
}
push @data, [@dataX];
push @data, [@dataY];
#print Dumper \@data;

my $OrderedFile = "./Results/RD.$Epsilon.$MinPts.$ARGV[2]";
open (OUT, ">$OrderedFile");
foreach my $x (1...scalar keys %SetOfNodes) {
	if (defined $OrderedNodes[$x-1][1]) {
		print OUT "$OrderedNodes[$x-1][0]\t $OrderedNodes[$x-1][1]\n";
	}
	else {
		print OUT "$OrderedNodes[$x-1][0]\t 10\n";
	}
}
close (OUT);

#system ("Rscript PlotR.R $OrderedFile RD.$Epsilon.$MinPts.$ARGV[2].pdf $Epsilon $MinPts");

print "Done.\n";

#########################  Functions  ############################

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