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
    print "\nUsage: SuperClustersID.pl RD.out_file_name_in_./Results Epsilon MinPts\n\n";
    exit;
}

my $inputFile = "./Results/$ARGV[0]";
my $Epsilon = $ARGV[1];
my $MinPts = $ARGV[2];
#my $Eta = $ARGV[3];
#my $CutoffEnd = $ARGV[4];

##########################  Reading from RD.out  #############################

my @InitialSet;

open(IN, "<$inputFile") || die "Can't open $inputFile: $!";
while (my $line = <IN>) {
	chomp $line;
	my @tabs2 = split(/\t/,$line);
	push @InitialSet, [$tabs2[0],$tabs2[1]];
}

###############################################################################

print "Preparing Clusters\n";
my %Clusters;
my @ClusterArray;


# Identify super clusters:

my $scs = 0; # super cluster start
for (my $i = 1; $i < scalar @InitialSet; $i++) {
	#print "i=$i\n";
	if ( $InitialSet[$i][1] == 10 ) {
		#print "RD($i)=10\n";

		if ( $InitialSet[$i-1][1] == 10 ) {
			$scs = $i;
		}
		else {
			$Clusters{SuperClusters}{$scs} = $i;
			$scs = $i;
		}		
	}
}
#print Dumper \%Clusters;

# Go inside each super cluster and start scanning from smallest r(x):
my $idSuperCluster = -1;

foreach my $superC ( sort { $Clusters{SuperClusters}{$a} <=> $Clusters{SuperClusters}{$b} } keys $Clusters{SuperClusters} ) {
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
	if ($Clusters{SuperClusters}{$superC} - $superC >= $MinPts) {
		$idSuperCluster++;
		my $ClusTot = 0;
		for (my $t = $superC+1; $t <= $Clusters{SuperClusters}{$superC}-1; $t++) { # start+1 b/c I want to get rid of the first tall peak, (-1 b/c end has included the start of the next)
			$ClusTot = $ClusTot + $InitialSet[$t][1];
		}
		my $ClusAvg = $ClusTot/($Clusters{SuperClusters}{$superC}-1 - $superC);
		push @ClusterArray, [$superC,$Clusters{SuperClusters}{$superC}-1,9.5,$ClusAvg,$idSuperCluster."."."0"."."."0"]; # Artificially recording the super cluster. Randomly picked 9.5. (-1 b/c end has included the start of the next)
		
	}

	# Start scanning from epsilon prime= minimum r(x) to maximum r(x):

	my $nsubc = 0; # number of sub-clusters
	my $nsubcPre = 0; # previous number of sub-clusters
	my $TempSubClusterRef;
	my $idLevel = 1;

	for (my $Cutoff = $InitialSet[$MinRDat][1] + 0.1; $Cutoff <=  $MaxRD + 0.5 ; $Cutoff += 0.1) {
		my $idSubCluster = 1;
		my $n = $superC;
		my $counter = 0;
		my $s = $superC;
		while ($n <= $Clusters{SuperClusters}{$superC}) {
			if ($InitialSet[$n][1] < $Cutoff) {
				$counter++;
				if ($counter >= $MinPts) {
					$Clusters{SubClusters}{$superC}{$s} = $n;
				}
			$n++;
			}
			else {
				$s = $n;
				$counter = 0;
				$n++;
			}
		}
		
		if (exists $Clusters{SubClusters}{$superC}) {  # if at leat one sub-cluster is found in the corresponding super cluster:
			$nsubc = scalar keys $Clusters{SubClusters}{$superC};
			#print "At cutoff=$Cutoff\tnsubc = $nsubc\t nsubc_previous=$nsubcPre\n";
			#print "\tnew cluster\n\t";
			#print Dumper $Clusters{SubClusters}{$superC};
			
						
			if ($nsubc != $nsubcPre) {
				my @OrderedS = sort { $Clusters{SubClusters}{$superC}{$a} <=> $Clusters{SubClusters}{$superC}{$b} } keys $Clusters{SubClusters}{$superC};
				my @OrderedE = @{$Clusters{SubClusters}{$superC}}{@OrderedS};
				foreach my $start (@OrderedS) {
					my $ClusTot = 0;
					for (my $t = $start+1; $t <= $Clusters{SubClusters}{$superC}{$start}; $t++) { # start+1 b/c I want to get rid of the first tall peak
						$ClusTot = $ClusTot + $InitialSet[$t][1];
					}
					my $ClusAvg = $ClusTot/($Clusters{SubClusters}{$superC}{$start}-$start);
					push @ClusterArray, [$start,$Clusters{SubClusters}{$superC}{$start},$Cutoff,$ClusAvg,$idSuperCluster.".".$idLevel.".".$idSubCluster];
					$idSubCluster++;
				}
				$nsubcPre = $nsubc;
				$idLevel++;
			}
			elsif ($nsubc == $nsubcPre && defined $TempSubClusterRef) {
				# print "\tinside the second loop\n";
				# print "\tTemporary cluster\n\t";
				# print Dumper $TempSubClusterRef;
				my $OverlapCheck = 1;
				foreach my $subC (keys $Clusters{SubClusters}{$superC}) {
					foreach my $subCpre (keys $TempSubClusterRef) {
						my ($ss, $se, $es, $ee); # s-start, e-end ; new cluster comes first
						$ss = $subC - $subCpre;
						$se = $subC - $TempSubClusterRef->{$subCpre};
						$es = $Clusters{SubClusters}{$superC}{$subC} - $subCpre;
						$ee = $Clusters{SubClusters}{$superC}{$subC} - $TempSubClusterRef->{$subCpre};
						#print "\t\tFor subC=$subC, subCpre=$subCpre: ($ss, $es, $se, $ee)";
						if (($ss*$es)> 0 && ($se*$ee) > 0) {
							$OverlapCheck = 0; # one of new sub clusters does not overlap with the previous sub cluster in consideration, check the next previous sub cluster
							#print "\tOverlapCheck = $OverlapCheck  continue checking\n";
						}
						else {
							$OverlapCheck = 1;
							#print "\tOverlapCheck = $OverlapCheck  stop checking, found one overlap\n";
							last; # one of new sub clusters overlaps with one of the previous sub clusters, stop checking
						}
					}
					if ($OverlapCheck == 0) {
						#print "\tOverlapCheck=0 passed subC=$subC\n";
						last; # one of new sub clusters does not overlap with any of the previous sub clusters, stop checking for other non-overlapping new clusters. One is enough!
					}
				}
				if ($OverlapCheck == 0) {
					my @OrderedS = sort { $Clusters{SubClusters}{$superC}{$a} <=> $Clusters{SubClusters}{$superC}{$b} } keys $Clusters{SubClusters}{$superC};
					my @OrderedE = @{$Clusters{SubClusters}{$superC}}{@OrderedS};
					foreach my $start (@OrderedS) {
						my $ClusTot = 0;
						for (my $t = $start+1; $t <= $Clusters{SubClusters}{$superC}{$start}; $t++) { # start+1 b/c I want to get rid of the first tall peak
							$ClusTot = $ClusTot + $InitialSet[$t][1];
						}
						my $ClusAvg = $ClusTot/($Clusters{SubClusters}{$superC}{$start} - $start);
						push @ClusterArray, [$start,$Clusters{SubClusters}{$superC}{$start},$Cutoff,$ClusAvg,$idSuperCluster.".".$idLevel.".".$idSubCluster];
						$idSubCluster++;
					}
					$idLevel++;
				}
			}
		}
		$TempSubClusterRef = $Clusters{SubClusters}{$superC}; # we store the current subclusters hash to detect overlapping and covering
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

my $OutFile1 = "./Results/$ARGV[0].SuperClustersID.plot";
open (OUT, ">$OutFile1");

for (my $i = 0; $i < scalar @ClusterArray; $i++) {
	print OUT "$ClusterArray[$i][0]\t$ClusterArray[$i][1]\t$ClusterArray[$i][2]\t$ClusterArray[$i][3]\t$ClusterArray[$i][4]\n";
}
close (OUT);

########################################

my $OutFile2 = "./Results/$ARGV[0].SuperClustersID.clusters";

open (OUT, ">$OutFile2");

print OUT "Cluster\tGene/Drug\tMutation/Gene\tDegree_Connectivity\tCloseness_Centrality\tGeodesic_From_Centroid\tRecurrence\tEpsilon_prime\tAvg_density\tCovering_clusters\n";

for (my $i = 0; $i < scalar @ClusterArray; $i++) {
	for (my $j = $ClusterArray[$i][0]; $j <= $ClusterArray[$i][1]; $j++) {
		my @characters = split(":",$InitialSet[$j][0]);
		my $Gene = $characters[0];
		my $Mutation = $characters[1];
		print OUT "$ClusterArray[$i][4]\t$Gene\t$Mutation\t0\t0\t1\t0\t$ClusterArray[$i][2]\t$ClusterArray[$i][3]\t$ClusterArray[$i][5]\n";
	}
}

close (OUT);

########################################
# Shiny File

my $OutFile3 = new FileHandle;
my $outFilename3 = "./Results/$ARGV[0].clusters.shiny.R";
$OutFile3->open( $outFilename3 , "w"  );

$OutFile3->print( "library(shiny)\n" );
$OutFile3->print( "ui <- basicPage(
  plotOutput(\"plot1\",click = \"plot_click\", hover = \"plot_hover\", brush = \"plot_brush\" ,  height = 900),
  verbatimTextOutput(\"info\")
)\n" );
$OutFile3->print( "y = read.table(\"./$ARGV[0]\")\n
z = read.table(\"./$ARGV[0].SuperClustersID.plot\")\n" );
$OutFile3->print( "RD<-y[[2]]\nID<-y[[1]]\nx0<-z[[1]]\ny0<-z[[3]]\nx1<-z[[2]]+1\ny1<-z[[3]]\nCluster<-z[[5]]\n" );
$OutFile3->print( "server <- function(input, output) {\n  output\$plot1 <- renderPlot({\n	barplot(RD,names.arg=ID,main=\"Reachability Plot:Epsilon=8 MinPts=4\",col=\"Red\", cex.names=0.6,border=NA, space=0, las=2, ylab=\"Reachabilty Distance (A)\")\n	segments (x0,y0,x1,y1)\n
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

system ("Rscript ClustersLines.R $inputFile ./Results/$ARGV[0].SuperClustersID.plot SuperClustersID.$ARGV[0].pdf $Epsilon $MinPts");
system ("Rscript HorizClustersLines.R $inputFile ./Results/$ARGV[0].SuperClustersID.plot SuperClustersID.$ARGV[0].Horiz.pdf $Epsilon $MinPts");
system ("Rscript EasyClustersLines.R $inputFile ./Results/$ARGV[0].SuperClustersID.plot SuperClustersID.$ARGV[0].Easy.pdf $Epsilon $MinPts");
print "Done.\n";

#########################################################################