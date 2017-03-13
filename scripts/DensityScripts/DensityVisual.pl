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
    print "\nUsage: DensityVisual.pl *.clusters_file_name_in_./Results PairwiseFile_name PDB_ID\n\n";
    exit;
}

my $ClustersFile = "./Results/$ARGV[0]";
my $PairwiseFile = $ARGV[1];
my $PDB = $ARGV[2];

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

my $this = {};
#my @spectrum = ("blue_yellow","rainbow","green_red","red_yellow","blue_green","cyan_red");
my $ChainColors = ['pink' , 'palegreen' , 'lightblue' , 'lightmagenta' , 'lightorange' , 'lightpink' , 'paleyellow' , 'lightteal' , 'aquamarine' , 'palecyan'];
my $SurfaceColors = ['violetpurple','violet','deeppurple', 'purple','lightmagenta', 'blue', 'lightblue', 'bluewhite', 'oxygen' , 'green', 'deepolive' , 'wheat' , 'lime' , 'brown' , 'orange' , 'yellow', 'salmon' , 'red',  'violetpurple' , 'limon' , 'sand' , 'raspberry' , 'slate', 'red','blue', 'green', 'yellow', 'deepolive' , 'deeppurple' , 'bluewhite' , 'lime' , 'purple' , 'dash' , 'orange' , 'brown' , 'salmon' , 'oxygen' , 'wheat' , 'violetpurple' , 'limon' , 'sand' , 'raspberry' , 'slate','red','blue', 'green', 'yellow', 'deepolive' , 'deeppurple' , 'bluewhite' , 'lime' , 'purple' , 'dash' , 'orange' , 'brown' , 'salmon' , 'oxygen' , 'wheat' , 'violetpurple' , 'limon' , 'sand' , 'raspberry' , 'slate'];
#my $spectrumPalette = shift @spectrum;
#my $ClusterSurfaceColor = shift $SurfaceColors;
my $currentColor = shift $SurfaceColors;


for (my $i = 0; $i < scalar @InputClusters; $i++) {
	my @IDs = split(/\./,$InputClusters[$i][0]);
	$this->{ClusterStructure}->{$IDs[0]}->{$IDs[1]} = $IDs[2];
}
#print Dumper \%ClusterStructure;

getMappingLocations($this);
print "Done getting mapping locations.\n";

#################  Clusters present in the given PDB structure  #####################

my $OutFile1 = "./Results/SuperClustersID.$ARGV[0].presentClusters";
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
my $outFilename = "./Results/$ARGV[0].$PDB.pml";
$fh->open( $outFilename , "w"  );
$fh->print( "reinitialize everything;\n" );
$fh->print( "load http://www.rcsb.org/pdb/files/".$PDB.".pdb;\n" );

$fh->print( "viewport 480,480;\n" );
$fh->print( "preset.publication(\"".$PDB."\");\n" );
$fh->print( "#show mesh;\n" );
$fh->print( "\n" );
$fh->print( "bg_color white;\n" );
$fh->print( "\n" );
foreach my $chain (keys $this->{chains}) {
	my $chain_color = shift $ChainColors;
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

#########################################   Functions   ###########################################

sub getMappingLocations {
	my ( $this ) = @_;
		my $fh = new FileHandle;
		print STDOUT "Getting mutation-mutation pairs\n";
		unless( $fh->open( $PairwiseFile ) ) { die "Could not open pairwise file\n" };
		while ( my $line = <$fh> ) {
			chomp( $line );
			if ( $line =~ $PDB ) {
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
	return 1;
}
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
			#$fh->print( "sele ".$currentColor."_".$GeneName."_".$MutName."_".$chain.", (resi ".$res." and chain ".$chain.");\n" ); #CHECK whether needed
			$fh->print( "color ".$currentColor.", (resi ".$res." and chain ".$chain.");\n" );
			#$fh->print( "show spheres, (resi ".$res." and chain ".$chain.");\n" );
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
