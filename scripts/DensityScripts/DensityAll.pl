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
if ($num_args != 4) {
	print "\nUsage: DensityAll.pl InputFileName_in_./Test epsilon MinPts PDB_ID/\n\n";
    exit;
}

my $PairwiseFileName = $ARGV[0];
my $Epsilon = $ARGV[1];
my $MinPts = $ARGV[2];
my $PDBName = $ARGV[3];


##################################################################################

print "\n Running OpticsWithR.pl....\n\n";
system ("perl OpticsWithR.pl $Epsilon $MinPts $PairwiseFileName");

print "\n Running SuperClustersID.pl....\n\n";
system ("perl SuperClustersID.pl RD.$Epsilon.$MinPts.$PairwiseFileName $Epsilon $MinPts");

print "\n Running DensityVisual.pl....\n\n";
system ("perl DensityVisual.pl RD.$Epsilon.$MinPts.$PairwiseFileName.SuperClustersID.clusters ./Test/$PairwiseFileName $PDBName");

print "\n Your Reachability plot with cluster marks, Clusters file, and the Pymol script has been produced.\n";
print "Done.\n";