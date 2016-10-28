package TGI::Mutpro::Main::Density;

use strict;
use warnings;

use LWP::Simple;
use FileHandle;
use Data::Dumper qw(Dumper);
use List::Util qw[min max];
use Getopt::Long;

our @ISA = qw( TGI::Mutpro::Main::Cluster );

my $EPSILONDEFAULT = 4;
my $MINPTSDEFAULT = 3;
my $CUTOFFSTARTDEFAULT = 2;
my $CUTOFFENDDEFAULT = 4;

sub new {
    my $class = shift;
    my $this = {};
    $this->{'pairwise_file'} = '3D_Proximity.pairwise';
    $this->{'epsilon'} = undef;
    $this->{'minpts'} = undef;
    $this->{'cutoffstart'} = undef;
    $this->{'cutoffend'} = undef;
    bless $this, $class;
    $this->process();
    return $this;
}

sub process {
	my $this = shift;
    my ( $help, $options );
    unless( @ARGV ) { die $this->help_text(); }
    $options = GetOptions (
        'pairwise-file=s' => \$this->{'pairwise_file'},
        'epsilon=f' => \$this->{'epsilon'},
        'minpts=f' => \$this->{'minpts'},
        'cut-off-start=f' => \$this->{'cut_off_start'},
        'cut-off-end=f' => \$this->{'cut_off_end'},
        'help' => \$help,
    );
    if ( $help ) { print STDERR help_text(); exit 0; }
    unless( $options ) { die $this->help_text(); }
    unless( $this->{'pairwise_file'} ) { warn 'You must provide a pairwise file! ', "\n"; die $this->help_text(); }
    unless( -e $this->{'pairwise_file'} ) { warn "The input pairwise file (".$this->{'pairwise_file'}.") does not exist! ", "\n"; die $this->help_text(); }
    if ( not defined $this->{'epsilon'} ) {
    	warn "HotSpot3D::Cluster warning: no Epsilon value given, setting to default Epsilon value = 4\n";
			$this->{'epsilon'} = $EPSILONDEFAULT;
    }
    if ( not defined $this->{'minpts'} ) {
    	warn "HotSpot3D::Cluster warning: no MinPts value given, setting to default MinPts value = 3\n";
			$this->{'minpts'} = $MINPTSDEFAULT;
    }
    if ( not defined $this->{'cut-off-start'} ) {
    	warn "HotSpot3D::Cluster warning: no cut_off_start given, setting to default value = 2\n";
			$this->{'cut-off-start'} = $CUTOFFSTARTDEFAULT;
    }
    if ( not defined $this->{'cut-off-end'} ) {
    	warn "HotSpot3D::Cluster warning: no cut_off_end given, setting to default value = 4\n";
			$this->{'cut-off-end'} = $CUTOFFENDDEFAULT;
    }
    print STDOUT "Preparing Density-based Clusters...\n";
##################
#     OPTICS     #
##################
	my %SetOfNodes; # Hash containing all the nodes
	my $file = "$this->{'pairwise_file'}";
	open( IN, "<$file" ) || die "Could not open the pairwise file $!\n";

	while ( my $line = <IN> ) {
		chomp $line;
		my @tabs = split(/\t/,$line);
		my @char19 = split("",$tabs[19]);
		my $dis = $char19[0].$char19[1].$char19[2].$char19[3].$char19[4];
		my $key1 = CombineWords($tabs[0],$tabs[4]);
		my $value1 = CombineWords($tabs[9],$tabs[13]);

		$SetOfNodes{$key1}{distances}{$value1} = $dis;
		$SetOfNodes{$value1}{distances}{$key1} = $dis;
	}

	close (IN);

	foreach my $i ( keys %SetOfNodes ) {
		$SetOfNodes{$i}{processInfo} = "False";
	}

    my @OrderedNodes;

    ################# Main OPTICS function ####################

    foreach my $p ( keys %SetOfNodes ) {
        if ( $SetOfNodes{$p}{processInfo} =~ "False" ) {
            ########## Expand Cluster Order ###########
            my %neighbors; # is a hash with keys neigbor indices whose values are mutual separations
            my %OrderSeeds; # is a hash to add seeds
            %neighbors = %{GetNeighbors($p,\%SetOfNodes)};
            $SetOfNodes{$p}{processInfo} = "True"; # set as processed
            my $RD = undef;
            my $CD;
            $CD = GetCoreDistance(\%neighbors,$this);
            push @OrderedNodes, [$p,$RD,$CD]; # write to the file 
            if (defined $CD) {
                OrderSeedsUpdate(\%neighbors,$p,$CD, \%OrderSeeds, \%SetOfNodes);
                while (scalar keys %OrderSeeds != 0) {
                    my @SeedKeys = sort { $OrderSeeds{$a} <=> $OrderSeeds{$b} } keys %OrderSeeds;
                    my @SeedValues = @OrderSeeds{@SeedKeys};
                    my $CurrentObject =  $SeedKeys[0]; # CurrentObject is the object having the least RD in OrderSeeds
                    %neighbors = %{GetNeighbors($CurrentObject,\%SetOfNodes)};
                    $SetOfNodes{$CurrentObject}{processInfo} = "True"; # set as processed
                    $RD = $SeedValues[0];
                    $CD = GetCoreDistance(\%neighbors,$this);
                    push @OrderedNodes, [$CurrentObject,$RD,$CD]; # write to the file 
                    delete $OrderSeeds{$CurrentObject};
                    if (defined $CD) {
                        OrderSeedsUpdate(\%neighbors,$CurrentObject,$CD, \%OrderSeeds, \%SetOfNodes);
                    }
                } # loop through nodes in the OrderedNodes
            } # if the node is a core
        } # if the node is not processed
    } # get nodes from the pool

    ##############################################################

    my $OrderedFile = "RD.$this->{'epsilon'}.$this->{'minpts'}.$this->{'pairwise_file'}.out";
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
    print STDOUT "Done.\n";
}
###################
#  Sub Functions  #
###################
sub GetNeighbors {
    my ($Obj, $Set_ref)=@_;
    my %neighborHash;
    foreach my $i (keys %{$Set_ref->{$Obj}->{distances}}) {
            $neighborHash{$i} = "$Set_ref->{$Obj}->{distances}->{$i}";
    }
    return \%neighborHash;
}

sub GetCoreDistance {
    my ($neighbors_ref, $this)=@_;
    my @keys = sort { $neighbors_ref->{$a} <=> $neighbors_ref->{$b} } keys %{$neighbors_ref}; # sort keys according to distances
    my @vals = @{$neighbors_ref}{@keys};
    my $CoreDist;
    if (scalar keys %{$neighbors_ref} >= $this->{'minpts'}){
            $CoreDist = $vals[$this->{'minpts'}-1]; # MinPt^th-distance
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

sub help_text{
    my $this = shift;
        return <<HELP

Usage: hotspot3d density [options]

                             REQUIRED
--pairwise-file              3D pairwise data file

                             OPTIONAL
--epsilon                    Epsilon value, default: 4
--minpts                     MinPts, default: 3
--cut-off-start              Starting Epsilon-prime value to look for clusters, default: 2 
--cut-off-end                Ending Epsilon-prime value to look for clusters, default: 4

--help                       this message

HELP

}

1;
