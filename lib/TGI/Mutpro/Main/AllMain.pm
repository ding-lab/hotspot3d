package TGI::Mutpro::Main::AllMain;
#
#----------------------------------
# $Authors: Adam D Scott
# $Date: 2016-11-10 $
# $Revision: 2016-09-29 $
# $URL: $
# $Doc: Run all main steps $ 
#----------------------------------
#
use strict;
use warnings;
our $VERSION = '0.0';

use Carp;
use Cwd;
use Getopt::Long;
use IO::File;
use FileHandle;

use TGI::Mutpro::Main::Proximity;
use TGI::Mutpro::Main::Post;
use TGI::Mutpro::Main::Cluster;
use TGI::Mutpro::Main::Summary;
use TGI::Mutpro::Main::Significance;
#use TGI::Mutpro::Main::Visual;

my $PVALUEDEFAULT = 0.05;
my $DISTANCEDEFAULT = 10;
my $MAXDISTANCE = 100;
my $WEIGHT = "weight";
my $RECURRENCE = "recurrence";
my $UNIQUE = "unique";
my $AVERAGEDISTANCE = "average";
my $SHORTESTDISTANCE = "shortest";
my $NETWORK = "network";
my $DENSITY = "density";

my $SEARCH = "search";
my $POST = "post";
my $CLUSTER = "cluster";
my $SUMMARY = "summary";
my $SIGNIFICANCE = "sigclus";
my $VISUAL = "visual";

sub new {
	my $class = shift;
	my $this = {};
	$this->{'command'} = "";
	$this->{'start'} = $SEARCH;
    $this->{'maf_file'} = undef;
    $this->{'data_dir'} = undef;
    $this->{'drugport_file'} = undef;
    $this->{'output_prefix'} = '3D_Proximity';
    $this->{'skip_silent'} = undef;
    $this->{'missense_only'} = undef;
    $this->{'p_value_cutoff'} = undef;
    $this->{'3d_distance_cutoff'} = undef;
    $this->{'linear_cutoff'} = 0;
    $this->{'amino_acid_header'} = "amino_acid_change";
    $this->{'transcript_id_header'} = "transcript_name";
    $this->{'weight_header'} = $WEIGHT;
    $this->{'max_radius'} = 10;
    $this->{'clustering'} = undef;
    $this->{'vertex_type'} = $RECURRENCE;
    $this->{'distance_measure'} = $AVERAGEDISTANCE;
	$this->{'simulations'} = 1000000;
	bless $this, $class;
	$this->process();
	return $this;
}

sub process {
	my $this = shift;
	my ($help, $options);
	unless (@ARGV) { die $this->help_text(); }
	$this->{'command'} = \@ARGV;
	$options = GetOptions (
		'output-dir=s'	=> \$this->{'_OUTPUT_DIR'},
		'start=s' => \$this->{'start'} ,

        'maf-file=s' => \$this->{'maf_file'},
        'prep-dir=s'    => \$this->{'data_dir'},
        'drugport-file=s' => \$this->{'drugport_file'},
        'skip-silent' => \$this->{'skip_silent'},
        'missense-only' => \$this->{'missense_only'},
        'output-prefix=s' => \$this->{'output_prefix'},
        'p-value-cutoff=f' => \$this->{'p_value_cutoff'},
        '3d-distance-cutoff=f' => \$this->{'3d_distance_cutoff'},
        'linear-cutoff=i' => \$this->{'linear_cutoff'},
        'amino-acid-header=s' => \$this->{'amino_acid_header'},
        'transcript-id-header=s' => \$this->{'transcript_id_header'},
        'weight-header=s' => \$this->{'weight_header'},
        'max-radius=f' => \$this->{'max_radius'},
        'clustering=s' => \$this->{'clustering'},
        'vertex-type=s' => \$this->{'vertex_type'},
        'distance-measure=s' => \$this->{'distance_measure'},
        'simulations=i' => \$this->{'simulations'},

		'help' => \$help,
	);
	if ($help) { print STDERR help_text(); exit 0; }
	unless ($options) { die $this->help_text(); }

	$this->steps();
	return;
}

sub command {
	my $this = shift;
	if ( @_ ) {
		$this->{'command'} = \@_;
	}
	return $this->{'command'};
}

sub maffile {
	my $this = shift;
	if ( @_ ) { $this->{'maf_file'} = shift; }
	return $this->{'maf_file'};
}

sub prepdir {
	my $this = shift;
	if ( @_ ) { $this->{'data_dir'} = shift; }
	return $this->{'data_dir'};
}

sub skipsilent {
	my $this = shift;
	if ( @_ ) { $this->{'skip_silent'} = shift; }
	return $this->{'skip_silent'};
}

sub missenseonly {
	my $this = shift;
	if ( @_ ) { $this->{'missense_only'} = shift; }
	return $this->{'missense_only'};
}

sub pvaluecutoff {
	my $this = shift;
	if ( @_ ) { $this->{'p_value_cutoff'} = shift; }
	return $this->{'p_value_cutoff'};
}

sub max3ddis {
	my $this = shift;
	if ( @_ ) { $this->{'max_3d_dis'} = shift; }
	return $this->{'max_3d_dis'};
}

sub linearcutoff {
	my $this = shift;
	if ( @_ ) { $this->{'linear_cutoff'} = shift; }
	return $this->{'linear_cutoff'};
}

sub aminoacidheader {
	my $this = shift;
	if ( @_ ) { $this->{'amino_acid_header'} = shift; }
	return $this->{'amino_acid_header'};
}

sub transcriptidheader {
	my $this = shift;
	if ( @_ ) { $this->{'transcript_id_header'} = shift; }
	return $this->{'transcript_id_header'};
}

sub weightheader {
	my $this = shift;
	if ( @_ ) { $this->{'weight_header'} = shift; }
	return $this->{'weight_header'};
}

sub maxradius {
	my $this = shift;
	if ( @_ ) { $this->{'max_radius'} = shift; }
	return $this->{'max_radius'};
}

sub clustering {
	my $this = shift;
	if ( @_ ) { $this->{'clustering'} = shift; }
	return $this->{'clustering'};
}

sub vertextype {
	my $this = shift;
	if ( @_ ) { $this->{'vertex_type'} = shift; }
	return $this->{'vertex_type'};
}

sub distancemeasure {
	my $this = shift;
	if ( @_ ) { $this->{'distance_measure'} = shift; }
	return $this->{'distance_measure'};
}

sub simulations {
	my $this = shift;
	if ( @_ ) { $this->{'simulations'} = shift; }
	return $this->{'simulations'};
}

sub start {
	my $this = shift;
	if ( @_ ) { $this->{'start'} = shift; }
	return $this->{'start'};
}

sub search {
	my $this = shift;
	my $cmd = "hotspot3d search --maf-file ".$this->maffile()." --prep-dir ".$this->outputDir();
	if ( $this->drugportfile() ) { $cmd .= " --drugport-file ".$this->drugportfile(); }
	if ( $this->outputprefix() ) { $cmd .= " --output-prefix ".$this->outputprefix(); }
	if ( $this->skipsilent() ) { $cmd .= " --skip-silent "; }
	if ( $this->missenseonly() ) { $cmd .= " --missense-only "; }
	if ( $this->pvaluecutoff() ) { $cmd .= " --p-value-cutoff ".$this->pvaluecutoff(); }
	if ( $this->max3ddis() ) { $cmd .= " --3d-distance-cutoff ".$this->max3ddis(); }
	if ( $this->linearcutoff() ) { $cmd .= " --linear-cutoff ".$this->linearcutoff(); }
	if ( $this->transcriptidheader() ) { $cmd .= " --transcript-id-header ".$this->transcriptidheader(); }
	if ( $this->aminoacidheader() ) { $cmd .= " --amino-acid-header ".$this->aminoacidheader(); }
	print STDOUT "running: ".$cmd."\n";
	return system( $cmd );
}

sub post {
	my $this = shift;
	my $cmd = "hotspot3d post --maf-file ".$this->maffile();
	if ( $this->transcriptidheader() ) { $cmd .= " --transcript-id-header ".$this->transcriptidheader(); }
	if ( $this->aminoacidheader() ) { $cmd .= " --amino-acid-header ".$this->aminoacidheader(); }
	print STDOUT "running: ".$cmd."\n";
	return system( $cmd );
}

sub cluster {
	my $this = shift;
	my $cmd = "hotspot3d cluster --pairwise-file 3D_Proximity.pairwise ";
	if ( $this->outputprefix() ) { $cmd .= " --output-prefix ".$this->outputprefix(); }
	if ( $this->pvaluecutoff() ) { $cmd .= " --p-value-cutoff ".$this->pvaluecutoff(); }
	if ( $this->max3ddis() ) { $cmd .= " --3d-distance-cutoff ".$this->max3ddis(); }
	if ( $this->linearcutoff() ) { $cmd .= " --linear-cutoff ".$this->linearcutoff(); }
	if ( $this->clustering() ) { $cmd .= " --clustering ".$this->clustering(); }
	if ( $this->maxradius() ) { $cmd .= " --max-radius ".$this->maxradius(); }
	if ( $this->distancemeasure() ) { $cmd .= " --distance-measure ".$this->distancemeasure(); }
	if ( $this->vertextype() ) { $cmd .= " --vertex-type ".$this->vertextype(); }
	if ( $this->transcriptidheader() ) { $cmd .= " --transcript-id-header ".$this->transcriptidheader(); }
	if ( $this->aminoacidheader() ) { $cmd .= " --amino-acid-header ".$this->aminoacidheader(); }
	if ( $this->weightheader() ) { $cmd .= " --weight-header ".$this->weightheader(); }
	print STDOUT "running: ".$cmd."\n";
	return system( $cmd );
}

sub summary {
	my $this = shift;
#	my $cmd = "hotspot3d summary --clusters-file *.clusters";
#	print STDOUT "running: ".$cmd."\n";
#	return system( $cmd );
}

sub sigclus {
	my $this = shift;
#	my $cmd = "hotspot3d sigclus --prep-dir ".$this->outputDir()." --pairwise-file ";
#	print STDOUT "running: ".$cmd."\n";
#	return system( $cmd );
}

sub visual {
	my $this = shift;
#	my $cmd = "hotspot3d visual --output-dir ".$this->outputDir();
#	$cmd .= " --p-value-cutoff ".$this->pvaluecutoff();
#	$cmd .= " --3d-distance-cutoff ".$this->max3ddis();
#	$cmd .= " --linear-cutoff ".$this->minseqdis();
#	print STDOUT "running: ".$cmd."\n";
#	return system( $cmd );
}

sub steps {
	my $this = shift;
	my $ok = 0;
	if ( $this->{'start'} eq $SEARCH ) {
		$this->search();
		$this->post();
		$this->cluster();
		#$this->sigclus();
		#$this->summary();
		#$this->visual();
	} elsif ( $this->{'start'} eq $POST ) {
		$this->post();
		$this->cluster();
		#$this->sigclus();
		#$this->summary();
		#$this->visual();
	} elsif ( $this->{'start'} eq $CLUSTER ) {
		$this->cluster();
		#$this->sigclus();
		#$this->summary();
		#$this->visual();
	} elsif ( $this->{'start'} eq $SIGNIFICANCE ) {
		#$this->sigclus();
		#$this->summary();
		#$this->visual();
	} elsif ( $this->{'start'} eq $SUMMARY ) {
		#$this->summary();
		#$this->visual();
	} elsif ( $this->{'start'} eq $VISUAL ) {
		#$this->visual();
	} else {
		die "HotSpot3D::AllMain error: desired starting step unclear.\n".$this->help_text();
	}

	return;
}

sub help_text{
	my $this = shift;
		return <<HELP

Usage: hotspot3d main [options]

                             REQUIRED
--output-dir                 Output directory of proximity files
--maf-file                   .maf file used in proximity search step (used if vertex-type = recurrence)

                             OPTIONAL
--start                      Step to start on (search, post, cluster, sigclus, summary, visual), default: search
--drugport-file              DrugPort database parsing results file
--output-prefix              Prefix of output files, default: 3D_Proximity
--skip-silent                skip silent mutations, default: no
--missense-only              missense mutation only, default: no
--p-value-cutoff             P_value cutoff (<), default: 0.05 (if 3d-distance-cutoff also not set)
--3d-distance-cutoff         3D distance cutoff (<), default: 100 (if p-value-cutoff also not set)
--linear-cutoff              Linear distance cutoff (> peptides), default: 0
--amino-acid-header          .maf file column header for amino acid changes, default: amino_acid_change
--transcript-id-header       .maf file column header for transcript id's, default: transcript_name
--weight-header              .maf file column header for mutation weight, default: weight (used if vertex-type = weight)
--max-radius                 Maximum cluster radius (max network geodesic from centroid, <= Angstroms), default: 10
--clustering                 Cluster using network or density-based methods (network or density), default: network
--vertex-type                Graph vertex type for network-based clustering (recurrence, unique, or weight), default: recurrence
--distance-measure           Pair distance to use (shortest or average), default: average
--simulations                Number of simulations, default = 1000000

--help                       this message

HELP

}

1;
