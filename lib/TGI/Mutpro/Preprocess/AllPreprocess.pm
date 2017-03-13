package TGI::Mutpro::Preprocess::AllPreprocess;
#
#----------------------------------
# $Authors: Adam D Scott
# $Date: 2016-09-29 $
# $Revision: 2016-09-29 $
# $URL: $
# $Doc: Run all preprocessing steps (assuming bsub capable) $ 
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
use Parallel::ForkManager;

#use TGI::Mutpro::Preprocess::Uppro;
use TGI::Mutpro::Preprocess::Calroi;
use TGI::Mutpro::Preprocess::Statis;
use TGI::Mutpro::Preprocess::Anno;
use TGI::Mutpro::Preprocess::Trans;
use TGI::Mutpro::Preprocess::Cosmic;
use TGI::Mutpro::Preprocess::Prior;

my $MINDISTANCE = "minDistance";
my $AVGDISTANCE = "averageDistance";
my $CALROI = "calroi";
my $STATIS = "statis";
my $ANNO = "anno";
my $TRANS = "trans";
my $COSMIC = "cosmic";
my $PRIOR = "prior";
my $CPU = "cpu";
my $BSUB = "bsub";
my $NONE = "none";

sub new {
	my $class = shift;
	my $this = {};
	$this->{'command'} = "";
	$this->{'_OUTPUT_DIR'} = getcwd;
	$this->{'_BLAT'} = "blat";
	$this->{'GRCh'} = undef;
	$this->{'release'} = undef;
	$this->{'max_3d_dis'} = 100;
	$this->{'p_value_cutoff'} = 1;
	$this->{'min_seq_dis'} = 0;
	$this->{'start'} = $CALROI; 
	$this->{'parallel'} = $CPU;
	$this->{'status'} = 0;
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
		'blat=s'   => \$this->{'_BLAT'},
		'grch=i'	=> \$this->{'GRCh'},
		'release=i'	=> \$this->{'release'},
		'3d-distance-cutoff=i'	=> \$this->{'max_3d_dis'},
		'p-value-cutoff=i'	=> \$this->{'p_value_cutoff'},
		'linear-cutoff=i'   => \$this->{'min_seq_dis'},
		'start=s' => \$this->{'start'} ,
		'help' => \$help,
	);
	if ($help) { print STDERR help_text(); exit 0; }
	unless ($options) { die $this->help_text(); }

	$this->steps();
	return;
}

sub checkStatus {
	my ( $this , $step ) = @_;
	if ( $this->{'status'} ) {
		warn "HotSpot3D::AllPreprocess::checkStatus error: prior step (".$step.") failed!\n";
		exit $this->{'status'};
	}
	return $this->{'status'};
}

sub command {
	my $this = shift;
	if ( @_ ) {
		$this->{'command'} = \@_;
	}
	return $this->{'command'};
}

sub outputDir {
	my $this = shift;
	if ( @_ ) { $this->{'_OUTPUT_DIR'} = shift; }
	return $this->{'_OUTPUT_DIR'};
}

sub blat {
	my $this = shift;
	if ( @_ ) { $this->{'_BLAT'} = shift; }
	return $this->{'_BLAT'};
}

sub grch {
	my $this = shift;
	if ( @_ ) { $this->{'GRCh'} = shift; }
	return $this->{'GRCh'};
}

sub release {
	my $this = shift;
	if ( @_ ) { $this->{'release'} = shift; }
	return $this->{'release'};
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

sub minseqdis {
	my $this = shift;
	if ( @_ ) { $this->{'min_seq_dis'} = shift; }
	return $this->{'min_seq_dis'};
}

sub start {
	my $this = shift;
	if ( @_ ) { $this->{'start'} = shift; }
	return $this->{'start'};
}

sub status {
	my $this = shift;
	if ( @_ ) { $this->{'status'} = shift; }
	return $this->{'status'};
}

sub calroi {
	my $this = shift;
	$this->checkStatus( "uppro/calpro" );
	my $cmd = "hotspot3d calroi --output-dir ".$this->outputDir();
	print STDOUT "running: ".$cmd."\n";
	$this->status( system( $cmd ) );
	return;
}

sub statis {
	my $this = shift;
	$this->checkStatus( $CALROI );
	my $cmd = "hotspot3d statis --output-dir ".$this->outputDir();
	print STDOUT "running: ".$cmd."\n";
	$this->status( system( $cmd ) );
	return $this->status();
}

sub anno {
	my $this = shift;
	$this->checkStatus( $STATIS );
	my $cmd = "hotspot3d anno --output-dir ".$this->outputDir();
	print STDOUT "running: ".$cmd."\n";
	$this->status( system( $cmd ) );
	return $this->status();
}

sub trans {
	my $this = shift;
	$this->checkStatus( $ANNO );
	my $cmd = "hotspot3d trans --output-dir ".$this->outputDir();
	$cmd .= " --blat ".$this->blat();
	if ( $this->grch() ) {
		$cmd .= " --grch ".$this->grch();
	}
	if ( $this->release() ) {
		$cmd .= " --release ".$this->release();
	}
	print STDOUT "running: ".$cmd."\n";
	$this->status( system( $cmd ) );
	return $this->status();
}

sub cosmic {
	my $this = shift;
	$this->checkStatus( $TRANS );
	my $cmd = "hotspot3d cosmic --output-dir ".$this->outputDir();
	print STDOUT "running: ".$cmd."\n";
	$this->status( system( $cmd ) );
	return $this->status();
}

sub prior {
	my $this = shift;
	$this->checkStatus( $COSMIC );
	my $cmd = "hotspot3d prior --output-dir ".$this->outputDir();
	$cmd .= " --p-value-cutoff ".$this->pvaluecutoff();
	$cmd .= " --3d-distance-cutoff ".$this->max3ddis();
	$cmd .= " --linear-cutoff ".$this->minseqdis();
	print STDOUT "running: ".$cmd."\n";
	$this->status( system( $cmd ) );
	return $this->status();
}

sub steps {
	my $this = shift;
	my $ok = 0;
	#if ( $this->{'start'} eq $UPPRO ) {
	#	fork 
	#} elsif ( $this->{'start'} eq $CALROI ) {
	if ( $this->{'start'} eq $CALROI ) {
		$this->calroi();
		$this->statis();
		$this->anno();
		$this->trans();
		$this->cosmic();
		$this->prior();
	} elsif ( $this->{'start'} eq $STATIS ) {
		$this->statis();
		$this->anno();
		$this->trans();
		$this->cosmic();
		$this->prior();
	} elsif ( $this->{'start'} eq $ANNO ) {
		$this->anno();
		$this->trans();
		$this->cosmic();
		$this->prior();
	} elsif ( $this->{'start'} eq $TRANS ) {
		$this->trans();
		$this->cosmic();
		$this->prior();
	} elsif ( $this->{'start'} eq $COSMIC ) {
		$this->cosmic();
		$this->prior();
	} elsif ( $this->{'start'} eq $PRIOR ) {
		$this->prior();
	} else {
		die "HotSpot3D::AllPreprocess error: desired starting step unclear.\n".$this->help_text();
	}

	return;
}

sub help_text{
	my $this = shift;
		return <<HELP

Usage: hotspot3d prep [options]

                                   REQUIRED
--output-dir                       Output directory of proximity files

                                   OPTIONAL
--start                            What step to start on ( calroi , statis , anno , trans , cosmic , prior ), default is calroi
--blat                             Installation of blat to use for trans (defaults to your system default)
--grch                             Genome build (37 or 38), defaults to 38 or according to --release input
--release                          Ensembl release verion (55-87), defaults to 87 or to the latest release according to --grch input
                                       Note that releases 55-75 correspond to GRCh37 & 76-87 correspond to GRCh38
--p-value-cutoff                   p_value cutoff(<=) for prior, default is 0.05
--3d-distance-cutoff               3D distance cutoff (<= Angstroms) for prior, default is 20
--linear-cutoff                    Linear distance cutoff (> peptides) for prior, default is 0


--help                       this message

HELP

}

1;

__END__
 
=head1 NAME

TGI::Mutpro::Preprocess::Uppro - Create & update proximity files.

=head1 SYNOPSIS

  use TGI::Mutpro::Preprocess::Uppro;

=head1 DESCRIPTION

TGI::Mutpro::Preprocess::Uppro is to be used to create & update proximity files.
It is the first step of preprocessing procedure.


=head1 AUTHOR

Beifang Niu E<lt>beifang.cn@gmail.comE<gt>

=head1 SEE ALSO

https://github.com/ding-lab/hotspot3d

=head1 LICENSE

This library is free software with MIT licence; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

