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

#use TGI::Mutpro::Preprocess::Uppro;
use TGI::Mutpro::Preprocess::Calroi;
use TGI::Mutpro::Preprocess::Statis;
use TGI::Mutpro::Preprocess::Anno;
use TGI::Mutpro::Preprocess::Trans;
use TGI::Mutpro::Preprocess::Cosmic;
use TGI::Mutpro::Preprocess::Prior;

our @ISA = qw( TGI::Mutpro::Preprocess::Calroi TGI::Mutpro::Preprocess::Statis TGI::Mutpro::Preprocess::Anno TGI::Mutpro::Preprocess::Trans TGI::Mutpro::Preprocess::Cosmic TGI::Mutpro::Preprocess::Prior );
my $MINDISTANCE = "minDistance";
my $AVGDISTANCE = "averageDistance";
my $CALROI = "calroi";
my $STATIS = "statis";
my $ANNO = "anno";
my $TRANS = "trans";
my $COSMIC = "cosmic";
my $PRIOR = "prior";

sub new {
	my $class = shift;
	my $this = {};
	$this->{'_OUTPUT_DIR'} = getcwd;
	$this->{'_BLAT'} = "blat";
	$this->{'max_3d_dis'} = 100;
	$this->{'min_seq_dis'} = 0;
	$this->{'start'} = $CALROI; 
	bless $this, $class;
	$this->process();
	return $this;
}

sub process {
	my $this = shift;
	my ($help, $options);
	unless (@ARGV) { die $this->help_text(); }
	$options = GetOptions (
		'output-dir=s'	=> \$this->{'_OUTPUT_DIR'},
		'blat=s'   => \$this->{'_BLAT'},
		'3d-distance-cutoff=i'	=> \$this->{'max_3d_dis'},
		'p-value-cutoff=i'	=> \$this->{'p_value_cutoff'},
		'linear-distance-cutoff=i'   => \$this->{'min_seq_dis'},
		'start=s' => \$this->{'start'} ,
		'help' => \$help,
	);
	if ($help) { print STDERR help_text(); exit 0; }
	unless ($options) { die $this->help_text(); }

	$this->steps();
	return;
}

sub steps {
	my $this = shift;
	if ( $this->{'start'} eq $CALROI ) {
		$this->SUPER::TGI::Mutpro::Preprocess::Calroi::process();
	} elsif ( $this->{'start'} eq $STATIS ) {
		$this->SUPER::TGI::Mutpro::Preprocess::Statis::process();
	} elsif ( $this->{'start'} eq $ANNO ) {
		$this->SUPER::TGI::Mutpro::Preprocess::Anno::process();
	} elsif ( $this->{'start'} eq $TRANS ) {
		$this->SUPER::TGI::Mutpro::Preprocess::Trans::process();
	} elsif ( $this->{'start'} eq $COSMIC ) {
		$this->SUPER::TGI::Mutpro::Preprocess::Cosmic::process();
	} elsif ( $this->{'start'} eq $PRIOR ) {
		$this->SUPER::TGI::Mutpro::Preprocess::Prior::process();
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

