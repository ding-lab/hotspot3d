package TGI::Mutpro::Pair;
#
#----------------------------------
# $Authors: Adam Scott
# $Date: 2016-11-01
# $Revision:  $
# $URL: $
# $Doc: $ pair class
# 
#----------------------------------
#
use strict;
use warnings;

use Carp;
use Getopt::Long;

use List::MoreUtils qw( uniq );
use List::Util qw( min max );

use IO::File;
use FileHandle;

use Data::Dumper;

use TGI::Variant;

sub new {
    my $class = shift;
    my $this = {};
    $this->{'variant1'} = "";
    $this->{'variant2'} = "";
    $this->{'distances'} = [];
    bless $this, $class;
    return $this;
}

sub reset {
	my $this = shift;
	my $new = new TGI::Mutpro::Pair;
	$new->set();
	return $new;
}

sub print {
	my $this = shift;
	my $delim = "\t";
	if ( @_ ) { $delim = shift; }
	foreach my $distance ( @{$this->distances()} ) {
		($this->variant1())->print();
		my $pv1s = ($this->variant1())->proteinVariants();
		my $p1s = "";
		foreach my $pv1 ( @{$pv1s} ) {
			$p1s .= $pv1->transcript().":".$pv1->aminoAcidChange();
		}
		print $delim;
		($this->variant2())->print();
		my $pv2s = ($this->variant2())->proteinVariants();
		my $p2s = "";
		foreach my $pv2 ( @{$pv2s} ) {
			$p2s .= $pv2->transcript().":".$pv2->aminoAcidChange();
		}
		print $delim;
		#print $distance;
		#print $delim;
	}
}

sub set {
	my $this = shift;
	if ( @_ ) {	$this->variant1( shift ); }
	if ( @_ ) { $this->variant2( shift ); }
	if ( @_ ) { $this->distances( shift ); }
	return $this;
}

sub addDistance {
	my $this = shift;
	push @{$this->{'distances'}} , shift;
	return $this;
}

sub variant1 {
	my $this = shift;
	if ( @_ ) { $this->{'variant1'} = shift; }
	return $this->{'variant1'};
}

sub variant2 {
	my $this = shift;
	if ( @_ ) { $this->{'variant2'} = shift; }
	return $this->{'variant2'};
}

sub distances {
	my $this = shift;
	if ( @_ ) { $this->{'distances'} = shift; }
	return $this->{'distances'};
}

1;
