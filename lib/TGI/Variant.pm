package TGI::Variant;
#
#----------------------------------
# $Authors: Adam Scott
# $Date: 2016-11-01
# $Revision:  $
# $URL: $
# $Doc: $ variant class
# 
#----------------------------------
#
use strict;
use warnings;

use Data::Dumper;

use TGI::ProteinVariant;

sub new {
    my $class = shift;
    my $this = {};
    $this->{'gene'} = "";
    $this->{'chromosome'} = "";
    $this->{'start'} = "";
    $this->{'stop'} = "";
    $this->{'reference'} = "";
    $this->{'alternate'} = "";
    $this->{'proteinVariants'} = [];
    bless $this, $class;
    return $this;
}

sub reset {
    my $this = shift;
    $this->gene( "" );
    $this->chromosome( "" );
    $this->start( "" );
    $this->stop( "" );
    $this->reference( "" );
    $this->alternate( "" );
    $this->proteinVariants( [] );
    return $this;
}

sub print {
	my $this = shift;
	my $delim = "\t";
	if ( @_ ) { $delim = shift; }
	print $this->gene()." ".$this->hgvsg();
}

sub gene {
	my $this = shift;
	if ( @_ ) { $this->{'gene'} = shift; }
	return $this->{'gene'};
}

sub chromosome {
	my $this = shift;
	if ( @_ ) { $this->{'chromosome'} = shift };
	return $this->{'chromosome'};
}

sub start {
	my $this = shift;
	if ( @_ ) { $this->{'start'} = shift };
	return $this->{'start'};
}

sub stop {
	my $this = shift;
	if ( @_ ) { $this->{'stop'} = shift };
	return $this->{'stop'};
}

sub reference {
	my $this = shift;
	if ( @_ ) { $this->{'reference'} = shift };
	return $this->{'reference'};
}

sub alternate {
	my $this = shift;
	if ( @_ ) { $this->{'alternate'} = shift };
	return $this->{'alternate'};
}

sub hgvsg {
	my $this = shift;
	my $hgvsg = "";
	$hgvsg .= $this->chromosome().":g.";
	$hgvsg .= $this->start();
	if ( $this->stop() ) {
		$hgvsg .= "-".$this->stop();
	}
	if ( $this->reference() ) {
		$hgvsg .= $this->reference();
		if ( $this->alternate() ) {
			$hgvsg .= ">";
			$hgvsg .= $this->alternate();
		}
	}
	return $hgvsg;
}

sub set {
	my $this = shift;
    $this->gene( shift );
    $this->chromosome( shift );
    $this->start( shift );
    $this->stop( shift );
    $this->reference( shift );
    $this->alternate( shift );
	return $this;
}

sub addProteinVariant {
	my $this = shift;
	if ( @_ ) { push @{$this->{'proteinVariants'}} , shift; }
	return $this;
}

sub proteinVariant {
	my $this = shift;
	if ( @_ ) { 
		my $i = shift;
		if ( scalar @{$this->{'proteinVariants'}} > $i ) {
			return $this->{'proteinVariants'}->[$i];
		}
	}
	return $this->{'proteinVariants'}->[0];
}

sub proteinVariants {
	my $this = shift;
	return $this->{'proteinVariants'};
}

1;
