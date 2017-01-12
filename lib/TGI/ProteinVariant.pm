package TGI::ProteinVariant;
#
#----------------------------------
# $Authors: Adam Scott
# $Date: 2016-11-01
# $Revision:  $
# $URL: $
# $Doc: $ protein variant class
# 
#----------------------------------
#
use strict;
use warnings;

use Carp;
use Getopt::Long;

use IO::File;
use FileHandle;

use Data::Dumper;

sub new {
    my $class = shift;
    my $this = {};
    $this->{'transcript'} = "";
    $this->{'amino_acid_change'} = "";
    #$this->{'domain'} = undef;
    bless $this, $class;
    return $this;
}

sub reset {
    my $this = shift;
    $this->enst( "" );
    $this->aminoAcidChange( "" );
    return $this;
}

sub print {
	my $this = shift;
	my $delim = ":";
	if ( @_ ) { $delim = shift; }
	print join( $delim , ( $this->transcript() , $this->aminoAcidChange() ) );
}

sub set {
	my $this = shift;
	$this->transcript( shift );
	$this->aminoAcidChange( shift );
	return $this;
}

sub transcript {
	my $this = shift;
	if ( @_ ) { $this->{'transcript'} = shift; }
	return $this->{'transcript'};
}

sub aminoAcidChange {
	my $this = shift;
	if ( @_ ) { $this->{'amino_acid_change'} = shift; }
	return $this->{'amino_acid_change'};
}

1;
