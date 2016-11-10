package TGI::Data::CleanNumber;
#
#----------------------------------
# $Authors: Adam Scott
# $Date: 2016-11-07
# $Revision: v0.0 $
# $URL: $
# $Doc: $ assures clean numbers
# 
#----------------------------------
#
use strict;
use warnings;

#sub new {
#	my $proto = shift;
#    my $class = ref( $proto ) || $proto;
#	my $this = {};
#	$this->{'number'} = shift;
#    bless $this, $class;
#    return $this;
#}
#
#sub nullIsZero {
#	my $this = shift;
#	if ( $this->{'number'} =~ /N\/A/ ) { return 0; }
#	return &numOnly( $this->{'number'} );
#}
#
#sub numOnly {
#	my $this = shift;
#	$this->{'number'} =~ s/\D*(\d+)\D*/$1/g;
#	return $this->{'number'};
#}
#
#1;
#
sub new {
	my $proto = shift;
	my $class = ref( $proto ) || $proto;
	my $this = {};
    bless $this, $class;
    return $this;
}

sub nullIsZero {
	my $num = shift;
	if ( $num =~ /N\/A/ ) { return 0; }
	return &numOnly( $num );
}

sub numOnly {
	my $num = shift;
	$num =~ s/\D*(\d+)\D*/$1/g;
	return $num;
}

1;
