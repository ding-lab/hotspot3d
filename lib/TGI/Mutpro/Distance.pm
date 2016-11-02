package TGI::Mutpro::Distance;
#
#----------------------------------
# $Authors: Adam Scott
# $Date: 2016-11-01
# $Revision:  $
# $URL: $
# $Doc: $ distance class
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

my $ABSURD = 10000;

sub new {
    my $class = shift;
    my $this = {};
    $this->{'pdb'} = "";
    $this->{'chain1'} = "";
    $this->{'chain2'} = "";
    $this->{'position1'} = "";
    $this->{'position2'} = "";
    $this->{'distance'} = $ABSURD;
    bless $this, $class;
    return $this;
}

sub reset {
    my $this = shift;
    $this->pdb( "" );
    $this->chain1( "" );
    $this->chain2( "" );
    $this->position1( "" );
    $this->position2( "" );
    $this->distance( $ABSURD );
    return $this;
}

sub print {
	my $this = shift;
	my $delim = "\t";
	if ( @_ ) { $delim = shift; }
	print join( $delim , ( $this->chain1().":".$this->position1() , $this->chain2().":".$this->position2() , $this->pdb() , $this->distance() ) );
}

sub set {
	my $this = shift;
	$this->pdb( shift );
	$this->chain1( shift );
	$this->chain2( shift );
	$this->position1( shift );
	$this->position2( shift );
	$this->distance( shift );
	return $this;
}

sub pdb {
	my $this = shift;
	if ( @_ ) { $this->{'pdb'} = shift; }
	return $this->{'pdb'};
}

sub chain1 {
	my $this = shift;
	if ( @_ ) { $this->{'chain1'} = shift; }
	return $this->{'chain1'};
}

sub chain2 {
	my $this = shift;
	if ( @_ ) { $this->{'chain2'} = shift; }
	return $this->{'chain2'};
}

sub position1 {
	my $this = shift;
	if ( @_ ) { $this->{'position1'} = shift; }
	return $this->{'position1'};
}

sub position2 {
	my $this = shift;
	if ( @_ ) { $this->{'position2'} = shift; }
	return $this->{'position2'};
}

sub distance {
	my $this = shift;
	if ( @_ ) { $this->{'distance'} = shift; }
	return $this->{'distance'};
}

1;
