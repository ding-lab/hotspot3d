package TGI::Files::List;
#
#----------------------------------
# $Authors: Adam Scott
# $Date: 2016-10-26
# $Revision: v0.0 $
# $URL: $
# $Doc: $ file handling for .mafs
# 
#----------------------------------
#
use strict;
use warnings;

use Carp;

use Scalar::Util qw( openhandle );
use TGI::Files::File;
our @ISA = qw( TGI::Files::File );

sub new {
    my $class = shift;
    my $this = $class->SUPER::new( shift );
    $this->{'items'} = 0;
    bless $this, $class;
    return $this;
}

sub getList {
	my $this = shift;
	my $column = 0;
	if ( @_ ) { $column = shift; }
	print STDOUT "\nReading in ".$this->{'file_name'}."...\n";
	if ( not defined openhandle( $this->{'handle'} ) ) {
		$this->open();
	}
	seek( $this->{'handle'} , 0 , 0 );
	my %items;
	map {
		chomp;
		print $_."\n";
		my $item = (split /\t/)[$column];
		print $item."\n";
		$items{$item} += 1;
	} $this->getlines();
	$this->close();

	return \%items;
}

1;
