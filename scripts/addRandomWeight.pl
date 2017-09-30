#!/usr/bin/perl
#25 September 2017 - Adam D Scott - 

use strict;
use warnings;

use IO::File;
use FileHandle;

my $usage = 'perl addRandomWeight.pl <maf> <output> 
';

die $usage , unless @ARGV == 2;
my ( $maf , $output ) = @ARGV;

my $IN1 = FileHandle->new( $maf , "r" );
if ( not defined $IN1 ) { die "ADSERROR: Could not open/read ".$maf."\n"; }

my $IN2 = FileHandle->new( $output , "w" );
if ( not defined $IN2 ) { die "ADSERROR: Could not open/write ".$output."\n"; }

while ( my $line = <$IN1> ) {
	chomp( $line );
	$IN2->print( $line."\t" );
	if ( $line =~ /Hugo/ ) {
		$IN2->print( "RandomWeight" );
	} else {
		my $r = int( rand( 40 ) );
		$r -= 20;
		$IN2->print( $r )
	}
	$IN2->print( "\n" );
}
$IN1->close();
$IN2->close();
