package TGI::Mutpro::Preprocess::Drugport;
#
#----------------------------------
# $Authors: Beifang Niu 
# $Date: 2014-01-14 14:34:50 -0500 (Tue Jan 14 14:34:50 CST 2014) $
# $Revision:  $
# $URL: $
# $Doc: $ drugport database processing module
#----------------------------------
#
use strict;
use warnings;

use Carp;
use Cwd;
use Getopt::Long;
use LWP::Simple;
use IO::File;
use FileHandle;

sub new {
    my $class = shift;
    my $this = {};
    $this->{'output_file'} = 'drugport_parsing_output';
    $this->{'pdb_file_dir'} = undef; 
    $this->{'stat'} = undef;
    $this->{'page'} = "";
    bless $this, $class;
    $this->process();
    return $this;
}

sub process {
    my $this = shift;
    my ( $help, $options );
    unless( @ARGV ) { die $this->help_text(); };
    $options = GetOptions (
        'pdb-file-dir'  => \$this->{'pdb_file_dir'},
        'output-file=s' => \$this->{'output_file'},
        'help' => \$help,
    );
    if ( $help ) { print STDERR help_text(); exit 0; };
    unless( $options ) { die $this->help_text(); };
    unless( $this->{'output_file'} ) { warn 'You must provide a output file for drugport database ! ', "\n"; die $this->help_text(); };
    unless( $this->{'pdb_file_dir'} and (-e $this->{'pdb_file_dir'})) { warn " $_ is not exist ! \n"; die $this->help_text(); };
    #### processing ####
    # parse drugport database
    #
    # drug name and ids
    #
    my $appdrugs_url = "http://www.ebi.ac.uk/thornton-srv/databases/drugport/data/appdrugs_pdb.dat";
    $this->{'page'} = get( $appdrugs_url );
    unless( $this->{'page'} ) { die "can not access drugport database file $! \n"; }
    map{ @_ = split / /; if ( /^GENERIC_NAME/ ) { print $_[1]."\t" }; if ( /^DRUGNAME_ID/ ) { print $_[1]."\n" } } split /\n/, $this->{'page'};














}

sub help_text {
    my $this = shift;
        return <<HELP

Usage: hotspot3d drugport [options]

--pdb-file-dir          PDB file directory 
--output-file		Output file of drugport parsing

--help			this message

HELP

}

1;

