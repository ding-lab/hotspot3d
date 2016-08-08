package TGI::Mutpro::Preprocess::Calroi;
#
#----------------------------------
# $Authors: Beifang Niu 
# $Date: 2014-01-14 14:34:50 -0500 (Tue Jan 14 14:34:50 CST 2014) $
# $Revision:  $
# $URL: $
# $Doc: $ generate Region of Interest information for each Uniprot ID 
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
use TGI::Mutpro::Preprocess::Uniprot;

sub new {
    my $class = shift;
    my $this = {};
    $this->{_OUTPUT_DIR} = getcwd;
    $this->{_STAT} = undef;
    bless $this, $class;
    $this->process();
    return $this;
}

sub process {
    my $this = shift;
    my ( $help, $options );
    unless( @ARGV ) { die $this->help_text(); };
    $options = GetOptions (
        'output-dir=s' => \$this->{_OUTPUT_DIR},
        'help' => \$help,
    );
    if ( $help ) { print STDERR help_text(); exit 0; };
    unless( $options ) { die $this->help_text(); };
    unless( $this->{_OUTPUT_DIR} ) { warn 'You must provide a output directory ! ', "\n"; die $this->help_text(); };
    unless( -e $this->{_OUTPUT_DIR} ) { warn 'output directory is not exist  ! ', "\n"; die $this->help_text(); };
    #### processing ####
    # generate region of interest information ( ROI )
    my ( $proDir, $annoDir, $hugoUniproFile, );
    $proDir = "$this->{_OUTPUT_DIR}\/proximityFiles";
    $annoDir = "$proDir\/annotationFiles";
    $hugoUniproFile = "$this->{_OUTPUT_DIR}\/hugo.uniprot.pdb.csv";
    unless( -e $annoDir ) { mkdir( $annoDir ) || die "can not make ROI annotation directory !\n"; };
    my $fhlog = new FileHandle;
    unless( $fhlog->open(">$annoDir/ERRORS.LOG") ) { die "Could not open ERROR LOG file !\n" };
    my ( $line, @entireFile, $uniprotId, %allUniprotIds, $uniprotRef, $annotationRef, $start, $stop, $key, $desc, $entry, $pdb, );
    my $fhuid = new FileHandle;
    unless( $fhuid->open("<$hugoUniproFile") ) { die "Could not open uniprot id file !\n" };
    @entireFile = <$fhuid>;
    $fhuid->close();
    foreach $line (@entireFile) { 
        chomp $line;
        ( undef, $uniprotId, $pdb, ) = split /\s+/, $line;
        # Only use Uniprot IDs with PDB structures
        next if ( $pdb eq "N/A" || $uniprotId !~ /\w+/ );
        $allUniprotIds{$uniprotId} = 1; 
    }
    foreach $uniprotId ( sort keys %allUniprotIds ) {
        $uniprotRef = TGI::Mutpro::Preprocess::Uniprot->new($uniprotId);
        defined ($uniprotRef) || die "no object for '$uniprotId'";
        print STDOUT $uniprotId."\n";
        # The annotation is a ref to array made here:
        # 'push @domains, 
        # "$key\t($dmStart, $dmStop)\t$desc'";
        $annotationRef = $uniprotRef->domainsAfterPosition(1);
        my $fhoneuid = new FileHandle;
        unless( $fhoneuid->open("> $annoDir/$uniprotId.annotation.txt") ) { die "Could not open annotation file to write !\n" };
        foreach my $annotation ( @{$annotationRef} ) {
            if ( $annotation =~ /(\w+)\s+\((\d+)\,\s+(\d+)\)\s+(.*)\.?$/ ) { 
                $key = $1; $start = $2; $stop = $3; $desc = $4;
            } else {
                print $fhlog "ERROR: Could not parse domain description for '$uniprotId': '$entry'";
                next;
            }
            if ( $start > $stop ) { 
                print $fhlog "ERROR: Start ($start) > Stop ($stop) in '$entry'\n";
                next;
            }
            print $fhoneuid "$start\t$stop\t'$key'\t'$desc'\n";
        }
        $fhoneuid->close();
    }
    $fhlog->close();
}

sub help_text {
    my $this = shift;
        return <<HELP

Usage: hotspot3d calroi [options]

                             REQUIRED
--output-dir                 Output directory of proximity files

--help                       this message

HELP

}

1;

