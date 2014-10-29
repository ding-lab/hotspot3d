package TGI::Mutpro::Preprocess::Uppro;
#
#----------------------------------
# $Authors: Beifang Niu 
# $Date: 2014-01-14 14:34:50 -0500 (Tue Jan 14 14:34:50 CST 2014) $
# $Revision:  $
# $URL: $
# $Doc: $ create & update proximity files (the first step of preprocessing procedure)
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
use TGI::Mutpro::Preprocess::HugoGeneMethods;

sub new {
    my $class = shift;
    my $this = {};
    $this->{_OUTPUT_DIR} = getcwd;
    $this->{_MAX_3D_DIS} = 10;
    $this->{_MIN_SEQ_DIS} = 5;
    $this->{_STAT} = undef;
    ## change to use local pdb file
    $this->{_PDB_FILE_DIR} = undef;
    ## add drug port database file (08.04.2014)
    $this->{_DRUGPORT_FILE} = undef;

    bless $this, $class;
    $this->process();
    return $this;
}

sub process {
    my $this = shift;
    my ( $help, $options );
    unless( @ARGV ) { die $this->help_text(); }
    $options = GetOptions (
        'maX-3d-dis=i'   => \$this->{_MAX_3D_DIS},
        'min-seq-dis=i'  => \$this->{_MIN_SEQ_DIS},
        'output-dir=s'   => \$this->{_OUTPUT_DIR},
        'pdb-file-dir=s' => \$this->{_PDB_FILE_DIR},
        'drugport-file=s' => \$this->{_DRUGPORT_FILE},

        'help' => \$help,
    );
    if ( $help ) { print STDERR help_text(); exit 0; }
    unless( $options ) { die $this->help_text(); }
    unless( $this->{_OUTPUT_DIR} ) {  warn 'You must provide a output directory ! ', "\n"; die $this->help_text(); }
    unless( -e $this->{_OUTPUT_DIR} ) { warn 'output directory is not exist  ! ', "\n"; die $this->help_text(); }
    unless( $this->{_PDB_FILE_DIR} ) {  warn 'You must provide a PDB file directory ! ', "\n"; die $this->help_text(); }
    unless( -e $this->{_PDB_FILE_DIR} ) { warn 'PDB file directory is not exist  ! ', "\n"; die $this->help_text(); }
    unless( $this->{_DRUGPORT_FILE} ) {  warn 'You must provide a drugport database file ! ', "\n"; die $this->help_text(); }
    unless( -e $this->{_DRUGPORT_FILE} ) { warn 'DrugPort databse file is not exist  ! ', "\n"; die $this->help_text(); }

    #### processing: update program ####
    my ( $updateProgram, $proDir, $inproDir, $pdbCorDir, $logFile );
    $updateProgram = 'hotspot3d calpro';
    $proDir = "$this->{_OUTPUT_DIR}\/proximityFiles";
    $inproDir = "$proDir\/inProgress";
    $pdbCorDir = "$proDir\/pdbCoordinateFiles";
    $logFile = "$this->{_OUTPUT_DIR}\/hugo.uniprot.pdb.csv";
    unless( -e $proDir ) { mkdir($proDir) || die "can not make proximity directory !\n"; }
    unless( -e $logFile ) { if ( system("touch $logFile") != 0 ) { die "can not make hugo uniprot file !\n"; } }
    unless( -e $inproDir ) { mkdir($inproDir) || die "can not make progress directory !\n"; }
    unless( -e $pdbCorDir ) { mkdir($pdbCorDir) || die "can not make progress directory !\n"; }
    my ( $UniprotToStructureRef, $UniprotFileRef, %UniprotIdToUpdate, );
    # Get ref to hash listing all PDB structures that have been analyzed
    $UniprotToStructureRef = $this->currentStructures( $logFile );
    # Get ref to hash with all Uniprot files already written to disk
    $UniprotFileRef = $this->currentUniprotFiles( $proDir );
    # Write a new file with current Hugo symbols, Uniprot IDs, PDB structure IDs, and Hugo aliases
    # Keep track of Uniprot IDs that are new or that have new PDB structures
    # These are the ones that need to be updated (%UniprotIdToUpdate)
    my $fh = new FileHandle;
    unless( $fh->open( ">$logFile" ) ) { die "Could not open hugo uniprot file !\n" };
    my ( $hugoGeneRef, $hugoId, $aliasRef, $previousRef, $aliasList, $uniprotId, $uniprotRef, $pdbRef, );
    $hugoGeneRef = TGI::Mutpro::Preprocess::HugoGeneMethods::makeHugoGeneObjects();
    foreach $hugoId ( sort keys %{$hugoGeneRef} ) {
        print STDERR 'HUGO: ', "$hugoId\n";
        $aliasRef = $$hugoGeneRef{$hugoId}->getAllAliases();
        $previousRef =  $$hugoGeneRef{$hugoId}->getAllPreviousSymbols();
        $aliasList = "";
        map{ $aliasList .= "$_ "; } keys %{$aliasRef};
        map{ $aliasList .= "$_ "; } keys %{$previousRef};
        if ( $aliasList !~ /\w+/ ) { $aliasList = "N/A"; };
        $uniprotId = $$hugoGeneRef{$hugoId}->uniprot();
        if ( !defined $uniprotId ) { $fh->print( "$hugoId\tN/A\tN/A\t$aliasList\n" ); next; }
        $uniprotRef = TGI::Mutpro::Preprocess::Uniprot->new($uniprotId);
        $pdbRef = $uniprotRef->annotations("PDB");
        if ( !defined $pdbRef || scalar(@{$pdbRef}) == 0 ) { $fh->print( "$hugoId\t$uniprotId\tN/A\t$aliasList\n" ); next; }
        # If the Uniprot ID does not have an associated file, it has to be updated
        if ( !defined $$UniprotFileRef{$uniprotId} ) { $UniprotIdToUpdate{$uniprotId} = 1; }
        # update hugo uniprot file
        $fh->print( "$hugoId\t$uniprotId\t" );
	# If the Uniprot ID has a new structure, it has to be updated
        map{ my ( $pdbId ) = $_ =~ /^(\w+)\;/; if ( defined $pdbId ){ if( !defined $$UniprotToStructureRef{$uniprotId}{$pdbId} ){ $UniprotIdToUpdate{$uniprotId} = 1; }; $fh->print( "$pdbId " ); } } @{$pdbRef};
        $fh->print( "\t$aliasList\n" );
    }
    $fh->close();
    map {
        system( "touch $inproDir/$_.ProximityFile.csv" );
        my $submit_cmd = "bsub -oo $_.err.log -R 'select[type==LINUX64 && mem>8000] rusage[mem=8000]' -M 8000000 '$updateProgram --output-dir=$this->{_OUTPUT_DIR} --pdb-file-dir=$this->{_PDB_FILE_DIR} --drugport-file=$this->{_DRUGPORT_FILE} --uniprot-id=$_ --max-3d-dis=$this->{_MAX_3D_DIS} --min-seq-dis=$this->{_MIN_SEQ_DIS}'";
        print STDERR $submit_cmd."\n"; system( "$submit_cmd" );
    } keys %UniprotIdToUpdate;

    return 1;
}

# Check which Uniprot files are on disk Return ref to 
# hash with key = Uniprot Id 
sub currentUniprotFiles {
    my ( $this, $dir ) = @_;
    my ( %uniprotIds, $file, );
    opendir( DIR, $dir ) || die "Could not open '$dir': $!";
    map{ if ( $_ =~ /(\w+)\.ProximityFile\.csv/ ){  $uniprotIds{$1} = 1; }  } ( readdir DIR );
    closedir DIR;
    return \%uniprotIds;
}

# Return ref to hash with 
# key = uniprotId, $pdbId; value = 1
sub currentStructures {
    my ( $this, $logfilef ) = @_;
    my ( $uniprotId, $pdbId, %uniprotToStructure, $pdbList );
    my $fh = new FileHandle;
    unless( $fh->open($logfilef) ) { die "Could not open hugo uniprot file\n" };
    map{ chomp; ( undef, $uniprotId, $pdbList ) = split /\t/, $_; unless( $uniprotId eq "N/A" ) { map { $uniprotToStructure{$uniprotId}{$_} = 1; } split /s+/, $pdbList; } } $fh->getlines;
    $fh->close();
    return \%uniprotToStructure;
}

sub help_text{
    my $this = shift;
        return <<HELP

Usage: hotspot3d uppro [options]

--output-dir		Output directory of proximity files
--pdb-file-dir          PDB file directory 
--drugport-file         DrugPort database file    

--max-3d-dis            Maximum 3D distance in angstroms befor two amino acids
                        are considered 'close', default = 10
--min-seq-dis           Minimum linear distance in primary sequence. If two amino acids are <= 5 positions 
                        apart in the primary sequence, don't record them, default = 5 

--help			this message

HELP

}

1;

