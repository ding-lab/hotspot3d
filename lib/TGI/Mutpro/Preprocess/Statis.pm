package TGI::Mutpro::Preprocess::Statis;
#
#----------------------------------
# $Original authors: Beifang Niu & Adam D Scott
# #Modified by: Fernanda Martins Rodrigues @WashU (fernanda@wustl.edu; mrodrigues.fernanda@gmail.com)
# $Date: 2014-01-14 14:34:50 -0500 (Tue Jan 14 14:34:50 CST 2014) $
# $Revision: 2023-03-15 $
# $URL: $
# $Doc: $ statics related infor 
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
    $this->{_OUTPUT_DIR} = getcwd;
    $this->{_STAT} = undef;
    bless $this, $class;
    $this->process();
    return $this;
}

sub process {
    my $this = shift;
	$this->setOptions();
    #### processing ####
    # pvalue calculation program
	my $fh = $this->getInputFile( );
	my $proDir = $this->getInputDir( );
	my $pvaluesDir = $this->getOutputDir( $proDir );
	$this->calculatePValues( $pvaluesDir , $fh , $proDir );
	return 0;
}

sub setOptions {
	my $this = shift;
    my ( $help, $options );
    unless( @ARGV ) { die $this->help_text(); };
    $options = GetOptions (
        'output-dir=s' => \$this->{_OUTPUT_DIR},
        'help' => \$help,
    );
    if ( $help ) { warn help_text(); exit 0; };
    unless( $options ) { die $this->help_text(); };
    unless( $this->{_OUTPUT_DIR} ) { warn 'HotSpot3D::Statis::setOptions error: You must provide a output directory!', "\n"; die $this->help_text(); };
    unless( -e $this->{_OUTPUT_DIR} ) { warn 'HotSpot3D::Statis::setOptions error: The output directory does not exist!', "\n"; die $this->help_text(); };
	return;
}

sub getOutputDir {
	my ( $this , $proDir ) = @_;
    my $pvaluesDir = "$proDir/pvalues";
    unless( -e $pvaluesDir ) { mkdir($pvaluesDir) || die "HotSpot3D::Statis::getOutputDir error: can not make pvalues directory!\n"; };
	return $pvaluesDir;
}

sub getInputFile {
	my $this = shift;
    my $hugoUniprot = "$this->{_OUTPUT_DIR}/hugo.uniprot.alphafolddb.csv";
    unless( -e $hugoUniprot ) { die "HotSpot3D::Statis::getInputFile error: no hugo uniprot file!\n"; }; 
    my $fh = new FileHandle;
    unless( $fh->open("<$hugoUniprot") ) { die "HotSpot3D::Statis::getInputFile error: Could not open hugo uniprot file!\n" };
	return $fh;
}

sub getInputDir {
	my $this = shift;
    my $proDir = "$this->{_OUTPUT_DIR}/proximityFiles";
    unless( -e $proDir ) { die "HotSpot3D::Statis::getInputDir error: no proximity file directory!\n"; };
	return $proDir;
}

sub calculatePValues {
	my ( $this , $pvaluesDir , $fh , $proDir ) = @_;
    my @entireFile = <$fh>;
    $fh->close();
    my $u = 0;
    foreach my $line (@entireFile) {
        chomp $line;
		my ( $uniprotId, $alphafolddb, );
        (undef, $uniprotId, $alphafolddb) = split /\t/, $line;
        # Only use Uniprot IDs with AlphaFold DB structures
        next if ( $alphafolddb eq "N/A" || $uniprotId !~ /\w+/ );
		$this->calculatePValuesOfProtein( $uniprotId , $pvaluesDir , $proDir );
	}
	return;
}

sub calculatePValuesOfProtein {
	my ( $this , $uniprotId , $pvaluesDir , $proDir ) = @_;
	# proximity file
	my $proximityFile = "$proDir\/$uniprotId\.ProximityFile\.csv";
	next unless(-e $proximityFile);
	my $outputFile = "$pvaluesDir\/$uniprotId\.ProximityFile\.csv";
	# p_value calculating
	my $numberlines = $this->getPvalue( $uniprotId , $proximityFile , $outputFile );
	#delete file if null
	if ($numberlines == 0) { unlink( $outputFile ) or warn "HotSpot3D::Statis::calculatePValuesOfProtein warning: failed to delete $outputFile: $!"; }
	return;
}

# pvalue calculating
sub getPvalue {
    my ( $this , $uniprotId , $proximityfile , $outputf ) = @_;
    my $fh = new FileHandle;
    unless( $fh->open("<$proximityfile") ) { die "HotSpot3D::Statis::getPvalue error: Could not open hugo uniprot file '$proximityfile'!\n" };
    my %distances;
    # get distances list
    while (my $a = <$fh>) {
        next if ($a =~ /^WARNING:/);
        next if ($a =~ /UniProt_ID1/);
        chomp($a);
        my @t = split /\t/, $a;
        next if ($t[0] !~ /^\w+$/);
        next if ($t[5] !~ /^\w+$/);
        next if ($t[1] !~ /^\[[A-Z]\]$/);
        next if ($t[6] !~ /^\[[A-Z]\]$/);
        my $distance = $t[10];
        if ( $distance !~ /^-?\d+\.?\d*$/ ) {
            warn "Wrong distance : $distance \n";
            next;
        }
        $distances{$distance} = 1;
    }
    $fh->close();
    # sort and calculate p_value
    my ( %pvalues, $i, @t, $total, );
    @t = sort {$a<=>$b} keys %distances;
    $total = scalar( @t );
    $i = 0;
    map{ $pvalues{$_} = $i/$total; $i++; } @t;
    undef @t;
    undef %distances;

    $fh = new FileHandle;
    unless( $fh->open("<$proximityfile") ) { die "HotSpot3D::Statis::getPvalue error: Could not open hugo uniprot file '$proximityfile'!\n" };
    my $fho = new FileHandle;
    unless( $fho->open(">$outputf") ) { die "HotSpot3D::Statis::getPvalue error: Could not open file '$outputf' to write!\n" };
	print STDOUT $uniprotId." HotSpot3D::Statis::getPvalue - making p-value annotated output: ".$outputf."\n";
    my $numberlines = 0;
    # load p_values
	$fho->print( "UniProt_ID1\tChain1\tPosition1\tOffset1\tResidue_Name1\t" );
	$fho->print( "UniProt_ID2\tChain2\tPosition2\tOffset2\tResidue_Name2\t" );
	$fho->print( "Distance\tAlphafoldDB_ID\tP_Value\n" );
    while ( my $a = <$fh> ) {
        next if ($a =~ /^WARNING:/);
        next if ($a =~ /UniProt_ID1/);
        chomp($a);
        my @t = split /\t/, $a;
        next if ( $t[0] !~ /^\w+$/ );
        next if ( $t[5] !~ /^\w+$/ );
        next if ( $t[1] !~ /^\[[A-Z]\]$/ );
        next if ( $t[6] !~ /^\[[A-Z]\]$/ );
        my $distance = $t[10];
        if ( $distance !~ /^-?\d+\.?\d*$/ ) {
            warn "Wrong distance : $distance \n";
            next;
        }
        my $rounded = sprintf( "%.6f", $pvalues{$distance} );
        print $fho $a."\t".$rounded."\n";
        $numberlines++;
        undef @t;
    }
    $fh->close();
    $fho->close();
    # clear
    undef %pvalues;
    return $numberlines;
}

sub help_text{
    my $this = shift;
        return <<HELP

Usage: hotspot3d statis [options]

                             REQUIRED
--output-dir                 Output directory of proximity files

--help                       this message

HELP

}

1;

