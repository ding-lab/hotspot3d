package TGI::Mutpro::Preprocess::Anno;
#
#----------------------------------
# $Authors: Beifang Niu & Adam D Scott
# $Date: 2014-01-14 14:34:50 -0500 (Tue Jan 14 14:34:50 CST 2014) $
# $Revision: 1 $
# $URL: $
# $Doc: $ add annotation information for pairs  
#----------------------------------
#

use strict;
use warnings;

use Carp;
use Cwd;
use FileHandle;
use Getopt::Long;

use TGI::Data::CleanNumber;

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
	$this->setOptions( );
    #### processing ####
    # add ROI annotations after get pvalues
    ## do that after get pvalues
	my $entireFile = $this->getInputFile( );
	my $proximityDir = $this->getProximityDir( );
	my $pvaluesDir = $this->getPValuesDir( $proximityDir );
	my ( $annotationFileDir , $annotationDir ) = $this->getAnnotationDirs( $proximityDir );
	$this->makeAnnotations( $entireFile , $annotationFileDir , $annotationDir , $pvaluesDir , $proximityDir );
	return 0;
}

sub makeAnnotations {
	my ( $this , $entireFile , $annotationFileDir , $annotationDir , $pvaluesDir , $proximityDir ) = @_;
    my $annotations = {};
    foreach my $line ( @{$entireFile} ) {
        chomp $line;
		my ( $uniprotId , $pdb );
        ( undef, $uniprotId, $pdb, ) = split /\t/, $line;
        # Only use Uniprot IDs with PDB structures
        next if ( $pdb eq "N/A" || $uniprotId !~ /\w+/ );
		$this->getAnnotation( $uniprotId , $annotationFileDir , $annotations );
	}
	foreach my $line ( @{$entireFile} ) {
		chomp( $line );
		my ( $uniprotId , $pdb );
		( undef , $uniprotId , $pdb ) = split /\t/ , $line;
		print STDOUT $uniprotId."\n";
		next if ( $pdb eq "N/A" || $uniprotId !~ /\w+/ );
        $this->addAnnotation( $uniprotId , $pvaluesDir , $annotationDir , $annotations );
	}
	return;
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
    unless( $this->{_OUTPUT_DIR} ) { warn 'HotSpot3D::Anno::setOptions error: You must provide a output directory!', "\n"; die $this->help_text(); };
    unless( -e $this->{_OUTPUT_DIR} ) { warn 'HotSpot3D::Anno::setOptions error: output directory does not exist!', "\n"; die $this->help_text(); };
	return;
}

sub getPValuesDir {
	my ( $this , $proximityDir ) = @_;
    my $pvaluesDir = "$proximityDir\/pvalues";
    unless( -d $pvaluesDir ) { warn "HotSpot3D::Anno::getPValuesDir error: You must provide a valid p_values annotation directory!\n"; die help_text(); }
	return $pvaluesDir;
}

sub getProximityDir {
	my $this = shift;
    return "$this->{_OUTPUT_DIR}\/proximityFiles";
}

sub getAnnotationDirs {
	my ( $this , $proximityDir ) = @_;
    my $annotationFileDir = "$proximityDir\/annotationFiles";
    my $annotationDir = "$proximityDir\/annotations";
    unless( -d $annotationFileDir ) { warn "HotSpot3D::Anno::getAnnotationDirs error: You must provide a valid annotation file directory!\n"; die help_text(); }
    unless( -e $annotationDir ) { mkdir( $annotationDir ) || die "HotSpot3D::Anno::getAnnotationDirs error: can not make annotations directory!\n"; };
	return ( $annotationFileDir , $annotationDir );
}

sub getInputFile {
	my $this = shift;
    my $fhuid = new FileHandle;
    my $hugoUniprotf = "$this->{_OUTPUT_DIR}\/hugo.uniprot.pdb.csv";
    unless( $fhuid->open("< $hugoUniprotf") ) { die "HotSpot3D::Anno::getInputFile error: Could not open hugo uniprot id file!\n" };
    my @entireFile = $fhuid->getlines;
    $fhuid->close();
	return \@entireFile;
}

# get annotation information
sub getAnnotation {
    my ( $this , $uniprotId , $annotationFileDir , $annotations ) = @_;
	my $annotationFile = "$annotationFileDir\/$uniprotId\.annotation\.txt";
    my $fhano = new FileHandle;
    unless( $fhano->open("< $annotationFile") ) { die "HotSpot3D::Anno::getAnnotation error: Could not open annotation file for ".$uniprotId."!\n" };
	print STDOUT $uniprotId." HotSpot3D::Anno::getAnnotation - collecting annotations from ".$annotationFile."\n";
	my $annotation = {};
    while ( my $a = $fhano->getline ) {
        chomp($a);
        my ( $start, $end, $type, $anno, ) = split /\t/, $a;
        $type =~ s/'//g;
        #warn $type."\n";
        #$anno =~ s/^'|'$|.$//;
        $anno =~ s/^'//; $anno =~ s/'$//; $anno =~ s/.$//;
        #print $anno."\n";
        if ( $type eq "DISULFID" ) {
            $annotation->{$uniprotId}->{$start} = $anno;
            $annotation->{$uniprotId}->{$end} = $anno;
        } else {
			foreach my $b ($start..$end) {
				$annotation->{$uniprotId}->{$b} = $anno;
			}
		}
    }
    $fhano->close();
	return;
}

# Add ROI annotation
sub addAnnotation {
    my ( $this , $uniprotId , $pvaluesDir , $annotationDir , $annotations ) = @_;
	my $proximityFile = "$pvaluesDir\/$uniprotId\.ProximityFile\.csv";
	if ( not -e $proximityFile ) {
		warn $uniprotId." HotSpot3D::Anno::addAnnotation warning: skipping because the annotation file does not exist: ".$proximityFile."\n";
		return;
	}
	my $outputFile = "$annotationDir\/$uniprotId\.ProximityFile\.csv";
    # add annotation information
    my $fhin = new FileHandle;
    unless( $fhin->open("<$proximityFile") ) { die "HotSpot3D::Anno::addAnnotation error: Could not open proximity file: ".$proximityFile."\n" };
    my $fhout = new FileHandle;
    unless( $fhout->open(">$outputFile") ) { die "HotSpot3D::Anno::addAnnotation error: Could not open the file for annotation: ".$outputFile."\n" };
	print STDOUT $uniprotId." HotSpot3D::Anno::addAnnotation - writing feature annotated file: ".$outputFile."\n";
	my ( $coord1 , $coord2 , $offset1 , $offset2 );
	$fhout->print( "UniProt_ID1\tChain1\tPosition1\tOffset1\tResidue_Name1\tDomain1\t" );
	$fhout->print( "UniProt_ID2\tChain2\tPosition2\tOffset2\tResidue_Name2\tDomain2\t" );
	$fhout->print( "Distance\tPDB_ID\tP_Value\n" );
    while ( my $a = $fhin->getline ) {
        next if ($a =~ /^WARNING:/);
        next if ($a =~ /UniProt_ID1/);
        chomp($a);
        my @t = split /\t/, $a;
        next if ($t[0] !~ /^\w+$/);
        next if ($t[5] !~ /^\w+$/);
        next if ($t[1] !~ /^\[[A-Z]\]$/);
        next if ($t[6] !~ /^\[[A-Z]\]$/);
        my $distance = $t[10];
        if ( $distance !~ /^-?\d+\.?\d*$/ ) { warn "HotSpot3D::Anno::addAnnotation warning: Wrong distance : $distance \n"; next; }
        my ( $annoOneEnd, $annoTwoEnd, $uniprotCoorOneEnd, $uniprotCoorTwoEnd, );
        $annoOneEnd = $annoTwoEnd = "N\/A";
		$t[2] = TGI::Data::CleanNumber::nullIsZero( $t[2] );
		$t[3] = TGI::Data::CleanNumber::nullIsZero( $t[3] );
		$t[7] = TGI::Data::CleanNumber::nullIsZero( $t[7] );
		$t[8] = TGI::Data::CleanNumber::nullIsZero( $t[8] );
        $uniprotCoorOneEnd = $t[2] + $t[3];
        $uniprotCoorTwoEnd = $t[7] + $t[8];
        #warn $uniprotCoorOneEnd."\t".$uniprotCoorTwoEnd."\n";
        if ( defined $annotations->{$uniprotId}->{$uniprotCoorOneEnd} ) { $annoOneEnd = $annotations->{$uniprotId}->{$uniprotCoorOneEnd}; }
        if ( defined $annotations->{$uniprotId}->{$uniprotCoorTwoEnd} ) { $annoTwoEnd = $annotations->{$uniprotId}->{$uniprotCoorTwoEnd}; }
        # warn $annoOneEnd."\t".$annoTwoEnd."\n";
        foreach my $d (0..4) { print $fhout $t[$d]."\t"; }
        print $fhout $annoOneEnd."\t";
        foreach my $d (5..9) { print $fhout $t[$d]."\t"; }
        print $fhout $annoTwoEnd."\t";
        foreach my $d (10..11) { print $fhout $t[$d]."\t"; }
        print $fhout $t[12]."\n";
    }
    $fhin->close();
    $fhout->close();
	return;
}

sub help_text {
    my $this = shift;
        return <<HELP

Usage: hotspot3d anno [options]

                             REQUIRED
--output-dir                 Output directory of proximity files

--help                       this message

HELP

}

1;

