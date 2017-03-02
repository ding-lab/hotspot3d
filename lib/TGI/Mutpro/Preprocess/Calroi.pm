package TGI::Mutpro::Preprocess::Calroi;
#
#----------------------------------
# $Authors: Beifang Niu & Adam D Scott
# $Date: 2014-01-14 14:34:50 -0500 (Tue Jan 14 14:34:50 CST 2014) $
# $Revision: 1 $
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
	$this->setOptions();
    #### processing ####
    # generate region of interest information ( ROI )
	my $annoDir = $this->getOutputDir();
	my $fhuid = $this->getInputFile( );
	my $allUniprotIds = $this->getUniprotIds( $fhuid );
	$this->makeROIannotations( $allUniprotIds , $annoDir );
	return 0;
}

sub makeROIannotations {
	my ( $this , $allUniprotIds , $annoDir ) = @_;

    foreach my $uniprotId ( keys %{$allUniprotIds} ) {
		$this->makeROIannotationFile( $uniprotId , $annoDir );
    }
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
    unless( $this->{_OUTPUT_DIR} ) { warn 'HotSpot3D::Calroi::setOptions error: You must provide a output directory!', "\n"; die $this->help_text(); };
    unless( -e $this->{_OUTPUT_DIR} ) { warn 'HotSpot3D::Calroi::setOptions error: The output directory does not exist!', "\n"; die $this->help_text(); };
	return;
}

sub getInputFile {
	my ( $this , $annoDir ) = @_;

    my $fhuid = new FileHandle;
    my $hugoUniprotFile = "$this->{_OUTPUT_DIR}\/hugo.uniprot.pdb.csv";
    unless( $fhuid->open("<$hugoUniprotFile") ) { die "HotSpot3D::Calroi::getInputFile error: Could not open uniprot id file ".$hugoUniprotFile."!\n" };

	return $fhuid;
}

sub getOutputDir {
	my $this = shift;

    my $proDir = "$this->{_OUTPUT_DIR}\/proximityFiles";
    my $annoDir = "$proDir\/annotationFiles";
    unless( -e $annoDir ) { mkdir( $annoDir ) || die "HotSpot3D::Calroi::getOutputDir error: can not make ROI annotation directory!\n"; };

	return $annoDir;
}

sub getUniprotIds {
	my ( $this , $fhuid ) = @_;
	my $allUniprotIds;

    my ( $line , @entireFile , $uniprotId , $pdb );
    @entireFile = <$fhuid>;
    $fhuid->close();
    foreach $line (@entireFile) { 
        chomp $line;
        ( undef , $uniprotId , $pdb ) = split /\s+/ , $line;
        # Only use Uniprot IDs with PDB structures
        next if ( $pdb eq "N/A" || $uniprotId !~ /\w+/ );
        $allUniprotIds->{$uniprotId} = 1; 
    }
	return $allUniprotIds;
}

sub makeROIannotationFile {
	my ( $this , $uniprotId , $annoDir ) = @_;
	my $uniprotRef = TGI::Mutpro::Preprocess::Uniprot->new($uniprotId);
	defined ($uniprotRef) || die "HotSpot3D::Calroi::makeROIannotationFile error: no object for '$uniprotId'";
	# The annotation is a ref to array made here:
	# 'push @domains, 
	# "$key\t($dmStart, $dmStop)\t$desc'";
	my $annotationRef = $uniprotRef->domainsAfterPosition(1);
	my $file = $annoDir."/".$uniprotId.".annotation.txt";
	my $fhoneuid = new FileHandle;
	unless( $fhoneuid->open( $file , "w" ) ) {
		die "HotSpot3D::Calroi::makeROIannotationFile error: Could not open annotation file to write for ".$uniprotId." at ".$file."!\n";
	};
	print STDOUT $uniprotId." HotSpot3D::Calroi - Making annotation file ".$file."\n";
	$fhoneuid->print( "Feature_Start\tFeature_End\tFeature_Type\tFeature_Description\n" );
	foreach my $annotation ( @{$annotationRef} ) {
		my ( $key , $start , $stop , $desc );
		if ( $annotation =~ /(\w+)\s+\((\d+)\,\s+(\d+)\)\s+(.*)\.?$/ ) { 
			$key = $1; $start = $2; $stop = $3; $desc = $4;
		} else {
			warn "HotSpot3D::Calroi::makeROIannotationFile warning: Could not parse domain description for '$uniprotId'\n";
			next;
		}
		if ( $start > $stop ) { 
			warn "HotSpot3D::Calroi::makeROIannotationFile warning: Start ($start) > Stop ($stop) in '$uniprotId'\n";
			next;
		}
		$fhoneuid->print( join( "\t" , ( $start , $stop , $key , $desc ) )."\n" );
	}
	$fhoneuid->close();
	return;
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

