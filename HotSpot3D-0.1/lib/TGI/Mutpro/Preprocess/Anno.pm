package TGI::Mutpro::Preprocess::Anno;
#
#----------------------------------
# $Authors: Beifang Niu 
# $Date: 2014-01-14 14:34:50 -0500 (Tue Jan 14 14:34:50 CST 2014) $
# $Revision:  $
# $URL: $
# $Doc: $ add annotation information for pairs  
#----------------------------------
#

use strict;
use warnings;

use Carp;
use Cwd;
use Getopt::Long;

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
    # add ROI annotations after get pvalues
    ## do that after get pvalues
    my ( $hugoUniprotf, $proximityDir, $pvaluesDir, $annotationFileDir, $annotationsDir, );
    $hugoUniprotf = "$this->{_OUTPUT_DIR}\/hugo.uniprot.pdb.csv";
    $proximityDir = "$this->{_OUTPUT_DIR}\/proximityFiles";
    $pvaluesDir = "$proximityDir\/pvalues";
    $annotationFileDir = "$proximityDir\/annotationFiles";
    $annotationsDir = "$proximityDir\/annotations";
    unless( -d $pvaluesDir ) { warn "You must provide a valid p_values annotation directory ! \n"; die help_text(); }
    unless( -d $annotationFileDir ) { warn "You must provide a valid annotation file directory ! \n"; die help_text(); }
    unless( -e $annotationsDir ) { mkdir( $annotationsDir ) || die "can not make annotations directory !\n"; };
    my ( $uniprotId, $pdb, @entireFile, );
    my $fhuid = new FileHandle;
    unless( $fhuid->open("< $hugoUniprotf") ) { die "Could not open hugo uniprot id file !\n" };
    @entireFile = $fhuid->getlines;
    $fhuid->close();
    # get annotation information
    my %annotations;
    foreach my $line ( @entireFile ) {
        chomp $line;
        ( undef, $uniprotId, $pdb, ) = split /\t/, $line;
        # Only use Uniprot IDs with PDB structures
        next if ( $pdb eq "N/A" || $uniprotId !~ /\w+/ );
        print STDERR $uniprotId."\n";
        my $annotationFile = "$annotationFileDir\/$uniprotId\.annotation\.txt";
        # get annotation infor
        $this->getAnnotation($uniprotId, $annotationFile, \%annotations);
    }
    # generate new proximity files
    my $u = 0;
    foreach my $line ( @entireFile ) {
        chomp $line;
        ( undef, $uniprotId, $pdb, ) = split /\t/, $line;
        # Only use Uniprot IDs with PDB structures
        next if ( $pdb eq "N/A" || $uniprotId !~ /\w+/ );
        print STDERR $uniprotId."\n";
        # proximity file
        my $proximityFile = "$pvaluesDir\/$uniprotId\.ProximityFile\.csv";
        next unless( -e $proximityFile );
        my $outputFile = "$annotationsDir\/$uniprotId\.ProximityFile\.csv";
        # add annotation infor
        $this->addAnnotation( $proximityFile, \%annotations, $outputFile, $uniprotId );
        #delete file if null
   }
}

# get annotation information
sub getAnnotation {
    my ( $this, $uniprot, $annotationf, $annoref, ) = @_;
    my $fhano = new FileHandle;
    unless( $fhano->open("< $annotationf") ) { die "Could not open annotation file !\n" };
    while ( my $a = $fhano->getline ) {
        chomp($a);
        my ( $start, $end, $type, $anno, ) = split /\t/, $a;
        $type =~ s/'//g;
        #print STDERR $type."\n";
        #$anno =~ s/^'|'$|.$//;
        $anno =~ s/^'//; $anno =~ s/'$//; $anno =~ s/.$//;
        #print $anno."\n";
        if ( $type eq "DISULFID" ) {
            $annoref->{$uniprot}->{$start} = $anno;
            $annoref->{$uniprot}->{$end} = $anno;
        } else { foreach my $b ($start..$end) { $annoref->{$uniprot}->{$b} = $anno; } }
    }
    $fhano->close();
}

# Add ROI annotation
sub addAnnotation {
    my ( $this, $proximityfile, $annotationRef, $outputf, $uniprotId, ) = @_;
    # add annotation information
    my $fhin = new FileHandle;
    unless( $fhin->open("<$proximityfile") ) { die "Could not open proximity file !\n" };
    my $fhout = new FileHandle;
    unless( $fhout->open(">$outputf") ) { die "Could not open proximity file to add annotation information !\n" };
    while ( my $a = $fhin->getline ) {
        next if ($a =~ /^WARNING:/);
        chomp($a);
        my @t = split /\t/, $a;
        next if ($t[0] !~ /^\w+$/);
        next if ($t[5] !~ /^\w+$/);
        next if ($t[1] !~ /^\[[A-Z]\]$/);
        next if ($t[6] !~ /^\[[A-Z]\]$/);
        my $distance = $t[10];
        if ( $distance !~ /^-?\d+\.?\d*$/ ) { print "Wrong distance : $distance \n"; next; }
        my ( $annoOneEnd, $annoTwoEnd, $uniprotCoorOneEnd, $uniprotCoorTwoEnd, );
        $annoOneEnd = $annoTwoEnd = "N\/A";
        $uniprotCoorOneEnd = $t[2] + $t[3];
        $uniprotCoorTwoEnd = $t[7] + $t[8];
        #print STDERR $uniprotCoorOneEnd."\t".$uniprotCoorTwoEnd."\n";
        if ( defined $annotationRef->{$uniprotId}->{$uniprotCoorOneEnd} ) { $annoOneEnd = $annotationRef->{$uniprotId}->{$uniprotCoorOneEnd}; }
        if ( defined $annotationRef->{$uniprotId}->{$uniprotCoorTwoEnd} ) { $annoTwoEnd = $annotationRef->{$uniprotId}->{$uniprotCoorTwoEnd}; }
        # print STDERR $annoOneEnd."\t".$annoTwoEnd."\n";
        foreach my $d (0..4) { print $fhout $t[$d]."\t"; }
        print $fhout $annoOneEnd."\t";
        foreach my $d (5..9) { print $fhout $t[$d]."\t"; }
        print $fhout $annoTwoEnd."\t";
        foreach my $d (10..11) { print $fhout $t[$d]."\t"; }
        print $fhout $t[12]."\n";
    }
    $fhin->close();
    $fhout->close();
}

sub help_text {
    my $this = shift;
        return <<HELP

Usage: hotspot3d anno [options]

--output-dir		Output directory of proximity files

--help			this message

HELP

}

1;

