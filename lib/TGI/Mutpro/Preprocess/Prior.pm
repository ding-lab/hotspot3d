package TGI::Mutpro::Preprocess::Prior;
#
#----------------------------------
# $Authors: Beifang Niu and Adam D Scott
# $Date: 2014-01-14 14:34:50 -0500 (Tue Jan 14 14:34:50 CST 2014) $
# $Revision:  $
# $URL: $
# $Doc: $ do prioritization 
#----------------------------------
#
use strict;
use warnings;

use Carp;
use Getopt::Long;
use IO::File;
use FileHandle;

sub new {
    my $class = shift;
    my $this = {};
    $this->{_OUTPUT_DIR} = undef;
    $this->{_PVALUE_CUTOFF} = 0.05;
    $this->{_3D_CUTOFF} = 20;
    $this->{_1D_CUTOFF} = 0;
    $this->{_STAT} = undef;
    bless $this, $class;
    $this->process();
    return $this;
}

sub process {
    my $this = shift;
    my ( $help, $options );
    unless( @ARGV ) { die $this->help_text(); }
    $options = GetOptions (
        'output-dir=s' => \$this->{_OUTPUT_DIR},
        'p-value-cutoff=f' => \$this->{_PVALUE_CUTOFF},
        '3d-distance-cutoff=i' => \$this->{_3D_CUTOFF},
        'linear-cutoff=i' => \$this->{_1D_CUTOFF},
        'help' => \$help,
    );
    if ( $help ) { print STDERR help_text(); exit 0; }
    unless( $options ) { die $this->help_text(); }
    unless( $this->{_OUTPUT_DIR} ) { warn 'You must provide output directory ! ', "\n"; die help_text(); }
    unless( -d $this->{_OUTPUT_DIR} ) { warn 'You must provide a valid output directory ! ', "\n"; die help_text(); }
    my $fh   = new FileHandle;
    #### processing ####
    # do prioritization 
    my ( $UniprotIdFile, $proximityDir, $cosmicDir, $prioritizationDir, );
    $UniprotIdFile = "$this->{_OUTPUT_DIR}\/hugo.uniprot.pdb.transcript.csv";
    $proximityDir = "$this->{_OUTPUT_DIR}\/proximityFiles";
    $cosmicDir = "$proximityDir\/cosmicanno";
    $prioritizationDir = "$this->{_OUTPUT_DIR}\/prioritization";
    unless( -d $cosmicDir ) { warn "You must provide a valid COSMIC annotations directory ! \n"; die help_text(); }
    unless( -e $prioritizationDir ) { mkdir($prioritizationDir) || die "can not make prioritization result files directory\n"; }
    my $fhunipro = new FileHandle;
    unless( $fhunipro->open("<$UniprotIdFile") ) { die "Could not open Uniprot ID file !\n" };
    my $u = 0;
    while ( my $line = <$fhunipro> ) {
        chomp $line;
        my ( undef, $uniprotId, ) = split /\t/, $line;
        # Only use Uniprot IDs with PDB structures
        next if ( $uniprotId !~ /\w+/ );
        # proximity file
        my $cosmicFile = "$cosmicDir\/$uniprotId\.ProximityFile\.csv";
        next unless( -e $cosmicFile );
        my $outputFile = "$prioritizationDir\/$uniprotId\.ProximityFile\.csv";
        print STDOUT $uniprotId."\n";
        # add annotation infor
        $this->doPrior( $cosmicFile, $outputFile, $uniprotId, $this->{_PVALUE_CUTOFF}, $this->{_3D_CUTOFF}, $this->{_1D_CUTOFF} );
        #delete file if null
    }
    $fhunipro->close();
}

# prioritization based on 
# COSMIC annotation results
sub doPrior {
    my ( $this, $proximityfile, $outputf, $uniprotId, $pvalue, $threed, $lineard, ) = @_;
    print STDOUT "-------$proximityfile---\n";
    # read COSMIC annotation information
    my $fhproximity = new FileHandle;
    unless( $fhproximity->open("<$proximityfile") ) { die "Could not open COSMIC annotated proximity file !\n" };
    # hash for filtering same pairs
    # but keep distances and P_vales
    my %ss;
    while ( my $a = <$fhproximity> ) {
        next if ($a =~ /^WARNING:/);
        next if ($a =~ /UniProt_ID1/);
        chomp($a);
        my @t = split /\t/, $a;
        my ( $annoOneEnd, $annoTwoEnd, $uniprotCoorOneEnd, $uniprotCoorTwoEnd, );
        $annoOneEnd = $annoTwoEnd = "N\/A";
        next if ( ($t[2] =~ /N\/A/) or 
                  ($t[3] =~ /N\/A/) or 
                  ($t[9] =~ /N\/A/) or ($t[10] =~ /N\/A/));
        next unless ( ($t[2] =~ /\d+/) and
                      ($t[3] =~ /\d+/) and
                      ($t[9] =~ /\d+/) and ($t[10] =~ /\d+/) );
        my $oneEndContent = join("\t", @t[0..6]);
        my $twoEndContent = join("\t", @t[7..13]);
        my $proInfo       = join(" ", @t[14..16]);
        $uniprotCoorOneEnd = $t[2] + $t[3];
        $uniprotCoorTwoEnd = $t[9] + $t[10];
        my $tlineard = abs($uniprotCoorOneEnd - $uniprotCoorTwoEnd);
        next if ( ($t[0] eq $t[7]) and ( $tlineard <= $lineard) );
        next if ( $t[16] > $pvalue );
        next if ( $t[14] > $threed );
        # load infor into %ss hash
        if ((defined $ss{$oneEndContent}{$twoEndContent}) or (defined $ss{$twoEndContent}{$oneEndContent}) ) {
            if (defined $ss{$oneEndContent}{$twoEndContent}) {
                $ss{$oneEndContent}{$twoEndContent}{$proInfo} = 1;
            } else { $ss{$twoEndContent}{$oneEndContent}{$proInfo} = 1; }
        } else { $ss{$oneEndContent}{$twoEndContent}{$proInfo} = 1; }
    }
    $fhproximity->close();
    my $fhout = new FileHandle;
    unless( $fhout->open(">$outputf") ) { die "Could not open prioritization output file to write !\n" };
	print STDOUT "Creating ".$outputf."\n";
    # write prioritization result into file 
	$fhout->print( "UniProt_ID1\tChain1\tPosition1\tOffset1\tResidue_Name1\t" );
	$fhout->print( "Feature_Description1\tCOSMIC_Info1\t" );
	$fhout->print( "UniProt_ID2\tChain2\tPosition2\tOffset2\tResidue_Name2\t" );
	$fhout->print( "Feature_Description2\tCOSMIC_Info2\t" );
	$fhout->print( "Distance\tPDB_ID\tP_Value\n" );
    foreach my $a (keys %ss) {
        foreach my $b (keys %{$ss{$a}}) {
            print $fhout $a."\t".$b."\t";
            foreach my $c (keys %{$ss{$a}{$b}}) { print $fhout $c."|"; }
            print $fhout "\n";
        }
    }
    $fhout->close();
}

sub help_text{
    my $this = shift;
        return <<HELP

Usage: hotspot3d prior [options]

                             REQUIRED
--output-dir                 Output directory

                             OPTIONAL
--p-value-cutoff             p_value cutoff(<=), default is 0.05
--3d-distance-cutoff         3D distance cutoff (<= Angstroms), default is 20
--linear-cutoff              Linear distance cutoff (> peptides), default is 0

--help                       this message

HELP

}

1;

