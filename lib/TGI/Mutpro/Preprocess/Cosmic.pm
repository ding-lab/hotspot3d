package TGI::Mutpro::Preprocess::Cosmic;
#
#----------------------------------
# $Authors: Beifang Niu 
# $Date: 2014-01-14 14:34:50 -0500 (Tue Jan 14 14:34:50 CST 2014) $
# $Revision:  $
# $URL: $
# $Doc: $ cosmic database processing 
#----------------------------------
#
use strict;
use warnings;

use Carp;
use Cwd;
use FileHandle;
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
    # add COSMIC annotations after get ROI annotation information
    ## do that after get ROI annotation
    my ( $hugoUniprotFile, $proximityDir, $annotationsDir, $cosmicDir, $cosmicAnno, );
    $hugoUniprotFile = "$this->{_OUTPUT_DIR}\/hugo.uniprot.pdb.transcript.csv";
    $proximityDir = "$this->{_OUTPUT_DIR}\/proximityFiles";
    $annotationsDir = "$proximityDir\/annotations";
    $cosmicDir = "$proximityDir\/cosmicanno";
    $cosmicAnno = "$this->{_OUTPUT_DIR}\/cosmic\/cosmic_67_for_HotSpot3D_missense_only.tsv";
    unless( -d $annotationsDir ) { warn "You must provide a valid annotations directory ! \n"; die help_text(); }
    unless( -e $cosmicDir ) { mkdir( $cosmicDir ) || die "can not make COSMIC annotated files directory\n"; }
    my $transMaptoUnipro = $this->getTransMaptoUniprot( $hugoUniprotFile );
    my $cosmicHashRef = $this->getCosmicInfor( $transMaptoUnipro, $cosmicAnno );

    my $fh = new FileHandle;
    unless( $fh->open("<$hugoUniprotFile") ) { die "Could not open hugo uniprot file !\n" };
    my $u = 0;
    while ( my $line = <$fh> ) {
        chomp $line;
        my ( undef, $uniprotId, ) = split /\t/, $line;
        # Only use Uniprot IDs with PDB structures
        next if ( $uniprotId !~ /\w+/ );
        # proximity file
        my $annotationFile = "$annotationsDir\/$uniprotId\.ProximityFile\.csv";
        next unless( -e $annotationFile );
        my $outputFile = "$cosmicDir\/$uniprotId\.ProximityFile\.csv";
        print STDERR $uniprotId."\n";
        # add annotation infor
        $this->addCosmic( $annotationFile, $cosmicHashRef, $outputFile, $uniprotId );
        #delete file if null
    }
    $fh->close();
}

# add Cosmic annotation
sub addCosmic {
    my ( $this, $proximityfile, $cosmicRef, $outputf, $uniprotId, ) = @_;
    # add annotation information
    my $fhin = new FileHandle;
    unless( $fhin->open("<$proximityfile") ) { die "Could not open proximity file !\n" };
    my $fhout = new FileHandle;
    unless( $fhout->open(">$outputf") ) { die "Could not open output proximity file to write !\n" };
    #print $uniprotId."\n";
    while ( my $a = <$fhin> ) {
        next if ($a =~ /^WARNING:/);
        chomp($a);
        my @t = split /\t/, $a;
        my ( $annoOneEnd, $annoTwoEnd, $uniprotCoorOneEnd, $uniprotCoorTwoEnd, );
        $annoOneEnd = $annoTwoEnd = "N\/A";
        next if ( ($t[2] =~ /N\/A/) or 
                  ($t[3] =~ /N\/A/) or 
                  ($t[8] =~ /N\/A/) or ($t[9] =~ /N\/A/));
        $uniprotCoorOneEnd = $t[2] + $t[3];
        $uniprotCoorTwoEnd = $t[8] + $t[9];
        #print STDERR $uniprotCoorOneEnd."\t".$uniprotCoorTwoEnd."\n";
        if ( defined $cosmicRef->{$uniprotId}->{$uniprotCoorOneEnd} ) {
            $annoOneEnd = "";
            map{ $annoOneEnd .= $_.","; } keys %{$cosmicRef->{$uniprotId}->{$uniprotCoorOneEnd}};
            chop($annoOneEnd)
        }
        if ( defined $cosmicRef->{$uniprotId}->{$uniprotCoorTwoEnd} ) {
            $annoTwoEnd = "";
            map{ $annoTwoEnd .= $_.","; } keys %{$cosmicRef->{$uniprotId}->{$uniprotCoorTwoEnd}};
            chop($annoTwoEnd)
        }
        my $newLine = "";
        $newLine .= join("\t", @t[0..5]);
        $newLine .= "\t"; 
        $newLine .= $annoOneEnd."\t";
        $newLine .= join("\t", @t[6..11]);
        $newLine .= "\t";
        $newLine .= $annoTwoEnd."\t";
        $newLine .= join("\t", @t[12..14]);
        #print STDERR $newLine."\n";
        print $fhout $newLine."\n";
    }
    $fhin->close();
    $fhout->close();
}

# get mapping information
# of transcript id to uniprot id
sub getTransMaptoUniprot {
    my ( $this, $uniprotf, ) = @_;
    my %transHash;
    my $fhunipro = new FileHandle;
    unless( $fhunipro->open("< $uniprotf") ) { die "Could not open hugo uniprot file !\n" };
    while ( my $a = $fhunipro->getline ) {
        chomp( $a );
        my ( undef, $uniprotId, undef, undef, $transcripts, ) = split /\t/, $a;
        next if $transcripts =~ (/N\/A/);
        map{ 
            /(\w+)\[(.*?)]/;
            my $tmp_transcript_id = $1;
            $transHash{$tmp_transcript_id}{'UNIPROT'} = $uniprotId;
            map{  /(\d+)\|(\d+)-(\d+)\|(\d+)/; 
                $transHash{$tmp_transcript_id}{'POSITION'}{$2}{'TEND'} = $4; 
                $transHash{$tmp_transcript_id}{'POSITION'}{$2}{'UBEGIN'} = $1;
                $transHash{$tmp_transcript_id}{'POSITION'}{$2}{'UEND'} = $3;
            } split /\:/, $2;
        } split /,/, $transcripts;
    }
    $fhunipro->close();
    return \%transHash;
}

# get COSMIC annotation
sub getCosmicInfor {
    my ( $this, $uniprotRef, $cosmicf, ) = @_;
    my %cosmicHash;
    my $fhcosmic = new FileHandle;
    unless( $fhcosmic->open( "<$cosmicf" ) ) { die "Could not open COSMIC annotation file !\n" };
    while ( my $a = $fhcosmic->getline ) {
        chomp( $a );
        next if ( $a =~ /^gene_name/ ); 
        my ( $gene, $transcript, $type, $aac, $domain, $tissue, ) = split /\t/, $a;
        next unless ( defined $uniprotRef->{$transcript} );
        next unless ( $aac =~ /p\.\D(\d+)\D/ );
        my $uniprot = $uniprotRef->{$transcript}->{'UNIPROT'};
        my $tmp_hit_bool = 0; my $tmp_uniprot_position;
        foreach my $tmp_pos ( keys %{$uniprotRef->{$transcript}->{'POSITION'}} ){
            if ( ($1 >= $tmp_pos) and ($1 <= $uniprotRef->{$transcript}->{'POSITION'}->{$tmp_pos}->{'TEND'}) ) {
                $tmp_uniprot_position = $1 - $tmp_pos + $uniprotRef->{$transcript}->{'POSITION'}->{$tmp_pos}->{'UBEGIN'};
                $tmp_hit_bool = 1; 
                last;
            } 
        }
        next if ( $tmp_hit_bool == 0 );
        next if (defined $cosmicHash{$uniprot}{$tmp_uniprot_position}{$aac."|".$tissue});
        $cosmicHash{$uniprot}{$tmp_uniprot_position}{$aac."|".$tissue} = 1;
        #print STDERR $aac."\t".$1."\t".$tissue."\t".$uniprot."\n";
        #print STDERR $aac."\t".$1."\t".$uniprot."\n";
    }
    $fhcosmic->close();
    return \%cosmicHash;
}

sub help_text {
    my $this = shift;
        return <<HELP

Usage: hotspot3d cosmic [options]

                             REQUIRED
--output-dir                 Output directory of proximity files

--help                       this message

HELP

}

1;

