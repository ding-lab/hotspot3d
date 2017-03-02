package TGI::Mutpro::Preprocess::Cosmic;
#
#----------------------------------
# $Authors: Beifang Niu & Adam D Scott
# $Date: 2014-01-14 14:34:50 -0500 (Tue Jan 14 14:34:50 CST 2014) $
# $Revision: 1 $
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
    # add COSMIC annotations after get ROI annotation information
    ## do that after get ROI annotation
	my ( $fh , $transMapToUniprot ) = $this->getFileInputs( );
	my ( $annotationsDir , $cosmicDir ) = $this->getInputDirs( );
	my $cosmicHashRef = $this->getCOSMICInput( $transMapToUniprot );
	$this->makeCOSMICAnnotations( $fh , $annotationsDir , $cosmicDir , $cosmicHashRef );
	return 0;
}

sub makeCOSMICAnnotations {
	my ( $this , $fh , $annotationsDir , $cosmicDir , $cosmicHashRef ) = @_;
    while ( my $line = <$fh> ) {
        chomp $line;
        my ( undef, $uniprotId, ) = split /\t/, $line;
        # Only use Uniprot IDs with PDB structures
        next if ( $uniprotId !~ /\w+/ );
        # proximity file
        # add annotation infor
        $this->addCosmic( $annotationsDir , $cosmicDir , $cosmicHashRef , $uniprotId );
        #delete file if null
    }
    $fh->close();
	return;
}

sub getCOSMICInput {
	my ( $this , $transMapToUniprot ) = @_;
    my $cosmicAnno = "$this->{_OUTPUT_DIR}\/cosmic\/cosmic_67_for_HotSpot3D_missense_only.tsv";
    my $cosmicHashRef = $this->getCosmicInfor( $transMapToUniprot , $cosmicAnno );
	return $cosmicHashRef;
}

sub getInputDirs {
	my $this = shift;
    my $proximityDir = "$this->{_OUTPUT_DIR}\/proximityFiles";
    my $annotationsDir = "$proximityDir\/annotations";
    unless( -d $annotationsDir ) { warn "You must provide a valid annotations directory ! \n"; die help_text(); }
    my $cosmicDir = "$proximityDir\/cosmicanno";
    unless( -e $cosmicDir ) { mkdir( $cosmicDir ) || die "can not make COSMIC annotated files directory\n"; }
	return ( $annotationsDir , $cosmicDir );
}

sub getFileInputs {
	my $this = shift;
    my $hugoUniprotFile = "$this->{_OUTPUT_DIR}\/hugo.uniprot.pdb.transcript.csv";
    my $transMapToUniprot = $this->getTransMaptoUniprot( $hugoUniprotFile );
    my $fh = new FileHandle;
    unless( $fh->open("<$hugoUniprotFile") ) { die "Could not open hugo uniprot file !\n" };
	return ( $fh , $transMapToUniprot );
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
    unless( $this->{_OUTPUT_DIR} ) { warn 'HotSpot3D Cosmic Error: You must provide a output directory! ', "\n"; die $this->help_text(); };
    unless( -e $this->{_OUTPUT_DIR} ) { warn 'HotSpot3D Cosmic Error: Output directory does not exist! ', "\n"; die $this->help_text(); };
	return;
}

# add Cosmic annotation
sub addCosmic {
    my ( $this , $annotationsDir , $cosmicDir , $cosmicRef , $uniprotId ) = @_;
    # add annotation information
	my $annotationFile = "$annotationsDir\/$uniprotId\.ProximityFile\.csv";
	unless ( -e $annotationFile ) {
		warn "HotSpot3D::Cosmic::addCosmic warning: cannot add COSMIC, because there is no annotation file: ".$annotationFile."\n";
		return;
	}
    my $fhin = new FileHandle;
    unless( $fhin->open("<$annotationFile ") ) { die "Could not open proximity file !\n" };
    my $fhout = new FileHandle;
	my $outputFile = "$cosmicDir\/$uniprotId\.ProximityFile\.csv";
    unless( $fhout->open(">$outputFile") ) { die "Could not open output proximity file to write !\n" };
	print STDOUT $uniprotId." HotSpot3D::Cosmic::addCosmic - writing cosmic annotations to ".$outputFile."\n";
    #print $uniprotId."\n";
	my $skips = 0;
	my %structures;
	my $newlines = 0;
    while ( my $line = <$fhin> ) {
        if ($line =~ /^WARNING:/) {
			$skips += 1;
		}
        chomp($line);
        my @fields = split /\t/, $line;
		$structures{$fields[13]} += 1;
        my ( $annoOneEnd, $annoTwoEnd, $uniprotCoorOneEnd, $uniprotCoorTwoEnd, );
        $annoOneEnd = $annoTwoEnd = "N\/A";
        next if ( ($fields[2] =~ /N\/A/) or 
                  ($fields[3] =~ /N\/A/) or 
                  ($fields[8] =~ /N\/A/) or
				  ($fields[9] =~ /N\/A/));
		$fields[2] = TGI::Data::CleanNumber::numOnly( $fields[2] );
		$fields[3] = TGI::Data::CleanNumber::numOnly( $fields[3] );
		$fields[8] = TGI::Data::CleanNumber::numOnly( $fields[8] );
		$fields[9] = TGI::Data::CleanNumber::numOnly( $fields[9] );
        $uniprotCoorOneEnd = $fields[2] + $fields[3];
        $uniprotCoorTwoEnd = $fields[8] + $fields[9];
        #warn $uniprotCoorOneEnd."\t".$uniprotCoorTwoEnd."\n";
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
        $newLine .= join("\t", @fields[0..5]);
        $newLine .= "\t"; 
        $newLine .= $annoOneEnd."\t";
        $newLine .= join("\t", @fields[6..11]);
        $newLine .= "\t";
        $newLine .= $annoTwoEnd."\t";
        $newLine .= join("\t", @fields[12..14]);
        #warn $newLine."\n";
        print $fhout $newLine."\n";
		$newlines += 1;
    }
    $fhin->close();
	print STDOUT $uniprotId." Skipped ".$skips." structures\n";
	my $nStructures = scalar keys %structures;
	print STDOUT $uniprotId." Processed ".$nStructures." structures\n";
	print STDOUT $uniprotId." There are ".$newlines." annotated lines in ".$annotationFile."\n";
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
        #warn $aac."\t".$1."\t".$tissue."\t".$uniprot."\n";
        #warn $aac."\t".$1."\t".$uniprot."\n";
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

