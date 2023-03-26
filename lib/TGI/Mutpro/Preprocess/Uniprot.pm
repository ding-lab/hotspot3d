package TGI::Mutpro::Preprocess::Uniprot;
#
#----------------------------------
# $Original authors: Beifang Niu
# $Modified by: Fernanda Martins Rodrigues @WashU (fernanda@wustl.edu; mrodrigues.fernanda@gmail.com)
# $Date: 2014-01-14 14:34:50 -0500 (Tue Jan 14 14:34:50 CST 2014) $
# $Revision: 2023-03-26 $
# $URL: $
# $Doc: $ uniprot download and processing. This has been adapted on Mar 15 2023 to handle AlphaFold DB (v4) files 
#----------------------------------
#
use strict;
use warnings;

use LWP::Simple;
use Carp;
#
#  Stores the text version of a Uniprot entry downloaded
#  using ""http://www.uniprot.org/uniprot/$uniprotId.txt";
#
sub new {
    my ($class, $uniprotId) = @_;    
    bless {
		ID => $uniprotId,
		PAGE => ""
    }, $class;
}

sub uniprotId {
    my $self = shift;    
    if (@_) { $self->{ID} = shift; }
    return $self->{ID};
}

sub entireRecord {
    # Return: entire text record stored as a mulit-line page
    # Can be parsed by doing 
    #  foreach $line (split /\n/, $page) { chomp $line; etc.. }
    my $self = shift;
    if ( $self->{PAGE} eq "" ) { $self->retrieveRecordFromUniprot(); }
    return $self->{PAGE};
}

sub retrieveRecordFromUniprot {
    # Download current record from Uniprot web site.
    my $self = shift;
    my $uniprotId = $self->uniprotId();
    my $uniProtUrl = "http://www.uniprot.org/uniprot/$uniprotId.txt";
    my $page = get($uniProtUrl);
    if ( !defined $page ) { $page = "empty"; }
    $self->{PAGE} = $page;
}

sub generalAnnotation {
    # Returns: ref to array of general comments 
    # These are in the record preceeded by 'CC' and the
    # general annotion type is given in all caps e.g.
    # CC   -!- FUNCTION: Acts as a tumor suppressor in many tumor types; induces
    # CC       growth arrest or apoptosis depending on the physiological
    # CC       etc.
    # There are about 30 possible entries, but don't know how each one is
    # formated.  These are some:
    # CC   -!- FUNCTION: Function
    # CC   -!- COFACTOR: Cofactor
    # CC   -!- SUBUNIT:  Subunit structure
    # CC   -!- SUBCELLULAR LOCATION: Subcellular location
    # CC   -!- DOMAIN:   Domain
    # CC   -!- PTM:      Post-translational modification
    # CC   -!- DISEASE:  Involvement in disease
    # CC   -!- SIMILARITY:  Sequence similarities
    # CC   -!- INTERACTION:
    my $self = shift;
    my $annotationType = shift;
    my ( @annotations, $readingType, );
    $readingType = 0;
    foreach my $line ( split /\n/, $self->entireRecord() ) {
        chomp $line;
        if ( $line =~ /CC\s+\-\!\-\s+\w+/ ) { $readingType = ( $line =~ /CC\s+\-\!\-\s+$annotationType/i ) ? 1 : 0; }
	if ( $line =~ /CC\s+Copyrighted/ ) {  $readingType = 0; last; }
	if ( $readingType ) { push @annotations, $line; }
    }
    return \@annotations;
}

# MODIFIED ON 03-15-2023: commented this function out since it is not actually used and we are using Alpha Fold DB for this version of the tool
# sub parsePDBAnnotation {
# 	my ( $self , $pdbannotation ) = @ARGV;
# 	my $details = {};
# 	my ( $pdbID , $type , $resolution , $chainInfo ) = split( ";" , $pdbannotation );
# 	$pdbID =~ s/\s+//g;
# 	$type =~ s/\s+//g;
# 	$resolution =~ s/\s+//g;
# 	$chainInfo =~ s/\s+\.*//g;
# 	my ( $chains , $positions ) = split( "=" , $chainInfo );


# 	return $details;
# }

# MODIFIED ON 03-15-2023: added AlphaFold DB example line from uniprot file
sub annotations {
    # Returns: ref to array of annotations of given type
    # These are in the record as 'DR $annotationName; $id; $record'
    # e.g. 
    # DR   InterPro; IPR000008; C2_Ca-dep.
    # DR   Pfam; PF00168; C2; 2
    # DR   PDB; 1V27; NMR; -; A=807-934.
    # DR   UniGene; Hs.655271; -.
    # DR   AlphaFoldDB; Q3USB1; -.
    # DR   GO; GO:0005783; C:endoplasmic reticulum; IEA:UniProtKB-SubCell.
    #    The GO annotation is described in detail below
    #
    my $self = shift;
    my $annotationType = shift;
    my @annotations;
    (defined $annotationType) || confess "Did not get an annotation type as input parameter";
    foreach my $line ( split /\n/, $self->entireRecord() ) { 
        chomp $line; 
        if ( $line =~ /DR\s+$annotationType\;\s+(.*)/i ) { push @annotations, $1; }
    }
    return \@annotations;
}

sub sequenceFeatures {
    # Large section that includes domains
    # It has 'FT   $Tag  $start  $stop $description'
    # Can be used to get MUTAGEN (sites that have been mutagenized) CONFLICT VAR_SEQ, etc.
    my $self = shift;
    my $tag = shift;
    (defined $tag) || carp "No tag sent to sub sequenceFeature", return undef;
    my ( $line, @features );
    @features = ();
    foreach $line ( split /\n/, $self->entireRecord() ) {
        chomp $line;
        if ( $line =~ /FT\s+$tag/ ) { push @features, $line; }
    }

    return \@features;
}

sub domains {
    # Input: start and <optional stop> for region of interest
    # Return: ref to array of domains that overlap the given region
    # Skips the following entries: CONFLICT STRAND HELIX TURN VAR_SEQ INIT_MET VARIANT CHAIN
    # UPDATE ON 07/13/2021 due to changes in Uniprot file formatting;
    my $self = shift;
    my ($start, $stop) = @_;
    (defined $start) || carp "Start not defined for $self->{ID}", return undef;
    if ( !defined $stop ) { $stop = $start; }
    my ( @domains, $line, %skipList, $desc, $key, $dmStart, $dmStop, $countEvidenceLine, $tmpvar, $addEntry );
    $countEvidenceLine = 0;
    $addEntry = 0;
    my @skipThese = qw (CONFLICT STRAND HELIX TURN VAR_SEQ INIT_MET VARIANT CHAIN);
    map{ $skipList{$_} = 1; } @skipThese;
    foreach $line (split /\n/, $self->entireRecord()) {
        chomp $line;
        if ($line =~ /^FT/ ){
            if ( $line =~ /^FT\s+(\S+)\s+(\d+)\.+(\d+)$/ && $stop >= $2 && $start <= $3 && $countEvidenceLine == 0 && $addEntry == 0){#START AND STOP ARE SEPARATED BY PERIOD INSTEAD OF SPACE
                next if ( defined $skipList{$1} );
                $key = $1; $dmStart = $2; $dmStop = $3; $countEvidenceLine=1;
            }elsif ($line =~ /^FT\s+(\S+)\s+(\d+)$/ && $stop >= $2 && $start <= $2 && $countEvidenceLine == 0 && $addEntry == 0){#CASES WHERE ONLY START POS IS GIVEN
                next if ( defined $skipList{$1} );
                $key = $1; $dmStart = $2; $dmStop = $2; $countEvidenceLine=1;
            }elsif ($line =~ /^FT\s+\/(\S+)\=\"(.*)\"$/ && $countEvidenceLine == 1 && $addEntry == 0){
                $desc=$2; $addEntry=1;
            }elsif ( $line !~ /^FT\s+\/(\S+)\=\"(.*)\"$/ && $line =~ /^FT\s+\/(\S+)\=\"(.*)$/ && $countEvidenceLine == 1 && $addEntry == 0){
                $desc=$2;
                $countEvidenceLine=2;
            }elsif ( $line !~ /^FT\s+\/(\S+)\=\"(.*)$/ && $line =~ /^FT\s+(.*)\"$/ && $countEvidenceLine >= 2 && $addEntry == 0){
                $tmpvar=$desc;
                $desc=$tmpvar . $1;
                $addEntry=1;
            }elsif ( $line !~ /^FT\s+\/(\S+)\=\"(.*)$/ && $line =~ /^FT\s+(.*)$/ && $countEvidenceLine >= 2 && $addEntry == 0){
                $tmpvar=$desc;
                $desc=$tmpvar . $1;
                $countEvidenceLine=++$countEvidenceLine;
            }
        }elsif ($line =~ /FT\s+DISULFID\s+(\d+)\.+(\d+)/ && ( ($stop >= $1 && $start <= $1) || ($stop >= $2 && $start <= $2) )) {
            push @domains, "DISULFID\t($1, $1)\tDisulfide bond $1 <-> $2.";
            push @domains, "DISULFID\t($2, $2)\tDisulfide bond $1 <-> $2.";
        }elsif ($line =~ /FT\s+DISULFID\s+(\d+)$/ && ( ($stop >= $1 && $start <= $1) || ($stop >= $1 && $start <= $1) )) {
            push @domains, "DISULFID\t($1, $1)\tDisulfide bond $1 <-> $1.";
            push @domains, "SIGNAL\t($1, $1)\tSignal peptide.";
        }
        if ($addEntry == 1){
            push @domains, "$key\t($dmStart, $dmStop)\t$desc";
            $countEvidenceLine = 0; $addEntry = 0;
        }
    }
    return \@domains;
}

sub domainsForMultiplePositions {
    # Input: ref to hash of positions in Uniprot coordinate system
    # Return: ref to array of domains that contain at least one of
    #         the positions
    # Skips the following entries: CONFLICT STRAND HELIX TURN VAR_SEQ INIT_MET VARIANT CHAIN
    # UPDATE ON 07/13/2021 due to changes in Uniprot file formatting;
    my ($self, $positionRef) = @_;
    my ( $position, $addEntry, @domains, $line, %skipList, $desc, $key, $dmStart, $dmStop, $countEvidenceLine, $tmpvar, $done );
    my @skipThese = qw (CONFLICT STRAND HELIX TURN VAR_SEQ INIT_MET VARIANT CHAIN);
    map{ $skipList{$_} = 1; } @skipThese;
    $countEvidenceLine=0;
    $addEntry=0;
    foreach $line ( split /\n/, $self->entireRecord() ) {
        chomp $line;
        if ( $line =~ /^FT/ ){
            if ( $line =~ /^FT\s+(\S+)\s+(\d+)\.+(\d+)$/ && !defined $skipList{$1} && $countEvidenceLine == 0 && $addEntry == 0) {
                $key = $1; $dmStart = $2; $dmStop = $3; $countEvidenceLine=1;
            }elsif ( $line =~ /^FT\s+(\S+)\s+(\d+)$/ && !defined $skipList{$1} && $countEvidenceLine == 0 && $addEntry == 0){
                $key = $1; $dmStart = $2; $dmStop = $2; $countEvidenceLine=1;
            }elsif ( $line =~ /^FT\s+\/(\S+)\=\"(.*)\"$/ && $countEvidenceLine == 1 && $addEntry == 0){
                $desc=$2; $addEntry=1;
            }elsif ( $line !~ /^FT\s+\/(\S+)\=\"(.*)\"$/ && $line =~ /^FT\s+\/(\S+)\=\"(.*)$/ && $countEvidenceLine == 1 && $addEntry == 0){
                $desc=$2;
                $countEvidenceLine=2;
            }elsif ( $line !~ /^FT\s+\/(\S+)\=\"(.*)$/ && $line =~ /^FT\s+(.*)\"$/ && $countEvidenceLine >= 2 && $addEntry == 0){
                $tmpvar=$desc;
                $desc=$tmpvar . $1;
                $addEntry=1;
            }elsif ( $line !~ /^FT\s+\/(\S+)\=\"(.*)$/ && $line =~ /^FT\s+(.*)$/ && $countEvidenceLine >= 2 && $addEntry == 0){
                $tmpvar=$desc;
                $desc=$tmpvar . $1;
                $countEvidenceLine=++$countEvidenceLine;
            }
        }
        if ( $addEntry==1){
            $done = 0;
            # See if any of the positions are within this domain
            foreach $position (keys %{$positionRef}) {
                if ( $position >= $dmStart && $position <= $dmStop ) { $done = 1; }
            }
            if ( $done ) { push @domains, "$key\t($dmStart, $dmStop)\t$desc"; }
            $countEvidenceLine = 0; $addEntry = 0; $done = 0;
        }
    }
    return \@domains;
}

sub domainsAfterPosition {
    # This is for nonsense mutations.
    # List all domains that occur after the given position
    # Input: amino acid start position
    # Return: ref to array of domains that are after the position
    # Skips the following entries: CONFLICT STRAND HELIX TURN VAR_SEQ INIT_MET VARIANT CHAIN
    # UPDATE ON 07/13/2021 due to changes in Uniprot file formatting;
    my $self = shift;
    my $position = shift;
    (defined $position) || carp "\$position not defined for $self->{ID}", return undef;
    my ( @domains, $line, %skipList, $desc, $key, $dmStart, $dmStop, $countEvidenceLine, $tmpvar, $addEntry );
    my @skipThese = qw (CONFLICT STRAND HELIX TURN VAR_SEQ INIT_MET VARIANT CHAIN);
    map{ $skipList{$_} = 1; } @skipThese;
    $countEvidenceLine=0;
    $addEntry=0;
    foreach $line ( split /\n/, $self->entireRecord() ) {
	chomp $line;
	if ( $line =~ /^FT/ ){
            if ( $line =~ /^FT\s+(\S+)\s+(\d+)\.+(\d+)$/ && $3 > $position && $countEvidenceLine == 0 && $addEntry == 0) {
                next if ( defined $skipList{$1} );
                $key = $1; $dmStart = $2; $dmStop = $3; $countEvidenceLine=1;
            }elsif ( $line =~ /^FT\s+(\S+)\s+(\d+)$/ && $2 > $position && $countEvidenceLine == 0 && $addEntry == 0){
                next if ( defined $skipList{$1} );
                $key = $1; $dmStart = $2; $dmStop = $2; $countEvidenceLine=1;
            }elsif ($line =~ /^FT\s+\/(\S+)\=\"(.*)\"$/ && $countEvidenceLine == 1 && $addEntry == 0){
                $desc=$2; $addEntry=1;
            }elsif ( $line !~ /^FT\s+\/(\S+)\=\"(.*)\"$/ && $line =~ /^FT\s+\/(\S+)\=\"(.*)$/ && $countEvidenceLine == 1 && $addEntry == 0){
                $desc=$2;
                $countEvidenceLine=2;
            }elsif ( $line !~ /^FT\s+\/(\S+)\=\"(.*)$/ && $line =~ /^FT\s+(.*)\"/ && $countEvidenceLine >= 2 && $addEntry == 0){
                $tmpvar=$desc;
                $desc=$tmpvar . $1;
                $addEntry=1;
            }elsif ( $line !~ /^FT\s+\/(\S+)\=\"(.*)$/ && $line =~ /^FT\s+(.*)$/ && $countEvidenceLine >= 2 && $addEntry == 0){
                $tmpvar=$desc;
                $desc=$tmpvar . $1;
                $countEvidenceLine=++$countEvidenceLine;
            }
        }
        if ( $addEntry==1 ){
            push @domains, "$key\t($dmStart, $dmStop)\t$desc";
            $countEvidenceLine=0;
            $addEntry=0;
        }
    }
    return \@domains;
}

sub sequence {
    # Returns sequence as a string
    my $self = shift;
    my ( $line, $seq, $readingSeq, );
    $readingSeq = 0;
    $seq = "";
    my $record = $self->entireRecord();
    if ( !defined $record ) { return undef; }
    foreach  $line ( split /\n/, $self->entireRecord() ) {
	chomp $line;
	if ( $line =~ /\/\// ) { $readingSeq = 0; }
	if ( $line =~ /SQ\s+SEQUENCE/ ) { $readingSeq = 1; next; }
	if ( $readingSeq ) { $line =~ s/\s+//g; $seq .= $line; }
    }
    return $seq;
}

# Beifang Niu 05-06-2013
# added this function to get the Ensemble thranscript id
# and corresponding Protein id
# Fernanda Martins Rodrigues 03-26-2030
# fixed bug with matching line in newly formatted uniprot file
sub transProteinHash{
    # Returns hash:
    # $hash{transcriptid} = protein id
    my $self = shift;
    my %transProtein;
    my ( $line, $seq, $readingSeq,);
    $readingSeq = 0; 
    $seq = "";
    my $record = $self->entireRecord();
    if ( !defined $record ) { return undef; }
    foreach $line ( split /\n/, $self->entireRecord()) {
	chomp $line;
	if ( $line =~ /^DR\s+Ensembl;\s+(\w+.\d+);\s+(\w+.\d+);\s+(\w+.\d+)/ ) { $transProtein{$1} = $2; }
    }
    return \%transProtein;
}

return 1;
package TGI::Mutpro::Preprocess::Uniprot;
#
#----------------------------------
