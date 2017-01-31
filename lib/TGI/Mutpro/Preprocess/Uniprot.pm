package TGI::Mutpro::Preprocess::Uniprot;
#
#----------------------------------
# $Authors: Beifang Niu 
# $Date: 2014-01-14 14:34:50 -0500 (Tue Jan 14 14:34:50 CST 2014) $
# $Revision:  $
# $URL: $
# $Doc: $ uniprot download and processing 
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

sub parsePDBAnnotation {
	my ( $self , $pdbannotation ) = @ARGV;
	my $details = {};
	my ( $pdbID , $type , $resolution , $chainInfo ) = split( ";" , $pdbannotation );
	$pdbID =~ s/\s+//g;
	$type =~ s/\s+//g;
	$resolution =~ s/\s+//g;
	$chainInfo =~ s/\s+\.*//g;
	my ( $chains , $positions ) = split( "=" , $chainInfo );


	return $details;
}

sub annotations {
    # Returns: ref to array of annotations of given type
    # These are in the record as 'DR $annotationName; $id; $record'
    # e.g. 
    # DR   InterPro; IPR000008; C2_Ca-dep.
    # DR   Pfam; PF00168; C2; 2
    # DR   PDB; 1V27; NMR; -; A=807-934.
    # DR   UniGene; Hs.655271; -.
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
    my $self = shift;
    my ($start, $stop) = @_;
    (defined $start) || carp "Start not defined for $self->{ID}", return undef;
    if ( !defined $stop ) { $stop = $start; }
    my ( @domains, $line, %skipList, $desc, $key, $dmStart, $dmStop );
    my @skipThese = qw (CONFLICT STRAND HELIX TURN VAR_SEQ INIT_MET VARIANT CHAIN);
    map{ $skipList{$_} = 1; } @skipThese;
    foreach $line (split /\n/, $self->entireRecord()) {
	chomp $line;
	if ($line =~ /FT\s+(\S+)\s+(\d+)\s+(\d+)\s+(.*)/ && $stop >= $2 && $start <= $3 ) { 
	    $key = $1; $dmStart = $2; $dmStop = $3; $desc = $4;
	    next if ( defined $skipList{$key} );
	    push @domains, "$key\t($dmStart, $dmStop)\t$desc";
	}elsif ($line =~ /FT\s+DISULFID\s+(\d+)\s+(\d+)/ && ( ($stop >= $1 && $start <= $1) || ($stop >= $2 && $start <= $2) )) { 
	    push @domains, "DISULFID\t($1, $1)\tDisulfide bond $1 <-> $2.";
	    push @domains, "DISULFID\t($2, $2)\tDisulfide bond $1 <-> $2.";
	}elsif ( $line =~ /FT\s+SIGNAL\s+(\d+)\s+(\d+)/ && $stop >= $2 && $start <= $1 ) { 
	    push @domains, "SIGNAL\t($1, $2)\tSignal peptide.";
	}
    }
    return \@domains;
}

sub domainsForMultiplePositions {
    # Input: ref to hash of positions in Uniprot coordinate system
    # Return: ref to array of domains that contain at least one of
    #         the positions
    # Skips the following entries: CONFLICT STRAND HELIX TURN VAR_SEQ INIT_MET VARIANT CHAIN
    my ($self, $positionRef) = @_;
    my ( $position, $addEntry, @domains, $line, %skipList, $desc, $key, $dmStart, $dmStop );
    my @skipThese = qw (CONFLICT STRAND HELIX TURN VAR_SEQ INIT_MET VARIANT CHAIN);
    map{ $skipList{$_} = 1; } @skipThese;
    foreach $line ( split /\n/, $self->entireRecord() ) {
        chomp $line;
	if ( $line =~ /FT\s+(\S+)\s+(\d+)\s+(\d+)\s+(.*)/ && !defined $skipList{$1} ) {
            $key = $1; $dmStart = $2; $dmStop = $3; $desc = $4;
	    $addEntry = 0;
	    # See if any of the positions are within this domain
	    foreach $position (keys %{$positionRef}) {
		if ( $position >= $dmStart && $position <= $dmStop ) { $addEntry = 1; }
	    }
	    if ( $addEntry ) { push @domains, "$key\t($dmStart, $dmStop)\t$desc"; }
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
    my $self = shift;
    my $position = shift;
    (defined $position) || carp "\$position not defined for $self->{ID}", return undef;
    my ( @domains, $line, %skipList, $desc, $key, $dmStart, $dmStop );
    my @skipThese = qw (CONFLICT STRAND HELIX TURN VAR_SEQ INIT_MET VARIANT CHAIN);
    map{ $skipList{$_} = 1; } @skipThese;
    foreach $line ( split /\n/, $self->entireRecord() ) {
	chomp $line;
	if ( $line =~ /FT\s+(\S+)\s+(\d+)\s+(\d+)\s+(.*)/ && $3 > $position ) { 
	    $key = $1; $dmStart = $2; $dmStop = $3; $desc = $4;
	    next if ( defined $skipList{$key} );
	    push @domains, "$key\t($dmStart, $dmStop)\t$desc";
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
#
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
	if ( $line =~ /^DR\s+Ensembl;\s+(\w+);\s+(\w+);\s+(\w+)/ ) { $transProtein{$1} = $2; }
    }
    return \%transProtein;
}

#sub getCanonicalTranscript {
#	my $self = shift;
#	foreach my $line ( split /\n/ , $self->annotations( "Ensembl" ) ) {
#		if ( $line =~ /(ENST\d+);.*\[(\w+)-(\d+)\]/ ) {
#			if ( $3 == 1 ) {
#				return $1;
#			}
#		}
#	}
#	return "";
#}
#
#sub getGeneID {
#	my $self = shift;
#	foreach my $line ( split /\n/ , $self->entireRecord() ) {
#		chomp( $line );
#		if ( $line =~ /^ID\s+(\w+)\s+\w+;\s+\d+\sAA\./ ) {
#			return $1;
#		}
#	}
#	return "";
#}
#
#sub getPhosphosites {
#	my $this = shift;
#	my $description = shift;
#	return $this->getModifiedResidues( "Phospho" );
#}
#
#sub getModifiedResidues {
#	my $this = shift;
#	my $description = shift;
#	my $sites = {};
#	foreach my $feature ( split /\n/ , $this->sequenceFeatures() ) {
#		if ( $feature =~ /MOD_RES/ ) {
#			my ( $position , $type ) = $feature =~ m/FT\s+MOD_RES\s+(\d+)\s+\d+\s+(.*)/;
#			next unless ( $type =~ /$description/ );
#			my $detail = "";
#			$sites->{$position}->{type} = $type;
#			if ( $type =~ m/(.*); (.*)\./ ) {
#				$sites->{$position}->{type} = $1;
#				$sites->{$position}->{detail} = $2;
#			} elsif ( $type =~ m/(.*)\.$/ ) {
#				$sites->{$position}->{type} = $1;
#				$sites->{$position}->{detail} = "";
#			} elsif ( $type =~ m/(.*)\.(.*)/ ) {
#				$sites->{$position}->{type} = $1;
#				$sites->{$position}->{detail} = $2;
#			}
#		}
#	}
#	return $sites;
#}

return 1;
