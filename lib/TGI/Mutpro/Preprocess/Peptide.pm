package TGI::Mutpro::Preprocess::Peptide;
#
#----------------------------------
# $Authors: Beifang Niu 
# $Date: 2014-01-14 14:34:50 -0500 (Tue Jan 14 14:34:50 CST 2014) $
# $Revision:  $
# $URL: $
# $Doc: $ peptide class 
#----------------------------------
#
use strict;
use Carp;
use TGI::Mutpro::Preprocess::AminoAcid;
# Ordered amino acids
sub new {    
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self = {};
    $self->{NAME} = undef; # Chain id from PDB file
    $self->{AA} = {}; # Hash of AminoAcid objects
                      # key = residue number
                      # value = ref to AminoAcid object
    bless ($self, $class);
    return $self;
}

sub name {
    # Chain name (single letter)
    my $self = shift;
    if (@_) { $self->{NAME} = shift; }
    return $self->{NAME};
}

sub addPointToAminoAcid {
    # Add to AminoAcid object that is in this peptide
    # If there is no AminoAcid object for the given position, 
    # make one
    # Input: 3-letter code for amino acid, 
    # position in this peptide, 
    # (x,y,z) coordinates
    my $self = shift;
    my ($name, $position, $x, $y, $z) = @_;
    my ( $aaRef, $aa );
    # Make a new AminoAcid object if this position 
    # has not yet been added to peptide
    if (!defined ${$self->{AA}}{$position}) {
		$aa = new TGI::Mutpro::Preprocess::AminoAcid;
		$aa->name($name);
		$aa->position($position);
		$aa->chain($self->{NAME});
		${$self->{AA}}{$position} = \$aa;
    }
    # Add point
    $aaRef = $self->getAminoAcidObject($position);
    # If there are two different names for the given amino acid. 
    # It is ambiguous and should be deleted
    if (defined $aaRef && $$aaRef->name() ne $name) { $$aaRef->ambiguous(1); }
    $$aaRef->addPoint($x,$y,$z);
}
  
sub removeAmbiguousAminoAcids {
    # If any of the amino acid positions in this peptide have more than 
    # one amino acid name associated with them
    # they were marked as ambiguous e.g. position 12 is ambiguous
    #    ATOM      1  N  APRO A   12       3.278  21.202  20.087  0.83 56.23           N  
    #    ATOM      8  N  BSER A   12       3.302  21.148  20.087  0.17 56.57           N 
    my $self = shift;
    foreach my $position ( keys %{$self->{AA}} ) {
	my $aaRef = $self->getAminoAcidObject($position);
	if ( $$aaRef->ambiguous() ) { delete ${$self->{AA}}{$position}; }
    }

}  

sub addAminoAcid {
    # Input: amino acid residue number, 
    # ref to AminoAcid object
    my $self = shift;
    my ($position, $aaRef) = @_;
	if ( $aaRef->isAA( $aaRef->name() ) ) {
		${$self->{AA}}{$position} = $aaRef;
	}
}

sub getAminoAcidObject {
    # Input: position number
    # Return: ref to AminoAcid object at that position
    my $self = shift;
    my $position = shift;
    if ( defined ${$self->{AA}}{$position} ) {
        return ${$self->{AA}}{$position};
    } else { return undef; }
}

sub getAllAminoAcidObjects {
    # Return: ref to hash with key = residue number, 
    # value = ref to AminoAcid object
    my $self = shift;
    return \%{$self->{AA}};
}

sub aminoAcidPositionNumbers {
    # Return ref to array with the position 
    # numbers of all AminoAcid objects (sorted numerically)
    my $self = shift;
    my @positions = keys %{$self->{AA}};
    @positions = sort {$a<=>$b} @positions;
    return \@positions;
}

1;

