package TGI::Mutpro::Preprocess::AminoAcid;
#
#----------------------------------
# $Authors: Beifang Niu & Adam D Scott
# $Date: 2014-01-14 14:34:50 -0500 (Tue Jan 14 14:34:50 CST 2014) $
# $Revision:  $
# $URL: $
# $Doc: $ amino acid processing class  
#----------------------------------
#
use strict;
use warnings;
use Carp;
use TGI::Mutpro::Preprocess::Point;
#use PostData;
my $Debug = 0;
# Collection of points (atoms) in crystal structure
sub new {    
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self = {};

    $self->{CHAIN} = undef;    # ID of peptide chain
    $self->{POSITION} = undef; # Residue number in peptide chain
    $self->{AA} = undef;       # Name of amino acid (three letter code -- ARG, ASN, PRO, etc.)
    $self->{AVGPOINT} = undef; # Ref to point with average of X, Y, Z coordinates of all points
    $self->{POINTS} = ();      # Array of refs to Point objects (unordered Points -- i.e. atoms)
    $self->{AMBIGUOUS} = 0;    # Some PDB peptides have more than one residue at a given position. 
                               # Don't use these
	$self->{ACCEPTED} = { 	"ALA" => "A" , 
							"ARG" => "R" ,
							"ASN" => "N" ,
							"ASP" => "D" ,
							"CYS" => "C" ,
							"GLN" => "Q" ,
							"GLU" => "E" ,
							"GLY" => "G" ,
							"HIS" => "H" ,
							"ILE" => "I" ,
							"LEU" => "L" ,
							"LYS" => "K" ,
							"MET" => "M" ,
							"PHE" => "F" ,
							"PRO" => "P" ,
							"SER" => "S" ,
							"THR" => "T" ,
							"TRP" => "W" ,
							"TYR" => "Y" ,
							"VAL" => "V"
	};
    bless ($self, $class);

    return $self;
}

sub convertNameToSingle {
	my $this = shift;
	if ( @_ ) {
		my $aa = shift;
		if ( $this->isAA( $aa ) ) {
			return $this->{'ACCEPTED'}->{$aa};
		}
		warn "AminoAcid::convertNameToSingle - warnging: ".$aa." is not a three letter amino acid name\n";
		return $aa;
	}
	if ( $this->isAA( $this->name() ) ) {
		return $this->{'ACCEPTED'}->{$this->name()};
	}
	warn "AminoAcid::convertNameToSingle - warnging: ".$this->name()." is not a three letter amino acid name\n";
	return $this->name();
}

sub ambiguous {
    my $self = shift;
    if (@_) { $self->{AMBIGUOUS} = shift; }
    return $self->{AMBIGUOUS};
}
    
sub chain {
    my $self = shift;
    if (@_) { $self->{CHAIN} = shift; }
    return $self->{CHAIN};
}

sub position {
    # Residue number in peptide chain
    my $self = shift;
    if (@_) { $self->{POSITION} = shift; }
    return $self->{POSITION};
}

sub name {
    # Name of amino acid (three letter code)
    my $self = shift;
    if (@_) {
		#TODO: check if name is a real AA
		$self->{AA} = shift;
	}
    return $self->{AA};
}

sub addPoint {
    my $self = shift;
    my ($x, $y, $z) = @_;
    my $point = new TGI::Mutpro::Preprocess::Point;
    $point->xyz($x, $y, $z);
    push @{$self->{POINTS}}, \$point;
}

sub getPoints {
    # Return: array of refs to Point objects
    my $self = shift;
    return  @{$self->{POINTS}};
}

sub averagePoint {
    # Return: ref to Point object with average of X coordinates,
    # average of Y coordinates, 
    # average of Z coordinates of all points (atoms) 
    # belonging to this amino acid
    my $self = shift;
    if ( !defined $self->{AVGPOINT} ) {
		my ( $xTotal, $yTotal, $zTotal, $totalAtoms, $pointRef, $avgPoint, );
		foreach $pointRef ( $self->getPoints() ) {
			#PostData($pointRef); print "\n";
			my ($x,$y,$z) = $$pointRef->xyz();
			$xTotal += $x;
			$yTotal += $y;
			$zTotal += $z;
			$totalAtoms++;
		}
		$avgPoint = new TGI::Mutpro::Preprocess::Point;
		$avgPoint->xyz($xTotal/$totalAtoms, $yTotal/$totalAtoms, $zTotal/$totalAtoms);
		$self->{AVGPOINT} = \$avgPoint;
    }
    return $self->{AVGPOINT};
}

sub averageDistance {
    # Input: ref to AminoAcid object
    # Return: distance between this amino acid and 
    #         the input amino acid
    #         based on the average position of 
    #         all atoms in each amino acid
    my $self = shift;
    my $aaRef = shift;
    my ( $thisAvgPointRef, $thatAvgPointRef, $avgDistance, );
    $thisAvgPointRef = $self->averagePoint();
    $thatAvgPointRef = $$aaRef->averagePoint();
    $avgDistance = $$thisAvgPointRef->distance($thatAvgPointRef);
    return $avgDistance;
}

sub shortestDistance {
    # Input: ref to AminoAcid object
    # Return: shortest distance between any 
    #         two points in this amino acid
    #         and the input amino acid
    my $self = shift;
    my $aaRef = shift;
    my ( $distance, $shortestDistance, $thisPointRef, $thatPointRef  );
    $shortestDistance = 1e10;
    foreach $thisPointRef ( $self->getPoints() ) {
		foreach $thatPointRef ( $$aaRef->getPoints() ) {
			$distance = $$thisPointRef->distance($thatPointRef);
			if ( $distance < $shortestDistance ) { $shortestDistance = $distance; }
		}
    }
    ($shortestDistance < 1e10) || confess "Did not get any distance values";
    return $shortestDistance;
}

sub isHOH {
    my $this = shift;
	my $residue;
	if ( @_ ) { $residue = shift; } else { $residue = $this->name(); }
	if ( $residue eq "HOH" ) {
		return 1;
	}
	return 0;
}

sub isAA {
    my $this = shift;
	my $residue;
	if ( @_ ) { $residue = shift; } else { $residue = $this->name(); }
    if ( not exists $this->{ACCEPTED}->{$residue} ) {
        return 0;
    }
    return 1;
}

1;
