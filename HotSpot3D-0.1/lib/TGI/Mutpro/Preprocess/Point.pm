package TGI::Mutpro::Preprocess::Point;
#
#----------------------------------
# $Authors: Beifang Niu 
# $Date: 2014-01-14 14:34:50 -0500 (Tue Jan 14 14:34:50 CST 2014) $
# $Revision:  $
# $URL: $
# $Doc: $ do prioritization 
#----------------------------------
#
use strict;
use Carp;
# Point with X, Y, Z coordinates
# Used to calculate distance between atoms in PDB structure
sub new {    
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self = {};
    $self->{XYZ} = ();  # array of values in order x,y,z
                        # Values rounded to nearest thousandth
    bless ($self, $class);
    return $self;
}

sub xyz {
    # Get/Set array that holds X,Y,Z values
    # Returns an array
    my $self = shift;
    if (@_) { 
	foreach my $coord ( @_ ) {
            #print $coord."\n";
	    $coord = $self->round($coord);
	    push @{$self->{XYZ}}, $coord;
	}
    }
    return @{$self->{XYZ}}; 
}

sub distance {
    # Input: ref to a Point object
    # Return: distance between this point and input point rounded to nearest thousandths
    my $self = shift;
    my $pointRef = shift;
    my ($aX, $aY, $aZ) = $$pointRef->xyz();
    my ($bX, $bY, $bZ) = $self->xyz();
    my $distance = sqrt( ($aX-$bX)*($aX-$bX) + ($aY-$bY)*($aY-$bY) + ($aZ-$bZ)*($aZ-$bZ) );
    # Round off number to nearest 0.001
    $distance = $self->round($distance);
    return $distance;
}
  
sub samePoint {
    # Input: ref to a Point object
    # Return: 1 if the point has same coordinates as this one, 0 if not
    # This was changed from == comparison to 'eq' since seemingly equivalent
    # numbers were not comparing as expected.  Does Perl store as a float?
    my $self = shift;
    my $pointRef = shift;
    my ($aX, $aY, $aZ) = $$pointRef->xyz();
    my ($bX, $bY, $bZ) = $self->xyz();
    return ( $aX eq $bX && $aY eq $bY && $aZ eq $bZ );
}

sub round {    
    # Round off number to nearest 0.001
    my $self = shift;
    my $num = shift;
    $num += 0.0005;
    if ( $num =~ /(-?\d+\.\d{3})/ ) { $num = $1; } 
    return $num;
}
 
return 1;  
    
