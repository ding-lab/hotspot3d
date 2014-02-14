package TGI::Mutpro::Preprocess::Prioritization;
#
#----------------------------------
# $Authors: Beifang Niu 
# $Date: 2014-01-14 14:34:50 -0500 (Tue Jan 14 14:34:50 CST 2014) $
# $Revision:  $
# $URL: $
# $Doc: $ statics related infor 
#----------------------------------
#
use strict;
use warnings;

use LWP::Simple;
use Cwd;
use Carp;
#
#  do prioritization based on given parameters
#  the parameters: 
#       1D distance cutoff; 
#       3D distance cutoff;
#
## Initiate based on uniprotfile and annotation file
#
sub new
{
    my ( $class, $uniprotProximityFile, $annotationFile, $outputFile ) = @_;    
    bless
    {
	UFILE    => $uniprotProximityFile,
	ANNOFILE => $annotationFile,
        OUTPUT   => $outputFile,
    }, $class;
}
# prioritization 
sub doPrioritization {
    my $self = shift;
    my $uniproGeneHashRef = shift;
    my $linearDistanceCutoff = shift;
    my $structureDistanceCutoff = shift;
    my $annoFile   = $self->{ANNOFILE};
    my $uniproFile = $self->{UFILE};
    # annotation hash
    my $annoHashRef = $self->parseAnnotation();
    # content after prioritization
    my @outputProximityFile;
    # process proximity file line by line
    open(IN, "< $uniproFile") || confess "Could not open '$uniproFile': $!";
    my @entireFile = <IN>;
    close IN;
    foreach my $a (@entireFile) {
        next if ( $a =~ /^WARNING/);
	chomp $a;
        # split values
        my ($u1, $chain1, $loc1, $offset1, $resi1, 
            $u2, $chain2, $loc2, $offset2, $resi2, 
            $distance, $protein) = split /\t/, $a;
        next if (($offset1 =~ /N\/A/) or ($offset2 =~ /N\/A/));
        next unless ( ($loc1 =~ /^\d+$/) and ($loc2 =~ /^\d+$/) 
                                         and ($offset1 =~ /^\d+$/) 
                                         and ($offset2 =~ /^\d+$/) );
        my $uniprotLocation1 = $loc1 + $offset1;
        my $uniprotLocation2 = $loc2 + $offset2;
        # unvalid coordinate. i am not sure I need to make sure that with John
        my $tempLinearDistance = abs( $uniprotLocation1 - $uniprotLocation2 );
        # linear distance filtering
        next if ( ($u1 eq $u2) and ($tempLinearDistance < $linearDistanceCutoff) );
        # 3d distance filtering 
        next if ($distance > $structureDistanceCutoff);
        my $tempLine      = "$u1\t";
        my $tempGenesList = "";
        print "$u1\n";
        my @tempGenes = keys %{$uniproGeneHashRef->{$u1}};
        foreach my $g (@tempGenes) { $tempGenesList .= "$g,"; }
        chop($tempGenesList);
        $tempGenesList = "N\/A" unless($tempGenesList);
        my $tempAnnotation;
        if ( defined($annoHashRef->{$uniprotLocation1}) ) {
            $tempAnnotation = $annoHashRef->{$uniprotLocation1};
        } else { $tempAnnotation = "N\/A"; }     
        $tempLine .= "$tempGenesList\t$chain1\t$loc1\t$offset1\t$resi1\t$tempAnnotation\t$u2\t";
        $tempGenesList = "";
        @tempGenes = keys %{$uniproGeneHashRef->{$u2}};
        foreach my $g (@tempGenes) { $tempGenesList .= "$g,"; }
        chop($tempGenesList);
        $tempGenesList = "N\/A" unless($tempGenesList);
        if ( defined($annoHashRef->{$uniprotLocation2}) ) { 
            $tempAnnotation = $annoHashRef->{$uniprotLocation2}; 
        } else { $tempAnnotation = "N\/A"; }     
        $tempLine .= "$tempGenesList\t$chain2\t$loc2\t$offset2\t$resi2\t$tempAnnotation\t$distance\t$protein\n";
        push( @outputProximityFile, $tempLine );
    }
    # pour out the results
    return \@outputProximityFile;
}

# parse annotation file
sub parseAnnotation
{
    my $self = shift;    
    my $annoFile = $self->{ANNOFILE};
    my %annotationHash;
    return \%annotationHash unless(-e $annoFile);
    open(IN, "< $annoFile") || confess "Could not open '$annoFile': $!";
    my @entireFile = <IN>;
    close IN;
    foreach (@entireFile) {
        chomp;
        my ( $start, $end, $region, $annotation ) = split /\t/, $_;
        next if ($start !~ /^\d+$/) or ($end !~ /^\d+$/);
        # notice: DISULFID
        if ($region =~ /DISULFID/) {
            $annotationHash{$start} = $annotation;
            $annotationHash{$end} = $annotation;
        } else { foreach my $a ($start..$end) { $annotationHash{$a} = $annotation; } }
    }
    # return annotation
    return \%annotationHash;
}

1;

