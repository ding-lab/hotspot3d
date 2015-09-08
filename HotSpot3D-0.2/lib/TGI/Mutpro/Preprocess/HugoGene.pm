package TGI::Mutpro::Preprocess::HugoGene;
#
#----------------------------------
# $Authors: Beifang Niu 
# $Date: 2014-01-14 14:34:50 -0500 (Tue Jan 14 14:34:50 CST 2014) $
# $Revision:  $
# $URL: $
# $Doc: $ hogo gene processing another class 
#----------------------------------
#
use strict;
use Carp;
# Used to hold data from Hugo web site
sub new {
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self = {};
    # Unique values
    $self->{ID} = undef;             # Unique ID number (HGNC ID after removing 'HGNC:')
    $self->{SYMBOL} = undef;         # HUGO symbol (gene name)
    $self->{NAME} = undef;           # Gene description (string with white space)
    $self->{TYPE} = undef;           # Locus type ('pseudogene', 'gene with protein product, inferred', etc.)
    $self->{CHR} = undef;            # Chromosome assignment (includes band)
    $self->{UNIPROT} = undef;        # Uniprot name (mapped data supplied by UniProt)
    $self->{ENSEMBL} = undef;        # Ensembl name (mapped data supplied by Ensembl)
    $self->{OMIM} = undef;           # Online Mendelian Inheritance in Man name (mapped data supplied by NCBI)
    $self->{UCSC} = undef;           # UCSC ID (mapped data supplied by UCSC)
    $self->{START} = undef;          # Nucleotide start on chr (not supplied by Hugo)
    $self->{STOP} = undef;           # Nucleotide end on chr (not supplied by Hugo) 
    # Not necessarily unique
    $self->{ENTREZ} = ();            # Hash of Entrez (LocusId) name(s) 
    $self->{REFSEQ} = ();            # Hash of Ref Seq name(s) 
    $self->{VEGA} = ();              # Hash of Vega name(s) 
    $self->{PREVIOUS_SYMBOLS} = ();  # Hash of previous Hugo symbols 
                                     # (previous symbol may associate with more than one unique Hugo ID/symbol)
    $self->{ALIAS} = ();             # Hash of aliases 
                                     # (an alias may associate with more than one unique Hugo ID/symbol) 
    $self->{OLD_NAME} = ();          # Combination of 'Previous Names' and 'Name Aliases' (essentially descriptions)
    bless ($self, $class);
    return $self;
}

### Set/Get routines for all unique variables
sub id {
    # Don't want whitespace
    my $self = shift;
    if (@_)  { 
        $self->{ID} = shift; 
        $self->{ID} =~ s/\s+//g;
    }
    return $self->{ID};
}

sub symbol {
    # Don't want whitespace
    my $self = shift;
    if (@_) { 
        $self->{SYMBOL} = shift;
        $self->{SYMBOL} =~ s/\s+//g;
    }
    return $self->{SYMBOL};
}
    
sub name {
    # Whitespace expected
    my $self = shift;
    if (@_) { $self->{NAME} = shift; }
    return $self->{NAME};
}
  
sub type {
    # Whitespace expected
    my $self = shift;
    if (@_) { $self->{TYPE} = shift; }
    return $self->{TYPE};
}

sub chr {
    # Don't want whitespace
    my $self = shift;
    if (@_) { 
	$self->{CHR} = shift; 
	$self->{CHR} =~ s/\s+//g;  
    }
    return $self->{CHR};
}

sub uniprot { 
    # Don't want whitespace
    my $self = shift;
    if (@_) { 
        $self->{UNIPROT} = shift; 
        $self->{UNIPROT} =~ s/\s+//g;
    }
    return $self->{UNIPROT};
}
    
sub ensembl {
    # Don't want whitespace
    my $self = shift;
    if (@_)  { 
        $self->{ENSEMBL} = shift; 
        $self->{ENSEMBL} =~ s/\s+//g;
    }

    return $self->{ENSEMBL};
}
      
sub omim {
    # Don't want whitespace
    my $self = shift;
    if (@_) { 
        $self->{OMIM} = shift; 
        $self->{OMIM} =~ s/\s+//g;
    }
    return $self->{OMIM};
}   

sub ucsc {
    # Don't want whitespace
    my $self = shift;
    if (@_) { 
        $self->{UCSC} = shift;
        $self->{UCSC} =~ s/\s+//g;
    }
    return $self->{UCSC};
}

sub start {
    # Don't want whitespace
    my $self = shift;
    if (@_) { 
        $self->{START} = shift;
        $self->{START} =~ s/\s+//g;
    }
    return $self->{START};
}
     
sub stop {
    # Don't want whitespace
    my $self = shift;
    if (@_) { 
        $self->{STOP} = shift;
        $self->{STOP} =~ s/\s+//g;  
    }
    return $self->{STOP};
}

### Set routines for variables stored in hash

sub addEntrezName {
    my $self = shift;
    my $id = shift;
    $id =~ s/\s+//g;
    ${$self->{ENTREZ}}{$id} = 1;
}

sub addRefSeqName {
    my $self = shift;
    my $id = shift;
    $id =~ s/\s+//g;
    ${$self->{REFSEQ}}{$id} = 1;
}

sub addVegaName {
    my $self = shift;
    my $id = shift;
    $id =~ s/\s+//g;
    ${$self->{VEGA}}{$id} = 1;
}


sub addPreviousSymbol {
    my $self = shift;
    my $id = shift;
    $id =~ s/\s+//g;
    ${$self->{PREVIOUS_SYMBOLS}}{$id} = 1;
}

sub addAlias {
    my $self = shift;
    my $id = shift;
    $id =~ s/\s+//g;
    ${$self->{ALIAS}}{$id} = 1;
}

sub addOldDescription {
    my $self = shift;
    my $id = shift;
    ${$self->{OLD_NAME}}{$id} = 1;
}


## Get routines for variables stored in hash
# All return ref to hash 

sub getAllOldDescriptions {
    # Return ref to hash
    my $self = shift;
    return \%{$self->{OLD_NAME}};
}
    

sub getAllEntrezNames {
    # Return ref to hash
    my $self = shift;
    return \%{$self->{ENTREZ}};
}
 
sub getAllRefSeqNames {
    # Return ref to hash
    my $self = shift;
    return \%{$self->{REFSEQ}};
}
 
sub getAllVegaNames {
    # Return ref to hash
    my $self = shift;
    return \%{$self->{VEGA}};
}

sub getAllPreviousSymbols {
    # Return ref to hash
    my $self = shift;
    return \%{$self->{PREVIOUS_SYMBOLS}};
}


sub getAllAliases {
    # Return ref to hash
    my $self = shift;
    return \%{$self->{ALIAS}};
}

## Methods to see if this Hugo symbol refers to the same
## gene as identified 
sub isPreviousSymbol {
    # Input: Hugo name
    # Return: 1 if the input name matches one of the
    #         previous symbols of this object
    my ($self, $hugoName) = @_;
    $hugoName =~ s/\s+//g;
    foreach ( keys %{$self->getAllPreviousSymbols()} ) {
	if ( $_ eq $hugoName ) { return 1; }
    }

    return 0;
}


sub isAlias {
    # Input: Hugo name
    # Return: 1 if the input name matches one of the
    #         alias symbols of this object
    my ($self, $hugoName) = @_;
    $hugoName =~ s/\s+//g;
    foreach ( keys %{$self->getAllAliases()} ) {
	if ( $_ eq $hugoName ) { return 1; }
    }

    return 0;
}

sub hasVega {
    # Input: Vega gene ID
    # Return: 1 if this object has the Vega gene ID
    my ($self, $id) = @_;
    $id =~ s/\s+//g;
    foreach ( keys %{$self->getAllVegaNames()} )  {
	if ( $_ eq $id ) { return 1; }
    }

    return 0;
}

sub hasEntrez {
    # Input: Entrez gene ID
    # Return: 1 if this object has the Entrez gene ID
    my ($self, $id) = @_;
    $id =~ s/\s+//g;
    foreach ( keys %{$self->getAllEntrezNames()} ) {
	if ( $_ eq $id ) { return 1; }
    }

    return 0;

}

sub hasRefSeq {
    # Input: RefSeq gene ID
    # Return: 1 if this object has the RefSeq gene ID
    my ($self, $id) = @_;
    $id =~ s/\s+//g;
    foreach ( keys %{$self->getAllRefSeqNames()} ) {
	if ( $_ eq $id ) { return 1; }
    }

    return 0;
}

sub overlaps {
    # Input: start, stop, tolerance (optional)
    # Return: 1 if the start, stop overlap the start, stop of this gene
    #         within 'tolerance' base pairs
    # Returns undef if this gene does not have start, stop
    # ASSUMES: we are talking about the same chromosome.  Does not check 
    # chromosome.  
    my ( $self, $start, $stop, $tolerance ) = @_;
    if ( !defined $tolerance ) { $tolerance = 0; }
    if ( !defined $self->{START} || !defined $self->{STOP} ) {
	carp "WARNING: $self->{SYMBOL} does not have start, stop.  Can not determine overlap \n";
	return undef;
    }
    # Make sure start < stop for both
    my $thisStart = $self->start();
    my $thisStop = $self->stop();
    if ( $thisStart > $thisStop ) { 
        ($thisStart, $thisStop) = ($thisStop, $thisStart);
    }
    if ( $start > $stop ) { ($start, $stop) = ($stop, $start); }
    return ( $stop + $tolerance > $thisStart && $start - $tolerance < $thisStop );
}

sub sameMajorChrBand {
    # Input: chromosome (in any format)
    # Return: 1 if this object is on the same major chromosome band
    # i.e. It does not look at the number after the decimal point
    # 7p14.20 ~ 7p14 ~ 7p14.2  etc.
    # If one band is 'ter' or 'cen', it will not necessarily return the
    # correct answer since the 'ter' band can have a standard name
    my ( $self, $inputChr ) = @_;
    my ( $chr, $arm, $band, $inputArm, $inputBand );
    $chr = $self->chr();
    # This should not be called if $chr is not defined
    (defined $chr) || carp $self->symbol(), " does not have a defined chromosome";
    # Remove common junk
    $chr =~ s/chr//i; $inputChr =~ s/chr//i;
    $chr =~ s/\_//;  $inputChr =~ s/\_//;
    $chr =~ s/\s+//g; $inputChr =~ s/\s+//g;
    # Change to uppercase in case chromosomes are X or Y
    $chr = uc $chr;
    $inputChr = uc $inputChr;
    # Return 0 if the chromosomes are different
    my ($entireChr, $entireInputChr );
    if ( $chr =~ /(X|Y|\d+)(p|q|cen|ter)?/i ) { $entireChr = $1; }
    if ( $inputChr =~ /(X|Y|\d+)(p|q|cen|ter)?/i ) { $entireInputChr = $1; }
    if ( !defined $entireChr || !defined $entireInputChr )
    {
	carp "Unexpected format for \$inputChr: '$inputChr' sent to $self->{SYMBOL} or \$chr unexpected format";
	return 0;
    }	    
    if ( $entireChr ne $entireInputChr ) { return 0; }

    # One or both bands might be expressed as a range.

    #### maybe make a helper function to look at two bands......
    # so it can be called with potentially 4 different combinations....

    # They might match now
    if ( $chr eq $inputChr ) { return 1; }

    # If one of the chromosomes has cen, pter, or qter,
    # just go by chromosome number
    if ( $chr =~ /ter/i || $inputChr =~ /ter/i ||
	 $chr =~ /cen/i || $inputChr =~ /cen/i ) {
	if ( $chr =~ /(\w+)(p|q|cen)/i ) { $chr = $1; }
	if ( $inputChr =~ /(\w+)(p|q|cen)/i ) { $inputChr = $1; }
	return ( $inputChr eq $chr );
    }
    
    # Chromosomes can be in format '$chr[pq]\d+.\d+'
    # or can be given as a range 
    if ( $chr =~ /(\w+)([pq])(\d*)?/i ) {
	$chr = $1;
        $arm = $2;
        $band = $3;
    }
    if ( $inputChr =~ /(\w+)([pq])(\d*)?/i ) {
	$inputChr = $1;
        $inputArm = $2;
        $inputBand = $3;
    }
       
    # Return 0 if chromosomes are different
    if ( $inputChr ne $chr ) { return 0; }

    # If one of the chromosomes does not have an arm,
    # it is OK to return 1 since the chromosomes are the same
    if ( !defined $arm 
         || !defined $inputArm 
         || $arm eq "" 
         || $inputArm eq "" ) { return 1; }
    
    # If the arms are different, return 0
    if ( $arm ne $inputArm ) { return 0; }
    
    # The chromosomes and arms are the same, allow the
    # bands to differ by 1
    if ( defined $band && defined $inputBand && $band ne "" && $inputBand ne "") {
	return ( abs($band - $inputBand) <= 1 );
    }

    # The chromosomes are the same, one or more bands 
    # are not defined

    return 1;
}

sub sameChr {

    # Input: chromosome (in any format)
    # Return: 1 if this object is on the same chromosome

    my ( $self, $inputChr ) = @_;

    my $chr = $self->chr();

    # This should not be called if $chr is not defined
    (defined $chr)
    ||
    carp $self->symbol(), " does not have a defined chromosome";
  
    # Remove common junk
    $chr =~ s/chr//i; $inputChr =~ s/chr//i;
    $chr =~ s/\_//;  $inputChr =~ s/\_//;
    $chr =~ s/\s+//g; $inputChr =~ s/\s+//g;

    # Change to uppercase in case chromosomes are X or Y
    $chr = uc $chr;
    $inputChr = uc $inputChr;

    if ( $inputChr =~ /^(X|Y|\d+)/ ) {
	$inputChr = $1;
    } else {
	carp "Unexpected format for \$inputChr: '$inputChr' sent to $self->symbol{SYMBOL}";
    }
    if ( $chr =~ /^(X|Y|\d+)/ ) {
	$chr = $1;
    } else {
	carp "Unexpected format for chromosome \$chr: '$chr' in object $self->symbol{SYMBOL}, $self->symbol{CHR}";
    }
    
    return ( $chr eq $inputChr );
}

sub displayObject {
    # For debugging.  Print some of the contents to STDOUT
    my $self = shift;
    
    my ( $idRef );

    print $self->symbol();
    if ( defined $self->id() && $self->id() ne "" ) { 
        print "\tID: ", $self->id();
    } 
    if ( defined $self->ensembl() && $self->ensembl() ne "" ) { 
        print "\tEnsembl: ", $self->ensembl();
    } 
    $idRef = $self->getAllVegaNames();
    if ( defined $idRef && scalar(keys %{$idRef}) > 0 ) {
	print "\tVega: ";
	foreach ( keys %{$idRef} ) { print "$_ "; }
    }
    $idRef = $self->getAllEntrezNames();
    if ( defined $idRef && scalar(keys %{$idRef}) > 0 ) {
	print "\tEntrez: ";
	foreach ( keys %{$idRef} ) { print "$_ "; }
    }
    $idRef = $self->getAllPreviousSymbols();
    if ( defined $idRef && scalar(keys %{$idRef}) > 0 ) {
	print "\tPrevious symbols: ";
	foreach ( keys %{$idRef} ) { print "$_ "; }
    }

    $idRef = $self->getAllAliases();
    if ( defined $idRef && scalar(keys %{$idRef}) > 0 ) {
	print "\tAliases: ";
	foreach ( keys %{$idRef} ) { print "$_ "; }
    }
    if ( defined $self->chr() && $self->chr() ne "" ) { 
        print "\tChr: ", $self->chr();
    } 

}

return 1;

