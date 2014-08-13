package TGI::Mutpro::Preprocess::PdbStructure;
#
#----------------------------------
# $Authors: Beifang Niu 
# $Date: 2014-01-14 14:34:50 -0500 (Tue Jan 14 14:34:50 CST 2014) $
# $Revision:  $
# $URL: $
# $Doc: $ peptide class 
#----------------------------------
#
#
# Stores the text version of PDB structure downloaded using 
#   'http://www.rcsb.org/pdb/files/$pdbID.pdb'
#
use strict;
use warnings;

use Carp;
use LWP::Simple;

use TGI::Mutpro::Preprocess::Peptide;

sub new {
    my ($class, $pdbId, $pdbDir, ) = @_;    
    bless {
        ID => $pdbId,
        DIR => $pdbDir,
        PAGE => ""
    }, $class;
}

sub pdbId {
    my $self = shift;    
    if (@_) { $self->{ID} = shift; }
    return $self->{ID};
}

sub pdbDir {
    my $self = shift;    
    if (@_) { $self->{DIR} = shift; }
    return $self->{DIR};
}

sub retrieveStructureFromPdb {
    # Download current record from PDB web site.
    my $self = shift;
    my $pdbId = $self->pdbId();
    my $pdbUrl = "http://www.rcsb.org/pdb/files/$pdbId.pdb"; 
    $self->{PAGE} = get($pdbUrl);
}

sub localPdbFile {
    # This is when the PDB structure is from a local file
    my $self = shift;
    my $pdbId = $self->pdbId();
    my $pdbDir = $self->pdbDir();
    my ( $pdbFile, $page, @entireFile, $line, );
    $pdbFile = "$pdbDir/$pdbId.pdb";
    #open(IN, "< $pdbFile") || confess "Could not open '$pdbFile': $!";
    open(IN, "< $pdbFile") || warn "Could not open '$pdbFile': $!";
    @entireFile = <IN>;
    close IN;
    foreach (@entireFile) { $page .= $_; }
    $self->{PAGE} = $page;
}

#sub localPdbFile {
#    # This is when the PDB structure is 
#    # from a local file
#    my $self = shift;
#    my $pdbFile = shift;
#    my ( $page, @entireFile,, $line, );
#    open(IN, "< $pdbFile") || confess "Could not open '$pdbFile': $!";
#    @entireFile = <IN>;
#    close IN;
#    foreach (@entireFile) { $page .= $_; }
#    $self->{PAGE} = $page;
#}

sub entireRecord {
    # Return: entire text record stored as a mulit-line page
    # Can be parsed by doing 
    #  foreach $line (split /\n/, $page) 
    #  { chomp $line; etc.. }
    my $self = shift;
    #if ($self->{PAGE} eq "") { $self->retrieveStructureFromPdb(); }
    if ($self->{PAGE} eq "") { $self->localPdbFile(); }
    return $self->{PAGE};
}

sub offset {
    # Input: chain
    # Return: Offset needed to convert crystal coordinates to 
    # Uniprot coordinates
    #
    # Add 'offset' to crystal coordinate to get Uniprot coordinate
    # Subtract 'offset' from Uniprot coordinate to get cyrstal coordinate
    my $self = shift;
    my $chain = shift;
    my ($line, $offset,);
    $offset = undef;
    foreach $line (split /\n/, $self->entireRecord()) {
	chomp $line;        
	# match line like this:
        # DBREF  2AR5 A    2   113  UNP    Q14CQ9   Q14CQ9_HUMAN  1421   1532
        if ($line =~ /DBREF\s+\w+\s+$chain\s+(\-?\d+)\s+\d+\s+UNP\s+\S+\s+\S+\s+(\d+)/) {
	    # Don't use offset if it is ambiguous
	    if (defined $offset && $offset != $2-$1) { return undef; }
	    $offset = $2-$1;
        }
    }
    return $offset;
}

##
#  filtering the pdb files with too many NUMMDL modles
## Do not need this function when only one model
## is picked up
#
sub modelsfilter {
    # Return: the NUMMDL number if there is
    my $self = shift;
    my ($line, $models,);
    $models = 0;
    foreach $line (split /\n/, $self->entireRecord()) {
	chomp $line;        
	# match line like this:
        # NUMMDL    640
        if ($line =~ /NUMMDL\s+(\d+)/) { $models = $1; }
    }

    return $models;
}

sub chainStartStop {
    # Input: chain id
    # Return: start and stop of the chain in crystal coordinates
    # Note: this is not ideal since a chain can have two different entries.
    # It just returns the first start and last stop 
    # (without worrying about any missing sequence in middle)
    my $self = shift;
    my $chain = shift;
    my ( $line, $start, $stop, );
    #print "Warning: need to fix 'sub chainStartStop' 
    #deal with multiple entries for a chain\n";
    my $pdbId = $self->pdbId();
    foreach $line (split /\n/, $self->entireRecord()) {
	chomp $line;        
	# match line like this:
	# DBREF  2H32 A    1   126  UNP    P12018   VPREB_HUMAN     20    145
	#
	# DBREF  2FO0 A    1    14  UNP    P00519-2 ABL1_HUMAN       1     14 
	# DBREF  2FO0 A   57   531  UNP    P00519-2 ABL1_HUMAN      57    531 
        if ($line =~ /DBREF\s+$pdbId\s+$chain\s+(\-?\d+)\s+(\d+)\s+/) {
	    my $tempStart = $1; 
            my $tempStop =$2;
	    # Not ideal solution for multiple entries
	    if (!defined $start || $tempStart < $start) { $start = $tempStart; }
	    if (!defined $stop ||$tempStop > $stop) { $stop = $tempStop; }
	}
    }
    return ($start, $stop);
}

sub chainToUniprotId {
    # Return: ref to hash with key = chain; value = Uniprot ID of peptide 
    # If there is more than one 'DBREF' entry for a given chain, 
    # preferentially uses the one with 'UNP'
    #
    # If there is more than one with 'UNP' and the Uniprot ID is different,
    # then do not use return anything
    #
    # If there is more than one non-UNP entry, 
    # doesn't matter which one is used
    #
    # In some cases the Uniprot ID of the chain is not known. 
    # They will give a Genebank ID or something else
    #  
    my $self = shift;
    my ( $line,
         $chain,
         $uniprot,
         %chainToUniprot,
         %chainLength, 
         $chainId, 
         %ambiguous, );
    my $pdbId = $self->pdbId();
    foreach $line (split /\n/, $self->entireRecord()) {
	# match line like this:
	# DBREF  2H32 A    1   126  UNP    P12018   VPREB_HUMAN     20    145
	# 'UNP' is for Uniprot.  There are also entries with 'GB', 'PIR', 'EMBL', 'SWS' and possibly others
        if ($line =~ /DBREF\s+$pdbId\s+(\w+)\s+\-?\d+\s+\d+\s+UNP\s+(\S+)\s+\S+\s+\d+/) {
	    $chain = $1;
	    $uniprot = $2;
	    if (defined $chainToUniprot{$chain} && $chainToUniprot{$chain} ne $uniprot && $chainToUniprot{$chain} !~ /\_/) {
		# This can be ambiguous when the peptide chain is constructed in vitro from
		# two or more peptides
		# (e.g. 1LBK "Crystal structure of a recombinant glutathione transferase, created by replacing 
		#       the last seven residues of each subunit of the human class pi isoenzyme with the 
		#       additional C-terminal helix of human class")
		# DBREF  1LBK A    1   202  UNP    P09211   GTP_HUMAN        1    202             
		# DBREF  1LBK B    1   202  UNP    P09211   GTP_HUMAN        1    202             
		# DBREF  1LBK A  203   208  UNP    P08263   GTA1_HUMAN     208    213             
		# DBREF  1LBK B  203   208  UNP    P08263   GTA1_HUMAN     208    213     
		carp "Ambiguous chain to Uniprot in $pdbId: $chain $uniprot $chainToUniprot{$chain}";
		$ambiguous{$chain} = 1;
	    }
	    # If there is a non-Uniprot id associated with the chain, this will replace it
	    $chainToUniprot{$chain} = $uniprot; 
	}elsif ($line =~ /DBREF\s+$pdbId\s+(\w+)\s+\-?\d+\s+\d+\s+(\w+)\s+(\S+)\s+\S+\s+\d+/) {
	    # This is not a Uniprot ID.
            # Make the ID with an underscore '_' so it can be 
            # distinguished from a real Uniprot ID
	    $chain = $1;
	    my $chainId = $2."_".$3;
	    # Do not replace a real Uniprot ID with a non-Uniprot ID
	    if (!defined $chainToUniprot{$chain}) { $chainToUniprot{$chain} = $chainId; }
	}
    }
    # Remove ambiguous chains
    foreach $chain (keys %ambiguous) { delete $chainToUniprot{$chain}; }

    return \%chainToUniprot;
}

sub makePeptides {
    # Return: Ref to hash of Peptide objects
    # Creates Point, AminoAcid, and Peptide objects for each chain in the crystal
    # structure
    
    # The entries are in specific columns.  Whitespace is not relevant
    # 1 - 6 Record name "ATOM "
    # 7 - 11 Integer serial Atom serial number.
    # 13 - 16 Atom name Atom name.
    # 17 Character altLoc Alternate location indicator.
    # 18 - 20 Residue name resName Residue name.
    # 22 Character chainID Chain identifier.
    # 23 - 26 Integer resSeq Residue sequence number.
    # 27 AChar iCode Code for insertion of residues.
    # 31 - 38 Real(8.3) x Orthogonal coordinates for X in Angstroms
    # 39 - 46 Real(8.3) y Orthogonal coordinates for Y in Angstroms
    # 47 - 54 Real(8.3) z Orthogonal coordinates for Z in Angstroms
    # match lines like these:
    # 12345678901234567890123456789012345678901234567890123456789012345678901234567890
    # ATOM   7749  CD1 ILE A 999      51.224  54.016  84.367  1.00 83.16           C
    # ATOM   7750  N   ASN A1000      55.338  57.542  83.643  1.00 80.67           N

    my ( $self, $drugport_ref, ) = @_;
    my ( %peptideObjects, @entireFile, $line, @columns, $chain, $aminoAcid, $position, $x, $y, $z);
    foreach $line ( split /\n/, $self->entireRecord() ) {
	chomp $line;
	$aminoAcid = $position = $x = $y = $z = undef;
        ## Beifang Niu 05.14.2013
        # only grab one model if there are multiple models
        # exit when there is a line: 
        #
        # ATOM   1963  HD3 PRO A 123     -13.866   5.943  20.640  1.00  0.00           H
        # TER    1964      PRO A 123
        # ENDMDL 
        # MODEL        2
        #
        # Add code for drug information processing in PDB structure
        # TODO: choose only lines with HET_GROUP
        #
        last if ( $line =~ /^ENDMDL/ );
	if ($line !~ /^ATOM|^HETATM/) { next; }
	my @columns = split //, $line;
	# Get chain letter and make new Peptide object if necessary
	$chain = $columns[21];
	# If there is no chain, then give default value of 'Z'
	if ( $chain eq " " ) { $chain = "Z"; }
	if ( !defined $peptideObjects{$chain} ) {
	    $peptideObjects{$chain} = new TGI::Mutpro::Preprocess::Peptide;
	    $peptideObjects{$chain}->name($chain);
	}
	# Amino acid name
	foreach (17..19) { $aminoAcid .= $columns[$_]; }
	# Amino acid position
	#foreach (22..25) { $position .= $columns[$_]; } 
	foreach (22..29) { $position .= $columns[$_]; } 	   
	$position =~ s/\s+//g;
	# X,Y,Z coordinates of atom in space
	foreach (30..37) { $x .= $columns[$_]; }   
	foreach (38..45) { $y .= $columns[$_]; } 
	foreach (46..53) { $z .= $columns[$_]; } 
	# Add point to amino acid in peptide
	$peptideObjects{$chain}->addPointToAminoAcid($aminoAcid, $position, $x, $y, $z);
    }
    # Remove ambiguous amino acids from each peptide 
    # (these are peptide positions with more than one amino acid name)
    foreach $chain (keys %peptideObjects) { $peptideObjects{$chain}->removeAmbiguousAminoAcids(); }

    return \%peptideObjects;
}


sub residuesNearPosition {
    # Input: chain, position that defined an amino acid position
    # (in crystal coordinates) 
    #        $distance in angstroms
    #
    # Return: ref to hash with 
    # key = $chain $position within $distance of residue 
    # value = ref to  AminoAcid object
    # Positions are in crystal coordinates; 
    # not Uniprot coordinates
    my ( $self, 
         $residueChain, 
         $residuePosition, 
         $distance ) = @_;
    my ( $peptideRef,
         $inputAminoAcidRef,
         $allAaObjRef,
         $aaObjRef,
         $position,
         %closeAminoAcids,
         $chain, );
    my $debug = 0; 
    # Get Peptide objects to describe crystal structure
    $peptideRef = $self->makePeptides();
    # Get AminoAcid object for residue in chain '$residueChain', 
    # at position $position
    $inputAminoAcidRef = $$peptideRef{$residueChain}->getAminoAcidObject($residuePosition);
    # If the position is not in the peptide, return undef
    if (!defined $inputAminoAcidRef) { return undef; } 
    foreach $chain (keys %{$peptideRef}) {
	$allAaObjRef = $$peptideRef{$chain}->getAllAminoAcidObjects();
	foreach $position ( sort by_number keys %{$allAaObjRef} ) {
	    # Get this particular AminoAcid object
 	    $aaObjRef = $$peptideRef{$chain}->getAminoAcidObject($position);
	    if ( $debug ) { print "$chain$position  ", $$aaObjRef->minDistance($inputAminoAcidRef); }
	    if ( $$aaObjRef->minDistance($inputAminoAcidRef) <= $distance ) {
		# This amino acid residue is within $distance
		$closeAminoAcids{"$chain$position"} = $aaObjRef;
		if ( $debug ) { print "  <== added"; }
	    }
	    if ( $debug ) { print "\n"; }
	}
    }

    return \%closeAminoAcids;
}

sub chainSequence {
    # Input: chain
    # Return: ref to hash with key = position of residue in cyrstal, 
    # value = single-letter amino acid code

    # The entries are in specific columns.  Whitespace is not relevant
    # The Chain is in column 21 (starting with col 0)
    # The position is in columns 22 through 25
    # match lines like these:
    # 12345678901234567890123456789012345678901234567890123456789012345678901234567890
    # ATOM   7749  CD1 ILE A 999      51.224  54.016  84.367  1.00 83.16           C
    # ATOM   7750  N   ASN A1000      55.338  57.542  83.643  1.00 80.67           N
    # ATOM    126  CD1 ILE A -15      41.428  88.478 110.141  1.00100.92           C
    # ATOM    127  N   GLY A   8      86.873  61.296 112.806  1.00 87.60           N
    my $self = shift;
    my $chain = shift;
    my ( $line, $position, $residue, %sequence, );
    foreach $line ( split /\n/, $self->entireRecord() ) {
	chomp $line;
        ## Beifang Niu 05.14.2013
        # only grab one model if there are multiple models
        # exit when there is a line: 
        #
        # ATOM   1963  HD3 PRO A 123     -13.866   5.943  20.640  1.00  0.00           H
        # TER    1964      PRO A 123
        # ENDMDL 
        # MODEL        2
        #
        last if ( $line =~ /^ENDMDL/ );
	if ( $line !~ /^ATOM|^HETATM/ ) {
        #if ( $line =~ /^ATOM/ ) {
	    $position = "";
            $residue = "";
	    my @columns = split //, $line;
	    if ( $columns[21] eq $chain ) { 
		# Amino acid residue
		foreach (17..19) { $residue .= $columns[$_]; }
		$residue =~ s/\s+//g;
		$residue = convertAA($residue);
		# Position
		foreach (22..29) { $position .= $columns[$_]; }
		$position =~ s/\s+//g;
		if ( $position > 0 ) { $sequence{$position} = $residue; }
	    } # matches 'if ( $columns[21] eq $chain )'
	} # matches 'if ( $line =~ /^ATOM|^HETATM/ )'
    } # matches 'foreach $line'

    return \%sequence;
}


sub convertAA {
    # Convert single amino acid code to triplet
    # or triplet to single amino acid code
    my $residue = $_[0];
    $residue = uc $residue;
    my %oneToThree =  (
		       A => 'ALA',
		       R => 'ARG',
		       N => 'ASN',
		       D => 'ASP',
		       C => 'CYS',
		       Q => 'GLN',
		       E => 'GLU',
		       G => 'GLY',
		       H => 'HIS',
		       I => 'ILE',
		       L => 'LEU',
		       K => 'LYS',
		       M => 'MET',
		       F => 'PHE',
		       P => 'PRO',
		       S => 'SER',
		       T => 'THR',
		       W => 'TRP',
		       Y => 'TYR',
		       V => 'VAL',
		       );
    my %threeToOne;
    foreach (keys %oneToThree) { $threeToOne{$oneToThree{$_}} = $_; }
    if ( defined $oneToThree{$residue} ) {
	return $oneToThree{$residue};
    }elsif ( defined $threeToOne{$residue} ) {
	return $threeToOne{$residue};
    }elsif ( length( $residue ) <= 3 ) {
        length( $residue ) == 1 ? return "Z" : return $residue;
    } else { carp "Unrecognized format for amino acid '$residue'"; }
}

sub printStructure {
    my $self = shift;
    foreach my $line (split /\n/, $self->entireRecord()) { chomp $line; print "$line\n"; }
}

sub by_number { $a <=> $b; }

return 1;

