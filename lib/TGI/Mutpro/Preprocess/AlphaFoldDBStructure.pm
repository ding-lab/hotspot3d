package TGI::Mutpro::Preprocess::AlphaFoldDBStructure;
#
#----------------------------------
# $Authors: Fernanda Martins Rodrigues (mrodrigues.fernanda@gmail.com; fernanda@wustl.edu)
# $Date: 2014-01-14 14:34:50 -0500 (Tue Jan 14 14:34:50 CST 2014) $
# $Revision: 2023-03-15 $
# $URL: $
# $Doc: $ peptide class; script based on the original PDBStructure.pm file in original hotspot3d code (by Beifang Niu and Amila Weerasinghe)
#----------------------------------
#
#
# Stores the text version of AlphaFold DB (v4) structure downloaded using 
#   'https://alphafold.ebi.ac.uk/files/AF-$alphafolddbId-F1-model_v4.pdb'
# It can also parse AlphaFold DB .pdb files stored locally
# NOTE: Latest AlphaFoldDB version (v4) is hardcoded in the script as of right now (03-15-2023) - will update script to handle different versions in the future
#
# Notes from AlphaFoldBD downloads page:
# "For downloading all predictions for a given species, use the download links below. Note that this option is only available on the desktop version of the site.
# The uncompressed archive files (.tar) contain all the available compressed PDB and mmCIF files (.gz) for a reference proteome. 
# In the case of proteins longer than 2700 amino acids (aa), AlphaFold provides 1400aa long, overlapping fragments. For example, 
# Titin has predicted fragment structures named as Q8WZ42-F1 (residues 1–1400), Q8WZ42-F2 (residues 201–1600), etc. These fragments 
# are currently only available for the human proteome in these proteome archive files, not on the website."
#
# From the above, there are 208 proteins longer 2700 aa in AlphaFold DB v4. 121 of these map to unreviewed Uniprot IDs. The remaining 87 are proteins 
# like Titin, which are not very relevant at the moment.
# Will currently only consider proteins smaller than 2700 aa.
# TO DO: add support for proteins >2700 aa - will need to modify functions which read alphafold.pdb files locally, since these structures need to be 
# downloaded (not available via url, as described above).


use strict;
use warnings;

use Carp;
use LWP::Simple;
use Data::Dumper;

use TGI::Mutpro::Preprocess::Peptide;

sub new {
    my ($class, $alphafolddbId, $alphafolddbDir, ) = @_;    
    bless {
        ID => $alphafolddbId,
        DIR => $alphafolddbDir,
        PAGE => ""
    }, $class;
}

sub alphafolddbId {
    my $self = shift;    
    if (@_) { $self->{ID} = shift; }
    return $self->{ID};
}

sub alphafolddbDir {
    my $self = shift;    
    if (@_) { $self->{DIR} = shift; }
    return $self->{DIR};
}

sub retrieveStructureFromAlphaFoldDB {
    # Download current record from AlphaFoldDB web site.
    my $self = shift;
    my $alphafolddbId = $self->alphafolddbId();
    # NOTE: Latest AlphaFoldDB version (v4) is hardcoded in the script as of right now (03-15-2023) - will update script to handle different versions in the future
    my $alphafolddbUrl = "https://alphafold.ebi.ac.uk/files/AF-$alphafolddbId-F1-model_v4.pdb"; 
    $self->{PAGE} = get($alphafolddbUrl);
}

sub localAlphaFoldDBFile {
    # This is when the AlphaFoldDB structure is from a local file - this is recommended
    my $self = shift;
    my $alphafolddbId = $self->alphafolddbId();
    my $alphafolddbDir = $self->alphafolddbDir();
    my ( $alphafolddbFile, $page, @entireFile, $line );
    # NOTE: Latest AlphaFoldDB version (v4) is hardcoded in the script as of right now (03-15-2023) - will update script to handle different versions in the future
    # TO DO: add support for dealing with proteins >2700 aa (this will have more than one .pdb file [AF-$alphafolddbId-F2-model_v4.pdb] per uniprot id/alphafold id)
    $alphafolddbFile = "$alphafolddbDir/AF-$alphafolddbId-F1-model_v4.pdb";
    unless( -e $alphafolddbFile ) { $self->{PAGE} = ""; return }
    open(IN, "< $alphafolddbFile") || warn "Could not open '$alphafolddbFile': $!";
    @entireFile = <IN>;
    close IN;
    foreach (@entireFile) { $page .= $_; }
    $self->{PAGE} = $page;
}

sub entireRecord {
    # Return: entire text record stored as a mulit-line page
    # Can be parsed by doing 
    #  foreach $line (split /\n/, $page) 
    #  { chomp $line; etc.. }
    my $self = shift;
    if ($self->{PAGE} eq "") { $self->localAlphaFoldDBFile(); }
    if ($self->{PAGE} eq "") { $self->retrieveStructureFromAlphaFoldDB(); }
    return $self->{PAGE};
}

sub offset {
    # NOTE: should have no offset (offset=0), as uniprot and alphafold db coordinates should be the same; just leaving it here not to change the code too much for now
    # TO DO: this will likely change when dealing with proteins longer than 2700aa - needs to be implemented still
    # Input: an array of chains
    # Return: Offset needed to convert crystal coordinates to 
    # Uniprot coordinates for different chains and regions
    #
    # Add 'offset' to crystal coordinate to get Uniprot coordinate
    # Subtract 'offset' from Uniprot coordinate to get cyrstal coordinate
    my $self = shift;
    my $chainToUniprotIdRef = shift;
    my @chainArray = keys %{$chainToUniprotIdRef};
    my $alphafolddbId = $self->{ID};
    my $offsetHash = {}; # a hash to store chain->afStart->afStop->offset corresponding to 
    #each DBREF entry (allowed multiple DBREF for a given chain, although this is not the case for AlphaFoldDB). But do this only for the chains already been cleared 
    #and stored in $chainToUniprotIdRef. 
    # If SEQADV is also present this hash will be modified to return offsets and corresponding chains and ranges - THIS IS NOT APPLICABLE TO ALPHAFOLD DB - SEQADV LINES NOT PRESENT
    my $infoHash = {}; # a hash to store chain->afStart->afStop->uStart->uStop
    #my $seqadvHash = {}; # a hash to store chain->SEQADV_event->SEQADV_comment # THIS IS NOT APPLICABLE TO ALPHAFOLD DB - SEQADV LINES NOT PRESENT
    #my $skipChains = {}; # a hash to store chains with REMARK 999 in their SEQADV line #- THIS IS NOT APPLICABLE TO ALPHAFOLD DB - SEQADV LINES NOT PRESENT
    #my $hasSEQADV = 0; # a bool to identify whether there's a SEQADV line
    my ($line, $offset,);
    foreach $line (split /\n/, $self->entireRecord()) {
        chomp $line;         
        #if ( $line =~ /^DBREF/ or $line =~ /^SEQADV/ ) { #- THIS IS NOT APPLICABLE TO ALPHAFOLD DB - SEQADV LINES NOT PRESENT
        if ( $line =~ /^DBREF/ ) {
            foreach my $chain ( @chainArray ) {
                # match line like this: (possible to have multiple lines for a chain, although this is only the case in PDB, not alphafold)
                # DBREF  XXXX A    1  1484  UNP    Q9BXF3   CECR2_HUMAN      1   1484
                if ( $line =~ /DBREF\s+\w+\s+$chain\s+(\-?\d+)\s+(\d+)\s+UNP\s+\S+\s+\S+\s+(\d+)\s+(\d+)/ ) {
                    my $afStart = $1; # af=alphafolddb
                    my $afStop = $2;
                    my $uStart = $3; # u=uniprot
                    my $uStop = $4;
                    my $offset = $uStart - $afStart;
                    $offsetHash->{$chain}->{$afStart}->{$afStop} = $offset;
                    $infoHash->{$chain}->{$afStart}->{$afStop}->{$uStart} = $uStop;
                }
                # THIS IS NOT APPLICABLE TO ALPHAFOLD DB - SEQADV LINES NOT PRESENT
                # match line like this:
                # SEQADV 3VW8     A       UNP  P08581    ALA  1057 DELETION
                # elsif ( $line =~ /SEQADV\s+\w+\s+$chain\s+UNP\s+\w+\s+\w+\s+(\d+)\s+(\w+)/ ) {
                #     my $position = $1;
                #     my $type = $2;
                #     if ( $type eq "DELETION" or $type eq "INSERTION" ) {
                #         $hasSEQADV = 1;
                #         $seqadvHash->{$chain}->{$position} = $type;
                #     }
                # }
                # # match line like this:
                # # SEQADV 1ZEO GLY A  201  UNP  P37231              SEE REMARK 999
                # elsif ( $line =~ /SEQADV.*REMARK/ ) {
                #     $hasSEQADV = 1;
                #     $skipChains->{$chain} = 0;
                # }
            }
        }
    }
    # THIS IS NOT APPLICABLE TO ALPHAFOLD DB - SEQADV LINES NOT PRESENT
    # do some clean up and make new regions depending on events reported at SEQADV lines
    # if ( $hasSEQADV ) {
    #     print STDOUT "TGI::Mutpro::Preprocess::AlphaFoldDBStructure::offset() ALERT: SEQADV line present for alphafolddbId = $alphafolddbId\n";
    #     # add chains like below to skip chains hash (there are deletions but delta(alphafolddb) = delta(uniprot))
    #     # DBREF  3TAD A  866  1193  UNP    O75334   LIPA2_HUMAN    866   1193
    #     # SEQADV 3TAD     A       UNP  O75334    PRO   976 DELETION 
    #     # SEQADV 3TAD     A       UNP  O75334    SER   977 DELETION
    #     # ...                      
    #     foreach my $chain ( keys %{$infoHash} ) {
    #         foreach my $afStart ( keys %{$infoHash->{$chain}} ) {
    #             foreach my $afStop ( keys %{$infoHash->{$chain}->{$afStart}} ) {
    #                 foreach my $uStart ( keys %{$infoHash->{$chain}->{$afStart}->{$afStop}} ) {
    #                     my $uStop = $infoHash->{$chain}->{$afStart}->{$afStop}->{$uStart};
    #                     if ( exists $seqadvHash->{$chain} ) { # if there are reported DELETIONs or INSERTIONs in this chain
    #                         my $deltaAlphaFoldDB = $afStop - $afStart;
    #                         my $deltaUNIPROT = $uStop - $uStart;
    #                         if ( $deltaAlphaFoldDB == $deltaUNIPROT ) {
    #                             $skipChains->{$chain} = 0; # add to the set of skip chains
    #                         }
    #                     }
    #                 }
    #             }
    #         }
    #     }
    #     # remove chains with "REMARK 999" and delta(alphafolddb) = delta(uniprot) even if SEQADV DEL,INS lines are present
    #     foreach my $chain ( keys %{$skipChains} ) {
    #         if ( defined $chain ) { # and exists $chainToUniprotIdRef->{$chain} 
    #             print STDOUT "TGI::Mutpro::Preprocess::AlphaFoldDBStructure::offset() WARNING:: $alphafolddbId chain $chain is deleted bc SEQADV=REMARK 999 and wrong sequence matches\n";
    #             delete $chainToUniprotIdRef->{$chain};
    #             delete $offsetHash->{$chain};
    #             delete $seqadvHash->{$chain};
    #         }
    #     }
    #     # make new regions and define offsets for them
    #     foreach my $chain ( keys %{$seqadvHash} ) {
    #         foreach my $position ( keys %{$seqadvHash->{$chain}} ) {
    #             my $type = $seqadvHash->{$chain}->{$position};
    #             foreach my $afStart ( keys %{$offsetHash->{$chain}} ) {
    #                 foreach my $afStop ( keys %{$offsetHash->{$chain}->{$afStart}} ) {
    #                     if ( $afStart <= $position and $position <= $afStop ) { # DEL or INT is inside the corresponding range
    #                         my $offset = $offsetHash->{$chain}->{$afStart}->{$afStop};
    #                         if ( $type eq "DELETION" ) {
    #                             delete $offsetHash->{$chain}->{$afStart}->{$afStop}; # replace this region by two regions below
    #                             if ( $afStart <= $position - 1 ) { # first region
    #                                 my $newPend = $position - 1;
    #                                 my $newOffset = $offset;
    #                                 $offsetHash->{$chain}->{$afStart}->{$newPend} = $newOffset;
    #                             }
    #                             # second region
    #                             my $newPstart = $position;
    #                             my $newOffset = $offset + 1;
    #                             $offsetHash->{$chain}->{$newPstart}->{$afStop} = $newOffset;
    #                         }
    #                         elsif ( $type eq "INSERTION" ) {
    #                             delete $offsetHash->{$chain}->{$afStart}->{$afStop}; # replace this region by two regions below
    #                             if ( $afStart <= $position - 1 ) { # first region
    #                                 my $newPend = $position - 1;
    #                                 my $newOffset = $offset;
    #                                 $offsetHash->{$chain}->{$afStart}->{$newPend} = $newOffset;
    #                             }
    #                             if ( $afStop >= $position + 1 ) { # second region
    #                                 my $newPstart = $position + 1;
    #                                 my $newOffset = $offset - 1;
    #                                 $offsetHash->{$chain}->{$newPstart}->{$afStop} = $newOffset;
    #                             }
    #                         }
    #                     }
    #                 }
    #             }
    #         }
    #     }
    # }
    return $offsetHash;
}


## THIS IS NOT APPLICABLE TO ALPHAFOLD DB - NUMMDL LINES NOT PRESENT; this function is also not called anywhere else ...
#  filtering the pdb files with too many NUMMDL modles
## Do not need this function when only one model
## is picked up
#
# sub modelsfilter {
#     # Return: the NUMMDL number if there is
#     my $self = shift;
#     my ($line, $models,);
#     $models = 0;
#     foreach $line (split /\n/, $self->entireRecord()) {
# 	chomp $line;        
# 	# match line like this:
#         # NUMMDL    640
#         if ($line =~ /NUMMDL\s+(\d+)/) { $models = $1; }
#     }

#     return $models;
# }

sub chainStartStop {
    # Input: chain id
    # Return: start and stop of the chain in crystal coordinates
    # Note: this is not ideal in PDB since a chain can have two different entries - not a problem in alphafoldDB
    # It just returns the first start and last stop 
    # (without worrying about any missing sequence in middle)
    my $self = shift;
    my $chain = shift;
    my ( $line, $start, $stop, );
    #print "Warning: need to fix 'sub chainStartStop' 
    #deal with multiple entries for a chain\n"; - this comes from pdb, it is not really the case for AlphaFoldDB
    my $alphafolddbId = $self->alphafolddbId();
    foreach $line (split /\n/, $self->entireRecord()) {
		chomp $line;        
		# match line like this:
		# DBREF  XXXX A    1  1484  UNP    Q9BXF3   CECR2_HUMAN      1   1484
        if ($line =~ /DBREF\s+\w+\s+$chain\s+(\-?\d+)\s+(\d+)\s+UNP\s+$alphafolddbId+\s+/) {
			my $tempStart = $1; 
			my $tempStop =$2;
			# Not ideal solution for multiple entries
			if (!defined $start || $tempStart < $start) { $start = $tempStart; }
			if (!defined $stop ||$tempStop > $stop) { $stop = $tempStop; }
		}
    }
    return ($start, $stop);
}

sub recordExists {
	my $self = shift;
	if ( $self->entireRecord() eq "" ) {
		warn "HotSpot3D::AlphaFoldDBStructure::chainToUniprotId warning: ".$self->alphafolddbId()." not found\n";
		return 0;
	}
	return 1;
}

sub chainToUniprotId {
    # Return: ref to hash with key = chain; value = Uniprot ID of peptide 
    # If there is more than one 'DBREF' entry for a given chain, 
    # preferentially uses the one with 'UNP' - this is the only option in AlphaFoldDB
    #
    # If there is more than one with 'UNP' and the Uniprot ID is different,
    # then do not use return anything
    #  
    my $self = shift;
    my ( $line,
         $chain,
         $uniprot,
         %chainToUniprot,
         %chainLength, 
         $chainId, 
         %ambiguous, );
    my $alphafolddbId = $self->alphafolddbId();
	if ( $self->recordExists() == 0 ) {
		return {};
	}
    foreach $line (split /\n/, $self->entireRecord()) {
    	# match line like this:
    	# DBREF  XXXX A    1  1484  UNP    Q9BXF3   CECR2_HUMAN      1   1484
    	# 'UNP' is for Uniprot.  There are also entries with 'GB', 'PIR', 'EMBL', 'SWS' and possibly others
        if ($line =~ /DBREF\s+\w+\s+(\w+)\s+\-?\d+\s+\d+\s+UNP\s+$alphafolddbId\s+\S+\s+\d+/) {
    	    $chain = $1;
    	    $uniprot = $alphafolddbId;
            # The code in if statement below does not apply to alphafoldDB, but leaving it here for now just in case I may have missed something 
            # TO DO: comment this whole function out for future versions, once results are curated
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
                carp "Ambiguous chain to Uniprot in $alphafolddbId: $chain $uniprot $chainToUniprot{$chain}";
                $ambiguous{$chain} = 1;
            }
            # If there is a non-Uniprot id associated with the chain, this will replace it
            $chainToUniprot{$chain} = $uniprot; 
        # Again, AlphaFoldDB will only have DBREF lines associated with a UNP id, but leaving it here for now in case I have missed something - won't hurt results
    	}elsif ($line =~ /DBREF\s+\w+\s+(\w+)\s+\-?\d+\s+\d+\s+(\w+)\s+$alphafolddbId\s+\S+\s+\d+/) {
    	    # This is not a Uniprot ID.
            # Make the ID with an underscore '_' so it can be 
            # distinguished from a real Uniprot ID
    	    $chain = $1;
    	    my $chainId = $2."_".$alphafolddbId;
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
    # ATOM      1  N   MET A   1      46.953  -1.972 -15.722  1.00 31.86           N

    #my ( $self, $drugport_ref, ) = @_;
    my ( $self ) = @_;
    my ( %peptideObjects, @entireFile, $line, @columns, $chain, $aminoAcid, $position, $x, $y, $z);
	if ( $self->recordExists() == 0 ) {
		return {};
	}
    foreach $line ( split /\n/, $self->entireRecord() ) {
		chomp $line;
		$aminoAcid = $position = $x = $y = $z = undef;
        ## NOTE: AlphaFoldDB only have one model per id, but leaving it here in case things change in future versions
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
        ## NOTE: Only lines starting with ATOM are present in AlphaFoldDB
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
    foreach $chain (keys %peptideObjects) {
		$peptideObjects{$chain}->removeAmbiguousAminoAcids();
	}

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
    # ATOM      1  N   MET A   1      46.953  -1.972 -15.722  1.00 31.86           N
    my $self = shift;
    my $chain = shift;
    my ( $line, $position, $residue, %sequence, );
    foreach $line ( split /\n/, $self->entireRecord() ) {
        chomp $line;
        ## NOTE: AlphaFoldDB only have one model per id, but leaving it here in case things change in future versions
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
        ## NOTE: Only lines starting with ATOM are present in AlphaFoldDB
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
    my $convertedAA = undef;
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
        length( $residue ) == 1 ? return "Z" : return $convertedAA;
    } else { return $convertedAA; carp "Unrecognized format for amino acid '$residue'"; }
}

sub printStructure {
    my $self = shift;
    foreach my $line (split /\n/, $self->entireRecord()) { chomp $line; print "$line\n"; }
}

sub by_number { $a <=> $b; }

return 1;