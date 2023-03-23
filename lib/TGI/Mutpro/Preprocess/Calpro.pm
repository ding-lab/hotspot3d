package TGI::Mutpro::Preprocess::Calpro;
#
#----------------------------------
# $Authors: Beifang Niu and Amila Weerasinghe
# $Modified by: Fernanda Martins Rodrigues @WashU (fernanda@wustl.edu; mrodrigues.fernanda@gmail.com)
# $Date: 2014-01-14 14:34:50 -0500 (Tue Jan 14 14:34:50 CST 2014) $
# $Revision: 2023-03-15 $
# $URL: $
# $Doc: $ calculate proximity file for one uniprotid (used by first step)
#----------------------------------
#
#package Calpro;
#
#  Input: Uniprot id
#  Output: 1) File with structure-based proximity data
#          2) File with list of structures where coordinates 
#             can be mapped to Uniprot coordinates
# 
#----------------------------------
#

use strict;
use warnings;

use Carp;
use Cwd;
use IO::File;
use FileHandle;
use File::Copy;
use LWP::Simple;
use Getopt::Long;

use TGI::Mutpro::Preprocess::Uniprot;
# MODIFIED ON 03-15-2023: indicate alpha fold script instead of pdb
use TGI::Mutpro::Preprocess::AlphaFoldDBStructure;
use TGI::Mutpro::Preprocess::HugoGeneMethods;

my $SHORTESTDISTANCE = "shortest";
my $AVGDISTANCE = "average";

sub new {
	my $class = shift;
	my $this = {};
	$this->{'output_dir'} = getcwd;
	$this->{'max_3d_dis'} = 100;
	$this->{'min_seq_dis'} = 0;
	$this->{'uniprot_id'} = undef;
	$this->{'restrictedaa_pairs'} = 0;
	$this->{'uniprot_ref'} = undef;
	$this->{'stat'} = undef;
	# MODIFIED ON 03-15-2023: changed pdb to alphafold dir
	$this->{'alphafolddb_file_dir'} = undef;
	$this->{'distance_measure'} = $AVGDISTANCE;
## add drug port database file (08.04.2014)
#$this->{'drugport_file'} = undef;
	bless $this, $class;
	$this->process();
	return $this;
	}

	sub process {
	my $this = shift;
	my ( $help, $options );
	unless( @ARGV ) { die $this->help_text(); }
	$options = GetOptions (
		'3d-distance-cutoff=i'    => \$this->{'max_3d_dis'},
		'linear-cutoff=i'   => \$this->{'min_seq_dis'},
		'output-dir=s'    => \$this->{'output_dir'},
		'uniprot-id=s'    => \$this->{'uniprot_id'},
		# MODIFIED ON 03-15-2023: changed pdb to alphafold dir
		'alphafolddb-file-dir=s'  => \$this->{'alphafolddb_file_dir'},
		'measure=s'  => \$this->{'distance_measure'},
		#'drugport-file=s' => \$this->{'drugport_file'},

		'help' => \$help,
	);
	if ($help) { print STDERR help_text(); exit 0; }
	unless($options) { die $this->help_text(); }
	unless(defined $this->{'uniprot_id'}) { warn 'You must provide a Uniprot ID !', "\n"; die $this->help_text(); }
	unless( $this->{'output_dir'} and (-e $this->{'output_dir'} ) ) { warn 'You must provide a output directory ! ', "\n"; die $this->help_text(); }
	# MODIFIED ON 03-15-2023: changed pdb to alphafold dir
	unless( $this->{'alphafolddb_file_dir'} and (-e $this->{'alphafolddb_file_dir'}) ) { warn 'You must provide a AlphaFold DB file directory ! ', "\n"; die $this->help_text(); }
	if ( $this->{'distance_measure'} ne $SHORTESTDISTANCE and $this->{'distance_measure'} ne $AVGDISTANCE ) {
		warn "HotSpot3D::Calpro warning: measure not recognized, resetting to default = average\n";
		$this->{'distance_measure'} = $AVGDISTANCE;
	}
	#unless( $this->{'drugport_file'} and (-e $this->{'drugport_file'}) ) { warn 'You must provide a drugport database file ! ', "\n"; die $this->help_text(); }
	#### processing ####
	# This can be used if only want to write out pairs when one is in a Uniprot-annotated domain
	# Default is to write out every pair that is within $MaxDistance 
	# angstroms (provided they are > $PrimarySequenceDistance
	# amino acids away in the primary sequence)

	## processing drugport database file ( 08-05-2014 )
	#my $drugport_infor_hash_ref = $this->get_drugport_database_info( $this->{'drugport_file'} );
	$this->{'uniprot_ref'} = TGI::Mutpro::Preprocess::Uniprot->new( $this->{'uniprot_id'} );
	my $Dir = "$this->{'output_dir'}/proximityFiles";
	my $ProximityFile = "$Dir/inProgress/$this->{'uniprot_id'}.ProximityFile.csv";
	# MODIFIED ON 03-15-2023: changed alphafolddb information
	# get the linkage information between AlphaFoldDB ids if one AlphaFoldDB ids includes two or more Uniprot 
	# IDs, it should not be removed and also, sometimes, two Uniprot IDs locate in two
	# similar molecules
	my $hugoUniprot = "$this->{'output_dir'}/hugo.uniprot.alphafolddb.csv";
	my $linkage = $this->getLinkageInfo( $hugoUniprot );
	# Make proximity file to a temporary directory so current 
	# file can be used until this is finished 
	$this->writeProximityFile( $ProximityFile, $linkage );
	#$this->writeProximityFile( $ProximityFile, $linkage, $drugport_infor_hash_ref, );
	# Write a file that says if the amino acid sequence in Uniprot
	# is consistent with the sequence in the AlphaFold DB structure(s)
	my $AlphFoldDBCoordinateFile = "$Dir/alphafolddbCoordinateFiles/$this->{'uniprot_id'}.coord";
	$this->checkOffsets($ProximityFile, $AlphFoldDBCoordinateFile);
	# Now move the file from the $dir/inProgress subdirectory
	# to the root dir
	move($ProximityFile, "$Dir/$this->{'uniprot_id'}.ProximityFile.csv");

	return 1;
}

# MODIFIED ON 03-15-2023: changed function to parse alphafolddb instead of pdb
## Get proteins involved multiple uniprots ##
sub getLinkageInfo {
	my ( $this, $hugof, ) = @_;
	my $hugofh = new FileHandle;
	unless( $hugofh->open( $hugof ) ) { die "Could not open hugo uniprot file !\n" };
	my ( %multiUn, %temph, );
	foreach ( $hugofh->getlines ) {
		chomp; my ( undef, $uniprotId, $alphafolddb, ) = split /\t/;
		# Only use Uniprot IDs with AlphaFold DB structures
		next if ( $alphafolddb eq "N/A" || $uniprotId !~ /\w+/ );
		map{ if ( (defined $temph{$_}) and ($temph{$_} ne $uniprotId) ) { $multiUn{$_} = 1; } else { $temph{$_} = $uniprotId; } } split /\s+/, $alphafolddb;
	}
	return \%multiUn;
}

# MODIFIED ON 03-15-2023: changed function to parse alphafolddb instead of pdb
sub writeProximityFile {
	my ( $this, $file, $ulink, ) = @_;
	#my ( $this, $file, $ulink, $drugport_ref, ) = @_;
	# $ulink linkage infor of unirpot IDs
	my $fh = new FileHandle;
	unless( $fh->open( ">$file" ) ) { die "Could not open proximity file to write !\n" };
	print STDOUT "Creating ".$file."\n";
	my ( $uniprotDomainRef,
		$alphafolddbRef, %alphafolddbIds, 
		#%allOffsets,
		$structureRef, $peptideRef, $chainToUniprotIdRef, 
		$uniprotChain,
		$otherChainOffset, 
		$residuePositionsRef,
		$residuePosition,  
		$uniprotAminoAcidRef,
		$otherChainUniprotId,
		$allAaObjRef, 
		$skip, 
		$position,  
		$distanceBetweenResidues, 
		$aaObjRef,
		$otherDomainRef, 
		$uniprotChainOffset,
		$uniprotAaName, 
		$otherAaName, 
		$correctedPosition, 
		$uniprotIdToDomainsRef, 
		$hugoGeneRef, $hugoId, );
	# Get Uniprot-defined domains for $this->{'uniprot_id'}
	### 121105  Don't make this dependent on a specific type of annotation.
	#           Keep distances of all pairs even if the two amino acids are 
	#           not within annotated region
	#           Should ignore amino acids that would be found using the primary 
	#           sequence-based proximity analysis
	$uniprotDomainRef = uniprotDomains( $this->{'uniprot_id'} ) if ( $this->{'restrictedaa_pairs'} );
	# Get all AlphaFoldDB structures
	$alphafolddbRef = $this->{'uniprot_ref'}->annotations("AlphaFoldDB");
	%alphafolddbIds = ();
	if ( !defined $alphafolddbRef || scalar(@{$alphafolddbRef}) == 0 ) {
		carp "Did not get AlphaFoldDB IDs for $this->{'uniprot_id'} \n";
		return;
	}
	# Note: new added function
	# Filtering some alphafolddb files here to avoid unnecessary heavy load from CPU and Memory
	# Found this problem by testing all uniprot IDs initiation
	my $alphafolddbfilterRef = $this->filteringAlphaFoldDB( $alphafolddbRef, $ulink );
	map{
		print STDOUT "$_\n";
		$alphafolddbIds{$1} = 1 if ( $_ =~ /^(\w+)\;/ );
	} @{$alphafolddbfilterRef};
	# Download and parse AlphaFoldDB files
	$fh->print( "UniProt_ID1\tChain1\tPosition1\tOffset1\tResidue_Name1\t" );
	$fh->print( "UniProt_ID2\tChain2\tPosition2\tOffset2\tResidue_Name2\t" );
	$fh->print( "Distance\tAlphaFoldDB_ID\n" );
	foreach my $alphafolddbId (keys %alphafolddbIds) {
		my $allOffsets = {};
		$structureRef = TGI::Mutpro::Preprocess::AlphaFoldDBStructure->new( $alphafolddbId, $this->{'alphafolddb_file_dir'} );
		## don't need this part when only one model be picked up
		# Get Peptide objects to describe crystal structure
		$peptideRef = $structureRef->makePeptides();
		# Get chain representing given Uniprot ID.
		# Choose the chain that is the longest.
		$chainToUniprotIdRef = $structureRef->chainToUniprotId();
		$uniprotChain = undef; 
		my $chainLength = 0;
		foreach my $chain ( keys %{$chainToUniprotIdRef} ) { # choose the longest chain (to be used as the base chain)
			#print STDERR $chain."\n";
			next if ( $$chainToUniprotIdRef{$chain} ne $this->{'uniprot_id'} );
			my ( $chainStart, $chainStop, ) = $structureRef->chainStartStop( $chain );
			if ( $chainStop - $chainStart + 1 > $chainLength ) { 
				$uniprotChain = $chain;
				print STDOUT "$alphafolddbId\t$uniprotChain\t";
				$chainLength = $chainStop - $chainStart + 1; 
				print STDOUT $chainLength."\n";
			}
		}
		unless ( defined $uniprotChain ) {
			print $fh "WARNING: Did not get chain for '$this->{'uniprot_id'}' in '$alphafolddbId'.";
			print $fh "  Skipping structure\n";
			next; 
		}
		#  120905  Don't tie this to Uniprot annotation
		# Get all domains for all of the Uniprot IDs in this structure
		#print STDERR "RestrictAminoAcidPairs: $this->{'restrictedaa_pairs'}\n"; 
		if ( $this->{'restrictedaa_pairs'} ) {
			$uniprotIdToDomainsRef = getDomainsForAllUniprotIds( $chainToUniprotIdRef );
		}
		# Get all offsets needed to convert crystal coordinates 
		# to Uniprot coordinates
		# Add 'offset' to crystal coordinate to get Uniprot coordinate
		####
		# map{
		# 	my $offset = $structureRef->offset($_);
		# 	$offset = "N/A" if ( !defined $offset || $offset !~ /\-?\d+$/ );
		# 	#print STDERR "$_\t offset: $offset\n";
		# 	$allOffsets{$_} = $offset;
		# } keys %{$chainToUniprotIdRef};
		###
		$allOffsets = $structureRef->offset( $chainToUniprotIdRef );
		# Get offset needed to convert the crystal coordinates to Uniprot coordinates
		# for the chain corresponding to $UniprotId
		
		#print STDERR "uniprotChain:  $uniprotChain\t  uniprotChainOffset: $uniprotChainOffset\n";
		# Get position numbers (in crystal coordinates) of all residues in the peptide 
		# chain resolved in the structure
		unless ( defined $$peptideRef{$uniprotChain} ) {
			print STDERR "\$\$peptideRef{\$uniprotChain} not defined for \$uniprotChain = $uniprotChain.";
			print STDERR "  \$UniprotId = '$this->{'uniprot_id'}' in '$alphafolddbId'. Skipping structure\n";
			next; 
		}
		$residuePositionsRef = $$peptideRef{$uniprotChain}->aminoAcidPositionNumbers();
		# Go through every position and see if it is close to an annotated domain 
		# (but not in that domain)
		# or if it is close to another peptide chain in the crystal that is a different 
		# protein that '$UniprotId'
		###
		# Updated 120905.  Record any pair of amino acids that are within $MaxDistance angstroms,
		#                  but have amino acid positions > $PrimarySequenceDistance away
		foreach $residuePosition ( @{$residuePositionsRef} ) {
			#print STDERR "residuePosition: $residuePosition\n";
			# Get AminoAcid object for residue in chain '$uniprotChain', 
			# at position $position
			$uniprotAminoAcidRef = $$peptideRef{$uniprotChain}->getAminoAcidObject( $residuePosition );
			next if ( $$uniprotAminoAcidRef->isHOH() == 1 );
			$uniprotAaName = $$uniprotAminoAcidRef->name();
			# 	Updated	170510 : use a hash with chain and regions for retrieving the offset
			$uniprotChainOffset = getOffset( $allOffsets, $uniprotChain, $residuePosition );
			#      120905  Don't tie this to Uniprot annotation
			# See if this is in an annotated Uniprot domain.  
			# If so, write it out.
			if ( $this->{'restrictedaa_pairs'} ) {
				#print STDERR "RestrictAminoAcidPairs\n";
				if ( defined $$uniprotDomainRef{$residuePosition+$uniprotChainOffset} ) {
					foreach my $domain ( keys %{$$uniprotDomainRef{$residuePosition + $uniprotChainOffset}} ) {
					$fh->print( "$this->{'uniprot_id'}\t[$uniprotChain]\t$residuePosition\t" );
					$fh->print( "$uniprotChainOffset\t$uniprotAaName\t" );
					$fh->print( "$this->{'uniprot_id'}\t[$uniprotChain]\t$residuePosition\t" );
					$fh->print( "$uniprotChainOffset\t$uniprotAaName\t" );
					$fh->print( "$domain\t0\t$alphafolddbId\n" );
					}
				}
			}
			# Now compare this amino acid to all other amino acids in all of the peptide 
			# chain(s) in the crystal
			foreach my $chain ( keys %{$peptideRef} ) {
				# Get the UniprotId of this peptide chain
				$otherChainUniprotId = $$chainToUniprotIdRef{$chain};
				#  Updated 121105.  Following is not necessary if not 
				#  restricting pairs written out
				#  Get the domains of $otherChainUniprotId 
				if ( $this->{'restrictedaa_pairs'} ) {
					$otherDomainRef = $$uniprotIdToDomainsRef{ $otherChainUniprotId };
				}
				# Get the offset needed to convert crystal coordinates of this 
				# peptide chain to Uniprot coordinates 
				
				$allAaObjRef = $$peptideRef{$chain}->getAllAminoAcidObjects();
				## remove the positions with insertion code
				my @tmp_array_positions = grep{ /\d+$/ } keys %{$allAaObjRef}; 
				foreach $position ( sort {$a<=>$b} @tmp_array_positions ) {
					$otherChainOffset = getOffset( $allOffsets, $chain, $position );
					$aaObjRef = $$peptideRef{$chain}->getAminoAcidObject($position);
					next if ( $$aaObjRef->isHOH() == 1 );
					if ( (defined $otherChainOffset) and ($otherChainOffset eq "N/A") ) { 
						$correctedPosition = $position;
					} else { 
						if ( defined $otherChainOffset ) {
							$correctedPosition = $position + $otherChainOffset;
						} 
					};
					# Skip if the amino acid at '$position' of peptide chain '$chain' 
					# is not close to the amino acid at '$residuePosition' 
					# of peptide chain '$uniprotChain'
					if ( $this->{'distance_measure'} eq $SHORTESTDISTANCE ) {
						$distanceBetweenResidues = $$aaObjRef->shortestDistance($uniprotAminoAcidRef);
					} elsif ( $this->{'distance_measure'} eq $AVGDISTANCE ) {
						$distanceBetweenResidues = $$aaObjRef->averageDistance($uniprotAminoAcidRef);
					} else {
						$distanceBetweenResidues = $$aaObjRef->averageDistance($uniprotAminoAcidRef);
					}
					if ( $distanceBetweenResidues > $this->{'max_3d_dis'} ) { next; }
					# Also skip if the two amino acids are in the same chain and 
					# within <= $PrimarySequenceDistance
					# residues of each other
					# If two amino acids are close to each other in the primary sequence, 
					# don't record them.
					# They will be detected by the proximity analysis in other application
					next if ( $chain eq $uniprotChain && abs($position - $residuePosition) <= $this->{'min_seq_dis'} );
					$otherAaName = $$aaObjRef->name();
					unless( defined $otherChainUniprotId ){ $otherChainUniprotId = "N/A"; };
					unless( defined $otherChainOffset )   { $otherChainOffset    = "N/A"; };
					if ( ! $this->{'restrictedaa_pairs'} ) {
						$fh->print( "$this->{'uniprot_id'}\t[$uniprotChain]\t$residuePosition\t" );
						$fh->print( "$uniprotChainOffset\t$uniprotAaName\t" );
						$fh->print( "$otherChainUniprotId\t[$chain]\t$position\t" );
						$fh->print( "$otherChainOffset\t$otherAaName\t" );
						$fh->print( "$distanceBetweenResidues\t$alphafolddbId\n" );
					}
					##### This is just to restrict what is written out to pairs in which 
					##### one of the residues is in an annotated domain
					if ( $this->{'restrictedaa_pairs'} ) {
						# OK. The amino acid at '$residuePosition' of peptide chain '$uniprotChain' 
						# is close to the amino acid 
						# at $position of $chain.  If '$chain' and '$uniprotChain' 
						# represent different proteins, then
						# print it out.  This is an interaction site between different proteins
						# It doesn't matter if there is an annotated domain at $position 
						# of $chain, but print it out if there is  
						if ( $otherChainUniprotId ne $this->{'uniprot_id'} ) {
							# Initialize to no annoated domains at $position
							my @otherChainDomains = ();
										push @otherChainDomains, "-";
							if ( defined $$otherDomainRef{$correctedPosition} ) {
								@otherChainDomains = keys %{$$otherDomainRef{$correctedPosition}};
							}
							foreach ( @otherChainDomains ) {
								$fh->print( "$this->{'uniprot_id'}\t[$uniprotChain]\t$residuePosition\t" );
								$fh->print( "$uniprotChainOffset\t$uniprotAaName\t" );
								$fh->print( "$otherChainUniprotId\t[$chain]\t$position\t" );
								$fh->print( "$otherChainOffset\t$otherAaName\t$_\t" );
								$fh->print( "$distanceBetweenResidues\t$alphafolddbId\n");
							}
							next;
						}
						# If we are here, the two chains represent the same protein
						# (or the two chains are the same)
						# Only print out if '$position' is in an annotated domain 
						# AND '$residuePosition' is not in the same domain 
						# (since that domain has already been noted)
						# Skip if this is not an annotated domain in $UniprotId 
						next if ( !defined $$otherDomainRef{$correctedPosition} );
						foreach my $domain ( keys %{$$otherDomainRef{$correctedPosition}} ) { 
							# Skip if the domain is the same as 
							next if ( defined $$uniprotDomainRef{$residuePosition+$uniprotChainOffset} && defined $$uniprotDomainRef{$residuePosition+$uniprotChainOffset}{$domain} );
							$fh->print( "$this->{'uniprot_id'}\t[$uniprotChain]\t$residuePosition\t" );
							$fh->print( "$uniprotChainOffset\t$uniprotAaName\t" );
							$fh->print( "$otherChainUniprotId\t[$chain]\t$position\t" );
							$fh->print( "$otherChainOffset\t$otherAaName\t" );
							$fh->print( "$domain\t$distanceBetweenResidues\t$alphafolddbId\n" );
						}
					} #if restrictedaa_pairs
				} #foreach other residue
			} #foreach chain
		} #foreach residue position ref
	} #foreach alphafolddb id
	return 1;
}

sub getOffset {
	my ( $allOffsets, $chain ,$residuePosition ) = @_;
	my $offset = undef;
	foreach my $afStart ( keys %{$allOffsets->{$chain}} ) {
		foreach my $afStop ( keys %{$allOffsets->{$chain}->{$afStart}} ) {
			if ( $afStart <= $residuePosition and $residuePosition <= $afStop ) {
				$offset = $allOffsets->{$chain}->{$afStart}->{$afStop};
				last;
			}
		}
		if ( defined $offset ) {
			last;
		}
	}
	$offset = "N/A" if ( !defined $offset );
	return $offset;
}

sub getDomainsForAllUniprotIds {
	# Input: ref to hash with key = chain;
	#        value = Uniprot ID
	# Return: ref to hash with key = Uniprot ID; 
	#        value ref to hash of domains
	my ( $this, $chainToUniprotIdRef, ) = @_;
	my ( $uniprotId, %uniprotIdToDomains, $uniprotDomainRef, );
	foreach ( keys %{$chainToUniprotIdRef} ) {
		$uniprotId = $$chainToUniprotIdRef{$_};
		if ( !defined $uniprotId ) { next; }
		$uniprotDomainRef = $this->uniprotDomains($uniprotId);
		$uniprotIdToDomains{$uniprotId} = $uniprotDomainRef;
	}
	return \%uniprotIdToDomains;
}

sub uniprotDomains {
	# Input: Uniprot ID
	# Return: ref to hash '$uniprotDomains{$position}{"$key: $description"}'
	my ( $this, $uniprotId, ) = @_;
	my ( %uniprotDomains, $start, $stop, $key, $description, $position, );  
	my $uniprotRef = TGI::Mutpro::Preprocess::Uniprot->new( $uniprotId ); 
	# Check to see if $uniprotId returned a valid Uniprot record
	my $uniprotRecord = $uniprotRef->entireRecord();
	return \%uniprotDomains unless ( defined $uniprotRecord );
	my @recordLines = split /\n/, $uniprotRecord;
	if ( !defined $uniprotRecord || scalar(@recordLines) <= 10 ) {  return \%uniprotDomains; }
	# Get all Uniprot-defined domains.
	# Ref to array from 'push @domains, "$key\t($dmStart, $dmStop)\t$desc";'
	# Need length of protein
	my $proteinLength = length( $uniprotRef->sequence() );
	my $domainRef = $uniprotRef->domains( 1, $proteinLength );
	foreach my $entry ( @{$domainRef} ) {
		if ( $entry =~ /(\w+)\s+\((\d+)\,\s+(\d+)\)\s+(.*)\.?$/ ){
			 ( $key, $start, $stop, $description ) = ($1, $2, $3, $4);
		} else {
			print STDERR "WARNING: Could not parse domain description for '$uniprotId': '$entry'\n";
		}
		if ( $start > $stop) {
			print STDERR "WARNING: Error parsing domain for '$uniprotId'. Start ($start) > Stop ($stop) in '$entry'\n";
		}
		foreach $position ( $start..$stop ) {
			$uniprotDomains{$position}{"$key: $description"} = 1;
		}
	}
	return \%uniprotDomains;
}

sub checkOffsets {
	my ( $this, $proximityFile, $coordFile, ) = @_;
	my ( $line, 
		 $uniprotA,
		 $positionA, 
		 $offsetA, 
		 $aminoAcidA, 
		 $uniprotB, 
		 $positionB, 
		 $offsetB, 
		 $aminoAcidB, 
		 $alphafolddbId,
		 $uniprot, 
		 $uniprotRef, 
		 $uniprotSequenceRef, 
		 $position, 
		 %alphafolddbUniprotPosition,  );
	my $profh = new FileHandle;
	unless( $profh->open( "< $proximityFile" ) ) {  die "Could not open proximity file $proximityFile to read !\n"  };
	my @entireFile = <$profh>;
	$profh->close();
	my $coorfh = new FileHandle;
	unless( $coorfh->open( "> $coordFile" ) ) {  die "Could not open coordinate file $coordFile to write !\n" };
	$coorfh->print( "AlphaFoldDB_ID\tUniProt_ID\tErrors\tTotal_Position\tFraction_Errors\n" );
	print STDOUT "Creating ".$coordFile."\n";
	foreach $line ( @entireFile ) {
		chomp $line;
		next if ( $line =~ /WARNING/ );
		( $uniprotA, undef, $positionA, $offsetA, $aminoAcidA, $uniprotB, undef, $positionB, $offsetB, $aminoAcidB ) = split /\t/, $line;
		if ( $line =~ /(\S+)\s*$/ ) { $alphafolddbId = $1; }
			# print STDERR "Unexpected format for \$uniprotA ($uniprotA) in $line.  Skipping. \n"; }
		next if ( $uniprotA !~ /^\w+$/ ); 
		next if ( $uniprotB !~ /^\w+$/ || $offsetA !~ /^-?\d+$/ || $offsetB !~ /^-?\d+$/ || $positionA !~ /^-?\d+$/ || $positionB !~ /^-?\d+$/ );
		#if not an AA, convertAA returns original value if length aminoAcid <= 3, 'Z' if length == 1, undef otherwise
		$aminoAcidA = TGI::Mutpro::Preprocess::AlphaFoldDBStructure::convertAA( $aminoAcidA );
		$aminoAcidB = TGI::Mutpro::Preprocess::AlphaFoldDBStructure::convertAA( $aminoAcidB );
		next if ( !defined $aminoAcidA & !defined $aminoAcidB );
		#next unless ( TGI::Mutpro::Preprocess::AminoAcid::checkAA( $aminoAcidA ) and TGI::Mutpro::Preprocess::AminoAcid::checkAA( $aminoAcidB )
		if ( defined $alphafolddbUniprotPosition{$alphafolddbId}{$uniprotA}{$positionA+$offsetA} && $alphafolddbUniprotPosition{$alphafolddbId}{$uniprotA}{$positionA+$offsetA} ne $aminoAcidA ) {
			print $coorfh "Inconsistent amino acids for $uniprotA position $positionA+$offsetA in $alphafolddbId: '$alphafolddbUniprotPosition{$alphafolddbId}{$uniprotA}{$positionA+$offsetA}' and $aminoAcidA \n";
		}
		$alphafolddbUniprotPosition{$alphafolddbId}{$uniprotA}{$positionA+$offsetA} = $aminoAcidA;
		if ( defined $alphafolddbUniprotPosition{$alphafolddbId}{$uniprotB}{$positionB+$offsetB} && $alphafolddbUniprotPosition{$alphafolddbId}{$uniprotB}{$positionB+$offsetB} ne $aminoAcidB ) {
			print $coorfh "Inconsistent amino acids for $uniprotB position $positionB+$offsetB in $alphafolddbId: '$alphafolddbUniprotPosition{$alphafolddbId}{$uniprotB}{$positionB+$offsetB}' and $aminoAcidB \n";
		}
		$alphafolddbUniprotPosition{$alphafolddbId}{$uniprotB}{$positionB+$offsetB} = $aminoAcidB;
	} #foreach line of entireFile
	my %alphafolddbUniprotErrorCount;
	foreach $alphafolddbId ( keys %alphafolddbUniprotPosition ) {
		foreach $uniprot ( keys %{$alphafolddbUniprotPosition{$alphafolddbId}} ) {
			$uniprotRef = TGI::Mutpro::Preprocess::Uniprot->new($uniprot);
			next if ( !defined $uniprotRef );
			$uniprotSequenceRef = $this->getUniprotSeq( $uniprot );
			$alphafolddbUniprotErrorCount{$alphafolddbId}{$uniprot} = 0;
			foreach $position ( sort {$a<=>$b} keys %{$alphafolddbUniprotPosition{$alphafolddbId}{$uniprot}} ) {
				#check if position exists in sequence and it has correct amino acid
				if ( !defined $$uniprotSequenceRef{$position} || $$uniprotSequenceRef{$position} ne $alphafolddbUniprotPosition{$alphafolddbId}{$uniprot}{$position} ) {
					$alphafolddbUniprotErrorCount{$alphafolddbId}{$uniprot}++;
				}
			}
		}
	} #foreach alphafolddbId
	foreach $alphafolddbId ( keys %alphafolddbUniprotErrorCount ) {
		foreach $uniprot ( keys %{$alphafolddbUniprotErrorCount{$alphafolddbId}} ) {
			print $coorfh "$alphafolddbId \t $uniprot \t errors: $alphafolddbUniprotErrorCount{$alphafolddbId}{$uniprot} \t total: ";
			print $coorfh scalar(keys %{$alphafolddbUniprotPosition{$alphafolddbId}{$uniprot}}), "\t";
			my $fraction = $alphafolddbUniprotErrorCount{$alphafolddbId}{$uniprot}/scalar(keys %{$alphafolddbUniprotPosition{$alphafolddbId}{$uniprot}});
			if ( $fraction != 0 && $fraction != 1 ) {
				$fraction += 0.005;
				if ( $fraction =~ /(0\.\d{2})/ ) { $fraction = $1; }
			}
			print $coorfh "$fraction\n";
		}
	}
	$coorfh->close();

	return 1;
}
		
# Extract uniprot sequence
sub getUniprotSeq {
	my ( $this, $uniprot, ) = @_;
	my %seq;
	my $uniprotRef = TGI::Mutpro::Preprocess::Uniprot->new($uniprot);
	my $sequence = $uniprotRef->sequence();
	if ( !defined $sequence ) { return \%seq; }
	my @seqarray = split //, $sequence;
	map{  $seq{$_+1} = $seqarray[$_];  } (0..$#seqarray);

	return \%seq;
}

# MODIFIED ON 03-15-2023: changed function to parse alphafolddb instead of pdb
# These types not really available in AlphaFold since that databse is based on predictions - will likely remove this function once results are curated
# TO DO: remove function once we make sure we don't need it?
# Note: throw away some AlphaFoldDBs  
# Filtering some AlphaFoldDB files here to avoid unnecessary heavy load from CPU and Memory
# Found this problem by testing all uniprot IDs initiation
sub filteringAlphaFoldDB {
	my ( $this, $entrysRef, $ulinks, ) = @_;
	my ( @entrysafterFiltered, $NMRs, $Xrays, $neutron, $other, %tmph, $total, );
	$NMRs = $Xrays = $neutron = $other = $total = 0;
	foreach my $a ( @{$entrysRef} ) {
		my ( $alphafolddbtype ) = $a =~ /^\w+;\s+(.*?);\s+/;
		SWITCH:{
			$alphafolddbtype eq 'X-ray'   && do { $Xrays++;   last SWITCH; };
			$alphafolddbtype eq 'NMR'     && do { $NMRs++;    last SWITCH; };
			$alphafolddbtype eq 'Neutron' && do { $neutron++; last SWITCH; };
			$other++;
		}
	}
	## new added in order to retrieve more AlphaFoldDBs
	$total = $NMRs + $Xrays + $neutron + $other;
	map{ push(@entrysafterFiltered, $_); } @{$entrysRef}; 
	return \@entrysafterFiltered;
	#if ( $total < 50 ) {  map{ push(@entrysafterFiltered, $_); } @{$entrysRef}; return \@entrysafterFiltered; }
	# only Neutron and X-ray
	if ( ($Xrays > 0) || ($neutron > 0) ) {
		foreach my $a ( @{$entrysRef} ) {
			my ( $alphafolddbtype, $tresolution, $chaind, ) = $a =~ /^\w+;\s+(.*?);\s+(.*?);\s+(.*?)\./;
			next unless ( $tresolution =~ /A$/ );
			my ( $resolution ) = $tresolution =~ /(.*?)\s+A/;
			$tmph{$chaind}{$resolution} = $a;
		}
		foreach my $d ( keys %tmph ) {
			my $mark = 0;
			foreach my $c ( sort {$a <=> $b} keys %{$tmph{$d}} ) {
				my ( $alphafolddb ) = $tmph{$d}{$c} =~ /^(\w+);\s+/;
				# load filtered alphafolddbs
				if ( $mark == 0 ) {
					push( @entrysafterFiltered, $tmph{$d}{$c} );
					$mark++;
				}elsif ( defined $ulinks->{$alphafolddb} ) { push(@entrysafterFiltered, $tmph{$d}{$c}); }
			}
		}
	} else { map{ push(@entrysafterFiltered, $_); } @{$entrysRef}; }

	return \@entrysafterFiltered;
}

sub help_text {
    my $this = shift;
        return <<HELP
Usage: hotspot3d calpro [options]
                             		REQUIRED
--output-dir                 		Output directory of proximity files
--alphafolddb-file-dir              AlphaFoldDB file directory 
--uniprot-id                		Uniprot ID
                             		OPTIONAL
--3d-distance-cutoff         		Maximum 3D distance (<= Angstroms), default: 100
--linear-cutoff              		Minimum linear distance (> peptides), default: 0
--measure                    		Distance measurement between residues (shortest or average), default: average
--help                       		this message
HELP

}

1;
