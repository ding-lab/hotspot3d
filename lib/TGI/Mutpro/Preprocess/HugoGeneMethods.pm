package TGI::Mutpro::Preprocess::HugoGeneMethods;
#
#----------------------------------
# $Authors: Beifang Niu & Adam D Scott
# $Date: 2014-01-14 14:34:50 -0500 (Tue Jan 14 14:34:50 CST 2014) $
# $Revision: $
# $URL: $
# $Doc: $ hugo gene processing class 
#----------------------------------
#
use strict;
use warnings;
our $VERSION = '0.3';

use Carp;
use LWP::Simple;
use TGI::Mutpro::Preprocess::HugoGene;
my $HugoUrl = "https://www.genenames.org/cgi-bin/hgnc_downloads?".
    "title=HGNC+output+data&hgnc_dbtag=on&".
    "col=gd_hgnc_id&".
    "col=gd_app_sym&".
    "col=gd_app_name&".
    "col=gd_status&".
    "col=gd_locus_type&".
    "col=gd_prev_sym&".
    "col=gd_prev_name&".
    "col=gd_aliases&".
    "col=gd_name_aliases&".
    "col=gd_pub_chrom_map&".
    "col=gd_date2app_or_res&".
    "col=gd_date_mod&".
    "col=gd_date_sym_change&".
    "col=gd_date_name_change&".
    "col=gd_pub_acc_ids&".
    "col=gd_enz_ids&".
    "col=gd_pub_eg_id&".
    "col=gd_mgd_id&".
    "col=gd_other_ids&".
    "col=gd_other_ids_list&".
    "col=gd_pubmed_ids&".
    "col=gd_pub_refseq_ids&".
    "col=gd_record_type&".
    "col=gd_primary_ids&".
    "col=gd_secondary_ids&".
    "col=gd_vega_ids&".
    "col=gd_lsdb_links&".
    "col=md_gdb_id&".
    "col=md_eg_id&".
    "col=md_mim_id&".
    "col=md_refseq_id&".
    "col=md_prot_id&".
    "col=md_ensembl_id&".
    "col=md_ucsc_id&".
    "status=Approved&".
    "status=Entry+Withdrawn&".
    "status_opt=2&".
    "level=pri_sec&=on&where=&order_by=gd_app_sym_sort&limit=&".
    "format=text&".   
    "submit=submit&.cgifields=&.cgifields=level&.cgifields=chr&.cgifields=status&.cgifields=hgnc_dbtag";
# 0. HGNC ID
# 1. Approved Symbol
# 2. Approved Name
# 3. Status
# 4. Locus Type
# 5. Previous Symbols
# 6. Previous Names
# 7. Aliases
# 8. Name Aliases
# 9. Chromosome
# 10. Date Approved
# 11. Date Modified
# 12. Date Symbol Changed
# 13. Date Name Changed
# 14. Accession Numbers
# 15. Enzyme IDs
# 16. Entrez Gene ID
# 17. MGD ID (mouse genome database)
# 18. Specialist Database Links
# 19. Specialist Database IDs
# 20. Pubmed IDs
# 21. RefSeq IDs
# 22. Record Type
# 23. Primary IDs
# 24. Secondary IDs
# 25. VEGA IDs  
# 26. Locus Specific Databases
# 27. GDB ID (mapped data)
# 28. Entrez Gene ID (mapped data supplied by NCBI)
# 29. OMIM ID (mapped data supplied by NCBI)
# 30. RefSeq (mapped data supplied by NCBI)
# 31. UniProt ID (mapped data supplied by UniProt)
# 32. Ensembl ID (mapped data supplied by Ensembl)
# 33. UCSC ID (mapped data supplied by UCSC)

sub makeHugoGeneObjects {
    # Return: ref to HugoGene objects with key = Hugo name
    #
    # 0. HGNC ID
    # 1. Approved Symbol
    # 2. Approved Name
    # 4. Locus Type   ('pseudogene', 'gene with protein product, function unknown',
    #                  'pseudogene, transcribed', etc.)
    # 5. Previous Symbols
    # 6. Previous Names
    # 7. Aliases
    # 8. Name Aliases
    # 9. Chromosome
    # 16. Entrez Gene ID
    # 28. Entrez Gene ID (mapped data supplied by NCBI)
    # 21. RefSeq IDs
    # 30. RefSeq (mapped data supplied by NCBI)
    # 25. VEGA IDs 
    # 29. OMIM ID (mapped data supplied by NCBI)
    # 31. UniProt ID (mapped data supplied by UniProt)
    # 32. Ensembl ID (mapped data supplied by Ensembl)
    # 33. UCSC ID (mapped data supplied by UCSC) 
	my $list = shift;
    my %hugo_objects;
    my $page = get($HugoUrl);
	if ( not $page ) {
		die "HotSpot3D::HugoGeneMethods::makeHugoGeneObjects ERROR: no data from url request.\n";
	}
    foreach my $line (split /\n/, $page) {
		if ($line =~ /withdrawn/i || $line =~ /HGNC ID\s+Approved\s+Symbol\s+/) { next; }
		my @entries = split /\t/, $line;
		my $hugo = $entries[1];
		$hugo =~ s/\s+//g;
		if ( not exists $list->{$hugo} ) {
			next unless ( scalar keys %{$list} == 0 );
		}
		$hugo_objects{$hugo} = new TGI::Mutpro::Preprocess::HugoGene;
		$hugo_objects{$hugo}->symbol($hugo);
		(defined $entries[0] && $entries[0] ne "" ) || confess "No ID for '$hugo'";
		# Remove 'HGNC:' prefix in ID
		if ($entries[0] =~ /HGNC\:(\d+)/) {
			$hugo_objects{$hugo}->id($1); 
		} else { 
			confess "Unexpected format for Hugo id '$entries[0]'"; 
		}
		if (defined $entries[2]  && $entries[2]  ne "") { $hugo_objects{$hugo}->name($entries[2]); }
		if (defined $entries[4]  && $entries[4]  ne "") { $hugo_objects{$hugo}->type($entries[4]); }
		if (defined $entries[9]  && $entries[9]  ne "") { $hugo_objects{$hugo}->chr($entries[9]); }
		if (defined $entries[29] && $entries[29] ne "") { $hugo_objects{$hugo}->omim($entries[29]); }
		if (defined $entries[31] && $entries[31] ne "") { $hugo_objects{$hugo}->uniprot($entries[31]); }
		if (defined $entries[32] && $entries[32] ne "") { $hugo_objects{$hugo}->ensembl($entries[32]); }	
		if (defined $entries[33] && $entries[33] ne "") { $hugo_objects{$hugo}->ucsc($entries[33]); }
		if (defined $entries[16] && $entries[16] ne "") { $hugo_objects{$hugo}->addEntrezName($entries[16]); }
		if (defined $entries[28] && $entries[28] ne "") { $hugo_objects{$hugo}->addEntrezName($entries[28]); }
		if (defined $entries[21] && $entries[21] ne "") { $hugo_objects{$hugo}->addRefSeqName($entries[21]); }
		if (defined $entries[30] && $entries[30] ne "") { $hugo_objects{$hugo}->addRefSeqName($entries[30]); }
		my ($id);
		if (defined $entries[5] && $entries[5] ne "") {
			foreach $id (split /\,/, $entries[5]) { 
				$id =~ s/\s+//g;
				$hugo_objects{$hugo}->addPreviousSymbol($id); 
			} 
		}
		if (defined $entries[7] && $entries[7] ne "") {
			foreach $id (split /\,/, $entries[7]) { 
				$id =~ s/\s+//g;
				$hugo_objects{$hugo}->addAlias($id); 
			} 
		}
		if (defined $entries[6] && $entries[6] ne "") {
			foreach $id (split /\,/, $entries[6]) { 
				$hugo_objects{$hugo}->addOldDescription($id); 
			} 
		}
		if (defined $entries[8] && $entries[8] ne "") {
			foreach $id (split /\,/, $entries[6]) { 
				$hugo_objects{$hugo}->addOldDescription($id);
			} 
		}
		if (defined $entries[25] && $entries[25] ne "") {
			foreach $id (split /\,/, $entries[25]) { 
			$id =~ s/\s+//g;
			$hugo_objects{$hugo}->addVegaName($id); 
			} 
		}
    }
    return \%hugo_objects;
}

sub downloadPreviousHugoSymbols {
    # Return ref to hash with key = previous Hugo symbol; value = array of current Hugo name(s)
    #        ref to hash with key = current Hugo name; value = array of previous name(s)
    #  return (\%previousToHugo, \%hugoToPrevious);
    # The previous symbol(s) are in $entries[6]
    my (%previousToHugo, %hugoToPrevious, $previous,);
    my $page = get($HugoUrl);
    foreach my $line (split /\n/, $page) {
	if ( $line =~ /withdrawn/i || $line =~ /HGNC ID\s+Approved\s+Symbol\s+/ ) { next; }
	my @entries = split /\t/, $line;
	if ( !defined $entries[5] || $entries[5] eq "" ) { next; }
	foreach $previous ( split /\s+/, $entries[5] ) {
	    $previous =~ s/\,//g;
	    push @{$previousToHugo{$previous}}, $entries[1];
	    push @{$hugoToPrevious{$entries[1]}}, $previous;
	}
    }
    return (\%previousToHugo, \%hugoToPrevious);
}

sub previousHugoSymbols {
    # Input: optional ref to hash of HugoGene objects
    # Return ref to hash with key = previous Hugo symbol; value = array of current Hugo name(s)
    #        ref to hash with key = current Hugo name; value = array of previous name(s)
    #  return (\%previousToHugo, \%hugoToPrevious);
    # 
    my $hugoObjRef = $_[0];
    if (!defined $hugoObjRef ) { 
	my ( $previousToHugoRef, $hugoToPreviousRef ) = downloadPreviousHugoSymbols();
	return ($previousToHugoRef, $hugoToPreviousRef);
    }
    my ($hugo, $previous, $previousRef, %previousToHugo, %hugoToPrevious,  );
    foreach $hugo (keys %{$hugoObjRef} ) {
	$previousRef = $$hugoObjRef{$hugo}->getAllPreviousSymbols();
	foreach $previous ( keys %{$previousRef} ) {
	    push @{$previousToHugo{$previous}}, $hugo;
	    push @{$hugoToPrevious{$hugo}}, $previous;
	}
    }
    return (\%previousToHugo, \%hugoToPrevious);
}

sub downloadHugoAlias {
    # Return ref to hash with key = hugo alias, value = array of Hugo name(s)
    # and another ref to hash with key = hugo name, value = array of aliases
    #   return (\%aliasToHugo, \%hugoToAlias);
    # The aliases are in $entries[7], the Hugo name is in $entries[1]
    my ( %aliasToHugo, %hugoToAlias, $alias );
    my $page = get($HugoUrl);
    foreach my $line (split /\n/, $page) {
	if ( $line =~ /withdrawn/i || $line =~ /HGNC ID\s+Approved\s+Symbol\s+/ ) { next; }
	my @entries = split /\t/, $line;
	# The Hugo name is the second array element and the Uniprot name is the 29th array element
	if ( !defined $entries[7] || $entries[7] eq "" ) { next; }
	foreach $alias ( split /\s+/, $entries[7] ) {
	    $alias =~ s/\,//g;
	    push @{$aliasToHugo{$alias}}, $entries[1];
	    push @{$hugoToAlias{$entries[1]}}, $alias;
	}
    }
    return (\%aliasToHugo, \%hugoToAlias);
}

sub hugoAlias {
    # Input: optional ref to hash of HugoGene objects
    # Return ref to hash with key = hugo alias, value = array of Hugo name(s)
    # and another ref to hash with key = hugo name, value = array of aliases
    #   return (\%aliasToHugo, \%hugoToAlias);
    my $hugoObjRef = $_[0];
    if ( !defined $hugoObjRef ) { 
	my ( $aliasToHugoRef, $hugoToAliasRef ) = downloadHugoAlias();
	return ($aliasToHugoRef, $hugoToAliasRef);
    }
    my ( $hugo, $alias, $aliasRef, %aliasToHugo, %hugoToAlias, );
    foreach $hugo ( keys %{$hugoObjRef} ) {
	$aliasRef = $$hugoObjRef{$hugo}->getAllAliases();
	foreach $alias ( keys %{$aliasRef} ) {
	    push @{$aliasToHugo{$alias}}, $hugo;
	    push @{$hugoToAlias{$hugo}}, $alias;
	}
    }
    return (\%aliasToHugo, \%hugoToAlias);
}

sub isUpdatedHugoGeneRecord {
    # Input: ref to HugoGene object with valid symbol
    #        ref to HugoGeneDb object with symbol that is not currently valid
    #        Maximum number of mismatches to external gene names
    #        Minimum number of matches to external gene names (if not an alias or previous symbol)
    # Return: 1 if the first object is the updated record of the second
    #   
    # Criteria are the following:
    # If chromosomes are defined, they must match (although the HugoGeneDb object (based on
    # information in database) could have more than one chromosome assignment.
    # See if they have external IDs (Ensembl, Vega, etc.)
    # Record number of external ID matches and mismatches
    # If the symbol of the second object is a previous symbol or alias of first,
    # then return 1 if the number of ID mismatches is <= $maxMisMatches
    # If the only thing to go on is the external IDs, then allow given number of mismatches
    # and require $minMatches matches 
    # Kept track of each type of ID in case it makes sense to weigh one more than
    # others.
    my ( $currentRef, $invalidRef, $maxMismatches, $minMatches ) = @_;
    my $debug = 1;
    my ( $entrezMismatch, $ensemblMismatch, $omimMismatch, $vegaMismatch,
	 $refSeqMismatch, $uniprotMismatch, $ucscMismatch, $mismatchTotal, 
	 $entrezMatch, $ensemblMatch, $omimMatch, $vegaMatch,
	 $refSeqMatch, $uniprotMatch, $ucscMatch, $matchTotal,
	 $invalidSymbol, 
	 $hugoIdMatch, $hugoIdMismatch, );
    $invalidSymbol = $$invalidRef->symbol();
    $uniprotMismatch = $ensemblMismatch = $omimMismatch = $ucscMismatch = $vegaMismatch = $refSeqMismatch = $entrezMismatch = $entrezMatch = $ensemblMatch = $omimMatch = $vegaMatch = $refSeqMatch = $uniprotMatch = $omimMatch = $ucscMatch = $hugoIdMatch = $hugoIdMismatch = 0;
    # If they are on different chromosomes, they are not the same gene records
    #### This needs to change so it will see if there is a match to one of (maybe) multiple
    # chromosome assignments
    # i.e. 
    if ( $$invalidRef->hasMultipleChr() ) { 
         # then do something different
    } else {
	# There is a single chromosome assignment for $$invalidRef
	if ( defined $$invalidRef->chr() && defined $$currentRef->chr() &&
	     $$invalidRef->chr() ne "" && $$currentRef->chr() ne "" &&
	     !$$currentRef->sameChr($$invalidRef->chr()) )
         { 
	    print "$invalidSymbol and ", $$currentRef->symbol(), " are on differenc chromosomes \n";
	    return 0; 
	}
    }
    # The chromosomes are the same.  See if they are at same position 
    # If the nucleotide positions are defined for each, require the positions overlap
    # within a tolerance of 100000
    my $tolerance = 100000;
    if ( defined $$invalidRef->start() && defined $$invalidRef->stop() &&
	 defined $$currentRef->start() && defined $$currentRef->stop() )
     {
	my $overlap = $$currentRef->overlaps($$invalidRef->start(), $$invalidRef->stop(), $tolerance);
	if ( defined $overlap && !$overlap ) { return 0; }
    }
    # Hugo numeric id
    if ( defined $$invalidRef->id() && $$invalidRef->id() ne "" &&
	 defined $$currentRef->id() && $$currentRef->id() ne "" )
     {
	if ( $$currentRef->id() ne $$invalidRef->id() ) { 
	    $hugoIdMismatch = 1; 
	} else { $hugoIdMatch = 1; }
    }
    # Uniprot
    if ( defined $$invalidRef->uniprot() && $$invalidRef->uniprot() ne "" &&
	 defined $$currentRef->uniprot() && $$currentRef->uniprot() ne "" ) {
	if ( $$currentRef->uniprot() ne $$invalidRef->uniprot() ) { 
	    $uniprotMismatch = 1; 
	} else { $uniprotMatch = 1; }
    }
    # Ensembl
    if ( defined $$invalidRef->ensembl() && $$invalidRef->ensembl() ne "" &&
	 defined $$currentRef->ensembl() && $$currentRef->ensembl() ne "" ) {
	 if ( $$currentRef->ensembl() ne $$invalidRef->ensembl() ) { 
	     $ensemblMismatch = 1; 
	 } else { 
	     $ensemblMatch = 1; 
	 }
     }	     
    # OMIM
    if ( defined $$invalidRef->omim() && $$invalidRef->omim() ne "" &&
	 defined $$currentRef->omim() && $$currentRef->omim() ne "" ) {
	if ( $$currentRef->omim() ne $$invalidRef->omim() ) { 
	    $omimMismatch = 1; 
	} else {
	    $omimMatch = 1;
	}
    }
    # UCSC
    if ( defined $$invalidRef->ucsc() && $$invalidRef->ucsc() ne "" &&
	 defined $$currentRef->ucsc() && $$currentRef->ucsc() ne "" ) {
	if ( $$currentRef->ucsc() ne $$invalidRef->ucsc() ) { 
	    $ucscMismatch = 1; 
	} else {
	    $ucscMatch = 1;
	}
    }
    # Vega (multiple IDs possible)
    if ( scalar(keys %{$$currentRef->getAllVegaNames()}) > 0 &&
	 scalar(keys %{$$invalidRef->getAllVegaNames()}) > 0 ) {
	$vegaMismatch = 1;
	foreach my $entrez ( keys %{$$invalidRef->getAllVegaNames()} ) {
	    if ( $$currentRef->hasVega($entrez) ) { $vegaMismatch = 0; }
	}
	if ( !$vegaMismatch ) { $vegaMatch = 1; }
    }
    # RefSeq (multiple IDs possible)
    if ( scalar(keys %{$$currentRef->getAllRefSeqNames()}) > 0 &&
	 scalar(keys %{$$invalidRef->getAllRefSeqNames()}) > 0 ) {
	$refSeqMismatch = 1;
	foreach my $entrez ( keys %{$$invalidRef->getAllRefSeqNames()} ) {
	    if ($$currentRef->hasRefSeq($entrez)) { $refSeqMismatch = 0; }
	}
	if (!$refSeqMismatch) { $refSeqMatch = 1; }
    }
    # Entrez (multiple IDs possible)
    if ( scalar(keys %{$$currentRef->getAllEntrezNames()}) > 0 &&
	 scalar(keys %{$$invalidRef->getAllEntrezNames()}) > 0 ) {
	$entrezMismatch = 1;
	foreach my $entrez ( keys %{$$invalidRef->getAllEntrezNames()} ) {
	    if ($$currentRef->hasEntrez($entrez)) { $entrezMismatch = 0; }
	}
	if (!$entrezMismatch) { $entrezMatch = 1; }
    }
    # Get total for number of ID mismatches and matches
    $mismatchTotal = $uniprotMismatch + $ensemblMismatch + $omimMismatch + $ucscMismatch + $vegaMismatch + $refSeqMismatch + $entrezMismatch;
    $matchTotal = $entrezMatch + $ensemblMatch + $omimMatch + $vegaMatch + $refSeqMatch + $uniprotMatch + $omimMatch + $ucscMatch;
    ###    PREVIOUS SYMBOL
    # If the symbol of the invalid gene is listed as a previous symbol of $currentRef
    # then return 1 unless there are too many mismatches with other identifiers
    # But if the Hugo IDs match, it is OK even if there are 'too many' mismatches of other IDs
    if ($$currentRef->isPreviousSymbol($invalidSymbol)) {
	if ($mismatchTotal <= $maxMismatches || $hugoIdMatch) { 
	    return 1; 
	} else {
	    if ($debug) {
		print "WARNING: Invalid symbol $invalidSymbol is a previous symbol for ", $$currentRef->symbol(), 
		" but there are $mismatchTotal inconsistent external IDs and $matchTotal consistent IDs \n";
		$$invalidRef->displayObject(); print "\n";
		$$currentRef->displayObject(); print "\n";
	    }
	    return 0;
	}
    }
    ###    ALIAS
    # If the symbol of the invalid gene is listed as an alias of $currentRef
    # then return 1 unless there are too many mismatches with other identifiers
    # But if the Hugo IDs match, it is OK even if there are 'too many' mismatches of other IDs
    if ($$currentRef->isAlias($invalidSymbol)) { 
	if ($mismatchTotal <= $maxMismatches || $hugoIdMatch) { 
	    return 1; 
	} else {
	    if ($debug) {
		print "WARNING: Invalid symbol $invalidSymbol is an alias symbol for ", $$currentRef->symbol(), 
		" but there are $mismatchTotal inconsistent external IDs  and $matchTotal consistent IDs \n";
		$$invalidRef->displayObject(); print "\n";
		$$currentRef->displayObject(); print "\n";
	    }
	    return 0;
	}
    }
    ###   MATCH BY EXTERNAL IDENTIFIERS
    # 
    if ($matchTotal > $minMatches) {
	if ($mismatchTotal == 0) { 
	    return 1; 
	} else {
	    if ($debug) {
		print "WARNING: Invalid symbol $invalidSymbol has $matchTotal external IDs in common with ", 
		$currentRef->symbol(), 
		" but there are $mismatchTotal inconsistent external IDs \n";
		$$invalidRef->displayObject(); print "\n";
		$$currentRef->displayObject(); print "\n";
	    }
	    return 0;
	}
    }
    ###    Hugo ID number
    #
    if ($hugoIdMatch) {
	if ($mismatchTotal <= $maxMismatches) { 
	    return 1; 
	} else {
	    if ($debug) {
		print "WARNING: Invalid symbol $invalidSymbol has same Id as ", $$currentRef->symbol(), 
		" (", $$invalidRef->id(), 
		") but there are $mismatchTotal inconsistent external IDs and $matchTotal consistent IDs\n";
		$$invalidRef->displayObject(); print "\n";
		$$currentRef->displayObject(); print "\n";
	    }
	    return 0;
	}
    }	

}

sub ensemblToHugoId {
    # Return ref to hash of hash with key = Ensembl ID, Hugo gene name, value = 1
    # Get all current Hugo gene records from the Hugo web site.
    # Look in each record to get the unique Ensembl gene ID
    my ($hugoGeneRef, $hugoId, %ensemblToHugo, $ensemblId,);
    $hugoGeneRef = TGI::Mutpro::Preprocess::HugoGeneMethods::makeHugoGeneObjects();
    foreach $hugoId (keys %{$hugoGeneRef}) {
	$ensemblId = $$hugoGeneRef{$hugoId}->ensembl();
	(defined $ensemblId && $ensemblId ne "") || next;
	${$ensemblToHugo{$ensemblId}}{$hugoId} = 1;
    }
    return \%ensemblToHugo;
}

sub hugoToUniprotId {
    # Return hash with key = Hugo ID, value = Uniprot ID
    my ($hugoGeneRef, $hugoId, $uniprotId, %hugoToUniprot );
    $hugoGeneRef = TGI::Mutpro::Preprocess::HugoGeneMethods::makeHugoGeneObjects();
    foreach $hugoId (keys %{$hugoGeneRef}) {
	$uniprotId = $$hugoGeneRef{$hugoId}->uniprot();
	(!defined $hugoToUniprot{$hugoId}) ||
	    carp "Multiple Uniprot IDs for '$hugoId': '$uniprotId' and '$hugoToUniprot{$hugoId}'";
	$hugoToUniprot{$hugoId} = $uniprotId;
    }
    return \%hugoToUniprot;
}
