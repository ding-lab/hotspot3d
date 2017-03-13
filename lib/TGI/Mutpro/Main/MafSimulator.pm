package MafSimulator;
#
#----------------------------------
# $Author: R Jay Mashl
# $Date: 2016-10-10 10:42:04 -0500 $
# $Version: 0.1 $
# $URL: $
# $Doc: $ simulated MAF generator for HotSpot3D
# $Doc: $ subs getCoverage, determineCoverage, determineRatio based on preliminary code from Adam D Scott
#----------------------------------
#
use strict;
use warnings;
no  warnings 'uninitialized';

use FileHandle;
use TGI::Mutpro::Preprocess::Uniprot;

# Replacement for rand()
use Math::Random::Secure qw(irand);

sub getCoverage {
# Based on code from Adam D Scott
    my $opts         = shift;
    my $transcripts  = $opts->{'transcripts'};
    my $maf          = $opts->{'maf'};
    my $statsFn      = $opts->{'statsFn'};
    my $sitesFn      = $opts->{'sitesFn'};
    my $debug        = $opts->{'debug'};

    my $geneListFile = $opts->{'geneListFile'};
    my $targetGene   = $opts->{'geneList'};
    $targetGene      = join(",", split( /\s+/, `cat $geneListFile` )) if $geneListFile;
    my %Gene_h       = map{ ($_,0) } split( /,/, $targetGene );  # bool found status

    my $variants = {};

    my ( $TX_FH, $MAF_FH, $STATS_FH, $SITES_FH );
    $TX_FH    = FileHandle->new( $transcripts , "r" );
    die "ERROR: Could not open/read transcripts file ".$transcripts."\n" if ( ! defined $TX_FH );
    $MAF_FH   = FileHandle->new( $maf ,         "r" );
    die "ERROR: Could not open/read maf ".$maf."\n"                      if ( ! defined $MAF_FH );
    $STATS_FH = FileHandle->new( $statsFn ,     "w" );
    die "ERROR: Could not open/write stats file ".$statsFn."\n"          if ( ! defined $STATS_FH );
    $SITES_FH = FileHandle->new( $sitesFn ,     "w" );
    die "ERROR: Could not open/write sites file ".$sitesFn."\n"          if ( ! defined $SITES_FH );

    # Process variant MAF
    while ( <$MAF_FH> ) {
        chomp;
        next if( /^Hugo/i );
        my ( $gene , $chr , $start , $stop , $ref , $alt , $transcript , $aachange ) = (split( /\t/ , $_ ))[0,1,2,3,4,7,9,12];
        print "Current gene = $gene\n" if $debug;
        if( $targetGene ) {
            next if(! exists $Gene_h{$gene} );
            $Gene_h{$gene} = 1;
        }
        if( ! exists $variants->{$gene} ) {
        print "   not seen before\n" if $debug;
            $variants->{$gene}{'vars'}     = ();
            $variants->{$gene}{'pdbDepth'} = ();
            $variants->{$gene}{'numvars'}  = 0;
        } else {
            print "   seen before\n" if $debug;
        }
        my $position =  $aachange;
        $position    =~ s/\D+(\d+)\D*/$1/;
        print "    aachange, pos = $aachange , $position\n" if $debug;
        # Rules:
        # 1. Repeated variants (with different annotations, e.g., ACAD8 11:134132450) count as one
        # 2. Different mutations at same position (e.g., ABCB6  p.R375Q, p.R375W) are treated as independent hotspot.
        #     W A R N I N G: The hotspot clusters output for these simulations lists "position-based cluster" instead.
        #     However, this module produces histograms for both variant- and position-based clusters.
        if( ! exists $variants->{$gene}{'vars'}{$aachange} ) {
            $variants->{$gene}{'vars'}{$aachange}      =  { 'pos' => $position, 'isCovered' => 0 };
            $variants->{$gene}{'pdbDepth'}{$position}  =  0;
            $variants->{$gene}{'numvars'}              += 1;
        }
    }
    if( $targetGene ) { foreach( keys %Gene_h ) { print "# Warning: gene $_ was not found in the MAF\n" if( $Gene_h{$_} == 0 ); } }
    $MAF_FH->close();

    # Process transcripts
    my $statsHdr = join( "\t", "#Gene","UniProt","TxLen","TxFractionCovered","NumMafVars","NumMafVarsCovered", "AvgPdbDepthAtMafVarPositions",
                         "ListOfMutsWithCvg", "ResiOfMutsWithCvg", "ListOfMutWithoutCvg", "ResiOfMutsWithoutCvg" ) . "\n";
    my $bStatsHdr = 0;
    my $sitesHdr = join( "\t", "#Gene","NumSitesAvail","Sites" ) . "\n";
    my $bSitesHdr = 0;
    print ( $statsHdr );
    $STATS_FH->print( $statsHdr );

    while ( <$TX_FH> ) {
        chomp;
        my ( $gene, $uniprot ) = (split( /\t/ , $_ ))[0,1];
        if( $targetGene ) { next if( !exists $Gene_h{$gene} || !$Gene_h{$gene} ); }
        my ( $hasPdbCoverage , $Length  ) = determineCoverage( $uniprot , $gene , $variants );
        if( ! $Length ) {
            print join( "\n", "####", "#### Warning: no transcript found for id=$uniprot", "####" ),"\n";
        } else {   # transcript exists
            my ( $txFractionCovered, $numMafVarsCovered, $avgPdbDepth, $siteList, $coverageLists ) = determineRatio( $uniprot , $Length , $hasPdbCoverage , $gene , $variants );
            my $numMafVars     = ( $variants->{$gene}{'numvars'} || 0 );
            $avgPdbDepth       = sprintf( "%.3f", $avgPdbDepth );
            $txFractionCovered = sprintf( "%.6f", $txFractionCovered );
            my $printLine      = join( "\t",  $gene, $uniprot, $Length, $txFractionCovered, $numMafVars, $numMafVarsCovered, $avgPdbDepth,
                                       $coverageLists->[0], $coverageLists->[1], $coverageLists->[2], $coverageLists->[3] ) . "\n";
            print ( $printLine );
            $STATS_FH->print( $printLine );
            if( $numMafVarsCovered > 0 ) {
                $SITES_FH->print( $sitesHdr )  if ( !$bSitesHdr );
                $SITES_FH->print( $siteList );
                $bSitesHdr = 1;
            }
        }
    }
    $TX_FH->close();
    $STATS_FH->close();
    $SITES_FH->close();
    print "Command getCoverage completed.\n\n";
    return "";
}

sub determineCoverage {
	my ( $id , $gene , $variants ) = @_;
	my $protein        = new TGI::Mutpro::Preprocess::Uniprot( $id );
    my $tmp_seq        = $protein->sequence();
	my @sequence       = split( // , $protein->sequence() );
	my $Length         = scalar @sequence;      # max. sequence length from transcript, may be 0
    my @hasPdbCoverage = (0) x ($Length + 1);   # sets maximum size from transcript; residue numbering starts at 1

    sleep(1);
    foreach( @{ $protein->annotations("PDB") } ) {
        my ( $pdb , $type , $res , $this_resi_range ) = split( /; / );
#        next if ( $type eq "Model" );
        my @parts = split( /, / , $this_resi_range );
        foreach my $part ( @parts ) {
            my ( $chain , $positions ) = split( /=/ , $part      );
            my ( $start , $stop      ) = split( /-/ , $positions );
            $stop =~ s/\.//g;
            my $length = $stop - $start;
#            print( join( "\t" , ( $id , $gene, $pdb , $type , $res , $part , $chain , $positions , $start , $stop , $length ) )."\n" );
                for my $pos ( $start..$stop ) {
                    $hasPdbCoverage[$pos] = 1;   # mark as available
                    $variants->{$gene}{'pdbDepth'}{$pos} += 1 if( exists $variants->{$gene}{'pdbDepth'}{$pos} );
                } #foreach position covered
        } #foreach uneven chain
    } #foreach PDB
	return ( \@hasPdbCoverage , $Length );
}

sub determineRatio {
	my ( $id , $Length , $hasPdbCoverage , $gene , $variants ) = @_;
    my %txCovSitesList = ();   # list of sites having PDB coverage (all such considered available for randomization)
    my ( $sumPdbDepth, $numMafVarsCovered, $numMafVarPositions, $avgPdbDepth, $txFractionCovered ) = (0) x 5;
    my $listing = "";

    for( my $pos=1 ; $pos <= $Length ; $pos++)  {
        if( int( $hasPdbCoverage->[$pos] ) == 1 ) {
            $txCovSitesList{$pos} = 1;
            if( exists $variants->{$gene}{'pdbDepth'}{$pos} ) {
                $numMafVarPositions++;
                foreach( keys %{$variants->{$gene}{'vars'}} ) {  # search to mark covered positions
                    # print "Comparing 'pos' to pos:  ",$variants->{$gene}{'vars'}{$_}{'pos'}," :: ",$pos,"\n";
                    if( $variants->{$gene}{'vars'}{$_}{'pos'} eq $pos ) {   # note: 'pos' is not necessarily numeric e.g., p.*277in_frame_del
                        $variants->{$gene}{'vars'}{$_}{'isCovered'} = 1;
                        $numMafVarsCovered++;
                    }
                }
                $sumPdbDepth += $variants->{$gene}{'pdbDepth'}{$pos};
            }
        }
    }
    $txFractionCovered = scalar( keys %txCovSitesList ) / $Length;
    $avgPdbDepth       = ( $sumPdbDepth / $numMafVarPositions ) if( $numMafVarPositions > 0);
    # Print list of sites available for randomization (and sort for sanity)
    $listing           = sprintf( "%s\n", join("\t",  $gene, scalar( keys %txCovSitesList ), join(" ", (sort {$a <=> $b} keys %txCovSitesList) ), ) );
    # Make lists of un/covered variants
    my %coverages = ();
    foreach( keys %{$variants->{$gene}{'vars'}} ) {
        if( $variants->{$gene}{'vars'}{$_}{'isCovered'} == 1) {
            $coverages{'with'}{'vars'}{$_} = 1;
            $coverages{'with'}{'pos'}{ $variants->{$gene}{'vars'}{$_}{'pos'} } = 1;
        } else {
            $coverages{'without'}{'vars'}{$_} = 1;
            $coverages{'without'}{'pos'}{ $variants->{$gene}{'vars'}{$_}{'pos'} } = 1;
        }
    }
    my $coverageLists = [];
    $coverageLists->[0] = join( ",", keys %{$coverages{'with'}{'vars'}} );
    $coverageLists->[1] = join( ",", keys %{$coverages{'with'}{'pos'}} );
    $coverageLists->[2] = join( ",", keys %{$coverages{'without'}{'vars'}} );
    $coverageLists->[3] = join( ",", keys %{$coverages{'without'}{'pos'}} );
    for my $i (0..3) { $coverageLists->[$i] = "NA" if( ! length( $coverageLists->[$i] ) ); }  # negative response if no keys
	return ( $txFractionCovered, $numMafVarsCovered, $avgPdbDepth, $listing, $coverageLists );
}



sub getRand {
    my $range = shift;
#    return int( rand($range) ); # using standard lib
    return irand($range);
}

sub randomize {
    my $opts              = shift;
    my $statsFn           = $opts->{'statsFn'};
    my $sitesFn           = $opts->{'sitesFn'};
    my $distributionsFn   = $opts->{'distributionsFn'};
    my $numShuffles       = $opts->{'numShuffles'};
    my $entryNumberOffset = $opts->{'entryNumberOffset'};
    my $randomUnique      = $opts->{'randomUnique'};
    my $debug             = $opts->{'debug'};

    my $geneListFile      = $opts->{'geneListFile'};
    my $targetGene        = $opts->{'geneList'};
    $targetGene           = join(",", split( /\s+/, `cat $geneListFile` )) if $geneListFile;
    my %Gene_h            = map{ ($_,0) } split( /,/, $targetGene );  # bool found status

    # Load coverage data
    my $STATSINP   = FileHandle->new( $statsFn, "r" );
    die "ERROR: Could not open/read stats file $statsFn\n" if ( ! defined $STATSINP );
    my $validGenes = {};
    while( <$STATSINP> ) {
        chomp;
        next if( /^#/ );
        my ( $gene, $id, $txLength, $txFraction, $numMafVars, $numMafVarsCovered, $avgPdbDepth ) = split /\t/;
        if ( !$numMafVarsCovered ) {
            print "Notice: No PDB coverage available for MAF variants of gene $gene\n";
            next;
        }
        $validGenes->{$gene} = $numMafVarsCovered;
    }
    $STATSINP->close();

    # Load site availability
    my $RANGEINP    = FileHandle->new( $sitesFn, "r" );
    die "ERROR: Could not open/read sites file $sitesFn\n" if ( ! defined $RANGEINP );
    my $validRanges = {};
    while( <$RANGEINP> ) {
        chomp;
        next if( /^#/ );
        my ( $gene, $numPos, $sitesString ) = split /\t/;
        my $sites = [];
        @{$sites} = split /\s+/, $sitesString;
        die "Error: unexpected number of sites in range for gene $gene\n\n" if ( $numPos != scalar @{$sites} );
        $validRanges->{$gene} = $sites;
    }
    $RANGEINP->close();

    # Generate lists of random residues
    my $DISTOUT  = FileHandle->new( $distributionsFn, "a");
    die "ERROR: Could not open/write/append distributions file $distributionsFn\n\n" if ( ! defined $DISTOUT );
    my $header   = join("\t", "#Gene","EntryNumber","ResidueList" )."\n";
    my $noneMsg  = "# No matching genes found.\n";
    my ($bHeader, $bNoneMsg, $nGenes, $nLines) = (0) x 4;
    foreach my $gene (sort {$a cmp $b} keys %{$validRanges} ) {
        my $numSites = $validGenes->{$gene};
        my $idxMax   = scalar @{$validRanges->{$gene}};
        if( $targetGene ) {
            next if(! exists $Gene_h{$gene} );
            $Gene_h{$gene} = 1;
        }
        if( !$bHeader ) { $DISTOUT->print( $header ); $bHeader = 1; }
        $DISTOUT->print( "#gene = $gene, idxMax = $idxMax, numSites = $numSites\n" ) if $debug;
        print( "gene = $gene, idxMax = $idxMax, numSites = $numSites\n" ) if $debug;
        for( my $i = 0 ; $i < $numShuffles ; $i++ ) {
            if( $randomUnique ) {
                my %selected = ();
                while( scalar( keys %selected ) < $numSites ) { $selected{ $validRanges->{$gene}[ getRand($idxMax) ] } = 1; $bNoneMsg = 1; }
                $DISTOUT->print( join( "\t", $gene, ($i + $entryNumberOffset), join( ",", sort {$a <=> $b} keys %selected ), ), "\n" );
            } else {
                my @selected = ();
                foreach ( 1..$numSites ) { push @selected, $validRanges->{$gene}[ getRand($idxMax) ]; $bNoneMsg = 1; }
                $DISTOUT->print( join( "\t", $gene, ($i + $entryNumberOffset), join( ",", sort {$a <=> $b} @selected ), ), "\n" );
            }
            $nLines++;
        }
        $nGenes++;
    }
    if( $targetGene ) { foreach( keys %Gene_h ) { print "Warning: gene $_ was not found\n" if( $Gene_h{$_} == 0 ); } }
    print ($noneMsg) if( !$bNoneMsg );
    $DISTOUT->close();
    print "Number of new distributions generated: $nLines\n";
    print "Number of genes......................: $nGenes\n";
    print "Command randomize completed.\n\n";
    return "";
}


sub generateMafs {
    my $opts              = shift;
    my $distributionsFn   = $opts->{'distributionsFn'};
    my $srcMaf            = $opts->{'srcMaf'};
    my $mafStemOutFn      = $opts->{'mafStemOutFn'};
    my $entryNumberOffset = $opts->{'entryNumberOffset'};

    my $geneListFile      = $opts->{'geneListFile'};
    my $targetGene        = $opts->{'geneList'};
    $targetGene           = join(",", split( /\s+/, `cat $geneListFile` )) if $geneListFile;
    my %targets           = map{ ($_,0) } split( /,/, $targetGene );

    my %Gene_h  = ();

    # Read distributions
    my $DISTRIB = FileHandle->new( $distributionsFn, "r" );
    die "ERROR: Could not open/read distributions file ".$distributionsFn."\n" if ( ! defined $DISTRIB );
    my $distrib = {};
    while( <$DISTRIB> ) {
        chomp;
        next if( /^#/ );
        my ( $gene, $idx, $sitesString ) = split /\t/;
        if( $idx >= $entryNumberOffset ) {
            $distrib->{$gene}{$idx} = $sitesString;
            if( !$targetGene ) { $Gene_h{$gene} = 0; }
        }
    }
    $DISTRIB->close();

    # Apply/reset genes of interest
    if( $targetGene ) {
        foreach( keys %targets ) { if( ! exists $distrib->{$_} ) { print "Warning: no mafs will ge generated for gene $_ (no distributions available).\n"; } else { $Gene_h{$_} = 0; } }
    }

    # Read source maf
    my $SRCMAF = FileHandle->new( $srcMaf, "r" );
    die "ERROR: Could not open/read source maf ".$srcMaf."\n" if ( ! defined $SRCMAF );
    my $header;
    my %scrmafSeen = ();
    while (<$SRCMAF>) {
        chomp;
        next if( /^#/ );
        if( /^Hugo/i ) { $header = $_; next; }
        my @F = split /\t/;
        my ( $gene , $chr , $start , $stop , $ref , $alt , $transcript , $aachange ) = @F[0,1,2,3,4,7,9,12];
        next if( exists $scrmafSeen{$gene} );
        $scrmafSeen{$gene} = 1;
        next if( ! exists $distrib->{$gene} );
        if( $targetGene ) { next if( ! exists $targets{$gene} ); }

        # Use the first entry for a given gene as a template; the working assumption here is that
        # only the gene and amino acid position should matter. The Variant_Classification field also at least
        # needs to be a valid phrase.
        my ($aaFrom, $resi, $aaTo) = $aachange =~ /(\D+)(\d+)(\D*)/;
        my $nMaf = 0;
        foreach my $idx ( keys %{$distrib->{$gene}} ) {
            my $OUTMAF = FileHandle->new( "$mafStemOutFn.$gene.$idx", "w" );
            die "ERROR: Could not open/write output maf $mafStemOutFn.$gene.$idx\n" if ( ! defined $OUTMAF );
            $OUTMAF->print( "$header\n" );
            foreach( split( /,/, $distrib->{$gene}{$idx}) ) {
#               $F[4]  = "MutationSite";  # do not redefine
                $F[12] = sprintf( "%s%u%s",  $aaFrom, $_, $aaTo );
                $F[13] = "NULL";
                $OUTMAF->print( join("\t", @F),"\n" );
            }
            $OUTMAF->close();
            $nMaf++;
        }
        print "Created $nMaf new simulated maf(s) for gene $gene\n";
        $Gene_h{$gene} = 1;
    }
    $SRCMAF->close();

    # Summary report
    foreach( keys %Gene_h ) { print "Warning: possible maf/distribution mismatch: gene $_ not present in source maf.\n" if( $Gene_h{$_} == 0 ); }
    print "Command generateMafs completed.\n\n";
    return "";
}


sub getSizeHisto {
    my $opts            = shift;
    my $listFile        = $opts->{'clustersListFile'};
    my $sizeHistoFn     = $opts->{'sizeHistoFn'};
    my $distributionsFn = $opts->{'distributionsFn'};
    my $clustVarCntFn   = $opts->{'clusterVarCounts'};
    my $debug           = $opts->{'debug'};
    my $numReported     = { numFiles => 0, numFilesEmpty => 0, numPos => 0, numVar => 0 }; # global level
    my $numAvailable    = { numPos => 0, numVar => 0 };
    my $posSizeHisto    = {};
    my $varSizeHisto    = {};
    my $dist            = {};

    my $DISTIN = FileHandle->new( $distributionsFn, "r" );
    die "ERROR: Could not open/read distributions file $distributionsFn\n" if ( ! defined $DISTIN );
    while( my $inFile = <$DISTIN> ) {
        chomp $inFile;
        next if( $inFile =~ /^#/ );
        my ( $gene, $simID, $residueSimList ) = split /\t/, $inFile;  # ACADS   0  1,2,3,4,5,5,6,7
        my @tmp = split /,/, $residueSimList;
        foreach( @tmp ) {
            $dist->{$simID}{pos}{$_}  = 1;
            $dist->{$simID}{var}{$_} += 1;
        }
        $dist->{$simID}{numPos}  = scalar( keys %{$dist->{$simID}{pos}} ); # total number of positions in trial
        $dist->{$simID}{numVar}  = scalar @tmp;                            # total number of variants in trial
    }
    $DISTIN->close();

    # Count (sub)clusters based on label
    my $CFL = FileHandle->new( $listFile, "r" );
    die "ERROR: Could not open/read list file ".$listFile."\n" if ( ! defined $CFL );

    my $CLUSTVARCNT;
    if( $clustVarCntFn ) {
        $CLUSTVARCNT = FileHandle->new( $clustVarCntFn, "w" );
        die "ERROR: Could not open/write list file ".$clustVarCntFn."\n" if ( ! defined $CLUSTVARCNT );
        $CLUSTVARCNT->print( join("\t", "#SimID", "numVarsInClusters", "numPosInClusters") . "\n" );
    }

    while( my $inFile = <$CFL> ) {
        chomp $inFile;
        my ( $simID, $hsFile )   = split /\t/, $inFile;
        my $posClustInfo         = {}; # position-based clusters
        my $varClustInfo         = {}; # variant-based clusters
        my $bFileHadData         = 0;
        $numAvailable->{numPos} += $dist->{$simID}{numPos};
        $numAvailable->{numVar} += $dist->{$simID}{numVar};
        $numReported->{thisPos}  = 0;  # for current cluster file
        $numReported->{thisVar}  = 0;  # for current cluster file
        # Read clusters file; currently hotspot lists only positions/variants that are in clusters
        print( "Reading $hsFile ...\n" );
        foreach ( `cat $hsFile` ) {
            chomp;
            next if /^Cluster/;
            my ( $clusterLabel, $gene, $aachange, $geodesic ) = (split /\t/, $_)[0,1,2,5];
            $aachange =~ s/\D+(\d+)\D*/$1/;
            $bFileHadData = 1;
            tally( $posClustInfo, $clusterLabel, 1 );
            $numReported->{numPos}  += 1;
            $numReported->{thisPos} += 1;
            tally( $varClustInfo, $clusterLabel, $dist->{$simID}{var}{$aachange} );
            $numReported->{numVar}  += $dist->{$simID}{var}{$aachange};
            $numReported->{thisVar} += $dist->{$simID}{var}{$aachange};
            print "(simID, position, depth): $simID $aachange ".$dist->{$simID}{var}{$aachange}."\n" if $debug;
        }
        $numReported->{numFiles} += 1;
        if( $clustVarCntFn ) {
            $CLUSTVARCNT->print( join("\t", $simID, $numReported->{thisVar}, $numReported->{thisPos}) . "\n" );
        }

        # Update size histo
        if( !$bFileHadData ) {
            $numReported->{numFilesEmpty} += 1;
        } else {
            foreach( keys %{$posClustInfo} ) {
                tally( $posSizeHisto, $posClustInfo->{$_}, 1 );
                tally( $varSizeHisto, $varClustInfo->{$_}, 1 );
            }
        }
    }
    if( $clustVarCntFn ) { $CLUSTVARCNT->close(); }
    $CFL->close();
    doNormalizeAndPrint( $posSizeHisto, $varSizeHisto, $sizeHistoFn, $numReported, $numAvailable );
}

sub mergeHistos {
    my $opts            = shift;
    my $listFile        = $opts->{'histosListFile'};
    my $sizeHistoFn     = $opts->{'sizeHistoFn'};
    my $debug           = $opts->{'debug'};
    my $numFiles        = 0;
    my $posSizeHisto    = {};
    my $varSizeHisto    = {};
    my $numReported     = { numFiles => 0, numFilesEmpty => 0, numPos => 0, numVar => 0 };
    my $numAvailable    = { numPos => 0, numVar => 0 };

    my $CFL = FileHandle->new( $listFile, "r" );
    die "ERROR: Could not open/read list file ".$listFile."\n" if ( ! defined $CFL );
    while( my $inFile = <$CFL> ) {
        chomp $inFile;
        print( "Reading $inFile ...\n" );
        foreach ( `cat $inFile` ) {
            chomp;
            if( /^#/ ) {
                my $tmp = (split /\s+/)[1];
                if( /numFiles:/                   ) { $numReported->{numFiles}      += $tmp; }
                if( /numFilesEmpty:/              ) { $numReported->{numFilesEmpty} += $tmp; }
                if( /numberOfVariantsReported:/   ) { $numReported->{numVar}        += $tmp; }
                if( /numberOfPositionsReported:/  ) { $numReported->{numPos}        += $tmp; }
                if( /numberOfVariantsAvailable:/  ) { $numAvailable->{numVar}       += $tmp; }
                if( /numberOfPositionsAvailable:/ ) { $numAvailable->{numPos}       += $tmp; }
                next;
            }
            my ( $size, $varCount, $varProb, $posCount, $posProb ) = split /\t/;
            tally( $posSizeHisto, $size, $posCount );
            tally( $varSizeHisto, $size, $varCount );
        }
        $numFiles++;
        print( "$numFiles histo files were read. Computing...\n" );
    }
    doNormalizeAndPrint( $posSizeHisto, $varSizeHisto, $sizeHistoFn, $numReported, $numAvailable );
    print( "Done.\n" );
}

sub tally {
    my ( $obj, $bin, $incr ) = @_;
    $incr += 0;
    if( !exists $obj->{$bin} ) { $obj->{$bin} = $incr; } else { $obj->{$bin} += $incr; }
}

sub doNormalizeAndPrint {    # TODO: separate these operations
    my ( $posSizeHisto, $varSizeHisto, $sizeHistoFn, $numReported, $numAvailable ) = @_;

    my $HISTO = FileHandle->new( $sizeHistoFn, "w" );
    die "ERROR: Could not open/write histogram ".$sizeHistoFn."\n" if ( ! defined $HISTO );
    my @lines = ();
    push @lines, "#numFiles: " . $numReported->{numFiles} . "\n";
    push @lines, "#numFilesEmpty: " . $numReported->{numFilesEmpty} . "\n";
    push @lines, "#fractionOfFilesEmpty: " . $numReported->{numFilesEmpty} / $numReported->{numFiles} . "\n";
    push @lines, "#numberOfVariantsAvailable: " . $numAvailable->{numVar} . "\n";
    push @lines, "#numberOfVariantsReported: " . $numReported->{numVar} . "\n";
    push @lines, "#numberOfPositionsAvailable: " . $numAvailable->{numPos} . "\n";
    push @lines, "#numberOfPositionsReported: " . $numReported->{numPos} . "\n";
    push @lines, join( "\t", "#SubClusterSize", "Variant-based(Count,Prob)", "Position-based(Count,Prob)" ) . "\n";

    my %tableRows = ();
    my $posSum = 0;
    foreach( keys %{$posSizeHisto} ) { $posSum += $posSizeHisto->{$_}; }
    foreach( sort{ $a <=> $b } keys %{$posSizeHisto} ) {
        $tableRows{$_}{pos}{count} = $posSizeHisto->{$_};
        $tableRows{$_}{pos}{prob}  = $posSizeHisto->{$_} / $posSum;
    }
    my $varSum = 0;
    foreach( keys %{$varSizeHisto} ) { $varSum += $varSizeHisto->{$_}; }
    foreach( sort{ $a <=> $b } keys %{$varSizeHisto} ) {
        $tableRows{$_}{var}{count} = $varSizeHisto->{$_};
        $tableRows{$_}{var}{prob}  = $varSizeHisto->{$_} / $varSum;
    }
    # Supply missing columns
    foreach( keys %{$posSizeHisto} ) { if( !exists $tableRows{$_}{var}{count}) { $tableRows{$_}{var}{count} = $tableRows{$_}{var}{prob} = 0; } }
    foreach( keys %{$varSizeHisto} ) { if( !exists $tableRows{$_}{pos}{count}) { $tableRows{$_}{pos}{count} = $tableRows{$_}{pos}{prob} = 0; } }

    foreach( sort{ $a <=> $b } keys %tableRows ) {
        push @lines, join( "\t", $_,
                           $tableRows{$_}{var}{count}, $tableRows{$_}{var}{prob},
                           $tableRows{$_}{pos}{count}, $tableRows{$_}{pos}{prob} ) . "\n";
    }
    # Print
    foreach( @lines ) { print $_; $HISTO->print($_); }
    $HISTO->close();
}

1;
