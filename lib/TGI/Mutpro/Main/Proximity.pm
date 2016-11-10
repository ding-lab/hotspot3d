package TGI::Mutpro::Main::Proximity;
#
#----------------------------------
# $Authors: Beifang Niu and Adam D Scott
# $Date: 2014-01-14 14:34:50 -0500 (Tue Jan 14 14:34:50 CST 2014) $
# $Revision: 0.2 $
# $URL: $
# $Doc: $ proximity pairs searching (main function)
#----------------------------------
#
use strict;
use warnings;

use Carp;
use Getopt::Long;
use IO::File;
use FileHandle;

use TGI::Mutpro::Preprocess::AminoAcid;

my $PVALUEDEFAULT = 0.05;
my $DISTANCEDEFAULT = 10;
my $MAXDISTANCE = 100;

sub new {
    my $class = shift;
    my $this = {};
    $this->{'maf'} = undef;
    $this->{'skip_silent'} = undef;
    $this->{'missense_only'} = undef;
    $this->{'data_dir'} = undef;
    $this->{'drugport_file'} = undef;
    $this->{'output_prefix'} = '3D_Proximity';
    $this->{'pvalue_cutoff'} = undef;
    $this->{'3d_distance_cutoff'} = undef;
    $this->{'linear_cutoff'} = 0;
    $this->{'stat'} = undef;
    $this->{'acceptable_types'} = undef;
    $this->{'amino_acid_header'} = "amino_acid_change";
    $this->{'transcript_id_header'} = "transcript_name";
    bless $this, $class;
    $this->process();
    return $this;
}

sub process {
    my $this = shift;
    my ( $help, $options );
    unless( @ARGV ) { die $this->help_text(); }
    $options = GetOptions (
        'maf-file=s' => \$this->{'maf'},
        'prep-dir=s'    => \$this->{'data_dir'},
        'output-prefix=s' => \$this->{'output_prefix'},
        'drugport-file=s' => \$this->{'drugport_file'},
        'skip-silent' => \$this->{'skip_silent'},
        'missense-only' => \$this->{'missense_only'},
        'p-value-cutoff=f' => \$this->{'pvalue_cutoff'},
        '3d-distance-cutoff=i' => \$this->{'3d_distance_cutoff'},
        'linear-cutoff=i' => \$this->{'linear_cutoff'},
        'amino-acid-header=s' => \$this->{'amino_acid_header'},
        'transcript-id-header=s' => \$this->{'transcript_id_header'},
        'help' => \$help,
    );
    if ( $help ) { print STDERR help_text(); exit 0; }
    unless( $options ) { die $this->help_text(); }
	if ( not defined $this->{'pvalue_cutoff'} ) {
        if ( not defined $this->{'3d_distance_cutoff'} and not defined $this->{'pvalue_cutoff'} ) {
            warn "HotSpot3D Search Warning: no pair distance limit given, setting to default p-value cutoff = 0.05\n";
            $this->{'pvalue_cutoff'} = $PVALUEDEFAULT;
            $this->{'3d_distance_cutoff'} = $MAXDISTANCE;
        } else {
            $this->{'pvalue_cutoff'} = 1;
        }
    } else {
        if ( not defined $this->{'3d_distance_cutoff'} ) {
            $this->{'3d_distance_cutoff'} = $MAXDISTANCE;
        }
    }
    unless( $this->{'data_dir'} ) { warn 'You must provide a output directory ! ', "\n"; die help_text(); }
    unless( -d $this->{'data_dir'} ) { warn 'You must provide a valid data directory ! ', "\n"; die help_text(); }
    unless( $this->{'maf'} and (-e $this->{'maf'}) ) { warn 'You must provide a MAF format file ! ', "\n"; die $this->help_text(); }
    if ( $this->{'drugport_file'} ) { unless( -e $this->{'maf'} ) { warn 'Drugport parsing results file does not exist ! ', "\n"; die $this->help_text(); } }
    my $uniprot_file = "$this->{'data_dir'}\/hugo.uniprot.pdb.transcript.csv";
    unless( -e $uniprot_file ) { warn 'Uniprot parsing file does not exist ! ', "\n"; die $this->help_text(); }
    my $prior_dir = "$this->{'data_dir'}\/prioritization";
    unless( -d $prior_dir ) { die "the Prioritization directory $prior_dir does not exist ! \n"; };
    my @t = qw( num_muts
                num_missense
                num_silent
                num_with_uniprot
                num_unexpect_format
                num_expect_format
                num_trans 
                num_uniprot_involved
                num_trans_with_uniprot
                num_uniprot_with_trans
                num_aa_posmatch
                num_nearmatch
                num_aa_nearmatch
                num_novel
                num_nt_novel
                proximity_close_eachother
                );
    map{ $this->{'stat'}{$_} = 0; } @t;
	$this->setAcceptableMutationTypes();
    # parse uniprot file 
    my $trans_to_uniprot = $this->getTransMaptoUniprot( $uniprot_file );
    $this->{'stat'}{'num_trans_with_uniprot'} = keys %$trans_to_uniprot;
    # parse drugport resuls file
    my $drugport_hash_ref = $this->getDrugportInfo( $this->{'drugport_file'} );
    my $maf_hash_ref = $this->parseMaf( $this->{'maf'}, $trans_to_uniprot );
    $this->{'stat'}{'num_uniprot_involved'} = keys %$maf_hash_ref;
    my ( $pairoutref, $cosmicref, $roiref, $drugport_results_ref, $drugport_nonresults_ref ) = $this->proximitySearching( $maf_hash_ref, $prior_dir, $drugport_hash_ref );
    print STDOUT "searching done...\n";
    my %filterHash; map {
        chomp; @_ = split /\t/;
        my $geneOne = join("\t", @_[0..8]);
        my $geneTwo = join("\t", @_[9..17]);
        print STDOUT $geneOne."\t".$geneTwo."\n";
        if ( defined $filterHash{$geneOne}{$geneTwo} ) { $filterHash{$geneOne}{$geneTwo} .= $_[19]; 
        } elsif ( defined $filterHash{$geneTwo}{$geneOne} ) { $filterHash{$geneTwo}{$geneOne} .= $_[19];
        } else { $filterHash{$geneOne}{$geneTwo} .= $_[18] . "\t" . $_[19]; }
    } @$pairoutref;
    # pour out proximity pairs
    my %sortedHash;
    map { my $e = $_;  map {
            my $f = $_;
            my @t = split /\t/, $filterHash{$e}{$f};
            my $lDistance = $t[0];
            my %ss = map{ ($_, 1) } split /\|/, $t[1];
            my ( %dd, $miniP, $pvaluePart ); $pvaluePart = "";
            foreach my $d (keys %ss) { my @t0 = split / /, $d; $dd{$t0[2]}{$d} = 1; }
            my @t1 = sort {$a<=>$b} keys %dd;
            $miniP = $t1[0];
            foreach my $c ( @t1 ) { foreach my $g ( keys %{$dd{$c}} ) { $pvaluePart .= $g . "|"; } }
            $sortedHash{$miniP}{"$e\t$f\t$lDistance\t$pvaluePart"} = 1;
        } keys %{$filterHash{$e}};
    } keys %filterHash;
    # output pour out
    my $fh   = new FileHandle;
    die "Could not create pairwise close output file\n" unless( $fh->open(">$this->{'output_prefix'}.pairwise") );
	print STDOUT "Creating ".$this->{'output_prefix'}.".pairwise\n";
    map { map { $fh->print( $_."\n" ); $this->{'stat'}{'proximity_close_eachother'}++; } keys %{$sortedHash{$_}} } sort {$a<=>$b} keys %sortedHash;
    $fh->close();
    # pour out mutations close to cosmic
    die "Could not create cosmic close output file\n" unless( $fh->open( ">$this->{'output_prefix'}.cosmic" ) );
	print STDOUT "Creating ".$this->{'output_prefix'}.".cosmic\n";
    map { $fh->print( $_."\n" )  } @$cosmicref; $fh->close();
    # pour out mutations close to ROI
    die "Could not create Region of Interest(ROI) close output file\n" unless( $fh->open( ">$this->{'output_prefix'}.roi" ) );
	print STDOUT "Creating ".$this->{'output_prefix'}.".roi\n";
    map { $fh->print( $_."\n" )  } @$roiref; $fh->close();
    # pour out mutations close to drugs from drugport 
    die "Could not create drugport compound close output file\n" unless( $fh->open( ">$this->{'output_prefix'}.drugs.target" ) );
	print STDOUT "Creating ".$this->{'output_prefix'}.".drugs.target\n";
    map { $fh->print( $_."\n" )  } @$drugport_results_ref; $fh->close();
    # pour out mutations close to drugs from drugport (nontarget) 
    die "Could not create drugport compound close output file (nontarget)\n" unless( $fh->open( ">$this->{'output_prefix'}.drugs.nontarget" ) );
	print STDOUT "Creating ".$this->{'output_prefix'}.".drugs.nontarget\n";
    map { $fh->print( $_."\n" )  } @$drugport_nonresults_ref; $fh->close();
    # post processing like collapsed clean results
    $this->drug_proximity_postprocessing( $this->{'output_prefix'}, $this->{'drugport_file'} );
    print STDOUT "\n\n##################################################\n";
    print STDOUT "total mutations: " . $this->{'stat'}{'num_muts'} . "\n";
    print STDOUT "expected mutations: " . $this->{'stat'}{'num_expect_format'} . "\n";
    print STDOUT "unexpected format mutations: " . $this->{'stat'}{'num_unexpect_format'} . "\n";
    print STDOUT "mutations with matched uniprot: ". $this->{'stat'}{'num_with_uniprot'} . "\n";
    print STDOUT "total transcripts with valid uniprot sequences : " . $this->{'stat'}{'num_trans_with_uniprot'} . "\n";
    print STDOUT "total transcripts in maf : " . $this->{'stat'}{'num_trans'} . "\n";
    print STDOUT "total mutations to be analyzed:  " . $this->{'stat'}{'num_trans_with_uniprot'} . "\n";
    print STDOUT "total uniprots involved: " . $this->{'stat'}{'num_uniprot_involved'} . "\n";
    print STDOUT "\n\n##################################################\n";

    return 1;
}
# parse maf file 
sub parseMaf {
    my ( $this, $maf, $tuHashref )  = @_;
    my $fh = new FileHandle;
    unless( $fh->open($maf) ) { die "Could not open MAF format mutation file\n" };
    my $i = 0; my %fth;
    while ( my $ll = $fh->getline ) {
       if ( $ll =~ m/^Hugo_Symbol/ ) { chomp( $ll );
            %fth = map {($_, $i++)} split( /\t/, $ll );
            last;
       }
    }
    unless (    defined($fth{"Hugo_Symbol"}) 
            and defined($fth{"Chromosome"}) 
            and defined($fth{"Start_Position"})                                
            and defined($fth{"End_Position"}) 
            and defined($fth{"Reference_Allele"})                                
            and defined($fth{"Tumor_Seq_Allele1"}) 
            and defined($fth{"Tumor_Seq_Allele2"}) 
            and defined($fth{"Variant_Classification"}) 
            and defined($fth{$this->{"transcript_id_header"}})
            and defined($fth{$this->{"amino_acid_header"}}) ) {
        die "not a valid MAF annotation file with transcript and amino acid change !\n";
    }
    my @cols = ( $fth{"Hugo_Symbol"}, 
                 $fth{"Chromosome"}, 
                 $fth{"Start_Position"},                        
                 $fth{"End_Position"}, 
                 $fth{"Reference_Allele"}, 
                 $fth{"Tumor_Seq_Allele1"},                        
                 $fth{"Tumor_Seq_Allele2"}, 
                 $fth{"Variant_Classification"}, 
                 $fth{$this->{"transcript_id_header"}}, 
                 $fth{$this->{"amino_acid_header"}} );
    my ( %mafHash, %transHash );
    # reading file content
    while ( my $ll = $fh->getline ) {
        chomp( $ll );
        my ( $gene, $chr, $start, $end, $ref, $vart1, $vart2, $type, $trans, $aac ) = (split /\t/, $ll)[@cols];
        my $tc = join( "\t", $gene, $chr, $start, $end, $aac );
        $this->{'stat'}{'num_muts'}++;
        $transHash{ $trans } = 1;
		next if ( $this->unacceptable( $type ) );
        unless ( $aac =~ /p\.\w\D*\d+/ or $aac =~ /p\.\D*\d+in_frame_ins/i ) {
            print STDERR "Unexpected format for mutation ".$gene.":g.".$chr.":".$start.$ref.$vart2." of type ".$type.": '$aac'\n";
            $this->{'stat'}{'num_unexpect_format'}++;
            next;
        }
        my ( $residue, $position );
        if ( $aac =~ /p\.(\w)\D*(\d+)/ ) { $residue = $1; $position = $2; 
        } else { $position = $aac =~ /p\.\D*(\d+)in_frame_ins/i };
        next unless( (defined $position) and ($position =~ /^\d+$/) );
        $this->{'stat'}{'num_expect_format'}++;
        next unless( defined $tuHashref->{$trans} );
        my $tmp_uniprot_id = $tuHashref->{$trans}->{'UNIPROT'};
        my $tmp_hit_bool = 0; my $tmp_uniprot_position;
        foreach my $tmp_pos ( keys %{$tuHashref->{$trans}->{'POSITION'}} ){
            if ( ($position >= $tmp_pos) and ($position <= $tuHashref->{$trans}->{'POSITION'}->{$tmp_pos}->{'TEND'}) ) {
                $tmp_uniprot_position = $position - $tmp_pos + $tuHashref->{$trans}->{'POSITION'}->{$tmp_pos}->{'UBEGIN'};
                $tmp_hit_bool = 1; 
                last;
            } 
        }
        next if ( $tmp_hit_bool == 0 );
        $mafHash{ $tmp_uniprot_id }{ $tmp_uniprot_position }{ $tc } = 1;
        $this->{'stat'}{'num_with_uniprot'}++;
    }
    $this->{'stat'}{'num_trans'} = keys %transHash;
    $fh->close();

    return \%mafHash;
}
# get mapping information 
# of transcript id to uniprot id
sub getTransMaptoUniprot {
    my ( $this, $uniprotf ) = @_;
    my $fh = new FileHandle;
    unless( $fh->open($uniprotf) ) { die "Could not open uniprot transcript mapping file\n" };
    my %transHash;
    while ( my $a = $fh->getline ) {
        chomp($a);
        my (undef, $uniprotId, undef, undef, $transcripts) = split /\t/, $a;
        next if $transcripts =~ (/N\/A/);
        map { 
            /(\w+)\[(.*?)]/;
            my $tmp_transcript_id = $1;
            $transHash{$tmp_transcript_id}{'UNIPROT'} = $uniprotId;
            map {  /(\d+)\|(\d+)-(\d+)\|(\d+)/; 
                $transHash{$tmp_transcript_id}{'POSITION'}{$2}{'TEND'} = $4; 
                $transHash{$tmp_transcript_id}{'POSITION'}{$2}{'UBEGIN'} = $1;
                $transHash{$tmp_transcript_id}{'POSITION'}{$2}{'UEND'} = $3;
            } split /\:/, $2;
        } split /,/, $transcripts;
    }
    $fh->close();
    return \%transHash;
}

# get drugport database information
sub getDrugportInfo {
    my ( $this, $drugport_f ) = @_;
    my $fh = new FileHandle;
	my %drugport_hash;
    if ( $fh->open( $drugport_f ) ) {
		while ( my $a = $fh->getline ) {
			chomp($a);
			my ( $het, $target_pdb, $not_target_include_compound ) = (split /\t/, $a)[2,4,8];
			unless ( $target_pdb =~ /NULL/ ) { 
				map{ 
					my ($pdb, $chain, $loc) = $_ =~ /(\w+)\|(\w+)\|(\w+)/; 
					$pdb =~ s/ //g; $chain =~ s/ //g; $loc =~ s/ //g;
					unless ( $pdb and $chain and $loc ) { print STDOUT $a."\n"; }
					$drugport_hash{'TARGET'}{uc($pdb)}{$chain}{$loc} = $het;
				} split /,/,$target_pdb; 
			}
			unless ( $not_target_include_compound =~ /NULL/ ) { 
				map{ 
					my ($pdb, $chain, $loc) = $_ =~ /(\w+)\|(\w)\|(\w+)/; $pdb =~ s/ //g; $chain =~ s/ //g; $loc =~ s/ //g;
					unless ( $pdb and $chain and $loc ) { print STDOUT $a."\n"; }
					$drugport_hash{'NONTARGET'}{uc($pdb)}{$chain}{$loc} = $het;
				} split /,/, $not_target_include_compound; 
			}
		}
		$fh->close();
	} 
	return \%drugport_hash;
}

# proximity searching 
sub proximitySearching {
    my ( $this, $mafHashref, $proximityOutPrefix, $drugportref ) = @_;
    my ( @pairResults, @cosmicclose, @roiclose, @drugport_target_results, @drugport_nontarget_results, );
    my $fh = new FileHandle;
	my $AA = new TGI::Mutpro::Preprocess::AminoAcid;
    foreach my $a ( keys %$mafHashref ) {
        my $uniprotf = "$proximityOutPrefix\/$a.ProximityFile.csv";
        next unless( -e $uniprotf and $fh->open($uniprotf) ); 
        while ( my $b = <$fh> ) {
			next if ( $b =~ /UniProt_ID/g );
            chomp( $b ); my @ta = split /\t/, $b;
            my ( $uid1, $chain1, $pdbcor1, $offset1, $residue1, $domain1, $cosmic1, 
                 $uid2, $chain2, $pdbcor2, $offset2, $residue2, $domain2, $cosmic2, 
                 $proximityinfor ) = @ta;
			my $uniprotcor1 = $pdbcor1 + $offset1;
            my $uniprotcor2 = $pdbcor2 + $offset2;
            my $lineardis = undef;
            if ( $uid1 eq $uid2 ) { $lineardis = abs($uniprotcor1 - $uniprotcor2) } else { $lineardis = "N\/A"; }
            #print $a."\t".$uid2."\t".$uniprotcor1."\t".$uniprotcor2."\t".$lineardis."\n";
            if ( defined $mafHashref->{$a}->{$uniprotcor1} ) {
                if ( defined $mafHashref->{$uid2}->{$uniprotcor2} ) {
					#warn "check AA - ".$residue1." - ".$residue2;
					if ( !( $AA->isAA( $residue1 ) ) || !( $AA->isAA( $residue2 ) ) ) { 
						#warn " - bad AA pair"."\n";
						next;
					}
                    ## close each other
                    foreach my $c ( keys %{$mafHashref->{$a}->{$uniprotcor1}} ) {
                        foreach my $d ( keys %{$mafHashref->{$uid2}->{$uniprotcor2}} ) {
                            if ( $this->cutFiltering( $this->{'pvalue_cutoff'}, $this->{'linear_cutoff'}, $this->{'3d_distance_cutoff'}, $lineardis, $proximityinfor) ) {
                                push( @pairResults, join("\t", $c, @ta[1,2,5,6], $d, @ta[8,9,12,13], $lineardis, $proximityinfor) );
                            }
                        }
                    } 
                } else { 
                    ## close to COSMIC/Domain | to do
                    map { 
                        my $t_item = join("\t", $_, @ta[1,2,5,6,8,9,12,13], $lineardis, $proximityinfor); 
                        if ( $domain2 !~ /N\/A/ ) {
                            if ( $this->cutFiltering( $this->{'pvalue_cutoff'}, $this->{'linear_cutoff'}, $this->{'3d_distance_cutoff'}, $lineardis, $proximityinfor) ) {
                                push(@roiclose, $t_item);
                            }
                        };
                        if ( $cosmic2 !~ /N\/A/ ) {
                            if ( $this->cutFiltering( $this->{'pvalue_cutoff'}, $this->{'linear_cutoff'}, $this->{'3d_distance_cutoff'}, $lineardis, $proximityinfor) ) {
                                push(@cosmicclose, $t_item); 
                            }
                        }; 
                    } keys %{$mafHashref->{$a}->{$uniprotcor1}};
                }
            } else {
                if ( defined $mafHashref->{$uid2}->{$uniprotcor2} ) {
                    map {  
                        my $t_item = join("\t", $_, @ta[8,9,12,13,1,2,5,6], $lineardis, $proximityinfor); 
                        if ( $domain1 !~ /N\/A/ ) { 
                            if ( $this->cutFiltering( $this->{'pvalue_cutoff'}, $this->{'linear_cutoff'}, $this->{'3d_distance_cutoff'}, $lineardis, $proximityinfor) ) {
                                push(@roiclose, $t_item); 
                            }
                        }; 
                        if ( $cosmic1 !~ /N\/A/ ) {
                            if ( $this->cutFiltering( $this->{'pvalue_cutoff'}, $this->{'linear_cutoff'}, $this->{'3d_distance_cutoff'}, $lineardis, $proximityinfor) ) {
                                push(@cosmicclose, $t_item); 
                            }
                        } 
                    } keys %{$mafHashref->{$uid2}->{$uniprotcor2}};
                }
            }
            # drugport searching
			if ( $drugportref ) {
				if ( $AA->isHOH( $residue1 ) || $AA->isHOH( $residue2 ) ) {
					#warn "bad AA pair: ".$residue1." - ".$residue2."\n";
					next;
				}
				my %pdbs_hash = map{ my @t0 = split / /, $_; ($t0[1], 1) } split /\|/, $proximityinfor;
				my ( $real_chain1 ) = $chain1 =~ /\[(\w)\]/; my ( $real_chain2 ) = $chain2 =~ /\[(\w)\]/;
				my ( $e, $iter ); 
				map { 
					$iter = $_;
					if ( defined $drugportref->{'TARGET'}->{$iter}->{$real_chain1}->{$pdbcor1} ) {
						if ( defined $mafHashref->{$uid2}->{$uniprotcor2} ) {
							map { 
								if ( $this->cutFiltering( $this->{'pvalue_cutoff'}, $this->{'linear_cutoff'}, $this->{'3d_distance_cutoff'}, $lineardis, $proximityinfor) ) {
									push( @drugport_target_results, join("\t", $iter, $chain1, $pdbcor1, $drugportref->{'TARGET'}->{$iter}->{$real_chain1}->{$pdbcor1}, $_, @ta[8,9,11,12,13], $lineardis, $proximityinfor) );
								}
							} keys %{$mafHashref->{$uid2}->{$uniprotcor2}};
						}
					}
					if ( defined $drugportref->{'TARGET'}->{$iter}->{$real_chain2}->{$pdbcor2} ) {
						if ( defined $mafHashref->{$uid1}->{$uniprotcor1} ) {
							map { 
								if ( $this->cutFiltering( $this->{'pvalue_cutoff'}, $this->{'linear_cutoff'}, $this->{'3d_distance_cutoff'}, $lineardis, $proximityinfor) ) {
									push( @drugport_target_results, join("\t", $iter, $chain2, $pdbcor2, $drugportref->{'TARGET'}->{$iter}->{$real_chain2}->{$pdbcor2}, $_, @ta[1,2,4,5,6], $lineardis, $proximityinfor) );
								} 
							} keys %{$mafHashref->{$uid1}->{$uniprotcor1}};
						}
					}
					if ( defined $drugportref->{'NONTARGET'}->{$iter}->{$real_chain1}->{$pdbcor1} ) {
						if ( defined $mafHashref->{$uid2}->{$uniprotcor2} ) {
							map {
								if ( $this->cutFiltering( $this->{'pvalue_cutoff'}, $this->{'linear_cutoff'}, $this->{'3d_distance_cutoff'}, $lineardis, $proximityinfor) ) {
									push( @drugport_nontarget_results, join("\t", $iter, $chain1, $pdbcor1, $drugportref->{'NONTARGET'}->{$iter}->{$real_chain1}->{$pdbcor1}, $_, @ta[8,9,11,12,13], $lineardis, $proximityinfor) ); 
								}
							} keys %{$mafHashref->{$uid2}->{$uniprotcor2}};
						}
					}
					if ( defined $drugportref->{'NONTARGET'}->{$iter}->{$real_chain2}->{$pdbcor2} ) {
						if ( defined $mafHashref->{$uid1}->{$uniprotcor1} ) {
							map {
								if ( $this->cutFiltering( $this->{'pvalue_cutoff'}, $this->{'linear_cutoff'}, $this->{'3d_distance_cutoff'}, $lineardis, $proximityinfor) ) {
									push( @drugport_nontarget_results, join("\t", $iter, $chain2, $pdbcor2, $drugportref->{'NONTARGET'}->{$iter}->{$real_chain2}->{$pdbcor2}, $_, @ta[1,2,4,5,6], $lineardis, $proximityinfor) );
								}
							} keys %{$mafHashref->{$uid1}->{$uniprotcor1}};
						}
					}
				} keys %pdbs_hash;
			}
        }
        $fh->close();
    }

    return ( \@pairResults, \@cosmicclose, \@roiclose, \@drugport_target_results, \@drugport_nontarget_results );
}

# post processing of drug results
sub drug_proximity_postprocessing {
    my ( $this, $output_prefix, $drugport_parsing_results ) = @_;
	if ( not defined $drugport_parsing_results ) {
		warn "HotSpot3D Search Warning: Skipping drugport proximity, because no results file given.\n";
		return;
	}
    my $sub_fh_target = new FileHandle;
	my $sub_fh_drugport_parsing = new FileHandle;
	my $sub_fh_nontarget = new FileHandle; 
    my $sub_fh_output = new FileHandle;
	unless( $sub_fh_drugport_parsing->open( "$drugport_parsing_results" ) ) {
		warn "Could not open drugprot parsing output file\n";
		return 0;
	}
    die "Could not open drug target output file\n" unless( $sub_fh_target->open( "$output_prefix.drugs.target" ) );
    die "Could not create clean drug output file\n" unless( $sub_fh_output->open( ">$output_prefix.drugs.target.clean" ) );
	print STDOUT "Creating ".$this->{'output_prefix'}.".drugs.target.clean\n";
    $sub_fh_output->print( join( "\t", "Drug", "Drugport_ID", "PDB_ID", "Drug_Chain", "Compound_Location", "Res_Name", "Gene", "Chromosome", "Start", "Stop", "Amino_Acid_Change", "Res_Chain", "Mutation_Location_In_PDB", "Res_Name", "Domain_Annotation", "Cosmic_Annotation", "Linear_Distance_Between_Drug_and_Mutation", "3D_Distance_Information\n" ) );
    my %ss; map { 
        chomp; my @t = split /\t/; unless( $t[4] =~ /NULL/ ) { 
            map { 
                my ($pdb, $chain, $loc) = $_ =~ /(\w+)\|(\w+)\|(\w+)/; 
                $pdb =~ s/ //g; $chain =~ s/ //g; $loc =~ s/ //g; 
                $ss{uc($pdb)}{$chain}{$loc}{'DRUG'}{$t[0]} = 1; 
                $ss{uc($pdb)}{$chain}{$loc}{'DRUGID'}{$t[1]} = 1; 
            } split /,/,$t[4]; 
        } 
    } <$sub_fh_drugport_parsing>; 
    map { 
        chomp; my @t = split /\t/, $_; my ($chain) = $t[1] =~ /\[(\w+)\]/; my $d3info = "";
        map { 
            my @m = split / /, $_; if ( $t[0] eq $m[1] ) { $d3info = $_; } 
		} split /\|/, $t[15];
		my $drug = join( "\|", keys %{$ss{$t[0]}{$chain}{$t[2]}{'DRUG'}} );
		my $drugid = join( "\|", keys %{$ss{$t[0]}{$chain}{$t[2]}{'DRUGID'}} );
		$sub_fh_output->print( join( "\t", $drug, $drugid, @t[0..14], $d3info )."\n" );
	} <$sub_fh_target>;
	$sub_fh_output->close;
	$sub_fh_target->close;
	$sub_fh_drugport_parsing->close;
	die "Could not open drug nontarget output file\n" unless( $sub_fh_nontarget->open( "$output_prefix.drugs.nontarget" ) );
    die "Could not open drugprot parsing output file\n" unless( $sub_fh_drugport_parsing->open( "$drugport_parsing_results" ) );
    my $sub_fh_nontarget_output = new FileHandle;
    die "Could not create clean nontarget drug output file\n" unless( $sub_fh_nontarget_output->open( ">$output_prefix.drugs.nontarget.clean" ) );
	print STDOUT "Creating ".$this->{'output_prefix'}.".drugs.nontarget.clean\n";
    $sub_fh_nontarget_output->print( join( "\t", "Drug", "Drugport_ID", "PDB_ID", "Chain", "Compound_Location", "Res_Name", "Gene", "Chromosome", "Start", "Stop", "Amino_Acid_Change", "Chain", "Mutation_Location_In_PDB", "Res_Name", "Domain_Annotaiton", "Cosmic_Annotation", "Linear_Distance_Betweeen_Drug_and_Mutation", "3D_Distance_Information\n" ) );
    undef %ss; map { 
        chomp; my @t = split /\t/; unless ( $t[8] =~ /NULL/ ) {
            map { 
                my ($pdb, $chain, $loc) = $_ =~ /(\w+)\|(\w+)\|(\w+)/; $pdb =~ s/ //g; $chain =~ s/ //g; $loc =~ s/ //g; $ss{uc($pdb)}{$chain}{$loc}{'DRUG'}{$t[0]} = 1; $ss{uc($pdb)}{$chain}{$loc}{'DRUGID'}{$t[1]} = 1; 
            } split /,/,$t[8];} 
    } <$sub_fh_drugport_parsing>;
    map { 
        chomp; my @t = split /\t/, $_; my ($chain) = $t[1] =~ /\[(\w+)\]/; my $d3info = "";  
        map {  
            my @m = split / /, $_; if ( $t[0] eq $m[1] ) { $d3info = $_; } 
        } split /\|/, $t[15]; 
        my $drug = join( "\|", keys %{$ss{$t[0]}{$chain}{$t[2]}{'DRUG'}} ); 
        my $drugid = join( "\|", keys %{$ss{$t[0]}{$chain}{$t[2]}{'DRUGID'}} );
        $sub_fh_nontarget_output->print( join( "\t", $drug, $drugid, @t[0..14], $d3info )."\n" );
    } <$sub_fh_nontarget>;
    $sub_fh_nontarget_output->close;
    $sub_fh_nontarget->close;
    $sub_fh_drugport_parsing->close;

    return 1;
}

# get drugport database information
sub cutFiltering {
    my ( $this, $pvalue_cut, $linear_cut, $threed_cut, $linear_dis, $info_proximity ) = @_;
	my @infos = split( /\|/ , $info_proximity );
    my ( $dis_3d, $pvalue ) = (split / /, $infos[0])[0,2];
    if ( $linear_dis =~ /N\/A/ ) {
        return 1 if ( ($dis_3d <= $threed_cut) and ($pvalue <= $pvalue_cut) );
    } else {
        return 1 if ( ($dis_3d <= $threed_cut) and ($linear_dis >= $linear_cut) and ($pvalue <= $pvalue_cut) );
    }
    return undef;
}

sub setAcceptableMutationTypes {
	my $this = shift;
	@{$this->{'acceptable_types'}} = ( "Missense_Mutation" );
	if ( $this->{'missense_only'} ) {
		return 1;
	}
	push @{$this->{'acceptable_types'}} , "In_Frame_Ins";
	push @{$this->{'acceptable_types'}} , "In_Frame_Del";
	if ( not $this->{'skip_silent'} ) {
		push @{$this->{'accepatble_types'}} , "Silent";
	}
	return 1;
}

sub unacceptable {
	my ( $this , $type ) = @_;
	if ( grep{ $_ eq $type } @{$this->{'acceptable_types'}} ) {
		return 0;
	}
	return 1;
}

sub help_text{
    my $this = shift;
        return <<HELP

Usage: hotspot3d search [options]

                             REQUIRED
--maf-file                   Input MAF file used to search proximity results
--prep-dir                   HotSpot3D preprocessing results directory

                             OPTIONAL
--drugport-file              DrugPort database parsing results file
--output-prefix              Prefix of output files, default: 3D_Proximity 
--skip-silent                skip silent mutations, default: no
--missense-only              missense mutation only, default: no
--p-value-cutoff             p_value cutoff(<=), default: 0.05
--3d-distance-cutoff         3D distance cutoff (<=), default: 10
--linear-cutoff              Linear distance cutoff (>= peptides), default: 0 
--transcript-id-header       MAF file column header for transcript id's, default: transcript_name
--amino-acid-header          MAF file column header for amino acid changes, default: amino_acid_change 

--help                       this message

HELP

}

1;
