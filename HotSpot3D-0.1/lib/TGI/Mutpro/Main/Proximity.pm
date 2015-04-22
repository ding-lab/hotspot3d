package TGI::Mutpro::Main::Proximity;
#
#----------------------------------
# $Authors: Beifang Niu 
# $Date: 2014-01-14 14:34:50 -0500 (Tue Jan 14 14:34:50 CST 2014) $
# $Revision:  $
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

sub new {
    my $class = shift;
    my $this = {};
    $this->{'maf'} = undef;
    $this->{'skip_silent'} = undef;
    $this->{'missense_only'} = undef;
    $this->{'data_dir'} = undef;
    $this->{'drugport_file'} = undef;
    $this->{'output_prefix'} = '3D_Proximity';
    $this->{'pvalue_cutoff'} = 0.05;
    $this->{'3d_cufoff'} = 10;
    $this->{'1d_cutoff'} = 20;
    $this->{'stat'} = undef;
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
        'data-dir=s'    => \$this->{'data_dir'},
        'output-prefix=s' => \$this->{'output_prefix'},
        'drugport-file=s' => \$this->{'drugport_file'},
        'skip-silent' => \$this->{'skip_silent'},
        'missense-only' => \$this->{'missense_only'},
        'p-value=f' => \$this->{'pvalue_cutoff'},
        '3d-dis=i' => \$this->{'3d_cufoff'},
        'linear-dis=i' => \$this->{'1d_cutoff'},
        'help' => \$help,
    );
    if ( $help ) { print STDERR help_text(); exit 0; }
    unless( $options ) { die $this->help_text(); }
    unless( $this->{'data_dir'} ) { warn 'You must provide a output directory ! ', "\n"; die help_text(); }
    unless( -d $this->{'data_dir'} ) { warn 'You must provide a valid data directory ! ', "\n"; die help_text(); }
    unless( $this->{'maf'} and (-e $this->{'maf'}) ) { warn 'You must provide a MAF format file ! ', "\n"; die $this->help_text(); }
    if ( $this->{'drugport_file'} ) { unless( -e $this->{'maf'} ) { warn 'Drugport parsing results file is not exist ! ', "\n"; die $this->help_text(); } }
    my $uniprot_file = "$this->{'data_dir'}\/hugo.uniprot.pdb.transcript.csv";
    unless( -e $uniprot_file ) { warn 'Uniprot parsing file is not exist ! ', "\n"; die $this->help_text(); }
    my $prior_dir = "$this->{'data_dir'}\/prioritization";
    unless( -d $prior_dir ) { die "the Prioritization directory $prior_dir is not exist ! \n"; };
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
    # parse uniprot file 
    my $trans_to_uniprot = $this->getTransMaptoUniprot( $uniprot_file );
    $this->{'stat'}{'num_trans_with_uniprot'} = keys %$trans_to_uniprot;
    # parse drugport resuls file
    my $drugport_hash_ref = $this->getDrugportInfo( $this->{'drugport_file'} );
    my $maf_hash_ref = $this->parseMaf( $this->{'maf'}, $trans_to_uniprot );
    $this->{'stat'}{'num_uniprot_involved'} = keys %$maf_hash_ref;
    my ( $pairoutref, $cosmicref, $roiref, $drugport_results_ref, $drugport_nonresults_ref ) = $this->proximitySearching( $maf_hash_ref, $prior_dir, $drugport_hash_ref );
    print STDERR "searching done...\n";
    my %filterHash; map {
        chomp; @_ = split /\t/;
        my $geneOne = join("\t", @_[0..8]);
        my $geneTwo = join("\t", @_[9..17]);
        print STDERR $geneOne."\t".$geneTwo."\n";
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
    map { map { $fh->print( $_."\n" ); $this->{'stat'}{'proximity_close_eachother'}++; } keys %{$sortedHash{$_}} } sort {$a<=>$b} keys %sortedHash;
    $fh->close();
    # pour out mutations close to cosmic
    die "Could not create cosmic close output file\n" unless( $fh->open( ">$this->{'output_prefix'}.cosmic" ) );
    map { $fh->print( $_."\n" )  } @$cosmicref; $fh->close();
    # pour out mutations close to ROI
    die "Could not create Region of Interest(ROI) close output file\n" unless( $fh->open( ">$this->{'output_prefix'}.roi" ) );
    map { $fh->print( $_."\n" )  } @$roiref; $fh->close();
    # pour out mutations close to drugs from drugport 
    die "Could not create drugport compound close output file\n" unless( $fh->open( ">$this->{'output_prefix'}.drugs.target" ) );
    map { $fh->print( $_."\n" )  } @$drugport_results_ref; $fh->close();
    # pour out mutations close to drugs from drugport (nontarget) 
    die "Could not create drugport compound close output file (nontarget)\n" unless( $fh->open( ">$this->{'output_prefix'}.drugs.nontarget" ) );
    map { $fh->print( $_."\n" )  } @$drugport_nonresults_ref; $fh->close();
    print STDERR "\n\n##################################################\n";
    print STDERR "total mutations: " . $this->{'stat'}{'num_muts'} . "\n";
    print STDERR "expected mutations: " . $this->{'stat'}{'num_expect_format'} . "\n";
    print STDERR "unexpected format mutations: " . $this->{'stat'}{'num_unexpect_format'} . "\n";
    print STDERR "mutations with matched uniprot: ". $this->{'stat'}{'num_with_uniprot'} . "\n";
    print STDERR "total transcripts with valid uniprot sequences : " . $this->{'stat'}{'num_trans_with_uniprot'} . "\n";
    print STDERR "total transcripts in maf : " . $this->{'stat'}{'num_trans'} . "\n";
    print STDERR "total mutations to be analyzed:  " . $this->{'stat'}{'num_trans_with_uniprot'} . "\n";
    print STDERR "total uniprots involved: " . $this->{'stat'}{'num_uniprot_involved'} . "\n";
    print STDERR "\n\n##################################################\n";

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
            and defined($fth{"transcript_name"}) 
            and defined($fth{"amino_acid_change"}) ) {
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
                 $fth{"transcript_name"}, 
                 $fth{"amino_acid_change"} );
    my ( %mafHash, %transHash );
    # reading file content
    while ( my $ll = $fh->getline ) {
        chomp( $ll );
        my ( $gene, $chr, $start, $end, $ref, $vart1, $vart2, $type, $trans, $aac ) = (split /\t/, $ll)[@cols];
        my $tc = join( "\t", $gene, $chr, $start, $end, $aac );
        $this->{'stat'}{'num_muts'}++;
        $transHash{ $trans } = 1;
        next if ( ($this->{'skip_silent'}) and ($type eq "Silent") );
        next if ( ($this->{'missense_only'}) and ($type ne "Missense_Mutation") );
        unless ( $aac =~ /p\.\w\D*\d+/ or $aac =~ /p\.\D*\d+in_frame_ins/i ) {
            print STDERR "Unexpected format for mutation: '$aac'\n";
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
    unless( $fh->open( $drugport_f ) ) { die "Could not open drugport parsing results file ! \n" };
    my %drugport_hash;
    while ( my $a = $fh->getline ) {
        chomp($a);
        my ( $het, $target_pdb, $not_target_include_compound ) = (split /\t/, $a)[2,4,8];
        unless ( $target_pdb =~ /NULL/ ) { 
            map{ 
                my ($pdb, $chain, $loc) = $_ =~ /(\w+)\|(\w+)\|(\w+)/; 
                $pdb =~ s/ //g; $chain =~ s/ //g; $loc =~ s/ //g;
                unless ( $pdb and $chain and $loc ) { print $a."\n"; }
                $drugport_hash{'TARGET'}{uc($pdb)}{$chain}{$loc} = $het;
            } split /,/,$target_pdb; 
        }
        unless ( $not_target_include_compound =~ /NULL/ ) { 
            map{ 
                my ($pdb, $chain, $loc) = $_ =~ /(\w+)\|(\w)\|(\w+)/; $pdb =~ s/ //g; $chain =~ s/ //g; $loc =~ s/ //g;
                unless ( $pdb and $chain and $loc ) { print $a."\n"; }
                $drugport_hash{'NONTARGET'}{uc($pdb)}{$chain}{$loc} = $het;
            } split /,/, $not_target_include_compound; 
        }
    }
    $fh->close();
    return \%drugport_hash;
}

# proximity searching 
sub proximitySearching {
    my ( $this, $mafHashref, $proximityOutPrefix, $drugportref ) = @_;
    my ( @pairResults, @cosmicclose, @roiclose, @drugport_target_results, @drugport_nontarget_results, );
    my $fh = new FileHandle;
    foreach my $a ( keys %$mafHashref ) {
        my $uniprotf = "$proximityOutPrefix\/$a.ProximityFile.csv";
        next unless( -e $uniprotf and $fh->open($uniprotf) ); 
        while ( my $b = <$fh> ) {
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
                    ## close each other
                    foreach my $c ( keys %{$mafHashref->{$a}->{$uniprotcor1}} ) {
                        foreach my $d ( keys %{$mafHashref->{$uid2}->{$uniprotcor2}} ) {
                            if ( $this->cutFiltering( $this->{'pvalue_cutoff'}, $this->{'1d_cutoff'}, $this->{'3d_cufoff'}, $lineardis, $proximityinfor) ) {
                                push( @pairResults, join("\t", $c, @ta[1,2,5,6], $d, @ta[8,9,12,13], $lineardis, $proximityinfor) );
                            }
                        }
                    } 
                } else { 
                    ## close to COSMIC/Domain | to do
                    map { 
                        my $t_item = join("\t", $_, @ta[1,2,5,6,8,9,12,13], $lineardis, $proximityinfor); 
                        if ( $domain2 !~ /N\/A/ ) {
                            if ( $this->cutFiltering( $this->{'pvalue_cutoff'}, $this->{'1d_cutoff'}, $this->{'3d_cufoff'}, $lineardis, $proximityinfor) ) {
                                push(@roiclose, $t_item);
                            }
                        };
                        if ( $cosmic2 !~ /N\/A/ ) {
                            if ( $this->cutFiltering( $this->{'pvalue_cutoff'}, $this->{'1d_cutoff'}, $this->{'3d_cufoff'}, $lineardis, $proximityinfor) ) {
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
                            if ( $this->cutFiltering( $this->{'pvalue_cutoff'}, $this->{'1d_cutoff'}, $this->{'3d_cufoff'}, $lineardis, $proximityinfor) ) {
                                push(@roiclose, $t_item); 
                            }
                        }; 
                        if ( $cosmic1 !~ /N\/A/ ) {
                            if ( $this->cutFiltering( $this->{'pvalue_cutoff'}, $this->{'1d_cutoff'}, $this->{'3d_cufoff'}, $lineardis, $proximityinfor) ) {
                                push(@cosmicclose, $t_item); 
                            }
                        } 
                    } keys %{$mafHashref->{$uid2}->{$uniprotcor2}};
                }
            }
            # drugport searching
            my %pdbs_hash = map{ my @t0 = split / /, $_; ($t0[1], 1) } split /\|/, $proximityinfor;
            my ( $real_chain1 ) = $chain1 =~ /\[(\w)\]/; my ( $real_chain2 ) = $chain2 =~ /\[(\w)\]/;
            my ( $e, $iter ); 
            map { 
                $iter = $_;
                if ( defined $drugportref->{'TARGET'}->{$iter}->{$real_chain1}->{$pdbcor1} ) {
                    if ( defined $mafHashref->{$uid2}->{$uniprotcor2} ) {
                        map { 
                            if ( $this->cutFiltering( $this->{'pvalue_cutoff'}, $this->{'1d_cutoff'}, $this->{'3d_cufoff'}, $lineardis, $proximityinfor) ) {
                                push( @drugport_target_results, join("\t", $iter, $chain1, $pdbcor1, $drugportref->{'TARGET'}->{$iter}->{$real_chain1}->{$pdbcor1}, $_, @ta[8,9,11,12,13], $lineardis, $proximityinfor) );
                            }
                        } keys %{$mafHashref->{$uid2}->{$uniprotcor2}};
                    }
                }
                if ( defined $drugportref->{'TARGET'}->{$iter}->{$real_chain2}->{$pdbcor2} ) {
                    if ( defined $mafHashref->{$uid1}->{$uniprotcor1} ) {
                        map { 
                            if ( $this->cutFiltering( $this->{'pvalue_cutoff'}, $this->{'1d_cutoff'}, $this->{'3d_cufoff'}, $lineardis, $proximityinfor) ) {
                                push( @drugport_target_results, join("\t", $iter, $chain2, $pdbcor2, $drugportref->{'TARGET'}->{$iter}->{$real_chain2}->{$pdbcor2}, $_, @ta[1,2,4,5,6], $lineardis, $proximityinfor) );
                            } 
                        } keys %{$mafHashref->{$uid1}->{$uniprotcor1}};
                    }
                }
                if ( defined $drugportref->{'NONTARGET'}->{$iter}->{$real_chain1}->{$pdbcor1} ) {
                    if ( defined $mafHashref->{$uid2}->{$uniprotcor2} ) {
                        map {
                            if ( $this->cutFiltering( $this->{'pvalue_cutoff'}, $this->{'1d_cutoff'}, $this->{'3d_cufoff'}, $lineardis, $proximityinfor) ) {
                                push( @drugport_nontarget_results, join("\t", $iter, $chain1, $pdbcor1, $drugportref->{'NONTARGET'}->{$iter}->{$real_chain1}->{$pdbcor1}, $_, @ta[8,9,11,12,13], $lineardis, $proximityinfor) ); 
                            }
                        } keys %{$mafHashref->{$uid2}->{$uniprotcor2}};
                    }
                }
                if ( defined $drugportref->{'NONTARGET'}->{$iter}->{$real_chain2}->{$pdbcor2} ) {
                    if ( defined $mafHashref->{$uid1}->{$uniprotcor1} ) {
                        map {
                            if ( $this->cutFiltering( $this->{'pvalue_cutoff'}, $this->{'1d_cutoff'}, $this->{'3d_cufoff'}, $lineardis, $proximityinfor) ) {
                                push( @drugport_nontarget_results, join("\t", $iter, $chain2, $pdbcor2, $drugportref->{'NONTARGET'}->{$iter}->{$real_chain2}->{$pdbcor2}, $_, @ta[1,2,4,5,6], $lineardis, $proximityinfor) );
                            }
                        } keys %{$mafHashref->{$uid1}->{$uniprotcor1}};
                    }
                }
            } keys %pdbs_hash;
        }
        $fh->close();
    }

    return ( \@pairResults, \@cosmicclose, \@roiclose, \@drugport_target_results, \@drugport_nontarget_results );
}

# get drugport database information
sub cutFiltering {
    my ( $this, $pvalue_cut, $linear_cut, $threed_cut, $linear_dis, $info_proximity ) = @_;
    my ( $dis_3d, $pvalue ) = (split / /, $info_proximity)[0,2];
    if ( $linear_dis =~ /N\/A/ ) {
        return 1 if ( ($dis_3d <= $threed_cut) and ($pvalue <= $pvalue_cut) );
    } else {
        return 1 if ( ($dis_3d <= $threed_cut) and ($linear_dis >= $linear_cut) and ($pvalue <= $pvalue_cut) );
    }
    return undef;
}

sub help_text{
    my $this = shift;
        return <<HELP

Usage: hotspot3d search [options]

--maf-file              Input MAF file
                        In addition to the standard version 2.3 MAF headers, there needs to be 3 columns appended.
                        These column headers in the MAF must have these names in the header in order for the tool to
                        find them: 
                                transcript_name - the transcript name, such as NM_000028 
                                amino_acid_change - the amino acid change, such as p.R290H 
--data-dir		HotSpot3D preprocessing results directory

--drugport-file         Drugport database parsing results file ( optional )
--output-prefix         Prefix of output files, default: 3D_Proximity 
--skip-silent           skip silent mutations, default: no
--missense-only         missense mutation only, default: no
--p-value               p_value cutoff(<=), default: 0.05
--3d-dis                3D distance cutoff (<=), default: 10
--linear-dis            linear distance cutoff (>=), default: 20 

--help			this message

HELP

}

1;

