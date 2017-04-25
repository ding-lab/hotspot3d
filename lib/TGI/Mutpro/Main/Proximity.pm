package TGI::Mutpro::Main::Proximity;
#
#----------------------------------
# $Authors: Beifang Niu and Adam D Scott
# $Date: 2014-01-14 14:34:50 -0500 (Tue Jan 14 14:34:50 CST 2014) $
# $Revision: 0.3 $
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
	$this->{'maf_file'} = undef;
	$this->{'skip_silent'} = undef;
	$this->{'missense_only'} = undef;
	$this->{'data_dir'} = undef;
	$this->{'uniprot_file'} = undef;
	$this->{'site_file'} = undef;
	$this->{'drugport_file'} = undef;
	$this->{'output_prefix'} = '3D_Proximity';
	$this->{'p_value_cutoff'} = undef;
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
	my ( $prior_dir ) = $this->setOptions();
	$this->initializeStats();
	$this->setAcceptableMutationTypes();
	
	my $trans_to_uniprot = $this->getTranscriptsToUniprot();

	my $drugport_hash_ref = $this->getDrugportInfo();

	my $site_hash_ref = $this->getSites( $trans_to_uniprot );

	my $maf_hash_ref = $this->parseMaf( $this->{'maf_file'}, $trans_to_uniprot );

	print STDOUT "searching...\n";
	my ( $pairoutref , $cosmicref , $roiref , 
		 $drugport_results_ref , $drugport_nonresults_ref ,
		 $siteSiteRef , $mutationSiteRef ) = $this->proximitySearching( $maf_hash_ref , $prior_dir , $drugport_hash_ref , $site_hash_ref );
	print STDOUT "searching done...\n";

	my $sortedHash = $this->gatherOutput( $pairoutref );

	$this->writeMutationMutationPairwise( $sortedHash );
	$this->writeMutationCOSMIC( $cosmicref );
	$this->writeMutationROI( $roiref );
	$this->writeMutationDrugTarget( $drugport_results_ref );
	$this->writeMutationDrugNontarget( $drugport_nonresults_ref );
	$this->writeSiteSite( $siteSiteRef );
	$this->writeMutationSite( $mutationSiteRef );
	$this->drug_proximity_postprocessing( $this->{'output_prefix'}, $this->{'drugport_file'} );

	$this->printSummaryStats();

	return 1;
}

sub printSummaryStats {
	my $this = shift;
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
	return;
}

sub writeMutationSite {
	my ( $this , $mutationSite ) = @_;
	# pour out mutations close to sites
	my $fh   = new FileHandle;
	die "Could not create mutation-site output file\n" unless( $fh->open( ">$this->{'output_prefix'}.musites" ) );
	print STDOUT "Creating ".$this->{'output_prefix'}.".musites\n";
	my $header = $this->mutationSiteHeader( );
	$fh->print( $header."\n" );
	map { $fh->print( $_."\n" )  } @$mutationSite; $fh->close();
	return;
}

sub mutationSiteHeader {
	my $this = @_;
	my @header = ( "Gene1" , "Chromosome1" , "Start1" , "Stop1" , "Mutation1" ,
				   "Chain1" , "Position1" , "Feature1" , "COSMIC1" ,
				   "Gene2" , "Transcript2" , "TranscriptPosition2" , "Site2" , 
				   "Chain2" , "Position2" , "Feature2" , "COSMIC2" ,
				   "LinearDistance" , "DistanceInfo" );
	return join( "\t" , @header );
}

sub writeSiteSite {
	my ( $this , $siteSite ) = @_;
	# pour out sites close to sites
	my $fh   = new FileHandle;
	die "Could not create site-site output file\n" unless( $fh->open( ">$this->{'output_prefix'}.sites" ) );
	print STDOUT "Creating ".$this->{'output_prefix'}.".sites\n";
	my $header = $this->siteSiteHeader( );
	$fh->print( $header."\n" );
	map { $fh->print( $_."\n" )  } @$siteSite; $fh->close();
	return;
}

sub siteSiteHeader {
	my $this = @_;
	my @header = ( "Gene1" , "Transcript1" , "TranscriptPosition1" , "Site1" , 
				   "Chain1" , "Position1" , "Feature1" , "COSMIC1" ,
				   "Gene2" , "Transcript2" , "TranscriptPosition2" , "Site2" , 
				   "Chain2" , "Position2" , "Feature2" , "COSMIC2" ,
				   "LinearDistance" , "DistanceInfo" );
	return join( "\t" , @header );
}

sub writeMutationDrugNontarget {
	my ( $this , $drugport_nonresults_ref ) = @_;
	# pour out mutations close to drugs from drugport (nontarget) 
	my $fh   = new FileHandle;
	die "Could not create drugport compound close output file (nontarget)\n" unless( $fh->open( ">$this->{'output_prefix'}.drugs.nontarget" ) );
	print STDOUT "Creating ".$this->{'output_prefix'}.".drugs.nontarget\n";
	map { $fh->print( $_."\n" )  } @$drugport_nonresults_ref; $fh->close();
	return;
}

sub writeMutationDrugTarget {
	my ( $this , $drugport_results_ref ) = @_;
	# pour out mutations close to drugs from drugport 
	my $fh   = new FileHandle;
	die "Could not create drugport compound close output file\n" unless( $fh->open( ">$this->{'output_prefix'}.drugs.target" ) );
	print STDOUT "Creating ".$this->{'output_prefix'}.".drugs.target\n";
	map { $fh->print( $_."\n" )  } @$drugport_results_ref; $fh->close();
	return;
}

sub writeMutationROI {
	my ( $this , $roiref ) = @_;
	# pour out mutations close to ROI
	my $fh   = new FileHandle;
	die "Could not create Region of Interest(ROI) close output file\n" unless( $fh->open( ">$this->{'output_prefix'}.roi" ) );
	print STDOUT "Creating ".$this->{'output_prefix'}.".roi\n";
	my $header = $this->ROIHeader( );
	$fh->print( $header."\n" );
	map { $fh->print( $_."\n" )  } @$roiref; $fh->close();
	return;
}

sub ROIHeader {
	my $this = @_;
	my @header = ( "Gene1" , "Chromosome1" , "Start1" , "Stop1" , "Mutation1" ,
				   "Chain1" , "Position1" , "Feature1" , "COSMIC1" ,
				   "Chain2" , "Position2" , "Feature2" , "COSMIC2" ,
				   "LinearDistance" , "DistanceInfo" );
	return join( "\t" , @header );
}

sub writeMutationCOSMIC {
	my ( $this , $cosmicref ) = @_;
	# pour out mutations close to cosmic
	my $fh   = new FileHandle;
	die "Could not create cosmic close output file\n" unless( $fh->open( ">$this->{'output_prefix'}.cosmic" ) );
	print STDOUT "Creating ".$this->{'output_prefix'}.".cosmic\n";
	my $header = $this->COSMICHeader( );
	$fh->print( $header."\n" );
	map { $fh->print( $_."\n" )  } @$cosmicref; $fh->close();
	return;
}

sub COSMICHeader {
	my $this = @_;
	my @header = ( "Gene1" , "Chromosome1" , "Start1" , "Stop1" , "Mutation1" ,
				   "Chain1" , "Position1" , "Feature1" , "COSMIC1" ,
				   "Chain2" , "Position2" , "Feature2" , "COSMIC2" ,
				   "LinearDistance" , "DistanceInfo" );
	return join( "\t" , @header );
}

sub writeMutationMutationPairwise {
	my ( $this , $sortedHash ) = @_;
	# output pour out
	my $fh   = new FileHandle;
	die "Could not create pairwise close output file\n" unless( $fh->open(">$this->{'output_prefix'}.pairwise") );
	print STDOUT "Creating ".$this->{'output_prefix'}.".pairwise\n";
	my $header = $this->pairwiseHeader( );
	$fh->print( $header."\n" );
	map {
		map {
			$fh->print( $_."\n" );
			$this->{'stat'}{'proximity_close_eachother'}++;
		} keys %{$sortedHash->{$_}};
	} sort {$a<=>$b} keys %{$sortedHash};
	$fh->close();
	return;
}

sub pairwiseHeader {
	my $this = @_;
	my @header = ( "Gene1" , "Chromosome1" , "Start1" , "Stop1" , "Mutation1" ,
				   "Chain1" , "Position1" , "Feature1" , "COSMIC1" ,
				   "Gene2" , "Chromosome2" , "Start2" , "Stop2" , "Mutation2" ,
				   "Chain2" , "Position2" , "Feature2" , "COSMIC2" ,
				   "LinearDistance" , "DistanceInfo" );
	return join( "\t" , @header );
}

sub gatherOutput {
	my ( $this , $pairoutref ) = @_;
	my $filterHash = {};
	$this->getMutationMutationPairs( $pairoutref , $filterHash );
	# pour out proximity pairs
	my %sortedHash;
	map {
		my $mutation1Info = $_;
		map {
			my $mutation2Info = $_;
			my @t = split /\t/, $filterHash->{$mutation1Info}->{$mutation2Info};
			my $lDistance = $t[0];
			my %ss = map{ ($_, 1) } split /\|/, $t[1];
			my ( %dd, $miniP, $pvaluePart ); $pvaluePart = "";
			foreach my $d (keys %ss) {
				my @t0 = split / /, $d; $dd{$t0[2]}{$d} = 1;
			}
			my @t1 = sort {$a<=>$b} keys %dd;
			$miniP = $t1[0];
			foreach my $c ( @t1 ) {
				foreach my $g ( keys %{$dd{$c}} ) {
					$pvaluePart .= $g . "|";
				}
			}
			$sortedHash{$miniP}{"$mutation1Info\t$mutation2Info\t$lDistance\t$pvaluePart"} = 1;
		} keys %{$filterHash->{$mutation1Info}};
	} keys %{$filterHash};
	return \%sortedHash;
}

sub getMutationMutationPairs {
	my ( $this , $pairoutref , $filterHash ) = @_;
	map {
		chomp; @_ = split /\t/;
		my $mutation1Info = join("\t", @_[0..8]);
		my $mutation2Info = join("\t", @_[9..17]);
		print STDOUT $mutation1Info."\t".$mutation2Info."\n";
		if ( defined $filterHash->{$mutation1Info}->{$mutation2Info} ) {
			$filterHash->{$mutation1Info}->{$mutation2Info} .= $_[19]; 
		} elsif ( defined $filterHash->{$mutation2Info}->{$mutation1Info} ) {
			$filterHash->{$mutation2Info}->{$mutation1Info} .= $_[19];
		} else {
			$filterHash->{$mutation1Info}->{$mutation2Info} .= $_[18] . "\t" . $_[19];
		}
	} @$pairoutref;
	return;
}

sub getTranscriptsToUniprot {
	my $this = shift;
	# parse uniprot file 
	my $trans_to_uniprot = $this->getTransMaptoUniprot( $this->{'uniprot_file'} );
	$this->{'stat'}{'num_trans_with_uniprot'} = keys %$trans_to_uniprot;
	return $trans_to_uniprot;
}

sub initializeStats {
	my $this = shift;
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
	return;
}

sub setOptions {
	my ( $this ) = shift;
	my ( $help, $options );
	unless( @ARGV ) { die $this->help_text(); }
	$options = GetOptions (
		'maf-file=s' => \$this->{'maf_file'},
		'prep-dir=s'	=> \$this->{'data_dir'},
		'output-prefix=s' => \$this->{'output_prefix'},
		'drugport-file=s' => \$this->{'drugport_file'},
		'site-file=s' => \$this->{'site_file'},
		'skip-silent' => \$this->{'skip_silent'},
		'missense-only' => \$this->{'missense_only'},
		'p-value-cutoff=f' => \$this->{'p_value_cutoff'},
		'3d-distance-cutoff=i' => \$this->{'3d_distance_cutoff'},
		'linear-cutoff=i' => \$this->{'linear_cutoff'},
		'amino-acid-header=s' => \$this->{'amino_acid_header'},
		'transcript-id-header=s' => \$this->{'transcript_id_header'},
		'help' => \$help,
	);
	if ( $help ) { print STDERR help_text(); exit 0; }
	unless( $options ) { die $this->help_text(); }
	if ( not defined $this->{'p_value_cutoff'} ) {
		if ( not defined $this->{'3d_distance_cutoff'} and not defined $this->{'p_value_cutoff'} ) {
			warn "HotSpot3D Search Warning: no pair distance limit given, setting to default p-value cutoff = 0.05\n";
			$this->{'p_value_cutoff'} = $PVALUEDEFAULT;
			$this->{'3d_distance_cutoff'} = $MAXDISTANCE;
		} else {
			$this->{'p_value_cutoff'} = 1;
		}
	} else {
		if ( not defined $this->{'3d_distance_cutoff'} ) {
			$this->{'3d_distance_cutoff'} = $MAXDISTANCE;
		}
	}
	unless( $this->{'data_dir'} ) { warn 'You must provide a output directory ! ', "\n"; die help_text(); }
	unless( -d $this->{'data_dir'} ) { warn 'You must provide a valid data directory ! ', "\n"; die help_text(); }
	unless( $this->{'maf_file'} and (-e $this->{'maf_file'}) ) { warn 'You must provide a MAF format file ! ', "\n"; die $this->help_text(); }
	if ( $this->{'drugport_file'} ) { unless( -e $this->{'maf_file'} ) { warn 'Drugport parsing results file does not exist ! ', "\n"; die $this->help_text(); } }
	if ( $this->{'site_file'} ) { unless( -e $this->{'site_file'} ) { warn "site-file provided, but it does not exist!\n"; die $this->help_text(); } }
	$this->{'uniprot_file'} = "$this->{'data_dir'}\/hugo.uniprot.pdb.transcript.csv";
	unless( -e $this->{'uniprot_file'} ) { warn 'Uniprot parsing file does not exist ! ', "\n"; die $this->help_text(); }
	my $prior_dir = "$this->{'data_dir'}\/prioritization";
	unless( -d $prior_dir ) { die "the Prioritization directory $prior_dir does not exist ! \n"; };
	return ( $prior_dir );
}

# parse maf file 
sub parseMaf {
	my ( $this, $maf, $trans_to_uniprot )  = @_;
	my $fh = new FileHandle;
	unless( $fh->open($maf) ) { die "Could not open MAF format mutation file\n" };
	my $i = 0; my %fth;
	while ( my $ll = $fh->getline ) {
	   if ( $ll =~ m/^Hugo_Symbol/ ) { chomp( $ll );
			%fth = map {($_, $i++)} split( /\t/, $ll );
			last;
	   }
	}
	unless (	defined($fth{"Hugo_Symbol"}) 
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
	while ( my $line = $fh->getline ) {
		chomp( $line );
		my ( $gene, $chr, $start, $end, $ref, $vart1, $vart2, $type, $trans, $aac ) = (split /\t/, $line)[@cols];
		my $mutationKey = join( "\t", $gene, $chr, $start, $end, $aac );
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
		next unless( defined $trans_to_uniprot->{$trans} );
		my $tmp_uniprot_id = $trans_to_uniprot->{$trans}->{'UNIPROT'};
		my ( $tmp_hit_bool , $tmp_uniprot_position ) = $this->getPositionMatch( $trans_to_uniprot , $trans , $position );
		next if ( $tmp_hit_bool == 0 );
		$mafHash{ $tmp_uniprot_id }{ $tmp_uniprot_position }{ $mutationKey } = 1;
		$this->{'stat'}{'num_with_uniprot'}++;
	}
	$fh->close();

	$this->{'stat'}{'num_trans'} = keys %transHash;
	$this->{'stat'}{'num_uniprot_involved'} = keys %mafHash;

	return \%mafHash;
}

sub getPositionMatch {
	my ( $this , $trans_to_uniprot , $trans , $position ) = @_;
	my $tmp_uniprot_position = -1;
	my $tmp_hit_bool = 0;
	foreach my $tmp_pos ( keys %{$trans_to_uniprot->{$trans}->{'POSITION'}} ){
		if ( ( $position >= $tmp_pos ) 
		and  ( $position <= $trans_to_uniprot->{$trans}->{'POSITION'}->{$tmp_pos}->{'TEND'} ) ) { #pos >= tpos and pos <= tend
			#upos = pos - tpos + ubegin
			$tmp_uniprot_position = $position - $tmp_pos + $trans_to_uniprot->{$trans}->{'POSITION'}->{$tmp_pos}->{'UBEGIN'};
			$tmp_hit_bool = 1; 
			last;
		} 
	}
	return ( $tmp_hit_bool , $tmp_uniprot_position );
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

sub getSites {
	my ( $this , $trans_to_uniprot ) = @_;
	my $sites = {};
	my $fh = new FileHandle;
	unless( $fh->open( $this->{'site_file'} ) ) { die "Could not open sites file\n" };
	my $i = 0; my %fth;
	while ( my $ll = $fh->getline ) {
	   if ( $ll =~ m/^Hugo_Symbol/ ) { chomp( $ll );
			%fth = map {($_, $i++)} split( /\t/, $ll );
			last;
	   }
	}
	unless (	defined( $fth{ "Hugo_Symbol" } ) 
			and defined( $fth{ "Position" } )								
			and defined( $fth{ $this->{ "transcript_id_header" } } )
			and defined( $fth{ "Feature" } )
	) {
		die "not a valid site file with gene, transcript, and position!\n";
	}
	my @cols = ( $fth{ "Hugo_Symbol" } ,
				 $fth{ "Position" } ,
				 $fth{ $this->{ "transcript_id_header" } } ,
				 $fth{ "Feature" }
			   );
	# reading file content
	while ( my $line = $fh->getline ) {
		chomp( $line );
		my ( $gene , $position , $transcript , $feature ) = (split /\t/, $line)[@cols];
		next if ( not exists $trans_to_uniprot->{$transcript} );
		my ( $tmp_hit_bool , $tmp_uniprot_position ) = $this->getPositionMatch( $trans_to_uniprot , $transcript , $position );
		next if ( $tmp_hit_bool == 0 );
		my $mutationKey = join( "\t", $gene , $transcript , $position , $tmp_uniprot_position , $feature );
		my $uniprotID = $trans_to_uniprot->{$transcript}->{'UNIPROT'};
		$sites->{ $uniprotID }{ $tmp_uniprot_position } = $mutationKey;
	}
	$fh->close();
	return $sites;
}

# get drugport database information
sub getDrugportInfo {
	my ( $this ) = @_;
	my $drugport_f = $this->{'drugport_file'};
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
	my ( $this, $mafHashref, $proximityOutPrefix, $drugportref , $sites ) = @_;
	my ( @pairResults , @cosmicclose , @roiclose , 
		 @drugport_target_results , @drugport_nontarget_results , 
		 @mutationSiteResults , @siteSiteResults );
	my $fh = new FileHandle;
	my $AA = new TGI::Mutpro::Preprocess::AminoAcid;
	foreach my $uniprotID ( keys %$mafHashref ) {
		my $uniprotf = "$proximityOutPrefix\/$uniprotID.ProximityFile.csv";
		next unless( -e $uniprotf and $fh->open($uniprotf) ); 
		while ( my $b = <$fh> ) {
			next if ( $b =~ /UniProt_ID/g );
			chomp( $b ); my @line = split /\t/, $b;
			my ( $uid1, $chain1, $pdbcor1, $offset1, $residue1, $domain1, $cosmic1, 
				 $uid2, $chain2, $pdbcor2, $offset2, $residue2, $domain2, $cosmic2, 
				 $proximityinfor ) = @line;
			my $uniprotcor1 = $pdbcor1 + $offset1;
			my $uniprotcor2 = $pdbcor2 + $offset2;
			my $lineardis = undef;
			if ( $uid1 eq $uid2 ) { $lineardis = abs($uniprotcor1 - $uniprotcor2) } else { $lineardis = "N\/A"; }
			#print $uniprotID."\t".$uid2."\t".$uniprotcor1."\t".$uniprotcor2."\t".$lineardis."\n";
			my $muPart1 = join( "\t" , @line[1,2,5,6] ); #chain, position, domain, cosmic
			my $muPart2 = join( "\t" , @line[8,9,12,13] );
			if ( defined $mafHashref->{$uniprotID}->{$uniprotcor1} ) {
				if ( defined $mafHashref->{$uid2}->{$uniprotcor2} ) {
					#warn "check AA - ".$residue1." - ".$residue2;
					if ( !( $AA->isAA( $residue1 ) ) || !( $AA->isAA( $residue2 ) ) ) { 
						#warn " - bad AA pair"."\n";
						next;
					}
					## close each other
					foreach my $mutationKey1 ( keys %{$mafHashref->{$uniprotID}->{$uniprotcor1}} ) {
						foreach my $mutationKey2 ( keys %{$mafHashref->{$uid2}->{$uniprotcor2}} ) {
							if ( $this->cutFiltering( $lineardis, $proximityinfor) ) {
								push( @pairResults , join( "\t" , $mutationKey1 , $muPart1 , $mutationKey2 , $muPart2 , $lineardis , $proximityinfor ) );
							}
						}
					} 
				} else { 
					## close to COSMIC/Domain | to do
					map { 
						my $t_item = join("\t", $_, $muPart1 , $muPart2 , $lineardis, $proximityinfor); 
						if ( $domain2 !~ /N\/A/ ) {
							if ( $this->cutFiltering( $lineardis, $proximityinfor) ) {
								push(@roiclose, $t_item);
							}
						};
						if ( $cosmic2 !~ /N\/A/ ) {
							if ( $this->cutFiltering( $lineardis, $proximityinfor) ) {
								push(@cosmicclose, $t_item); 
							}
						}; 
					} keys %{$mafHashref->{$uniprotID}->{$uniprotcor1}};
				}
			} else {
				if ( defined $mafHashref->{$uid2}->{$uniprotcor2} ) {
					map {  
						my $t_item = join("\t", $_, $muPart2 , $muPart1 , $lineardis, $proximityinfor); 
						if ( $domain1 !~ /N\/A/ ) { 
							if ( $this->cutFiltering( $lineardis, $proximityinfor) ) {
								push(@roiclose, $t_item); 
							}
						}; 
						if ( $cosmic1 !~ /N\/A/ ) {
							if ( $this->cutFiltering( $lineardis, $proximityinfor) ) {
								push(@cosmicclose, $t_item); 
							}
						} 
					} keys %{$mafHashref->{$uid2}->{$uniprotcor2}};
				}
			}
			if ( not $AA->isHOH( $residue1 ) and not $AA->isHOH( $residue2 ) ) {
				# drugport searching
				if ( $drugportref ) {
					#warn "bad AA pair: ".$residue1." - ".$residue2."\n";
					my %pdbs_hash = map{ my @t0 = split / /, $_; ($t0[1], 1) } split /\|/, $proximityinfor;
					my ( $real_chain1 ) = $chain1 =~ /\[(\w)\]/; my ( $real_chain2 ) = $chain2 =~ /\[(\w)\]/;
					my ( $e, $pdb ); 
					map { 
						$pdb = $_;
						if ( defined $drugportref->{'TARGET'}->{$pdb}->{$real_chain1}->{$pdbcor1} ) {
							if ( defined $mafHashref->{$uid2}->{$uniprotcor2} ) {
								map { 
									if ( $this->cutFiltering( $lineardis, $proximityinfor) ) {
										push( @drugport_target_results, join("\t", $pdb, $chain1, $pdbcor1, $drugportref->{'TARGET'}->{$pdb}->{$real_chain1}->{$pdbcor1}, $_, @line[8,9,11,12,13], $lineardis, $proximityinfor) );
									}
								} keys %{$mafHashref->{$uid2}->{$uniprotcor2}};
							}
						}
						if ( defined $drugportref->{'TARGET'}->{$pdb}->{$real_chain2}->{$pdbcor2} ) {
							if ( defined $mafHashref->{$uid1}->{$uniprotcor1} ) {
								map { 
									if ( $this->cutFiltering( $lineardis, $proximityinfor) ) {
										push( @drugport_target_results, join("\t", $pdb, $chain2, $pdbcor2, $drugportref->{'TARGET'}->{$pdb}->{$real_chain2}->{$pdbcor2}, $_, @line[1,2,4,5,6], $lineardis, $proximityinfor) );
									} 
								} keys %{$mafHashref->{$uid1}->{$uniprotcor1}};
							}
						}
						if ( defined $drugportref->{'NONTARGET'}->{$pdb}->{$real_chain1}->{$pdbcor1} ) {
							if ( defined $mafHashref->{$uid2}->{$uniprotcor2} ) {
								map {
									if ( $this->cutFiltering( $lineardis, $proximityinfor) ) {
										push( @drugport_nontarget_results, join("\t", $pdb, $chain1, $pdbcor1, $drugportref->{'NONTARGET'}->{$pdb}->{$real_chain1}->{$pdbcor1}, $_, @line[8,9,11,12,13], $lineardis, $proximityinfor) ); 
									}
								} keys %{$mafHashref->{$uid2}->{$uniprotcor2}};
							}
						}
						if ( defined $drugportref->{'NONTARGET'}->{$pdb}->{$real_chain2}->{$pdbcor2} ) {
							if ( defined $mafHashref->{$uid1}->{$uniprotcor1} ) {
								map {
									if ( $this->cutFiltering( $lineardis, $proximityinfor) ) {
										push( @drugport_nontarget_results, join("\t", $pdb, $chain2, $pdbcor2, $drugportref->{'NONTARGET'}->{$pdb}->{$real_chain2}->{$pdbcor2}, $_, @line[1,2,4,5,6], $lineardis, $proximityinfor) );
									}
								} keys %{$mafHashref->{$uid1}->{$uniprotcor1}};
							}
						}
					} keys %pdbs_hash;
				} #if drugport
				if ( $sites ) {
					if ( $this->siteExists( $sites , $uid1 , $uniprotcor1 ) ) {
						#print "a site1: ".$sites->{$uid1}->{$uniprotcor1}."\n"; 
						my ( $gene , $transcript , $position , $uposition , $type ) = split( /\t/ , $sites->{$uid1}->{$uniprotcor1} );
						my $res = "p.".$AA->convertNameToSingle( $residue1 ).$uposition;
						my $sitePart1 = join( "\t" , $gene , $transcript , $position , $res , @line[1,2] , $type , $line[6] );
						if ( $this->siteExists( $sites , $uid2 , $uniprotcor2 ) ) { #site-site
							( $gene , $transcript , $position , $uposition , $type ) = split( /\t/ , $sites->{$uid2}->{$uniprotcor2} );
							$res = "p.".$AA->convertNameToSingle( $residue2 ).$uposition;
							my $sitePart2 = join( "\t" , $gene , $transcript , $position , $res , @line[8,9] , $type , $line[13] );
							push( @siteSiteResults , join( "\t" , ( $sitePart1 , $sitePart2 , $lineardis , $proximityinfor ) ) );
						}
						if ( defined $mafHashref->{$uid2}->{$uniprotcor2} ) {
							my $muKeys = join( "|" , keys( %{$mafHashref->{$uid2}->{$uniprotcor2}} ) ); 
							print "mutation: ".join( "\t" , ( $muKeys, $muPart2 , $sitePart1 , $lineardis , $proximityinfor ) );
							foreach my $mutationKey ( keys %{$mafHashref->{$uid2}->{$uniprotcor2}} ) {
								if ( $this->cutFiltering( $lineardis , $proximityinfor ) ) {
									push( @mutationSiteResults , join( "\t" , ( $mutationKey , $muPart2 , $sitePart1 , $lineardis , $proximityinfor ) ) );
								}
							}
						}
					} elsif ( $this->siteExists( $sites , $uid2 , $uniprotcor2 ) ) {
						if ( defined $mafHashref->{$uid1}->{$uniprotcor1} ) {
							#print "b site2: ".$sites->{$uid2}->{$uniprotcor2}."\n"; 
							my ( $gene , $transcript , $position , $uposition , $type ) = split( /\t/ , $sites->{$uid2}->{$uniprotcor2} );
							my $res = "p.".$AA->convertNameToSingle( $residue2 ).$uposition;
							my $sitePart = join( "\t" , $gene , $transcript , $position , $res , @line[8,9] , $type , $line[13] );
							foreach my $mutationKey ( keys %{$mafHashref->{$uid2}->{$uniprotcor2}} ) {
								if ( $this->cutFiltering( $lineardis , $proximityinfor ) ) {
									push( @mutationSiteResults , join( "\t" , ( $mutationKey , $muPart2 , $sitePart , $lineardis , $proximityinfor ) ) );
								}
							}
						}
					}
				} #if sites
			} #else { } #solvent accessibile
		}
		$fh->close();
	}

	return ( \@pairResults , \@cosmicclose , \@roiclose , 
			 \@drugport_target_results , \@drugport_nontarget_results ,
			 \@siteSiteResults , \@mutationSiteResults );
}

sub siteExists {
	my ( $this , $sites , $uid , $uposition ) = @_;
	if ( exists $sites->{$uid} ) {
		if ( exists $sites->{$uid}->{$uposition} ) {
			return 1;
		}
	}
	return 0;
}

# post processing of drug results
sub drug_proximity_postprocessing {
	my ( $this, $output_prefix, $drugport_parsing_results ) = @_;
	# post processing like collapsed clean results
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
	my ( $this, $linear_dis, $info_proximity ) = @_;
	my @infos = split( /\|/ , $info_proximity );
	my ( $dis_3d, $pvalue ) = (split / /, $infos[0])[0,2];
	if ( $linear_dis =~ /N\/A/ ) {
		return 1 if ( ( $dis_3d <= $this->{'3d_distance_cutoff'} ) 
				 and ( $pvalue <= $this->{'p_value_cutoff'} ) );
	} else {
		return 1 if ( ( $dis_3d <= $this->{'3d_distance_cutoff'} ) 
				 and ( $linear_dis >= $this->{'linear_cutoff'} ) 
				 and ( $pvalue <= $this->{'p_value_cutoff'} ) );
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
--site-file                  Protein site file
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
