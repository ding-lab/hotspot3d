package TGI::Mutpro::Preprocess::Trans;
#
#----------------------------------
# $Authors: Beifang Niu & Adam D Scott
# $Date: 2014-01-14 14:34:50 -0500 (Tue Jan 14 14:34:50 CST 2014) $
# $Revision: 3 $
# $URL: $
# $Doc: $ transcripts processing and added them in table 
#----------------------------------
#
use strict;
use warnings;

use Carp;
use Cwd;
use Getopt::Long;
use LWP::Simple;
use IO::File;
use FileHandle;
use File::Temp qw/ tempfile /;
use Archive::Extract;

use TGI::Mutpro::Preprocess::Uniprot;

my $LATESTGRCH = 38;
my $LATEST38RELEASE = 87; #TODO will need to be updated as Ensembl has new releases
my $EARLIEST38RELEASE = 76;
my $LATEST37RELEASE = 75;
my $EARLIEST37RELEASE = 55;

sub new {
    my $class = shift;
    my $this = {};
    $this->{_OUTPUT_DIR} = getcwd;
    $this->{_STAT} = undef;
    $this->{_BLAT} = "blat";
    $this->{GRCh} = undef;
    $this->{release} = undef;
    bless $this, $class;
    $this->process();
    return $this;
}

sub process {
    my $this = shift;
	$this->setOptions();
    # add transcript annotation for uniprot
	my $entireFile = $this->getInputFile( );
	my $fhout = $this->getOutputFileHandle( );
	my $peptides = $this->getPeptides( );
	$this->mapEnsemblTranscriptsToUniprot( $entireFile , $peptides , $fhout );
	return 0;
}

sub mapEnsemblTranscriptsToUniprot {
	my ( $this , $entireFile , $peptides , $fhout ) = @_;
	my $outputContent = "";
    foreach my $line ( @{$entireFile} ) {
        chomp $line;
		my ( $uniprotId , $pdb );
        ( undef, $uniprotId, $pdb, ) = split /\s+/, $line;
        # Only use Uniprot IDs with PDB structures
        next if ( $pdb eq "N/A" || $uniprotId !~ /\w+/ );
		my $transcriptContent = $this->mapEnsemblTranscriptsToThisUniprot( $uniprotId , $peptides );
		$outputContent .= $line."\t".$transcriptContent."\n";
    }
	print STDOUT "Writing to hupt\n";
    print $fhout $outputContent;
    $fhout->close();
	return;
}

sub mapEnsemblTranscriptsToThisUniprot {
	my ( $this , $uniprotId , $peptides ) = @_;
	my $uniprotRef = TGI::Mutpro::Preprocess::Uniprot->new($uniprotId);
	defined ($uniprotRef) || die "HotSpot3D::Trans::mapEnsemblTranscriptsToThisUniprot error: no object for '$uniprotId'";
	# The annotation is a ref to array made here:
	my $uniprotSequence = $uniprotRef->sequence();

	#print ">uniprotseq\n";
	#print $uniprotSequence."\n";

	my $enst2ensp = $uniprotRef->transProteinHash(); #non-versioned ENST ids
	next unless( defined $enst2ensp );
	my $transcriptContent = "";
	print STDOUT $uniprotId." HotSpot3D::Trans::mapEnsemblTranscriptsToThisUniprot - determining Ensembl->UniProt mapping\n";
	foreach my $enst ( keys %{$enst2ensp} ) { #non-versioned ENST id
		#print ">$enst\n";
		my $ensp = $enst2ensp->{$enst};
		my $proteinSequence = $this->getVersionedSequence( $peptides , $ensp ); #
		#print "$proteinSequence\n";
		next unless(defined $proteinSequence);
		if ( $uniprotSequence eq $proteinSequence ) {
			$transcriptContent .= $enst."[1|1-".length($proteinSequence)."|".length($proteinSequence)."],"; 
		} else {
			my ( undef, $tmp_uniprot_seq_file ) = tempfile();
			my $tmp_uniprot_fh = IO::File->new( $tmp_uniprot_seq_file, ">" ) or die "HotSpot3D::Trans::mapEnsemblTranscriptsToThisUniprot error: Temporary file could not be created. $!";
			$tmp_uniprot_fh->print( ">$uniprotId\n$uniprotSequence\n" );
			my ( undef, $tmp_transcript_protein_seq_file ) = tempfile();
			my $tmp_transcript_fh = IO::File->new( $tmp_transcript_protein_seq_file, ">" ) or die "HotSpot3D::Trans::mapEnsemblTranscriptsToThisUniprot error: Temporary file could not be created. $!";
			$tmp_transcript_fh->print( ">$enst\n$proteinSequence\n" );
			$tmp_uniprot_fh->close; $tmp_transcript_fh->close;
			my ( undef, $tmp_blat_output_file ) = tempfile();
			my $blat = "blat";
			system( "$this->{_BLAT} $tmp_uniprot_seq_file $tmp_transcript_protein_seq_file -t=prot -q=prot -out=blast $tmp_blat_output_file" );
			# parse blat output
			my $tmp_parse_cont = ""; 
			map{ $tmp_parse_cont .= $_.":"; } @{$this->parse_blat_output( $tmp_blat_output_file, $uniprotId, 0.90 )};
			if ( $tmp_parse_cont ne "" ) { chop( $tmp_parse_cont ); $transcriptContent .= $enst."[".$tmp_parse_cont."],"; };
			# clean files
			unlink $tmp_uniprot_seq_file;
			unlink $tmp_transcript_protein_seq_file;
			unlink $tmp_blat_output_file; 

		}
	}
	if ( $transcriptContent eq "" ) {
		print STDOUT $uniprotId." HotSpot3D::Trans::mapEnsemblTranscriptsToThisUniprot - no mappings\n";
		$transcriptContent = "N\/A";
	} else {
		print STDOUT $uniprotId." HotSpot3D::Trans::mapEnsemblTranscriptsToThisUniprot - complete\n";
		chop( $transcriptContent );
	}
	return $transcriptContent;
}

sub getPeptidesFile {
	my $this = shift;
    my $peptidesDir = "$this->{_OUTPUT_DIR}\/humanPeptides";
    unless( -e $peptidesDir ) { mkdir( $peptidesDir ) || die "HotSpot3D::Trans::getPeptidesFile error: can not make peptides directory!\n"; };
    my $peptidesFile = "$peptidesDir\/Homo_sapiens.GRCh".$this->{GRCh}.".".$this->{release}.".pep.all.fa";
	return $peptidesFile;
}
    
sub getPeptides {
	my ( $this ) = @_;
	my $peptidesFile = $this->getPeptidesFile( );
    ## get peptide file 
	my $url = $this->makeEnsemblFastaURL();
    my $downloadFile = $peptidesFile.".gz";
	my $decompressor = Archive::Extract->new( 'archive' => $downloadFile );
	if ( not -e $downloadFile ) {
		getstore( $url, $downloadFile );
	}
	$decompressor->extract( to => $peptidesFile );
    #system( "gzip -d $downloadFile" );
    # load peptide seqs
    return $this->loadPeptides( $peptidesFile );
}
    
sub getInputFile {
	my $this = shift;
	my $UniprotIdFile = "$this->{_OUTPUT_DIR}\/hugo.uniprot.pdb.csv";
    my $fhuid = new FileHandle;
    unless( $fhuid->open("<$UniprotIdFile") ) { die "HotSpot3D::Trans::getInputFile error: Could not open uniprot id file!\n" };
    my @entireFile = <$fhuid>;
    $fhuid->close();
	return \@entireFile;
}

sub getOutputFileHandle {
	my $this = shift;
	my $outputFile = "$this->{_OUTPUT_DIR}\/hugo.uniprot.pdb.transcript.csv";
    my $fhout = new FileHandle;
    unless( $fhout->open(">$outputFile") ) { die "HotSpot3D::Trans::getOutputFileHandle error: Could not open output file!\n" };
	print STDOUT "Creating hupt: ".$outputFile."\n";
	return $fhout;
}

sub setOptions {
	my ( $this ) = @_;
    my ( $help, $options );
    unless( @ARGV ) { die $this->help_text(); };
    $options = GetOptions (
        'output-dir=s' => \$this->{_OUTPUT_DIR},
        'blat=s' => \$this->{_BLAT},
        'grch=f' => \$this->{GRCh},
        'release=f' => \$this->{release},
        'help' => \$help,
    );
    if ( $help ) { warn help_text(); exit 0; };
    unless( $options ) { die $this->help_text(); };
    unless( $this->{_OUTPUT_DIR} ) { warn 'HotSpot3D::Trans::setOptions error: You must provide an output directory! ', "\n"; die $this->help_text(); };
    unless( -e $this->{_OUTPUT_DIR} ) { warn 'HotSpot3D::Trans::setOptions error: output directory does not exist! ', "\n"; die $this->help_text(); };
	if ( $this->checkBLAT() ) { warn 'HotSpot3D::Trans::setOptions error: blat not found - $this->{_BLAT}!' , "\n"; die $this->help_text(); }
	$this->setEnsemblVersions();
	return;
}

sub checkBLAT {
	my ( $this ) = @_;
	return system( "which $this->{_BLAT}" );
}

sub setEnsemblVersions {
	my ( $this ) = @_;
	if ( not defined $this->{GRCh} and not defined $this->{release} ) {
		warn "HotSpot3D::Trans::setEnsemblVersions warning: Ensembl GRCh & release versions not specified, defaulting to latest: GRCh".$LATESTGRCH.".".$LATEST38RELEASE."\n";
		$this->{GRCh} = $LATESTGRCH;
		$this->{release} = $LATEST38RELEASE;
	} elsif ( not defined $this->{GRCh} ) {
		if ( $this->{release} >= $EARLIEST38RELEASE and $this->{release} <= $LATEST38RELEASE ) {
			$this->{GRCh} = 38;
		} elsif ( $this->{release} <= $LATEST37RELEASE and $this->{release} >= $EARLIEST37RELEASE ) {
			$this->{GRCh} = 37;
		} elsif ( $this->{release} <= 54 ) {
			warn "HotSpot3D::Trans::setEnsemblVersions error: Ensembl releases <= 54 not supported\n";
			die $this->help_text();
		}
		warn "HotSpot3D::Trans::setEnsemblVersions warning: Ensembl GRCh version not specified. User specified release ".$this->{release}.", so using GRCh".$this->{GRCh}."\n";
	} elsif ( not defined $this->{release} ) {
		if ( $this->{GRCh} == 38 ) {
			warn "HotSpot3D::Trans::setEnsemblVersions warning: Ensembl release versions not specified, defaulting to latest: GRCh".$this->{GRCh}.".".$LATEST38RELEASE."\n";
			$this->{release} = $LATEST38RELEASE;
		} elsif ( $this->{GRCh} == 37 ) {
			warn "HotSpot3D::Trans::setEnsemblVersions warning: Ensembl release versions not specified, defaulting to latest: GRCh".$this->{GRCh}.".".$LATEST37RELEASE."\n";
			$this->{release} = 75;
		} else {
			warn "HotSpot3D::Trans::setEnsemblVersions error: Unrecognized GRCh version (user set ".$this->{GRCh}.", not 37 or 38)\n";
			die $this->help_text();
		}
	} else {
	}
}
    #### processing ####
sub getVersionedSequence {
	my ( $this , $peptideSeqs , $ensp ) = @_;
	foreach my $id ( keys %{$peptideSeqs} ) {
		if ( $id =~ /$ensp/ ) {
			return $peptideSeqs->{$id};
		}
	}
	return undef;
}

sub makeEnsemblFastaURL {
	my ( $this ) = @_;
	my $url = "ftp://ftp.ensembl.org/pub/release-";
	$url .= $this->{release};
	$url .= "/fasta/homo_sapiens/pep/Homo_sapiens.GRCh";
	$url .= $this->{GRCh};
	if ( $this->{GRCh} == 37 ) {
		$url .= ".";
		$url .= $this->{release};
	}
	$url .= ".pep.all.fa.gz";
	return $url;
}

# loading peptides
sub loadPeptides {
    my ( $this, $pepfile, ) = @_;
    my $fh = new FileHandle;
    unless( $fh->open("<$pepfile") ) { die "HotSpot3D::Trans::loadPeptides error: Could not open peptide file!\n" };
    my @entireFile = <$fh>;
    $fh->close();
    my ( %pep_hash, $content, $name, );
    $content = $name = "";
    foreach my $line (@entireFile) {
        if ($line =~ /^>/) {
			#print $line."\n";
            $pep_hash{$name} = $content if ($name ne "");
			if ( $this->{GRCh} == 38 ) {
				($name) = $line =~ /^>(\w+\.\d+) /;
			} else {
				($name) = $line =~ /^>(\w+) /;
			}
            $content = "";
        } else { chomp($line); $content .= $line; }
    }
    # last seq 
    if ($content ne "") { $pep_hash{$name} = $content; }
    return \%pep_hash;
}

# blat parsing
sub parse_blat_output {
    my ( $this, $blat_output, $uniprotid, $iden_cutoff, ) = @_;
    my ( $f, @top, %homos, %header, $index, $iden, );
    my ( $qstart, $qend, $qcont, $sstart, $send, $scont, ); 
    my ( @homoregions, );
    $f = $index = $iden = 0;
    my $fh = new FileHandle;
    unless( $fh->open( "<$blat_output" ) ) { die "HotSpot3D::Trans::parse_blat_output error: Could not open blat output file!\n" };
    foreach ( $fh->getlines ) { chomp; next unless ( /^>/ || $f == 1 ); last if ($f == 1 && /^>/ ); push @top, $_; $f = 1; };
    $fh->close;
    return \@homoregions unless( @top );
	print STDOUT $uniprotid." HotSpot3D::Trans::parse_blat_output - parsing BLAT output\n";
    foreach (@top) {
        #if ( /^>/ ) { ($pdb, $chain) = /^>(\w+)\_(\w) /; $header{'PDB'} = uc($pdb); $header{'CHAIN'} = uc($chain); }
        if ( /Identities/ ) { $index++; ($iden) = /\((\d+)\%\)/; $homos{$index}{'IDEN'} = $iden; }
        if ( /Query:/ ) { ( $qstart, $qcont, $qend ) = /Query:\s+(\d+)\s+(\S+)\s+(\d+)/; 
            unless (defined $homos{$index}{'QUESTART'}) { $homos{$index}{'QUESTART'} = $qstart; }
            $homos{$index}{'QUECONT'} .= $qcont;
            $homos{$index}{'QUEEND'} = $qend;
        }
        if ( /Sbjct:/ ) { ( $sstart, $scont, $send ) = /Sbjct:\s+(\d+)\s+(\S+)\s+(\d+)/; 
            unless (defined $homos{$index}{'SUBSTART'}) { $homos{$index}{'SUBSTART'} = $sstart; }
            $homos{$index}{'SUBCONT'} .= $scont;
            $homos{$index}{'SUBEND'} = $send;
        }
    }
    foreach my $d (keys %homos) { 
        next if ( $homos{$d}{'IDEN'} < $iden_cutoff * 100 );
        print STDOUT $homos{$d}{'IDEN'}."\n";
        #print $homos{$d}{'QUESTART'}." ".$homos{$d}{'QUECONT'}." ".$homos{$d}{'QUEEND'}."\n";      
        #print $homos{$d}{'SUBSTART'}." ".$homos{$d}{'SUBCONT'}." ".$homos{$d}{'SUBEND'}."\n";
        my ( @taq, @tas, $tqstart, $tsstart, $i, $j, $k, $open, );
        $tqstart = $i = $homos{$d}{'QUESTART'}; 
        $tsstart = $j = $homos{$d}{'SUBSTART'};
        @taq = split //, $homos{$d}{'QUECONT'};
        @tas = split //, $homos{$d}{'SUBCONT'};
        $k = $open = 0; 
        foreach my $e ( @taq ) {
            if ( $e ne '-' && $tas[$k] ne '-' ) {
                if ( $open == 1 ) { $tqstart = $i; $tsstart = $j; $open  = 0; }
                $i++; $j++; $k++;
            } else {
                if ( $e eq '-') {
                    if ( $open == 0 ) {
                        my $tt0 = $i - 1; my $tt1 = $j - 1;
                        push( @homoregions, "$tsstart|$tqstart-$tt1|$tt0" );
                        #print  "$tqstart,$tt0 -- $tsstart,$tt1\n";
                        $j++; $k++; $open = 1;
                    } else { $j++; $k++; }
                } 
                if ( $tas[$k] eq '-') {
                    if ( $open == 0 ) {
                        my $tt0 = $i - 1; my $tt1 = $j - 1;
                        push( @homoregions, "$tsstart|$tqstart-$tt1|$tt0" );
                        #print  "$tqstart,$tt0 -- $tsstart,$tt1\n";
                        $i++; $k++; $open = 1;
                    } else { $i++; $k++; }
                } 
            }
        }
        if ( $open == 0 ) {
            my $tt0 = $i - 1; my $tt1 = $j - 1;
            push( @homoregions, "$tsstart|$tqstart-$tt1|$tt0" );
            #print  "$tqstart,$tt0 -- $tsstart,$tt1\n";
        }
    }
    #
    return \@homoregions;
}

sub help_text {
    my $this = shift;
        return <<HELP

Usage: hotspot3d trans [options]

                             REQUIRED
--output-dir                 Output directory of proximity files

                             OPTIONAL
--blat                       Installation of blat to use (defaults to your system default)
--grch                       Genome build (37 or 38), defaults to 38 or according to --release input
--release                    Ensembl release verion (55-87), defaults to 87 or to the latest release according to --grch input
                                 Note that releases 55-75 correspond to GRCh37 & 76-87 correspond to GRCh38

--help                       this message

HELP

}

1;

