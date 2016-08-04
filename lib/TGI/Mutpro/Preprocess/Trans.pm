package TGI::Mutpro::Preprocess::Trans;
#
#----------------------------------
# $Authors: Beifang Niu 
# $Date: 2014-01-14 14:34:50 -0500 (Tue Jan 14 14:34:50 CST 2014) $
# $Revision:  $
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

use TGI::Mutpro::Preprocess::Uniprot;

sub new {
    my $class = shift;
    my $this = {};
    $this->{_OUTPUT_DIR} = getcwd;
    $this->{_STAT} = undef;
    $this->{_BLAT} = undef;
    bless $this, $class;
    $this->process();
    return $this;
}

sub process {
    my $this = shift;
    my ( $help, $options );
    unless( @ARGV ) { die $this->help_text(); };
    $options = GetOptions (
        'output-dir=s' => \$this->{_OUTPUT_DIR},
        'blat=s' => \$this->{_BLAT},
        'help' => \$help,
    );
    if ( $help ) { print STDERR help_text(); exit 0; };
    unless( $options ) { die $this->help_text(); };
    unless( $this->{_OUTPUT_DIR} ) { warn 'You must provide a output directory ! ', "\n"; die $this->help_text(); };
    unless( -e $this->{_OUTPUT_DIR} ) { warn 'output directory is not exist  ! ', "\n"; die $this->help_text(); };
    #### processing ####
    # add transcript annotation for uniprot
    my ( $peptidesDir, $UniprotIdFile, $peptidesFile, $outputFile, );
    $peptidesDir = "$this->{_OUTPUT_DIR}\/humanPeptides";
    $UniprotIdFile = "$this->{_OUTPUT_DIR}\/hugo.uniprot.pdb.csv";
    $peptidesFile = "$peptidesDir\/Homo_sapiens.GRCh37.74.pep.all.fa";
    $outputFile = "$this->{_OUTPUT_DIR}\/hugo.uniprot.pdb.transcript.csv";
    unless( -e $peptidesDir ) { mkdir( $peptidesDir ) || die "can not make peptides directory !\n"; };
    ## get peptide file 
    my $url = 'ftp://ftp.ensembl.org/pub/release-74/fasta/homo_sapiens/pep/Homo_sapiens.GRCh37.74.pep.all.fa.gz';
    my $downloadFile = "$peptidesDir\/Homo_sapiens.GRCh37.74.pep.all.fa.gz";
    getstore( $url, $downloadFile );
    system( "gzip -d $downloadFile" );
    # load peptide seqs
    my $peptideSeqsRef = $this->loadPeptides( $peptidesFile );
    my ( @entireFile, $uniprotId, %allUniprotIds, $uniprotRef, $annotationRef, $start, $stop, $key, $desc, $entry, $pdb, );
    my $fhuid = new FileHandle;
    unless( $fhuid->open("<$UniprotIdFile") ) { die "Could not open uniprot id file !\n" };
    @entireFile = <$fhuid>;
    $fhuid->close();
    my $outputContent = "";
    foreach my $line (@entireFile) {
        chomp $line;
        ( undef, $uniprotId, $pdb, ) = split /\s+/, $line;
        # Only use Uniprot IDs with PDB structures
        next if ( $pdb eq "N/A" || $uniprotId !~ /\w+/ );
        $allUniprotIds{$uniprotId} = 1;
        $uniprotRef = TGI::Mutpro::Preprocess::Uniprot->new($uniprotId);
        defined ($uniprotRef) || die "no object for '$uniprotId'";
        # The annotation is a ref to array made here:
        my $uniprotSequence = $uniprotRef->sequence();

        #print ">uniprotseq\n";
        #print $uniprotSequence."\n";

        my $transProtein = $uniprotRef->transProteinHash();
        next unless( defined $transProtein );
        my $transcriptContent = "";
        foreach my $transcript ( keys %{$transProtein} ) {
            #print ">$transcript\n";
            my $proteinid = $$transProtein{$transcript};
            my $proteinsequence = $$peptideSeqsRef{$proteinid};
            #print "$proteinsequence\n";
            next unless(defined $proteinsequence);
            if ( $uniprotSequence eq $proteinsequence ) {
                $transcriptContent .= $transcript."[1|1-".length($proteinsequence)."|".length($proteinsequence)."],"; 
            } else {
                my ( undef, $tmp_uniprot_seq_file ) = tempfile();
                my $tmp_uniprot_fh = IO::File->new( $tmp_uniprot_seq_file, ">" ) or die "Temporary file could not be created. $!";
                $tmp_uniprot_fh->print( ">$uniprotId\n$uniprotSequence\n" );
                my ( undef, $tmp_transcript_protein_seq_file ) = tempfile();
                my $tmp_transcript_fh = IO::File->new( $tmp_transcript_protein_seq_file, ">" ) or die "Temporary file could not be created. $!";
                $tmp_transcript_fh->print( ">$transcript\n$proteinsequence\n" );
                $tmp_uniprot_fh->close; $tmp_transcript_fh->close;
                my ( undef, $tmp_blat_output_file ) = tempfile();
				my $blat = "blat";
				if ( $this->{_BLAT} ) { $blat = $this->{_BLAT}; }
                system( "$blat $tmp_uniprot_seq_file $tmp_transcript_protein_seq_file -t=prot -q=prot -out=blast $tmp_blat_output_file" );
                # parse blat output
                my $tmp_parse_cont = ""; 
                map{ $tmp_parse_cont .= $_.":"; } @{$this->parse_blat_output( $tmp_blat_output_file, $uniprotId, 0.90 )};
                if ( $tmp_parse_cont ne "" ) { chop( $tmp_parse_cont ); $transcriptContent .= $transcript."[".$tmp_parse_cont."],"; };
                # clean files
                unlink $tmp_uniprot_seq_file; unlink $tmp_transcript_protein_seq_file; unlink $tmp_blat_output_file; 

            }
        }
        if ( $transcriptContent eq "" ) { $transcriptContent = "N\/A";
        } else { chop( $transcriptContent ); }
        print "$line\t$transcriptContent\n";
        $outputContent .= "$line\t$transcriptContent\n";
    }
    my $fhout = new FileHandle;
    unless( $fhout->open(">$outputFile") ) { die "Could not open output file !\n" };
    print $fhout $outputContent;
    $fhout->close();
}

# loading peptides
sub loadPeptides {
    my ( $this, $pepfile, ) = @_;
    my $fh = new FileHandle;
    unless( $fh->open("<$pepfile") ) { die "Could not open peptide file !\n" };
    my @entireFile = <$fh>;
    $fh->close();
    my ( %pep_hash, $content, $name, );
    $content = $name = "";
    foreach my $a (@entireFile) {
        if ($a =~ /^>/) {
            $pep_hash{$name} = $content if ($name ne "");
            ($name) = $a =~ /^>(\w+) /;
            $content = "";
        } else { chomp($a); $content .= $a; }
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
    unless( $fh->open( "<$blat_output" ) ) { die "Could not open blat output file !\n" };
    foreach ( $fh->getlines ) { chomp; next unless ( /^>/ || $f == 1 ); last if ($f == 1 && /^>/ ); push @top, $_; $f = 1; };
    $fh->close;
    return \@homoregions unless( @top );
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
        print $homos{$d}{'IDEN'}."\n";
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

--help                       this message

HELP

}

1;

