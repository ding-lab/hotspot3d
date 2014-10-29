package TGI::Mutpro::Preprocess::Homolog;
#
#----------------------------------
# $Authors: Beifang Niu 
# $Date: 2014-01-14 14:34:50 -0500 (Tue Jan 14 14:34:50 CST 2014) $
# $Revision:  $
# $URL: $
# $Doc: $ homology protein structure 
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

use TGI::Mutpro::Preprocess::Uniprot;

sub new {
    my $class = shift;
    my $this = {};

    $this->{_OUTPUT_DIR} = getcwd;
    $this->{_IDENTITY} = 0.3;
    $this->{_STAT} = undef;

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
        'identity=f' => \$this->{_IDENTITY},
        'help' => \$help,
    );
    if ( $help ) { print STDERR help_text(); exit 0; };
    unless( $options ) { die $this->help_text(); };
    unless( $this->{_OUTPUT_DIR} ) { 
        warn 'You must provide a output directory ! ', "\n";
        die $this->help_text();
    };
    unless( -e $this->{_OUTPUT_DIR} ) {
        warn 'output directory is not exist  ! ', "\n";
        die $this->help_text();
    };
    #### processing ####
    # get homology PDBs for uniprots without PDBs annotations
    my ( $pdbseqsDir, $UniprotIdFile, $pdbseqsFile, $outputFile, );
    $pdbseqsDir = "$this->{_OUTPUT_DIR}\/pdbsequences";
    $UniprotIdFile = "$this->{_OUTPUT_DIR}\/hugo.uniprot.pdb.csv";
    $pdbseqsFile = "$pdbseqsDir\/pdb_seqres.txt";
    unless( -e $pdbseqsDir ) { mkdir( $pdbseqsDir ) || die "can not make pdb sequences directory !\n"; };
    ## get pdbseqs file 
    my $url = 'ftp://ftp.rcsb.org/pub/pdb/derived_data/pdb_seqres.txt';
    getstore( $url, $pdbseqsFile );
    ## extract protein seqs
    #
    my %pdbseq_hash;
    my $pdbfh = new FileHandle;
    unless( $pdbfh->open("<$pdbseqsFile") ) { die "Could not open pdb sequences file !\n" };
    while (my $a = <$pdbfh>) {
        if ($a =~ /^>/) { 
            my ($name) = $a =~ /^>(\S+)\s/;
            my $b = <$pdbfh>;
            $pdbseq_hash{$name} = $b;
        }
    }
    $pdbfh->close();
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
        #
        #
        #
        #
        next if ($uniprotId !~ /\w+/ );
        next unless ( $pdb eq "N/A" );
        print STDERR "*** UNIPROT ID: $uniprotId\n";
        $allUniprotIds{$uniprotId} = 1;
        $uniprotRef = TGI::Mutpro::Preprocess::Uniprot->new($uniprotId);
        defined ($uniprotRef) || die "no object for '$uniprotId'";
        # The annotation is a ref to array made here:
        # 'push @domains, 
        # "$key\t($dmStart, $dmStop)\t$desc'";
        my $uniprotSequence = $uniprotRef->sequence();


        my $tfh = new FileHandle;
        unless( $tfh->open(">.temp_one_pseq") ) { die "Could not generate temp file !\n" };
        print $tfh ">$uniprotId\n$uniprotSequence\n";
        $tfh->close();

        my $t_blast_output = `blastall -p blastp -d $pdbseqsFile -i .temp_one_pseq -e 0.05`;

        #print $t_blast_output;
        #unlink(".temp_one_pseq");
        #
        #
        $this->parse_blastp_output( $t_blast_output, $uniprotId, $this->{_IDENTITY}, );

    }
}
#
# blastp parsing
sub parse_blastp_output {
    my ( $this, $blastp_output, $uniprotid, $iden_cutoff, ) = @_;

    my ( $f, @top, %homos, %header, $index, $pdb, $chain, $iden, );
    my ( $qstart, $qend, $qcont, $sstart, $send, $scont, ); 
    my ( @homoregions, );

    $f = 0; $index = 0; $iden = 0;
    foreach (split /\n/, $blastp_output) {
        next unless ( /^>/ || $f == 1 );
        last if ($f == 1 && /^>/ );
        push @top, $_; $f = 1;
    }

    return \@homoregions unless( @top );
    foreach (@top) {
        if ( /^>/ ) { 
            ($pdb, $chain) = /^>(\w+)\_(\w) /; 
            $header{'PDB'} = uc($pdb);
            $header{'CHAIN'} = uc($chain);
        }
        if ( /Identities/ ) { $index++; ($iden) = /\((\d+)\%\)/; $homos{$index}{'IDEN'} = $iden; }
        if ( /Query:/ ) {
            ( $qstart, $qcont, $qend ) = /Query:\s+(\d+)\s+(\S+)\s+(\d+)/; 
            unless (defined $homos{$index}{'QUESTART'}) { $homos{$index}{'QUESTART'} = $qstart; }
            $homos{$index}{'QUECONT'} .= $qcont;
            $homos{$index}{'QUEEND'} = $qend;
        }
        if ( /Sbjct:/ ) {
            ( $sstart, $scont, $send ) = /Sbjct:\s+(\d+)\s+(\S+)\s+(\d+)/; 
            unless (defined $homos{$index}{'SUBSTART'}) { $homos{$index}{'SUBSTART'} = $sstart; }
            $homos{$index}{'SUBCONT'} .= $scont;
            $homos{$index}{'SUBEND'} = $send;
        }
    }
    #map { print; print "\n"; } @top;
    #
    print $header{'PDB'}."_".$header{'CHAIN'}."\n";
    foreach my $d (keys %homos) { 
        next if ( $homos{$d}{'IDEN'} < $iden_cutoff * 100 );
        print "==".$d."\n";
        print $homos{$d}{'IDEN'}."\n";
        print $homos{$d}{'QUESTART'}." ".$homos{$d}{'QUECONT'}." ".$homos{$d}{'QUEEND'}."\n";      
        print $homos{$d}{'SUBSTART'}." ".$homos{$d}{'SUBCONT'}." ".$homos{$d}{'SUBEND'}."\n";
        my ( @taq, @tas, $tqstart, $tsstart, $i, $j, $k, $open, );
        $tqstart = $i = $homos{$d}{'QUESTART'}; 
        $tsstart = $j = $homos{$d}{'SUBSTART'};
        @taq = split //, $homos{$d}{'QUECONT'};
        @tas = split //, $homos{$d}{'SUBCONT'};
        $k = $open = 0; 
        foreach my $e ( @taq ) {
            if ( $e ne '-' && $tas[$k] ne '-' ) {
                if ($open == 1) {
                    $tqstart = $i;
                    $tsstart = $j;
                    $open  = 0;
                }
                $i++; $j++; $k++;
            } else {
                if ( $e eq '-') {
                    if ( $open == 0 ) {
                        my $tt0 = $i - 1;
                        my $tt1 = $j - 1;
                        push( @homoregions, "DBREF\t$header{'PDB'}\t$header{'CHAIN'}\t$tsstart\t$tt1\tUNP\t$uniprotid\tB2MG_HUMAN\t$tqstart\t$tt0\n" );
                        print "DBREF\t$header{'PDB'}\t$header{'CHAIN'}\t$tsstart\t$tt1\tUNP\t$uniprotid\tB2MG_HUMAN\t$tqstart\t$tt0\n";
                        #print  "$tqstart,$tt0 -- $tsstart,$tt1\n";
                        $j++; $k++; $open = 1;
                    } else { $j++; $k++; }
                } 
                if ( $tas[$k] eq '-') {
                    if ( $open == 0 ) {
                        my $tt0 = $i - 1;
                        my $tt1 = $j - 1;
                        push( @homoregions, "DBREF\t$header{'PDB'}\t$header{'CHAIN'}\t$tsstart\t$tt1\tUNP\t$uniprotid\tB2MG_HUMAN\t$tqstart\t$tt0\n" );
                        print "DBREF\t$header{'PDB'}\t$header{'CHAIN'}\t$tsstart\t$tt1\tUNP\t$uniprotid\tB2MG_HUMAN\t$tqstart\t$tt0\n";
                        #print  "$tqstart,$tt0 -- $tsstart,$tt1\n";
                        $i++; $k++; $open = 1;
                    } else { $i++; $k++; }
                } 
            }
        }

        if ($open == 0) {
            my $tt0 = $i - 1;
            my $tt1 = $j - 1;
            push( @homoregions, "DBREF\t$header{'PDB'}\t$header{'CHAIN'}\t$tsstart\t$tt1\tUNP\t$uniprotid\tB2MG_HUMAN\t$tqstart\t$tt0\n" );
            print "DBREF\t$header{'PDB'}\t$header{'CHAIN'}\t$tsstart\t$tt1\tUNP\t$uniprotid\tB2MG_HUMAN\t$tqstart\t$tt0\n";
            #print  "$tqstart,$tt0 -- $tsstart,$tt1\n";
        }
    }
    #
    #
    #
    return \@homoregions;

}

sub help_text {
    my $this = shift;
        return <<HELP

Usage: hotspot3d homo [options]

--output-dir		Output directory of proximity files
--identity              Identities cutoff, default is 30%

--help			this message

HELP

}

1;

