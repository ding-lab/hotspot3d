package TGI::Mutpro::Main::Post;
#
#----------------------------------
# $Authors: Beifang Niu 
# $Date: 2015-05-07 14:34:50 -0500 (Thu May  7 15:43:42 CDT 2015) $
# $Revision:  $
# $URL: $
# $Doc: $ proximity searching results postprocessing
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

sub new {
    my $class = shift;
    my $this = {};
    $this->{'maf'} = undef;
    $this->{'input_prefix'} = '3D_Proximity';
	$this->{'transcript_id_header'} = "transcript_name";
	$this->{'amino_acid_header'} = "amino_acid_change";
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
        'input-prefix=s' => \$this->{'input_prefix'},
        'transcript-id-header=s' => \$this->{'transcript_id_header'},
        'amino-acid-header=s' => \$this->{'amino_acid_header'},
        'help' => \$help,
    );
    if ( $help ) { print STDERR help_text(); exit 0; }
    unless( $options ) { die $this->help_text(); }
    my $input_pairwise_file = "$this->{'input_prefix'}.pairwise";
    unless( $this->{'maf'} and (-e $this->{'maf'}) ) { warn 'You must provide a MAF format file ! ', "\n"; die $this->help_text(); }
    unless( -e $input_pairwise_file ) { warn 'You must provide a 3D proximity pairwise result file ! ', "\n"; die $this->help_text(); }
    my $fh   = new FileHandle;
    die "Could not open pairwise file !\n" unless( $fh->open( $input_pairwise_file ) );
    my %ss = ();
    map {
        chomp; my @t = (split /\t/)[0,4,9,13,18,19];
        if ($t[4] =~ /N\/A/) {
            my @tt = split /\|/, $t[5];
            my $ld = $t[4];
            my $pvalue = (split / /, $tt[0])[-1];
            if ( defined $ss{"$t[0]\t$t[1]"}{"$t[2]\t$t[3]"} ) {
                my $p = $ss{"$t[0]\t$t[1]"}{"$t[2]\t$t[3]"}; 
                my ($pld, $p3dd) = split /\t/, $p;
                ## bug in one line
                #if ($p3dd < $pvalue) { 
                #    $ss{"$t[0]\t$t[1]"}{"$t[2]\t$t[3]"} = $ld."\t".$tt[0]; 
                #} 
                ## change to the following
                my $p3dd_pvalue = (split / /, $p3dd)[-1];
                if ($p3dd_pvalue > $pvalue) {
                    $ss{"$t[0]\t$t[1]"}{"$t[2]\t$t[3]"} = $ld."\t".$tt[0]; 
                }
            }elsif ( defined $ss{"$t[2]\t$t[3]"}{"$t[0]\t$t[1]"} ) { 
                my $p = $ss{"$t[2]\t$t[3]"}{"$t[0]\t$t[1]"}; 
                my ($pld, $p3dd) = split /\t/, $p; 
                ## bug in one line
                #if ($p3dd < $pvalue) { 
                #    $ss{"$t[2]\t$t[3]"}{"$t[0]\t$t[1]"} = $ld."\t".$tt[0]; 
                #}
                ## change to the following
                my $p3dd_pvalue = (split / /, $p3dd)[-1];
                if ($p3dd_pvalue > $pvalue) { 
                    $ss{"$t[2]\t$t[3]"}{"$t[0]\t$t[1]"} = $ld."\t".$tt[0]; 
                }
            }else { 
                $ss{"$t[0]\t$t[1]"}{"$t[2]\t$t[3]"} = $ld."\t".$tt[0]; 
            }
        }
    } <$fh>;
    $fh->close;
    my @complex_tmp1 = (); 
    my @complex_tmp2 = (); 
    map { my $a = $_; map { my $b = $_; my $v_tmp1 = join("\t", split(/ /, $ss{$a}{$b})); push( @complex_tmp1, "$a\t$b\t$v_tmp1" ) } keys %{$ss{$a}} } keys %ss;
    die "Could not open pairwise file !\n" unless( $fh->open( $input_pairwise_file ) );
    %ss = (); my %tt = ();  
    map {
        chomp; my @t = split /\t/, $_; 
        if ( $t[7] ne "N\/A" ) { $ss{$t[0]}{$t[4]} = $t[7] }; 
        if ( $t[8] ne "N\/A" ) { $tt{$t[0]}{$t[4]} = $t[8] };
        if ( $t[16] ne "N\/A" ) { $ss{$t[9]}{$t[13]} = $t[16] }; 
        if ( $t[17] ne "N\/A" ) { $tt{$t[9]}{$t[13]} = $t[17] }; 
    } <$fh>;
    $fh->close;

    map {
        my $domain0 = "N\/A"; my $cosmic0 = "N\/A"; my $domain1 = "N\/A"; my $cosmic1 = "N\/A"; 
        my @g = split /\t/, $_;  if ( defined $ss{$g[0]}{$g[1]} ) { $domain0 = $ss{$g[0]}{$g[1]}; };
        if ( defined $ss{$g[2]}{$g[3]} ) { $domain1 = $ss{$g[2]}{$g[3]}; };
        if ( defined $tt{$g[0]}{$g[1]} ) { $cosmic0 = $tt{$g[0]}{$g[1]}; };
        if ( defined $tt{$g[2]}{$g[3]} ) { $cosmic1 = $tt{$g[2]}{$g[3]}; };  
        my $v_tmp2 = join("\t", @g[0,1], $domain0, $cosmic0, @g[2,3], $domain1, $cosmic1, @g[4,5,6,7]);
        push( @complex_tmp2, $v_tmp2 ); 
    } @complex_tmp1;
    %ss = (); 
    my $fh1   = new FileHandle;
    die "Could not open maf file !\n" unless( $fh1->open( $this->{'maf'} ) );
	my $mafi = 0;
	my $mafhead = $fh1->getline(); chomp( $mafhead );
	my %mafcols = map {($_, $mafi++)} split( /\t/, $mafhead );
	unless (    defined($mafcols{"Hugo_Symbol"})
			and defined($mafcols{"Tumor_Sample_Barcode"})
			and defined($mafcols{$this->{'amino_acid_header'}}) ) {
		die "HotSpot3D Post Error: not a valid MAF annotation file with amino acid change!\n";
	}
	my @mafcols = (	$mafcols{"Hugo_Symbol"},
					$mafcols{"Tumor_Sample_Barcode"},
					$mafcols{$this->{'amino_acid_header'}} );
    map {
		chomp;
		my @t = split /\t/, $_;
		my ($sample) = $t[$mafcols[1]] =~ /(TCGA-\w\w-\w\w\w\w)/;
		if ( $t[$mafcols[0]] and $t[$mafcols[2]] and $sample ) {
			$ss{$t[$mafcols[0]]}{$t[$mafcols[2]]}{$sample}++;
		}
	} $fh1->getline;
    $fh1->close;
    @complex_tmp1 = ();
    map { 
        my ( $g0samples, $g1samples, $cocur );
        $g0samples = $g1samples = $cocur = "";
        my @g = split /\t/, $_;  
        foreach my $a ( keys %{$ss{$g[0]}{$g[1]}} ) { $g0samples .= "$a,$ss{$g[0]}{$g[1]}{$a}|"; if ( defined $ss{$g[4]}{$g[5]}{$a} ) { $cocur .= "$a,"; } };
        foreach my $b ( keys %{$ss{$g[4]}{$g[5]}} ) { $g1samples .= "$b,$ss{$g[4]}{$g[5]}{$b}|"; };
        chop($g0samples); chop($g1samples); chop($cocur); 
        if ( $cocur eq "" ) { $cocur="N\/A"; };
        my $v_tmp3 = join( "\t", @g[0..3], $g0samples, @g[4..7], $g1samples, $cocur, @g[8..11] );
        push( @complex_tmp1, $v_tmp3 );
    } @complex_tmp2;
    my %thash = ();
    map {
        my @t = split /\t/, $_;
        my $pdb = join("\t", @t[12..14]);
        my $g1 = join("\t", @t[0,2,3]);
        my $g2 = join("\t", @t[5,7,8]);
        $thash{$pdb}{$g1}{$g2}{'aac1'} .= $t[1].",";
        $thash{$pdb}{$g1}{$g2}{'aac2'} .= $t[6].",";
        $thash{$pdb}{$g1}{$g2}{'sample1'} .= $t[4].";";
        $thash{$pdb}{$g1}{$g2}{'sample2'} .= $t[9].";";
        $thash{$pdb}{$g1}{$g2}{'con'} .= $t[10].";";
        $thash{$pdb}{$g1}{$g2}{'ld'} .= $t[11].",";
    } @complex_tmp1;
    # output complex collapsed results 
    my $output_collapsed_file = "$this->{'input_prefix'}.pairwise.complex.collapsed";
    die "Could not create complex collapsed output file\n" unless( $fh1->open( ">$output_collapsed_file" ) );
    foreach my $pdb ( keys %thash ) {
        foreach my $g1 ( keys %{$thash{$pdb}} ) {
            foreach my $g2 ( keys %{$thash{$pdb}{$g1}}  ) {
                my @t1 = split /\t/, $g1;
                my @t2 = split /\t/, $g2;
                chop( $thash{$pdb}{$g1}{$g2}{'aac1'} );
                chop( $thash{$pdb}{$g1}{$g2}{'aac2'} );
                chop( $thash{$pdb}{$g1}{$g2}{'sample1'} );
                chop( $thash{$pdb}{$g1}{$g2}{'sample2'} );
                chop( $thash{$pdb}{$g1}{$g2}{'con'} );
                chop( $thash{$pdb}{$g1}{$g2}{'ld'} );
                $fh1->print( join( "\t", $t1[0], $thash{$pdb}{$g1}{$g2}{'aac1'}, @t1[1..2], $thash{$pdb}{$g1}{$g2}{'sample1'}, $t2[0], $thash{$pdb}{$g1}{$g2}{'aac2'}, @t2[1..2], $thash{$pdb}{$g1}{$g2}{'sample2'}, $thash{$pdb}{$g1}{$g2}{'con'}, $thash{$pdb}{$g1}{$g2}{'ld'}, $pdb ) );
                $fh1->print( "\n" );
            }
        }
    }
    $fh1->close;
    ## above is complex post-processing 
    ## single protein post-processing is the following
    #
    $fh   = new FileHandle;
    die "Could not open pairwise file !\n" unless( $fh->open( $input_pairwise_file ) );
    %ss = ();
    map {
        chomp; my @t = (split /\t/)[0,4,9,13,18,19];
        if ($t[4] !~ /N\/A/) {
            my @tt = split /\|/, $t[5];
            my $ld = $t[4];
            my $pvalue = (split / /, $tt[0])[-1];
            if ( defined $ss{"$t[0]\t$t[1]"}{"$t[2]\t$t[3]"} ) {
                my $p = $ss{"$t[0]\t$t[1]"}{"$t[2]\t$t[3]"}; 
                my ($pld, $p3dd) = split /\t/, $p;
                ## bug in one line
                #if ($p3dd < $pvalue) { 
                #    $ss{"$t[0]\t$t[1]"}{"$t[2]\t$t[3]"} = $ld."\t".$tt[0]; 
                #} 
                ## change to the following
                my $p3dd_pvalue = (split / /, $p3dd)[-1];
                if ($p3dd_pvalue > $pvalue) {
                    $ss{"$t[0]\t$t[1]"}{"$t[2]\t$t[3]"} = $ld."\t".$tt[0]; 
                }
            }elsif ( defined $ss{"$t[2]\t$t[3]"}{"$t[0]\t$t[1]"} ) { 
                my $p = $ss{"$t[2]\t$t[3]"}{"$t[0]\t$t[1]"}; 
                my ($pld, $p3dd) = split /\t/, $p; 
                ## bug in one line
                #if ($p3dd < $pvalue) { 
                #    $ss{"$t[2]\t$t[3]"}{"$t[0]\t$t[1]"} = $ld."\t".$tt[0]; 
                #}
                ## change to the following
                my $p3dd_pvalue = (split / /, $p3dd)[-1];
                if ($p3dd_pvalue > $pvalue) { 
                    $ss{"$t[2]\t$t[3]"}{"$t[0]\t$t[1]"} = $ld."\t".$tt[0]; 
                }
            }else { 
                $ss{"$t[0]\t$t[1]"}{"$t[2]\t$t[3]"} = $ld."\t".$tt[0]; 
            }
        }
    } <$fh>;
    $fh->close;
    @complex_tmp1 = (); 
    @complex_tmp2 = (); 
    map { my $a = $_; map { my $b = $_; my $v_tmp1 = join("\t", split(/ /, $ss{$a}{$b})); push( @complex_tmp1, "$a\t$b\t$v_tmp1" ) } keys %{$ss{$a}} } keys %ss;
    die "Could not open pairwise file !\n" unless( $fh->open( $input_pairwise_file ) );
    %ss = (); %tt = ();  
    map {
        chomp; my @t = split /\t/, $_; 
        if ( $t[7] ne "N\/A" ) { $ss{$t[0]}{$t[4]} = $t[7] }; 
        if ( $t[8] ne "N\/A" ) { $tt{$t[0]}{$t[4]} = $t[8] };
        if ( $t[16] ne "N\/A" ) { $ss{$t[9]}{$t[13]} = $t[16] }; 
        if ( $t[17] ne "N\/A" ) { $tt{$t[9]}{$t[13]} = $t[17] }; 
    } <$fh>;
    $fh->close;
    map {
        my $domain0 = "N\/A"; my $cosmic0 = "N\/A"; my $domain1 = "N\/A"; my $cosmic1 = "N\/A"; 
        my @g = split /\t/, $_;  if ( defined $ss{$g[0]}{$g[1]} ) { $domain0 = $ss{$g[0]}{$g[1]}; };
        if ( defined $ss{$g[2]}{$g[3]} ) { $domain1 = $ss{$g[2]}{$g[3]}; };
        if ( defined $tt{$g[0]}{$g[1]} ) { $cosmic0 = $tt{$g[0]}{$g[1]}; };
        if ( defined $tt{$g[2]}{$g[3]} ) { $cosmic1 = $tt{$g[2]}{$g[3]}; };  
        my $v_tmp2 = join("\t", @g[0,1], $domain0, $cosmic0, @g[2,3], $domain1, $cosmic1, @g[4,5,6,7]);
        push( @complex_tmp2, $v_tmp2 ); 
    } @complex_tmp1;
    %ss = (); 
    $fh1 = new FileHandle;
    die "Could not open maf file !\n" unless( $fh1->open( $this->{'maf'} ) );
    map {
		chomp;
		my @t = split /\t/, $_;
		my ($sample) = $t[$mafcols[1]] =~ /(TCGA-\w\w-\w\w\w\w)/;
		if ( $t[$mafcols[0]] and $t[$mafcols[2]] and $sample ) {
			$ss{$t[$mafcols[0]]}{$t[$mafcols[2]]}{$sample}++;
		}
	} <$fh1>;
    $fh1->close;
    @complex_tmp1 = ();
    map { 
        my ( $g0samples, $g1samples, $cocur );
        $g0samples = $g1samples = $cocur = "";
        my @g = split /\t/, $_;  
        foreach my $a ( keys %{$ss{$g[0]}{$g[1]}} ) { $g0samples .= "$a,$ss{$g[0]}{$g[1]}{$a}|"; if ( defined $ss{$g[4]}{$g[5]}{$a} ) { $cocur .= "$a,"; } };
        foreach my $b ( keys %{$ss{$g[4]}{$g[5]}} ) { $g1samples .= "$b,$ss{$g[4]}{$g[5]}{$b}|"; };
        chop($g0samples); chop($g1samples); chop($cocur); 
        if ( $cocur eq "" ) { $cocur="N\/A"; };
        my $v_tmp3 = join( "\t", @g[0..3], $g0samples, @g[4..7], $g1samples, $cocur, @g[8..11] );
        push( @complex_tmp1, $v_tmp3 );
    } @complex_tmp2;
    %thash = ();
    map {
        my @t = split /\t/, $_;
        my $pdb = join("\t", @t[12..14]);
        my $g1 = join("\t", @t[0,2,3]);
        my $g2 = join("\t", @t[5,7,8]);
        $thash{$pdb}{$g1}{$g2}{'aac1'} .= $t[1].",";
        $thash{$pdb}{$g1}{$g2}{'aac2'} .= $t[6].",";
        $thash{$pdb}{$g1}{$g2}{'sample1'} .= $t[4].";";
        $thash{$pdb}{$g1}{$g2}{'sample2'} .= $t[9].";";
        $thash{$pdb}{$g1}{$g2}{'con'} .= $t[10].";";
        $thash{$pdb}{$g1}{$g2}{'ld'} .= $t[11].",";
    } @complex_tmp1;
    # output complex collapsed results 
    $output_collapsed_file = "$this->{'input_prefix'}.pairwise.singleprotein.collapsed";
    die "Could not create singe protein collapsed output file\n" unless( $fh1->open( ">$output_collapsed_file" ) );
    foreach my $pdb ( keys %thash ) {
        foreach my $g1 ( keys %{$thash{$pdb}} ) {
            foreach my $g2 ( keys %{$thash{$pdb}{$g1}}  ) {
                my @t1 = split /\t/, $g1;
                my @t2 = split /\t/, $g2;
                chop( $thash{$pdb}{$g1}{$g2}{'aac1'} );
                chop( $thash{$pdb}{$g1}{$g2}{'aac2'} );
                chop( $thash{$pdb}{$g1}{$g2}{'sample1'} );
                chop( $thash{$pdb}{$g1}{$g2}{'sample2'} );
                chop( $thash{$pdb}{$g1}{$g2}{'con'} );
                chop( $thash{$pdb}{$g1}{$g2}{'ld'} );
                $fh1->print( join( "\t", $t1[0], $thash{$pdb}{$g1}{$g2}{'aac1'}, @t1[1..2], $thash{$pdb}{$g1}{$g2}{'sample1'}, $t2[0], $thash{$pdb}{$g1}{$g2}{'aac2'}, @t2[1..2], $thash{$pdb}{$g1}{$g2}{'sample2'}, $thash{$pdb}{$g1}{$g2}{'con'}, $thash{$pdb}{$g1}{$g2}{'ld'}, $pdb ) );
                $fh1->print( "\n" );
            }
        }
    }
    $fh1->close;

    return 1;
}

sub help_text{
    my $this = shift;
        return <<HELP

Usage: hotspot3d post [options]

						REQUIRED
--maf-file				Input MAF file (2.3 standard + columns for transcript id & amino acid change)
--input-prefix			The prefix of proximity searching output files

						OPTIONAL
--transcript-id-header	MAF file column header for transcript id's, default: transcript_name
--amino-acid-header		MAF file column header for amino acid changes, default: amino_acid_change

--help					this message

HELP

}

1;

