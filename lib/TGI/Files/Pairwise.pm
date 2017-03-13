package TGI::Mutpro::Files::Pairwise;
#
#----------------------------------
# $Authors: Adam Scott
# $Date: 2016-10-26
# $Revision: v0.0 $
# $URL: $
# $Doc: $ file handling for .pairwise
# 
#----------------------------------
#
use strict;
use warnings;

use Carp;

use TGI::Mutpro::Files::File;
our @ISA = qw( TGI::Mutpro::Files::File );

my $TRANSCRIPT = "transcript_name";
my $AMINOACIDCHANGE = "amino_acid_change";
my $WEIGHT = "weight";
my @REQUIRED = ( "Hugo_Symbol" , "Chromosome" , "Start_Position" ,
				"End_Position" , "Reference_Allele" ,
				"Tumor_Seq_Allele2" , "Tumor_Sample_Barcode" );

sub new {
    my $class = shift;
    my $this = $class->SUPER::new( shift );
    $this->{'samples'} = 0;
    $this->{'entries'} = 0;
    $this->{'variants'} = 0;
    $this->{'amino_acid_changes'} = 0;
	$this->{'required'} = \@REQUIRED;
	$this->{'transcript_id_header'} = $TRANSCRIPT;
	$this->{'amino_acid_header'} = $AMINOACIDCHANGE;
	$this->{'weight_header'} = $WEIGHT;
    bless $this, $class;
    return $this;
}

sub openPairwise {
    my $this = shift;
	if ( defined $this->{'file_name'} ) {
		$this->{'handle'}->open( $this->{'file_name'} , "r" );
		return 1;
	}
	return 0;
}

sub closePairwise {
	my $this = shift;
	$this->{'handle'}->close();
	return;
}

sub setTranscriptHeader {
	my $this = shift;
	if ( @_ ) {
		$this->{'transcript_id_header'} = shift;
	}
	return;
}

sub setAminoAcidChangeHeader {
	my $this = shift;
	if ( @_ ) {
		$this->{'amino_acid_header'} = shift;
	}
	return;
}

sub setWeightHeader {
	my $this = shift;
	if ( @_ ) {
		$this->{'weight_header'} = shift;
	}
	return;
}

sub requireTranscript {
	my $this = shift;
	if ( @_ ) {
		$this->setTranscriptHeader( shift );
	}
	push @{$this->{'required'}} , $this->{'transcript_id_header'};
}

sub requireAminoAcidChange {
	my $this = shift;
	if ( @_ ) {
		$this->setAminoAcidChangeHeader( shift );
	}
	push @{$this->{'required'}} , $this->{'amino_acid_header'};
}

sub requireWeight {
	my $this = shift;
	if ( @_ ) {
		$this->setWeightHeader( shift );
	}
	push @{$this->{'required'}} , $this->{'weight_header'};
}

sub setColumnIndices {
	my $this = shift;
	my $exception = "HotSpot3D::Pairwise error: required columns not present ";
	$exception .= "in this .maf (".$this->{'file_name'}."). Need the ";
	$exception .= "following columns:\n\t";
	$exception .= join( "\n\t" , @{$this->{'required'}} )."\n";
	my $header = $this->{'handle'}->getline(); chomp( $header );
	my $mafcols = $this->mapColumns( $header , $this->{'required'} , $exception );
	return $mafcols;
}

sub getlines {
	my $this = shift;
	return $this->{'handle'}->getlines;
}

sub readPairwise {
	my $this = shift;
	print STDOUT "\nReading in ".$this->{'file_name'}."...\n";
	seek( $this->{'handle'} , 0 , 0 );
	my @mafcols = @{$this->setColumnIndices()};
	my %mutations;
	map {
		chomp;
		my @line = split /\t/;
		my $variant = "";
		if ( $#line >= $mafcols[-1] && $#line >= $mafcols[-2] ) { #makes sure custom maf cols are in range
			my ( $gene , $chr , $start , $stop , $reference , $tumorAllele , $barID );
			my ( $transcript_name , $aachange , $weight );
			$weight = 1;
			if ( grep{ $_ eq $this->{'weight_header'} } @{$this->{'required'}} ) {
				( $gene , $chr , $start , $stop , $reference , $tumorAllele , $barID , $transcript_name , $aachange , $weight ) = @line[@mafcols];
			} elsif ( grep{ $_ eq $this->{'amino_acid_header'} } @{$this->{'required'}} ) {
				( $gene , $chr , $start , $stop , $reference , $tumorAllele , $barID , $transcript_name , $aachange ) = @line[@mafcols];
				$variant = join( "_" , ( $gene , $aachange , $chr , $start , $stop ) );
			} else {
				( $gene , $chr , $start , $stop , $reference , $tumorAllele , $barID ) = @line[@mafcols];
				$variant = join( "_" , ( $gene , $chr , $start , $stop ) );
			}
		}
	} $this->getlines();

	return;
}

sub getGenes {
	my $this = shift;
	print STDOUT "\nReading in ".$this->{'file_name'}." to get genes...\n";
	seek( $this->{'handle'} , 0 , 0 );
	my @mafcols = @{$this->setColumnIndices()};
	my %genes;
	map {
		chomp;
		my $gene = (split /\t/)[$mafcols[0]];
		$genes{$gene} += 1;
	} $this->getlines();

	return \%genes;
}

sub getSamples {
	my $this = shift;
	print STDOUT "\nReading in ".$this->{'file_name'}." to get genes...\n";
	seek( $this->{'handle'} , 0 , 0 );
	my @mafcols = @{$this->setColumnIndices()};
	my %samples;
	map {
		chomp;
		my $sample = (split /\t/)[$mafcols[6]];
		$samples{$sample} += 1;
	} $this->getlines();

	return \%samples;
}

1;

__DATA__
0	SMAD2
1	18
2	45368254
3	45368254
4	p.D450N
5	[A]
6	450
7	MH2. {ECO:0000255|PROSITE
8	p.D420N|lung,p.D450N|lung
9	SMAD4
10	18
11	48591918
12	48591918
13	p.R361C
14	[B]
15	361
16	MH2. {ECO:0000255|PROSITE
17	N/A
18	N/A
19	2.561 1U7V 0.004409|
