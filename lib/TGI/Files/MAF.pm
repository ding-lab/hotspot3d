package TGI::Files::MAF;
#
#----------------------------------
# $Authors: Adam Scott
# $Date: 2016-10-26
# $Revision: v0.0 $
# $URL: $
# $Doc: $ file handling for .mafs
# 
#----------------------------------
#
use strict;
use warnings;

use Carp;
use Scalar::Util qw( openhandle );
use List::MoreUtils qw( first_index );

use TGI::Files::File;
our @ISA = qw( TGI::Files::File );
use TGI::Files::List;
use TGI::Data::StringTemplate;

my $HUGOSYMBOL = "Hugo_Symbol";
#my $ENTREZGENEID = "Entrez_Gene_Id";
#my $CENTER = "Center";
#my $NCBIBUILD = "NCBI_Build";
my $CHROMOSOME = "Chromosome";
my $STARTPOSITION = "Start_Position";
my $ENDPOSITION = "End_Position";
#my $STRAND = "Strand";
my $REFERENCEALLELE = "Reference_Allele";
my $TUMORSEQALLELE2 = "Tumor_Seq_Allele2";
my $TUMORSAMPLEBARCODE = "Tumor_Sample_Barcode";

sub new {
    my $class = shift;
    my $this = $class->SUPER::new( shift );
    $this->{'samples'} = 0;
    $this->{'entries'} = 0;
    $this->{'variants'} = 0;
    bless $this, $class;
    return $this;
}

sub getHugo {
	my $this = shift;
	my $line = shift;
	return $this->getField( $HUGOSYMBOL , $line );
}

sub getChromosome {
	my $this = shift;
	my $line = shift;
	return $this->getField( $CHROMOSOME , $line );
}
sub getStart {
	my $this = shift;
	my $line = shift;
	return $this->getField( $STARTPOSITION , $line );
}
sub getStop {
	my $this = shift;
	my $line = shift;
	return $this->getField( $ENDPOSITION , $line );
}
sub getReference {
	my $this = shift;
	my $line = shift;
	return $this->getField( $REFERENCEALLELE , $line );
}
sub getVariant {
	my $this = shift;
	my $line = shift;
	return $this->getField( $TUMORSEQALLELE2 , $line );
}
sub getTumorSample {
	my $this = shift;
	my $line = shift;
	return $this->getField( $TUMORSAMPLEBARCODE , $line );
}

sub getGenes {
	my $this = shift;
	print STDOUT "\nReading in ".$this->{'file_name'}." to get genes...\n";
	my $closed = 0;
	if ( defined openhandle( $this->{'handle'} ) ) {
		$this->close();
		$closed = 1;
	}
	my $file = new TGI::Files::List( $this->{'file_name'} );
	$file->close();
	if ( $closed ) {
		$file->open();
	}
	my $genes = $file->getList( $this->getColumnIndex( $HUGOSYMBOL ) );

	return $genes;
}

sub getSamples {
	my $this = shift;
	print STDOUT "\nReading in ".$this->{'file_name'}." to get samples...\n";
	my $closed = 0;
	if ( defined openhandle( $this->{'handle'} ) ) {
		$this->close();
		$closed = 1;
	}
	my $file = new TGI::Files::List( $this->{'file_name'} );
	my $samples = $file->getList( $this->getColumnIndex( $TUMORSAMPLEBARCODE ) );
	$file->close();
	if ( $closed ) {
		$file->open();
	}

	return $samples;
}

1;
