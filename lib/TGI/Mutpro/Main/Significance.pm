package TGI::Mutpro::Main::Significance;
#
#----------------------------------
# $Authors: Adam D Scott and Carolyn Lou
# $Date: 2015-09-22
# $Revision:$2015-09-22
# $URL: $
# $Doc: $ cluster significance
#----------------------------------
#
use strict;
use warnings;

use Carp;
use Getopt::Long;
use IO::File;
use FileHandle;

use Scalar::Util qw( reftype );

sub new {
	my $class = shift;
	my $this = {};
	$this->{'prep-dir'} = undef;
	$this->{'pairwise'} = undef;
	$this->{'clusters'} = undef;
	$this->{'output'} = "hotspot3d.sigclus";
	$this->{'simulations'} = 1000000;
	bless $this, $class;
	$this->process();
	return $this;
}

sub process {
	my $this = shift;
	my ( $help, $options );
	unless( @ARGV ) { die $this->help_text(); }
	$options = GetOptions (
		'prep-dir=s' => \$this->{'prep_dir'},
		'pairwise=s' => \$this->{'pairwise'},
		'clusters=s' => \$this->{'clusters'},
		'output=s' => \$this->{'output_prefix'},
		'simulations=i' => \$this->{'simulations'},
		'help' => \$help,
	);
	#my $NSIMS = $this ->{'simulations'};

	if ( $help ) { print STDERR help_text(); exit 0; }
	unless( $options ) { die $this->help_text(); }

	my $hup_file = $this->{'prep_dir'}."/hugo.uniprot.pdb.transcript.csv";
	my $proximity_dir = $this->{'prep_dir'}."/proximityFiles/cosmicanno/";

#if hup file missing
	if(not defined $hup_file ) {
		warn 'You must provide a hup file !', "\n"; 
		die $this->help_text();
		if(not -e $hup_file){
			warn "The input hup file(".$hup_file.") does not exist! ", "\n";
			die $this->help_text();
		}
	}

#if proximity directory missing
	unless( $proximity_dir ) {
		warn 'You must provide a proximity file directory ! ', "\n";
		die help_text();
	}
	unless( -d $proximity_dir ) {
		warn 'You must provide a valid proximity file directory ! ', "\n";
		die help_text();
	}

#processing procedure
#retrieve cluster information
	my $cluster_f= $this->{'clusters'}; #debug
	my $clusters = $this->getClusters($cluster_f);

# parse hup file
	my ($genesOnPDB_hash_ref, $uniprot2HUGO_hash_ref, $HUGO2uniprot_hash_ref) = $this->processHUP( $hup_file ,$clusters);

#parse aaabbrv file
	my $AA_hash_ref = $this->processAA( );
#process everything
	$this ->getDistances($proximity_dir, $genesOnPDB_hash_ref, $HUGO2uniprot_hash_ref, $AA_hash_ref, $this->{'simulations'}, $clusters)
} 

#####
#	Subs
#####
# process HUP file
sub processHUP{
	my($this, $hup_f, $cluster) = @_;
	my $hupHandle = &getFile( $hup_f , "r" );
	unless ($hupHandle->open($hup_f)){die "Could not open hup file ".$hup_f."!\n"};
#	print "$cluster\n";#debug	
	my ( %genesOnPDB, %uniprot2HUGO, %HUGO2uniprot);
	while (my $line = $hupHandle->getline){
		chomp($line);
		my ( $hugo , $uniprot , $pdbs ) = ( split( "\t" , $line ) )[0,1,2];
	#	print "hup file pdbs: $pdbs\n"; #debug
		if( exists $cluster->{$hugo}){ #retrieve hup information only for specified clusters
			foreach my $pdb ( split(/\s/ , $pdbs)  ) {
				$genesOnPDB{$hugo}{$pdb} = 1;
				#print "hugos: $hugo\n";#debug
			}
			$uniprot2HUGO{$uniprot}{$hugo} = 1;
			$HUGO2uniprot{$hugo} = $uniprot;
			print "HUGO: $hugo\n";
		}
	}
	$hupHandle->close();
	return (\%genesOnPDB, \%uniprot2HUGO, \%HUGO2uniprot);
}

#read in amino acid information
sub processAA {
	my( $this ) = @_;
	my %AA;
	my ( $single , $triple );
	my @lines = <DATA>;
	foreach my $line ( @lines ) { #organized with single letter, three-letter, and full name of amino acids
		chomp( $line );
		( $single , $triple ) = ( split( /\t/ , $line ) )[0,1];
		$AA{$triple} = $single;
	}
	return \%AA;
}

sub getClusters {
	my ($this, $cluster_f) = @_;
	my $cluster;
	my %clusters;
	my $clustersHandle = &getFile( $cluster_f, "r" );
	while ( my $line = <$clustersHandle> ) {
		chomp( $line );
		my ( $clusterid , $gene , $mutation ) = (split( /\t/ , $line ))[0..2];
		$clusters{$gene}{$mutation} = $clusterid;
	}
	$clustersHandle->close();
	return \%clusters;
}

sub getPairwise {
#	print "getting pairwise now\n"; #debug
	my ( $this , $clusters , $pdb , $genesOnPDBref ) = @_;
	my $pairwiseHandle = &getFile( $this->{'pairwise'},"r" );
	my ( %mass2pvalues , %mapping , %withinClusterSumDistance );
	while ( my $line = <$pairwiseHandle> ) {
		chomp( $line );
		my @line = split( /\t/ , $line );
		my ( $gene1 , $mu1 , $chain1 , $res1 , $gene2 , $mu2 , $chain2 , $res2 , $info ) = @line[0,4,5,6,9,13,14,15,19];
		if ( exists $clusters->{$gene1}->{$mu1}
		&& exists $clusters->{$gene2}->{$mu2}
		&& $clusters->{$gene1}->{$mu1} == $clusters->{$gene2}->{$mu2} ) {
		#	print "pairwise info exists in clusters\n";#debug
			my @infos = split( /\|/ , $info );
			foreach my $dinfo ( @infos ) {
				my ( $dist , $pdb , $pval ) = split( /\s/ , $dinfo );
				if ( exists $genesOnPDBref->{$gene1}->{$pdb}#{$chain1}
				&& exists $genesOnPDBref->{$gene2}->{$pdb}#{$chain2}
#				&& $chain1 eq $chain2 ) { #specific check to mimic SpacePAC restriction
				){	
		#			print "pairwise exists in hup file\n"; #debug
					my $cluster = $clusters->{$gene1}->{$mu1};
					$mapping{$cluster}{$res1} = 1;
					$mapping{$cluster}{$res2} = 1;
					$withinClusterSumDistance{$cluster} += $dist;
					$mass2pvalues{$cluster} = $pval;
				} #if gene,pdb,chain1 exists
			} #foreach dinfo
		} #if genes/mus exists in clusters
	} #while pairwise file
	$pairwiseHandle->close();
#	print "got the pairwise info\n";
	return ( \%mapping , \%withinClusterSumDistance , \%mass2pvalues );
}

sub getDistances{
	my ($this, $proximity_d, $genesOnPDBref, $HUGO2uniprotref, $AAref, $NSIMS, $clusters) = @_;
#	print "getDistances started\n"; #debug
	my $hash_count = keys %$genesOnPDBref;#debug
#	print "number of genesonpdb keys: $hash_count\n";#debug
	foreach my $hugo(keys %$genesOnPDBref) {
#		print "hugo: $hugo\n";#debug
		my $uniprot = $HUGO2uniprotref->{$hugo};
		my %distances;
		my $proxFile = $proximity_d.$uniprot.".ProximityFile.csv";
#		print "proxfile: $proxFile\n"; #debug
		my $proxHandle = &getFile( $proxFile , "r" );
		while(my $line = <$proxHandle>) {
			chomp $line;
			my ( $up1 , $chain1 , $res1 , $aa1 , $up2 , $chain2 , $res2 , $aa2 , $dist , $pdb , $pval ) = (split( /\t/ , $line ))[0..2,4,7..9,11,14..16];
			if ( exists $genesOnPDBref->{$hugo}{$pdb}#{$chain1}
			#&& exists $genesOnPDBref->{$hugo}{$pdb}#{$chain2}
			#&& $chain1 eq $chain2 #specific check to mimic SpacePAC restriction
			&& exists $AAref->{$aa1}
			&& exists $AAref->{$aa2} ) {
				push @{$distances{$pdb}{$chain1}} , $dist;
				#print "WE GOT ONE\n";#debug
			}
		}
		$proxHandle->close();
		
		foreach my $pdb ( keys %distances ) {
			print "PDB: $pdb\n";#debug
			foreach my $chain ( keys %{$distances{$pdb}} ) {
#				print "chain: $chain\n";#debug
				my @distances = sort {$a <=> $b} @{$distances{$pdb}{$chain}};
				#my @distances = sort {$a <=> $b} keys %{$distances{$pdb}{$chain}};
				my $numDistances = scalar( @distances );

				#print STDOUT $tossedPairs{$pdb}{$chain}." pairs thrown out\n";
				
				my $fout = $this->{'output_prefix'}.".$hugo."."$pdb".".".$chain.".".$NSIMS."sims.avgDist_pvalue";
				my $OUT = &getFile( $fout , "w" );
#				print "Cluster\tGene\tPDB_ID\tChain\tNum_Residues\tNum_Pairs\tAvg_Distance\tEstimated_PValue\n";
				$OUT->print( "Cluster\tGene\tPDB_ID\tChain\tNum_Residues\tNum_Pairs\tAvg_Distance\tEstimated_PValue\n" );
				
				my ( $mapping , $withinClusterSumDistance , $mass2pvalues ) = $this->getPairwise( $clusters , $pdb , $genesOnPDBref );
	
				my %massDistribution;
				my %withinClusterAvgDistance;
				my ( $mass , $numPairs );
				foreach my $cluster ( keys %{$mapping} ) {
					$mass = scalar keys %{$mapping->{$cluster}};
					$numPairs = ( $mass ) * ( $mass - 1 ) / 2;
					if ( $mass > 2 ) {
						$massDistribution{$mass}{$cluster} = 1;
						$withinClusterAvgDistance{$cluster} = $withinClusterSumDistance->{$cluster} / $numPairs;
#						print "cluster: $cluster, avg distance: $withinClusterAvgDistance{$cluster}, sum distance: $withinClusterSumDistance->{$cluster}\n";#debug
						delete $mass2pvalues->{$cluster};
					} elsif ( $mass == 2 ) {
#print STDOUT "Cluster ".$cluster." has ".$mass." residues, so its p-value is ".$mass2pvalues->{$cluster}."\n"; #debug
						$OUT->print( join( "\t" , ( $cluster , $hugo , $pdb , $chain , $mass , $numPairs , $withinClusterSumDistance->{$cluster} , $mass2pvalues->{$cluster} ) )."\n" );
					} elsif ( $mass == 1 ) {
#						print STDOUT "Cluster ".$cluster." has ".$mass." residue, so it has no p-value\n";#debug
						$OUT->print( join( "\t" , ( $cluster , $hugo , $pdb , $chain , $mass , $numPairs , 0 , "NA" ) )."\n" );
					} #if mass block
				} #foreach cluster

				#my %permutationTestPValues;
#				print STDOUT $hugo."\t".$pdb."\t".$chain."\n";
				foreach my $numResidues ( keys %massDistribution ) {
					$numPairs = ( $numResidues ) * ( $numResidues - 1 ) / 2;
					my %below;
#					print STDOUT "Simulation\tSum_Distances\tAvg_Distance\n";
					for ( my $simulation = 0; $simulation < $NSIMS ; $simulation++ ) {
						my $sumDists = 0;
						for ( my $i = 0; $i < $numPairs ; $i++ ) {
							my $randIndex = int( rand( $numDistances ) );
							$sumDists += $distances[$randIndex];
						} #foreach random pick
						my $avgDist = $sumDists / $numPairs;
#						print STDOUT $simulation."\t".$sumDists."\t".$avgDist."\n";
						foreach my $cluster ( keys %{$massDistribution{$numResidues}} ) {
						#	print "cluster: $cluster, avg dist: $avgDist, withinetc: $withinClusterAvgDistance{$cluster}\n";
							if ( $avgDist < $withinClusterAvgDistance{$cluster} ) {
								$below{$cluster}++;
							}
						} #foreach cluster
					} #foreach simulation
					foreach my $cluster ( keys %{$massDistribution{$numResidues}} ) {
						my $permutationTestPValue;
						if(exists $below{$cluster}){
							$permutationTestPValue = $below{$cluster} / $NSIMS;

						#	print " below: $below{$cluster}\t";#debug
						#	print STDOUT "Estimated p-value: ".$permutationTestPValue."\n";
						} else{
							$permutationTestPValue=0;
							#print "why zero\n";#debug
						}
	
						$OUT->print( join( "\t" , ( $cluster , $hugo , $pdb , $chain , $numResidues , $numPairs , $withinClusterAvgDistance{$cluster} , $permutationTestPValue ) )."\n" );
					} #foreach cluster
				} #foreach numResidues
				$OUT->close();
			} #foreach chain
		} #foreach pdb
	} #foreach gene/protein from hup
	return 1;
}

sub getFile {
	my ( $file , $rw ) = @_;
#	print "$file\n"; #debug
	my $handle = FileHandle->new( $file, $rw );
	if ( not defined $handle ) {
		if ( $rw eq "w" ) {
			warn "ADSERROR: Could not open/write $file\n";
		} else {
			warn "ADSERROR: Could not open/read $file\n";
		}
		next;
	}
	return $handle;
}

sub help_text{
	my $this = shift;
	return <<HELP

Usage: hotspot3d sigclus [options]

	--prep-dir			Preprocessing directory 
	--pairwise			Pairwise file (pancan19.pairwise)
	--clusters			Cluster file (pancan19.intra.20..05.10.clusters)
	--output			Output file prefix (pancan19.intra.20..05.10)

	--simulations	Number of simulations, default = 1000000

	--help				This message

HELP

}

1;

__DATA__
A	ALA	Alanine
R	ARG	Arginine
N	ASN	Asparagine
D	ASP	Aspartic Acid
C	CYS	Cysteine
E	GLU	Glutamic Acid
Q	GLN	Glutamine
G	GLY	Glycine
H	HIS	Histidine
I	ILE	Isoleucine
L	LEU	Leucine
K	LYS	Lysine
M	MET	Methionine
F	PHE	Phenylalanine
P	PRO	Proline
S	SER	Serine
T	THR	Threonine
W	TRY	Tryptophan
Y	TYR	Tyrosine
V	VAL	Valine
*	S	Stop
