package TGI::Mutpro::Main::Significance;
#
#----------------------------------
# $Authors: Adam D Scott, Carolyn Lou, and Sohini Sengupta
# $Date: 2015-09-22
# $Revision:$2017-03-12
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
	$this->{'musites-file'} = undef;
	$this->{'sites-file'} = undef;
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
		'musites-file=s' => \$this->{'musites_file'},
		'sites-file=s' => \$this->{'sites_file'},
		'clusters=s' => \$this->{'clusters'},
		'output=s' => \$this->{'output_prefix'},
		'simulations=i' => \$this->{'simulations'},
		'help' => \$help,
	);
	#my $NSIMS = $this ->{'simulations'};

	if ( $help ) { print STDERR help_text(); exit 0; }
	unless( $options ) { die $this->help_text(); }

	my $hup_file = $this->{'prep_dir'}."/hugo.uniprot.pdb.transcript.csv";
#	my $hup_file = $this->{'prep_dir'}."/test_hugo";
	my $proximity_dir = $this->{'prep_dir'}."/proximityFiles/cosmicanno/";
	my $cluster_f= $this->{'clusters'}; #debug

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


#	print"working1\n";
#processing procedure

#parse aaabbrv file
	my $AA_hash_ref = $this->processAA( );

#retrieve cluster information
	my ($clusters, $clus_store, $recurrence) = $this->getClusters($cluster_f);

# parse hup file
	my ($genesOnPDB_hash_ref, $uniprot2HUGO_hash_ref, $HUGO2uniprot_hash_ref) = $this->processHUP( $hup_file ,$clusters);

#retrieve locations of mutations in clusters using pairfiles

#	print"working2\n";
        my (%res2Mut, %structures,%minMax);

	my ($withinClusterSumDistance,$clusterPvalue);
	#does mu-sites or sites-sites file exist
	if ( defined $this->{'sites_file'}) {
			if ( not -e $this->{'sites_file'} ) { 
				warn "The input site-site pair file (".$this->{'sites_file'}.") does not exist! ", "\n";
				die $this->help_text();
			}
			else{

				($withinClusterSumDistance,$clusterPvalue)=$this->readPairFile('sites_file',\%res2Mut,\%structures,\%minMax, $clusters,$recurrence);
				print"working\n";
			}
		}
		# else {
		#	warn "HotSpot3D::Cluster::setOptions warning: no sites-file included (cannot produce site-site clusters)!\n";

	if ( defined $this->{'musites_file'} ) {
			if ( not -e $this->{'musites_file'} ) { 
				warn "The input mutation-site pair file (".$this->{'musites_file'}.") does not exist! ", "\n";
				die $this->help_text();
			}
			else {
				($withinClusterSumDistance,$clusterPvalue)=$this->readPairFile('musites_file',\%res2Mut,\%structures,\%minMax, $clusters,$recurrence);
			}
		} 
	#	else {
	#		warn "HotSpot3D::Cluster::setOptions warning: no musites-file included (cannot produce mutation-site clusters)!\n";
	#	}


	if ( defined $this->{'pairwise'} ) {
			if ( not -e $this->{'pairwise'} ) { 
				warn "The input pairwise file (".$this->{'pairwise'}.") does not exist! ", "\n";
				die $this->help_text();
			}
			else {

				($withinClusterSumDistance,$clusterPvalue)=$this->readPairFile('pairwise',\%res2Mut, \%structures,\%minMax, $clusters, $recurrence);
			}
		} 
	#	else {
	#		warn "HotSpot3D::Cluster::setOptions warning: no musites-file included (cannot produce mutation-site clusters)!\n";
	#	}


#	my ($mass2pvalues,$res2Mut)= $this->getPairwise( $clusters,$bestStruc);



#retrieve best structure
	my $bestStruc=$this->getBestStructure($clus_store, \%minMax, $clusters, $recurrence, \%structures);

#process everything

	$this ->getDistances($proximity_dir, $HUGO2uniprot_hash_ref, $uniprot2HUGO_hash_ref, $AA_hash_ref, $this->{'simulations'}, $clusters,$bestStruc, $clus_store, $recurrence,$withinClusterSumDistance,$clusterPvalue)
} 

#####
#	Subs
#####

sub readPairFile{
	my ($this, $pair_file, $res2Mut, $structures,$minMax, $clusters, $recurrence)= @_;
	my $PairHandle = &getFile( $this->{$pair_file},"r" );
	print"readPair working\n";
	my ($withinClusterSumDistance,$clusterPvalue);
	while ( my $line = <$PairHandle> ) {
		chomp( $line );
		my @line = split( /\t/ , $line );
		my ( $gene1 , $mu1 , $chain1 , $res1 , $gene2 , $mu2 , $chain2 , $res2 , $info ); 
		if($pair_file eq 'musites_file'){
             		( $gene1 , $mu1 , $chain1 , $res1 , $gene2 , $mu2 , $chain2 , $res2 , $info ) = @line[0,4,5,6,9,11,13,14,18];
		}
		elsif($pair_file eq 'pairwise'){
              		( $gene1 , $mu1 , $chain1 , $res1 , $gene2 , $mu2 , $chain2 , $res2 , $info ) = @line[0,4,5,6,9,13,14,15,19];
		}
		elsif($pair_file eq 'sites_file'){
         		( $gene1 , $mu1 , $chain1 , $res1 , $gene2 , $mu2 , $chain2 , $res2 , $info ) = @line[0,2,4,5,8,10,12,13,17];
		}

		($withinClusterSumDistance,$clusterPvalue)=$this->getPairwise($gene1 , $mu1 , $chain1 , $res1 , $gene2 , $mu2 , $chain2 , $res2 , $info, $clusters, $recurrence);

		$chain1 =~ s/\[(.*)\]/$1/;
		$chain2 =~ s/\[(.*)\]/$1/;
		my @pdbInfos = split( '\|' , $info );
		foreach my $pdbInfo ( @pdbInfos ) {
			my $pdb = (split( ' ' , $pdbInfo ))[1];
			$structures->{$gene1}->{$mu1}->{$pdb}->{$chain1} = $res1;
			$structures->{$gene2}->{$mu2}->{$pdb}->{$chain2} = $res2;
			&checkMin( $minMax , $pdb , $gene1 , $chain1 , $res1 );
			&checkMin( $minMax , $pdb , $gene2 , $chain2 , $res2 );
			&checkMax( $minMax , $pdb , $gene1 , $chain1 , $res1 );
			&checkMax( $minMax , $pdb , $gene2 , $chain2 , $res2 );

		}
	}

        $PairHandle->close();
	return ($withinClusterSumDistance,$clusterPvalue);
}

sub checkMin {
	my ( $minMax , $pdb , $gene , $chain , $residue ) = @_;
	return unless ( $residue =~ /^\d+$/ );
	if ( exists $minMax->{$pdb}->{$gene}->{$chain}->{'min'} ) {
		if ( $minMax->{$pdb}->{$gene}->{$chain}->{'min'} <= $residue ) {
			return;
		}
	}
	$minMax->{$pdb}->{$gene}->{$chain}->{'min'} = $residue;
}

sub checkMax {
	my ( $minMax , $pdb , $gene , $chain , $residue ) = @_;
	return unless ( $residue =~ /^\d+$/ );
	if ( exists $minMax->{$pdb}->{$gene}->{$chain}->{'max'} ) {
		if ( $minMax->{$pdb}->{$gene}->{$chain}->{'max'} >= $residue ) {
			return;
		}
	}
	$minMax->{$pdb}->{$gene}->{$chain}->{'max'} = $residue;
}

sub getPairwise {
       # print "getting pairwise now\n"; #debug
        my ( $this , $gene1 , $mu1 , $chain1 , $res1 , $gene2 , $mu2 , $chain2 , $res2 , $info, $clusters, $bestStruc, $recurrence) = @_;
	my (%withinClusterSumDistance,%clusterPvalue);
 	if ( exists $clusters->{$gene1}->{$mu1} && exists $clusters->{$gene2}->{$mu2} && $clusters->{$gene1}->{$mu1} eq $clusters->{$gene2}->{$mu2} ) 
	{
		my @infos = split( /\|/ , $info );
		foreach my $dinfo ( @infos ) {
			my ( $dist , $pdb , $pval ) = split( /\s/ , $dinfo );
			my $cluster = $clusters->{$gene1}->{$mu1}; #store clusterid of mutation
#                       print "pairwise exists in hup file\n"; #debug
#			print"if statement in pairwise working\n";
	#		&storeRes2Mut($res2Mut,$pdb,$res1,$mu1, $chain1);
	#		&storeRes2Mut($res2Mut,$pdb,$res2,$mu2, $chain2);

	#		my $recsum1=&sumWeight($res2Mut,$pdb,$res1,$recurrence,$clus,$chain1); #get sum of weights for each mutation at a residue
	#		my $recsum2=&sumWeight($res2Mut,$pdb,$res2,$recurrence,$clus,$chain2);

			my $rec1=$recurrence->{$cluster}->{$mu1};
			my $rec2=$recurrence->{$cluster}->{$mu2};

			$withinClusterSumDistance{$cluster}{$pdb}+= $rec1*$rec2*$dist;
			$clusterPvalue{$cluster}{$pdb}=$pval;
			print "$mu1:$chain1\t$mu2:$chain2\t$pdb\n";
			print"Distance between pairs:$dist\n";
		}#foreach dinfo
	}			
#       print "got the pairwise info\n";
        return  (\%withinClusterSumDistance,\%clusterPvalue);
}


#function to add array of mutations at a specific residue site
sub storeRes2Mut{
        my ($res2Mut,$pdb,$res,$mut, $chain)= @_;
        if( exists $res2Mut->{$pdb}->{$res}->{$chain}){
                if(!grep $_ eq $mut, @{$res2Mut->{$pdb}->{$res}->{$chain}}){
                        push @{$res2Mut->{$pdb}->{$res}->{$chain}}, $mut; 
			print"$pdb\t$res\t$chain\t$mut";
                }
        }
        else{
                $res2Mut->{$pdb}->{$res}->{$chain}=();
                push @{$res2Mut->{$pdb}->{$res}->{$chain}}, $mut; 
        }

        return 1;
}


# process HUP file
sub processHUP{
	print "processing HUP\n";
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
				print "hugos: $hugo\n";#debug
			}
			$uniprot2HUGO{$uniprot}= $hugo;
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



sub getClusters { #get cluster id for gene and mutation
	print "getting clusters now\n";
	my ($this, $cluster_f,$structures) = @_;

	my $cluster;
	my %clusters;
	my %clus_store;
	my %recurrence;

	my $clustersHandle = &getFile( $cluster_f, "r" );
	while ( my $line = <$clustersHandle> ) {
		chomp( $line );
		if ($line!~/Cluster/){
			my ( $clusterid , $gene , $mutation, $weight, $resPos ) = (split( /\t/ , $line ))[0..2,6,8];
			my @genes;
			if (exists $clus_store{$clusterid}){
				@genes=@{$clus_store{$clusterid}};
				if(! grep $_ eq $gene, @genes){
					push(@genes, $gene);
				}
			
			}
			else{
				@genes=($gene);
			}
			if ($mutation=~/p.\S\d+$/){
				$mutation=$resPos;
			}

			$clusters{$gene}{$mutation} = $clusterid;
			$recurrence{$clusterid}{$mutation}=int($weight);
			$clus_store{$clusterid}=\@genes;

		}
	}
	$clustersHandle->close();
	return \%clusters, \%clus_store,\%recurrence;
}



sub getBestStructure{
	my ( $this,$clus_store,$minMax, $clusters,$recurrence, $structures) = @_;
	my %complex;
	my %bestStruc;
	my %represent;
	print"getBestStruc working\n";

	foreach my $gene (keys %{$clusters}){
		foreach my $mutation(keys %{$clusters->{$gene}}){
			my $clusterid=$clusters->{$gene}->{$mutation};
			my $weight=$recurrence->{$clusterid}->{$mutation};
			foreach my $pdb ( keys %{$structures->{$gene}->{$mutation}} ) {
				foreach my $chain ( keys %{$structures->{$gene}->{$mutation}->{$pdb}} ) {
					$represent{$pdb}{$gene}{$chain}{$clusterid}{$mutation.":".$structures->{$gene}->{$mutation}->{$pdb}->{$chain}} = $weight;
				}
		
			}
		}
	}

	foreach my $pdb ( sort keys %represent ) {
		foreach my $gene ( sort keys %{$represent{$pdb}}) {
			foreach my $chain ( sort keys %{$represent{$pdb}{$gene}} ) {
				$chain =~ m/\[(.*)\]/;
				foreach my $cluster ( sort keys %{$represent{$pdb}{$gene}{$chain}} ) {
					my ( $mutres , $mutation , $position );
					my @mutations;
					my %residues;
					my $recur = 0;
					foreach my $mutpos ( sort keys %{$represent{$pdb}{$gene}{$chain}{$cluster}} ) {
						( $mutation , $position ) = split( ":" , $mutpos );
						$recur += $represent{$pdb}{$gene}{$chain}{$cluster}{$mutpos};
						$mutres = join( "|" , ( $mutation , $position ) );
						push @mutations , $mutres;
						$residues{$position} = 1;
					}
					if ( exists $minMax->{$pdb}->{$gene}->{$chain} ) {
						my $min = $minMax->{$pdb}->{$gene}->{$chain}->{'min'};
						my $max = $minMax->{$pdb}->{$gene}->{$chain}->{'max'};
						my @outline = ( $pdb , $gene , $chain , $cluster , $min , $max , scalar @mutations , scalar keys %residues , $recur , join( ";" , @mutations ) );
						$complex{$cluster}{$pdb}{$gene}{$chain} = \@outline;
					}
				}
			}
		} #foreach pdb in represent => cluster
	} #foreach cluster in represent

	foreach my $cluster ( sort keys %complex ) {
		print"clusterID:$cluster\n";
		foreach my $pdb ( sort keys %{$complex{$cluster}} ) {
			print"PDB:$pdb\n";
			my ( @mutations , @geneChains );
			my ( $mutations , $residues , $recur );
			foreach my $gene ( sort keys %{$complex{$cluster}{$pdb}} ) {
				my @chains;
				foreach my $chain ( sort keys %{$complex{$cluster}{$pdb}{$gene}} ) {
					my @outline = @{$complex{$cluster}{$pdb}{$gene}{$chain}};
					$mutations += $outline[6];
					$residues += $outline[7];
					$recur += $outline[8];
					push @chains , $chain;
					push @mutations , $chains[-1]."\\";
					$mutations[-1] .= $outline[-1];
				} #foreach chain
				push @geneChains , $gene."|".join( "/" , @chains );
			} #foreach gene
			my @complexLine = ( $cluster , $pdb , join( ";" , @geneChains ) , $mutations , $residues , $recur , join( ";" , @mutations ) );

			my @pair=(int($recur),$pdb);
			print"$cluster\n";
			if (exists $clus_store->{$cluster} ){
				print"cluster id exists\n";
				if (exists $bestStruc{$cluster}){
						my $test=$bestStruc{$cluster}->[0];
						if(int($recur)> $test){
							$bestStruc{$cluster}=\@pair;
							print"new pair stored\n";
							print"$cluster\t$recur\t$pdb\n";
						}
				}
				else{
					$bestStruc{$cluster}=\@pair;
				}
			}


		} #foreach pdb
	} #foreach cluster

	return \%bestStruc;
}




sub getProx{
        my ($this, $proximity_d, $uniprotIDs, $AAref, $clusters, $best_pdb, $clus, $uniprot2HUGO) = @_;
	print"getProx started\n";
        my %distances;
        my $uniprotRef=$uniprotIDs->[0]; #pick 1 uniprot ID to get distances
        my $uni_length=scalar @{$uniprotIDs};
        my $proxFile = $proximity_d.$uniprotRef.".ProximityFile.csv";
#       print "proxfile: $proxFile\n"; #debug
        my $proxHandle = &getFile( $proxFile , "r" );
        my (%paircheck, %res_count);
        while(my $line = <$proxHandle>) {
                chomp $line;
                my ( $up1 , $chain1 , $res1 , $aa1 , $up2 , $chain2 , $res2 , $aa2 , $dist , $pdb , $pval ) = (split( /\t/ , $line ))[0..2,4,7..9,11,14..16];
                my $key1="$up1:$res1:$chain1";
                my $key2="$up2:$res2:$chain2";
		my $gene1=$uniprot2HUGO->{$up1};
		my $gene2=$uniprot2HUGO->{$up2};
       #        	print ("$aa1\t$aa2\t$pdb\t$best_pdb\t$key1\t$key2\n");  
                if (exists $AAref->{$aa1} && exists $AAref->{$aa2} && $pdb eq $best_pdb && !exists $paircheck{$pdb}{$key1}{$key2} && !exists $paircheck{$pdb}{$key2}{$key1}) {
	#			print"key has not been seen in hgetProx\n";
                	if(grep $_ eq $up2,@{$uniprotIDs}){ #if 2nd uniprot id is in uniprots that are in cluster
	#			print"uniprot is in cluster in getProx\n";
                                &addZeroDist($key1,\%res_count,$pdb,\%distances, $uni_length);
                                &addZeroDist($key2,\%res_count,$pdb,\%distances, $uni_length);

		#		print"distance getting pushed\n";
                                push @{$distances{$pdb}} , $dist;
                                $paircheck{$pdb}{$key1}{$key2}=1;

                        }
                }
        }
                $proxHandle->close();
        return \%distances;

}


sub sumWeight{
	print"sumWeight function working\n";
        my ($res2Mut,$pdb,$res,$recurrence,$clus,$chain)=@_;
	my $recsum=0;
        foreach my $mut(@{$res2Mut->{$pdb}->{$res}->{$chain}}){
                $recsum+=$recurrence->{$clus}->{$mut};
		print"$mut\n";    
	}

		print "$recsum\n";
        return $recsum
}


sub addZeroDist{
        my ($key,$res_count,$pdb, $distances,$uni_length)=@_;

        if (! grep $_ eq $key,@{$res_count->{$pdb}}){

                push @{$res_count->{$pdb}}, $key;

#                if ($uni_length==1){ #only push 0 distance if it's intramolecular cluster
                push @{$distances->{$pdb}}, 0;
#                }
        }

        return 1;
}

sub getDistances{
        my ($this, $proximity_d, $HUGO2uniprotref,$uniprot2HUGO, $AAref, $NSIMS, $clusters,$bestStruc, $clus_store, $recurrence,  $withinClusterSumDistance,$clusterPvalue) = @_;
        print "getDistances started\n"; #debug

        my $fout = $this->{'output_prefix'}.".".$NSIMS."sims.avgDist_pvalue";
        my $OUT = &getFile( $fout , "w" );
        $OUT->print( "Cluster\tGene\tPDB_ID\tNum_Residues\tNum_Mutations\tNum_Pairs\tAvg_Distance\tEstimated_PValue\n" );

        foreach my $clus ( keys %{$clus_store} ) {

                my $best_pdb=$bestStruc->{$clus}->[1]; #store best PDB for cluster
                my @uniprotIDs;

                foreach my $gene(@{$clus_store->{$clus}}){
                        my $uniprot = $HUGO2uniprotref->{$gene};
                        if (! grep $_ eq $uniprot,@uniprotIDs){
                                push @uniprotIDs, $uniprot;  #store all unique uniprotIds in cluster
                                print "$uniprot\n";

                        }
                }

		if (exists $bestStruc->{$clus}->[1]){		
			#get background distances for protein and the sum of all pairwise residue distances in cluster
			my $distances=$this->getProx($proximity_d, \@uniprotIDs, $AAref, $clusters, $best_pdb, $clus, $uniprot2HUGO);

			print "printing cluster: $clus\n";
			print "PDB: $best_pdb\n";#debug
			my @distances = sort {$a <=> $b} @{$distances->{$best_pdb}};
			#my @distances = sort {$a <=> $b} keys %{$distances{$pdb}{$chain}};
			my $numDistances = scalar( @distances );

			#print STDOUT $tossedPairs{$pdb}{$chain}." pairs thrown out\n";
			my $genes=join( ":" , @{$clus_store->{$clus}});
			my ( $mass , $numPairs, $numResidues );
			$mass=0;
			$numResidues=0;
			foreach my $mut(keys %{$recurrence->{$clus}}){
				my $rec=$recurrence->{$clus}->{$mut};
				print"$rec\t$mut\n";
				$mass+=$rec;
				$numResidues+=1;
			}

			$numPairs = ( $mass ) * ( $mass - 1 ) / 2;
	#               print"distances:@distances\n"; 
		       print"$mass\n";	
		       print"sumdistance:$withinClusterSumDistance->{$clus}->{$best_pdb}\n";
	               if ($numPairs!=0){                                  
			       #         delete $mass2pvalues->{$clus};
				if ($mass >2){
					my $withinClusterAvgDistance= $withinClusterSumDistance->{$clus}->{$best_pdb}/$numPairs;
					my %below;
					for ( my $simulation = 0; $simulation < $NSIMS ; $simulation++ ) {
						my $sumDists = 0;
						for ( my $i = 0; $i < $numPairs ; $i++ ) {
							my $randIndex = int( rand( $numDistances ) );
							$sumDists += $distances[$randIndex];
						} #foreach random pick
						my $avgDist = $sumDists / $numPairs;
						if ( $avgDist < $withinClusterAvgDistance) {
							$below{$clus}++;
							}
					} #foreach simulation
					my $permutationTestPValue;
					if(exists $below{$clus}){
						$permutationTestPValue = $below{$clus} / $NSIMS;

					} else{
						$permutationTestPValue=0;
						print "why zero\n";#debug
					}


					$OUT->print( join( "\t" , ( $clus , $genes , $best_pdb ,$numResidues,$mass , $numPairs , $withinClusterAvgDistance, $permutationTestPValue ) )."\n" );
				}

				elsif ( $mass == 2 ) {
					$OUT->print( join( "\t" , ( $clus , $genes , $best_pdb , $numResidues, $mass , $numPairs , $withinClusterSumDistance->{$clus}->{$best_pdb}, $clusterPvalue->{$clus}->{$best_pdb} ) )."\n" );
				} elsif ( $mass == 1 ) {
					$OUT->print( join( "\t" , ( $clus , $genes , $best_pdb , $numResidues, $mass , $numPairs , 0 , "NA" ) )."\n" );
				} #if mass block


			}
		}
        }#foreach cluster
                $OUT->close();
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
