package TGI::Mutpro::Preprocess::Uppro;
#
#----------------------------------
# $Original authors: Beifang Niu and Adam D Scott
# $Modified by: Fernanda Martins Rodrigues @WashU (fernanda@wustl.edu; mrodrigues.fernanda@gmail.com)
# $Date: 2014-01-14 14:34:50 -0500 (Tue Jan 14 14:34:50 CST 2014) $
# $Revision: 2023-03-15 $
# $URL: $
# $Doc: $ create & update proximity files (the first step of preprocessing procedure). This has been adapted on Mar 15 2023 to handle AlphaFold DB (v4) files
#----------------------------------
#
use strict;
use warnings;
our $VERSION = '0.4';

use Carp;
use Cwd;
use Getopt::Long;
use LWP::Simple;
use IO::File;
use FileHandle;
use List::MoreUtils qw( uniq );

use TGI::Mutpro::Preprocess::Uniprot;
use TGI::Mutpro::Preprocess::HugoGeneMethods;
use TGI::Files::File;
use TGI::Files::MAF;
use TGI::Files::List;

my $MINDISTANCE = "shortest";
my $AVGDISTANCE = "average";
my $BSUB = "bsub";
my $BSUBC1 = "compute1";
my $LOCAL = "local";
my $NONE = "none";

sub new {
    my $class = shift;
    my $this = {};
    my %sub_cmds = ( 
        'output_dir'  => getcwd,
        'max_3d_dis'  => 100,
        'min_seq_dis' => 0,
        'hold' => 0,
        'status' => undef,
        # MODIFIED ON 03-15-2023: changed pdb_file_dir to alphafolddb,
        'alphafolddb_file_dir' => undef,
        'genes' => undef,
        'distance_measure' => undef,
        #'drugport_file' => undef, 
        'cmd_list_submit_file' => "cmd_list_submit_file", 
		'parallel' => undef ,
		'max_processes' => undef ,
    );
    map{ $this->{$_} = $sub_cmds{$_} } keys %sub_cmds;
    bless $this, $class;
    $this->process();
    return $this;
}

sub process {
    my $this = shift;
    my ($help, $options);
    unless (@ARGV) { die $this->help_text(); }
    $options = GetOptions (
        '3d-distance-cutoff=i'    => \$this->{'max_3d_dis'},
        'linear-distance-cutoff=i'   => \$this->{'min_seq_dis'},
        'output-dir=s'    => \$this->{'output_dir'},
        # MODIFIED ON 03-15-2023: changed pdb_file_dir to alphafolddb,
        'alphafolddb-file-dir=s'  => \$this->{'alphafolddb_file_dir'},
        'gene-file=s'  => \$this->{'genes'},
        'measure=s'  => \$this->{'distance_measure'},
        #'drugport-file=s' => \$this->{'drugport_file'},
        'cmd-list-submit-file=s' => \$this->{'cmd_list_submit_file'},
        'hold'   => \$this->{'hold'},
		'parallel=s' => \$this->{'parallel'} ,
		'max-processes=s' => \$this->{'max_processes'} ,
        'help' => \$help,
    );
    if ($help) { print STDERR help_text(); exit 0; }
    unless ($options) { die $this->help_text(); }
    map{
		unless($this->{$_} and (-e $this->{$_})) {
			warn " $_ does not exist ! \n";
			die $this->help_text();
		}
	# MODIFIED ON 03-15-2023: changed pdb_file_dir to alphafolddb,
	} qw( output_dir alphafolddb_file_dir );
	if ( defined $this->{'parallel'} ) {
		if ( $this->{'parallel'} eq $LOCAL ) {
			if ( not defined $this->{'max_processes'} ) {
				die "HotSpot3D::Uppro error: local parallelization requested, but max-processes was not set\n    CAUTION: make sure you know your max CPU processes !\n";
			}
		} elsif ( $this->{'parallel'} ne $BSUB and $this->{'parallel'} ne $NONE and $this->{'parallel'} ne $BSUBC1 ) {
			my $error = "HotSpot3D::Uppro error: parallelization type not set, please use either: \n";
			$error .= "\tbsub (if you have an bsub/LSF server), \n";
			$error .= "\tcompute1 (if you want to use the LSF server on compute1 - only applicable to users members of the Ding Lab at WashU, \n";
			$error .= "\tlocal (if you want to parallelize on your CPU), or \n";
			$error .= "\tnone (if you want to run calpro for each protein serially) \n";
			die $error;
		}
	} elsif ( not defined $this->{'parallel'} ) {
		my $error = "HotSpot3D::Uppro error: parallelization type not set, please use either: \n";
		$error .= "\tbsub (if you have an bsub/LSF server), \n";
		$error .= "\tcompute1 (if you want to use the LSF server on compute1 - only applicable to users members of the Ding Lab at WashU, \n";
		$error .= "\tlocal (if you want to parallelize on your CPU), or \n";
		$error .= "\tnone (if you want to run calpro for each protein serially) \n";
		die $error;
	}
    my $pro_dir = "$this->{'output_dir'}\/proximityFiles";
    my $inpro_dir = "$pro_dir\/inProgress";
    # MODIFIED ON 03-15-2023: changed pdbcor_dir to alphafolddbcor_dir and logfile name from hugo.uniprot.pdb.csv to hugo.uniprot.alphafolddb.csv
    my $alphafolddbcor_dir = "$pro_dir\/alphafolddbCoordinateFiles";
    my $log_file = "$this->{'output_dir'}\/hugo.uniprot.alphafolddb.csv";
    my $log_dir = "$this->{'output_dir'}\/Logs\/";
    my $measure = $this->{'distance_measure'};
    unless (-e $pro_dir) { mkdir($pro_dir) || die "HotSpot3D Uppro Error: can not make $pro_dir !\n"; }
    unless (-e $log_file) { if (system("touch $log_file") != 0) { die "HotSpot3D Uppro Error: can not make hugo uniprot file !\n"; } }
    unless (-e $inpro_dir) { mkdir($inpro_dir) || die "HotSpot3D Uppro Error: can not make $inpro_dir !\n"; }
    # MODIFIED ON 03-15-2023: changed pdbcor_dir to alphafolddbcor_dir
    unless (-e $alphafolddbcor_dir) { mkdir($alphafolddbcor_dir) || die "HotSpot3D Uppro Error: can not make $alphafolddbcor_dir !\n"; }
    unless ( -e $log_dir ) { mkdir( $log_dir ) || die "HotSpot3D Uppro Error: can not make $log_dir !\n"; }
	unless ( $measure eq $MINDISTANCE or $measure eq $AVGDISTANCE ) {
		die "HotSpot3D::Uppro error: measure not set or recognized (use either shortest or average)\n";
	}
    my %uniprotid_toupdate;
    my $uniprot_to_structureref = $this->current_structures($log_file);
    my $uniprot_fileref = $this->currentuniprot_files($pro_dir);
    my $fh = new FileHandle;
    unless ($fh->open(">$log_file")) { die "Could not open hugo uniprot file !\n" };
	print STDOUT "Creating ".$log_file."\n";
	# MODIFIED ON 03-15-2023: changed pdb_ref to alphafolddb_ref
    my ($hugo_id, $alias_ref, $previous_ref, $alias_list, $uniprot_id, $uniprot_ref, $alphafolddb_ref);
	my $hugogene_ref;
	my ( $list , @fields );
	if ( $this->{'genes'} ) { 
		print STDOUT "Getting genes from ".$this->{'genes'}."\n";
		my $genesFH = new TGI::Files::List( $this->{'genes'} );
		$genesFH->open();
		$list = $genesFH->getList( 0 );
		$genesFH->close();
		my $ngenes = scalar keys %{$list};
		print STDOUT "Found ".$ngenes." (may include header)\n";
	}
	$hugogene_ref = TGI::Mutpro::Preprocess::HugoGeneMethods::makeHugoGeneObjects( $list );
    foreach $hugo_id (sort keys %{$hugogene_ref}) {
        print STDOUT 'HUGO: ', "$hugo_id\n";
        $alias_ref = $$hugogene_ref{$hugo_id}->getAllAliases();
        $previous_ref =  $$hugogene_ref{$hugo_id}->getAllPreviousSymbols();
        $alias_list = "";
        map { $alias_list .= "$_ "; } keys %{$alias_ref};
        map { $alias_list .= "$_ "; } keys %{$previous_ref};
        if ($alias_list !~ /\w+/) {$alias_list = "N/A"; };
        $uniprot_id = $$hugogene_ref{$hugo_id}->uniprot();
        if (!defined $uniprot_id) {
			$fh->print("$hugo_id\tN/A\tN/A\t$alias_list\n");
			next;
		}
        $uniprot_ref = TGI::Mutpro::Preprocess::Uniprot->new($uniprot_id);
        # MODIFIED ON 03-15-2023: changed pdb_ref to alphafolddb_ref
        $alphafolddb_ref = $uniprot_ref->annotations("AlphaFoldDB");
        if (!defined $alphafolddb_ref || scalar(@{$alphafolddb_ref}) == 0) {
			$fh->print( "$hugo_id\t$uniprot_id\tN/A\t$alias_list\n" );
			next;
		}
        if (!defined $$uniprot_fileref{$uniprot_id}) {
			$uniprotid_toupdate{$uniprot_id} = 1;
		}
        $fh->print("$hugo_id\t$uniprot_id\t");
        map {
			my ($alphafolddb_id) = $_ =~ /^(\w+)\;/;
			if (defined $alphafolddb_id) {
				if (!defined $$uniprot_to_structureref{$uniprot_id}{$alphafolddb_id}) {
					$uniprotid_toupdate{$uniprot_id} = 1;
				}
				$fh->print("$alphafolddb_id ");
			}
		} @{$alphafolddb_ref};
        $fh->print( "\t$alias_list\n" );
    }
    $fh->close();
	my $cmd_list_submit_file_fh = new FileHandle; #TGI::Files::File->new( $this->{'cmd_list_submit_file'} ) ;
        unless( $cmd_list_submit_file_fh -> open ( $this->{'cmd_list_submit_file'}, "w" ) )  { die "HotSpot3D Uppro Error: Could not open cmd file (".$this->{'cmd_list_submit_file'}.")"; }
	print STDOUT "Creating ".$this->{'cmd_list_submit_file'}."\n";
	my $command = "";
	if ( $this->{'parallel'} eq $BSUB ) {
		$command = $this->makeBSUBCalproCommand( $inpro_dir , $cmd_list_submit_file_fh , \%uniprotid_toupdate );
	# MODIFIED ON 03-15-2023: added BSUB for compute1 option - specific to dinglab folks 
	} elsif ( $this->{'parallel'} eq $BSUBC1 ) {
		$command = $this->makeBSUBCalproCommandCompute1( $inpro_dir , $cmd_list_submit_file_fh , \%uniprotid_toupdate );
	} elsif ( $this->{'parallel'} eq $LOCAL ) {
		$command = $this->makeLocalCalproCommand( $inpro_dir , $cmd_list_submit_file_fh , \%uniprotid_toupdate );
	} else {
		$command = $this->makeSerialCalproCommand( $inpro_dir , $cmd_list_submit_file_fh , \%uniprotid_toupdate );
	}
	if ( not $this->{'hold'} ) {
		system( $command );
	}
	$cmd_list_submit_file_fh->close();

    return 1;
}

sub makeCalproCommand {
	my ( $this , $protein ) = @_;
	my $command = "hotspot3d calpro ";
	$command .= "--output-dir=".$this->{'output_dir'}." ";
	$command .= "--alphafolddb-file-dir=".$this->{'alphafolddb_file_dir'}." ";
	$command .= "--uniprot-id=".$protein." ";
	$command .= "--3d-distance-cutoff=".$this->{'max_3d_dis'}." ";
	$command .= "--linear-cutoff=".$this->{'min_seq_dis'}." ";
	$command .= "--measure=".$this->{'distance_measure'};
	return $command;
}

sub makeBSUBCalproCommand {
	my ( $this , $inpro_dir , $cmd_list_submit_file_fh , $uniprotid_toupdate ) = @_;
    my $log_dir = "$this->{'output_dir'}\/Logs\/";
	map {
		system("touch $inpro_dir/$_.ProximityFile.csv");
		my $command = "";
		$command = "bsub -oo ".$log_dir.$_.".calpro.log -R 'select[type==LINUX64 && mem>16000] rusage[mem=16000]' -M 16000000 ";
		$command .= "'".$this->makeCalproCommand( $_ )."'";
		print STDOUT $command."\n"; 
		$cmd_list_submit_file_fh->print( $command."\n" );
	} keys %{$uniprotid_toupdate};
}

# MODIFIED ON 03-15-2023: Added code to create BSUB commands suitable for compute1 - specific to dinglab folks only
sub makeBSUBCalproCommandCompute1 {
	my ( $this , $inpro_dir , $cmd_list_submit_file_fh , $uniprotid_toupdate ) = @_;
    my $log_dir = "$this->{'output_dir'}\/Logs\/";
	map {
		system("touch $inpro_dir/$_.ProximityFile.csv");
		my $command = "";
		$command = "bsub -oo ".$log_dir.$_.".calpro.log -R 'select[mem>16000] rusage[mem=16000]' -M 16000000 -q general -G compute-dinglab -a 'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)' ";
		$command .= "'".$this->makeCalproCommand( $_ )."'";
		print STDOUT $command."\n"; 
		$cmd_list_submit_file_fh->print( $command."\n" );
	} keys %{$uniprotid_toupdate};
}

sub makeSerialCalproCommand {
	my ( $this , $inpro_dir , $cmd_list_submit_file_fh , $uniprotid_toupdate ) = @_;
    my $log_dir = "$this->{'output_dir'}\/Logs\/";
	my $serialCommand = "echo Serially Running Calpro Commands";
	map {
		system("touch $inpro_dir/$_.ProximityFile.csv");
		my $command = $this->makeCalproCommand( $_ );
		print STDOUT $command."\n"; 
		$cmd_list_submit_file_fh->print( $command."\n" );
		$serialCommand .= ";".$command;
	} keys %{$uniprotid_toupdate};
	return $serialCommand;
}

# MODIFIED ON 03-15-2023: updated parameters to reflect pdb to alphafold changes
sub makeLocalCalproCommand {
	my ( $this , $inpro_dir , $cmd_list_submit_file_fh , $proteinsToUpdate ) = @_;
    my $log_dir = "$this->{'output_dir'}\/Logs\/";
	my $parallel = "#!/usr/bin/perl\n#command script for parallelizing the calpro steps from HotSpot3D uppro\n\n";
	$parallel .= "use strict;\nuse warnings;\nuse Parallel::ForkManager;\n\n";
	$parallel .= "my \$usage = \'perl cmd_list_submit_file #n-processes#\n\'; die \$usage unless \@ARGV == 1;\n";
	$parallel .= "my ( \$Nprocesses ) = \@ARGV;\n\n";
	$parallel .= "my \$pm = Parallel::ForkManager->new( \$Nprocesses );\n";
	$parallel .= "DATA_LOOP:\n";
	$parallel .= "while ( my \$protein = <DATA> ) {\n";
	$parallel .= "\tchomp( \$protein );\n";
	$parallel .= "\tmy \$pid = \$pm->start and next DATA_LOOP;\n";
	$parallel .= "\tmy \$cmd = \'hotspot3d calpro ";
	$parallel .= "--output-dir=".$this->{'output_dir'}." ";
	$parallel .= "--alphafolddb-file-dir=".$this->{'alphafolddb_file_dir'}." ";
	$parallel .= "--uniprot-id=\'\.\$protein\.\' ";
	$parallel .= "--3d-distance-cutoff=".$this->{'max_3d_dis'}." ";
	$parallel .= "--linear-cutoff=".$this->{'min_seq_dis'}." ";
	$parallel .= "--measure=".$this->{'distance_measure'}." ";
	$parallel .= ">".$log_dir."\'\.\$protein\.\'.calpro.log\';\n";
	$parallel .= "\tsystem( \$cmd );\n";
	$parallel .= "\t\$pm->finish;\n";
	$parallel .= "}\n";
	$parallel .= "__DATA__\n";
	foreach my $protein ( keys %{$proteinsToUpdate} ) {
		$parallel .= $protein."\n";
	}
	$cmd_list_submit_file_fh->print( $parallel );
	my $command = "perl ".$this->{'cmd_list_submit_file'}." ".$this->{'max_processes'};
	print STDOUT $command."\n";
	return $command;
}

# Check which Uniprot files are on disk 
# Return ref to hash with 
# key = Uniprot Id 
sub currentuniprot_files {
    my ($this, $dir) = @_;
    my (%uniprot_ids, $file);
    opendir(DIR, $dir) || die "Could not open '$dir': $!";
    map {
		if ($_ =~ /(\w+)\.ProximityFile\.csv/) {
			$uniprot_ids{$1} = 1;
		}
	} (readdir DIR);
    closedir DIR;

    return \%uniprot_ids;
}

# MODIFIED ON 03-15-2023: changed to take alphafolddb_id instead of pdb_id
# Return ref to hash with 
# key = uniprot_id, $alphafolddb_id; 
# value = 1
sub current_structures {
    my ($this, $logfilef) = @_;
    my ($uniprot_id, $alphafolddb_id, %uniprot_tostructure, $alphafolddb_list);
    my $fh = new FileHandle;
    unless ($fh->open($logfilef)) { die "Could not open hugo uniprot file\n" };
    map {
		chomp;
		(undef, $uniprot_id, $alphafolddb_list) = split /\t/, $_;
		unless ($uniprot_id eq "N/A") {
			map {
				$uniprot_tostructure{$uniprot_id}{$_} = 1;
			} split /s+/, $alphafolddb_list;
		}
	} $fh->getlines;
    $fh->close();
    return \%uniprot_tostructure;
}
sub help_text{
    my $this = shift;
        return <<HELP
Usage: hotspot3d uppro [options]
                             REQUIRED
--output-dir                 Output directory of proximity files
--alphafolddb-file-dir       AlphaFold DB file directory 
--measure                    Distance measure between residues (shortest or average)
--parallel                   Parallelization to use (bsub, compute1, local, or none)
                                 use bsub if you want to use an LSF server
                                 use compute1 if you want to use the LSF server on compute1 (this is only applicable to users members of the Ding Lab at WashU)
                                 use local if you want to use your CPU (CAUTION: make sure you know your max CPU processes)
                                 use none if you want to run calpro for each protein serially
                             OPTIONAL
--max-processes              Set if using parallel type local (CAUTION: make sure you know your max CPU processes)
--gene-file                  File with HUGO gene names in the first column (like a .maf)
--3d-distance-cutoff         Maximum 3D distance (<= Angstroms), defaul: 100
--linear-distance-cutoff     Minimum linear distance (> peptides), default: 0
--cmd-list-submit-file       Batch jobs file to run calpro step in parallel, default: cmd_list_submit_file
--hold                       Do not submit batch jobs, just write cmd_list_submit_file, default: submits (takes no input)
--help                       this message
HELP
}
1;
__END__
 
=head1 NAME
TGI::Mutpro::Preprocess::Uppro - Create & update proximity files.
=head1 SYNOPSIS
  use TGI::Mutpro::Preprocess::Uppro;
=head1 DESCRIPTION
TGI::Mutpro::Preprocess::Uppro is to be used to create & update proximity files.
It is the first step of preprocessing procedure.
=head1 AUTHOR
Beifang Niu E<lt>beifang.cn@gmail.comE<gt>
For questions regarding HotSpot3d's AlphaFold DB version, please contact Fernanda Martins Rodrigues E<lt>mrodrigues.fernanda@gmail.comE<gt>
=head1 SEE ALSO
https://github.com/ding-lab/hotspot3d/
For version suitable for AlphaFold DB, please see: https://github.com/ding-lab/hotspot3d/blob/alphaFold_implementation
=head1 LICENSE
This library is free software with MIT licence; you can redistribute it and/or modify
it under the same terms as Perl itself.
=cut
