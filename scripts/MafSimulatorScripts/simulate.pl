#!/usr/bin/env perl
#
#----------------------------------
# $Author: R Jay Mashl
# $Date: 2016-10-10 10:42:04 -0500 $
# $Revision: 0.1 $
# $URL: $
# $Doc: $ driver for simulated MAF generator for HotSpot3D
#----------------------------------
#
use strict;
use warnings;
no  warnings 'uninitialized';

use Getopt::Long;
use TGI::Mutpro::Main::MafSimulator;

# Defaults
my $distributionsOutDefault = "distributions.out";
my $mafStemOutDefault       = "maf";
my $numShufflesDefault      = 5;
my $offsetDefault           = 0;
my $sitesOutDefault         = "sites.out";
my $sizeHistoOutDefault     = "sizes.histo";
my $statsOutDefault         = "stats.out";

# Commands and options
my $command;
my $opts = {};
$opts->{'clustersListFile'}  =  undef;
$opts->{'clusterVarCounts'}  =  undef;
$opts->{'debug'}             =  0;
$opts->{'distributionsFn'}   =  $distributionsOutDefault;
$opts->{'entryNumberOffset'} =  $offsetDefault;
$opts->{'geneList'}          =  undef;
$opts->{'geneListFile'}      =  undef;
$opts->{'histosListFile'}    =  undef;
$opts->{'maf'}               =  undef;
$opts->{'mafStemOutFn'}      =  $mafStemOutDefault;
$opts->{'numShuffles'}       =  $numShufflesDefault;
$opts->{'randomUnique'}      =  0;
$opts->{'sitesFn'}           =  $sitesOutDefault;
$opts->{'sizeHistoFn'}       =  $sizeHistoOutDefault;
$opts->{'srcMaf'}            =  undef;
$opts->{'statsFn'}           =  $statsOutDefault;
$opts->{'transcripts'}       =  undef;

GetOptions(
    'clustersListFile=s' => \$opts->{'clustersListFile'},
    'clusterVarCounts=s' => \$opts->{'clusterVarCounts'},
    'command=s'          => \$command,
    'debug'              => \$opts->{'debug'},
    'distributions=s'    => \$opts->{'distributionsFn'},
    'geneList=s'         => \$opts->{'geneList'},
    'geneListFile=s'     => \$opts->{'geneListFile'},
    'histosListFile=s'   => \$opts->{'histosListFile'},
    'maf=s'              => \$opts->{'maf'},
    'num-shuffles=i'     => \$opts->{'numShuffles'},
    'offset=i'           => \$opts->{'entryNumberOffset'},
    'out-maf-stem=s'     => \$opts->{'mafStemOutFn'},
    'random-unique'      => \$opts->{'randomUnique'},
    'site-stats=s'       => \$opts->{'statsFn'},
    'sites=s'            => \$opts->{'sitesFn'},
    'size-histo=s'       => \$opts->{'sizeHistoFn'},
    'source-maf=s'       => \$opts->{'srcMaf'},
    'transcripts=s'      => \$opts->{'transcripts'},
    ) or abort();

# Check for errors
if( !$command ) { abort(); }
if( $command eq "getCoverage" && ( !defined $opts->{'transcripts'} || !defined $opts->{'maf'} ) ) {
    print "\nError: Please specify both a transcripts file and a maf file.\n\n"; abort();
}
if( $command eq "randomize" && ( !defined $opts->{'sitesFn'} || !defined $opts->{'statsFn'} ) ) {
    print "\nError: Please specify both sites and stats files.\n\n"; abort();
}
if( $command eq "generateMafs" && ( !defined $opts->{'srcMaf'} || !defined $opts->{'distributionsFn'} ) ) {
    print "\nError: Please specify source maf (template) file and random distributions files\n\n"; abort();
}
if( $command eq "getSizeHisto" ) {
    if( !defined $opts->{'clustersListFile'} ) {
        print "\nError: Please specify clusters list file for HotSpot3D *.clusters files\n\n"; abort();
    } elsif( !defined $opts->{'distributionsFn'} ) {
        print "\nError: Please specify distributions file\n\n"; abort();
    }
}
if( defined $opts->{'geneList'} && defined $opts->{'geneListFile'} ) {
    print "\nError: Please either gene list or gene file, not both.\n\n"; abort();
}

my $result;

if( $command eq "getCoverage" ) { $result = MafSimulator::getCoverage( $opts ); }
elsif( $command eq "randomize" ) { $result = MafSimulator::randomize( $opts ); }
elsif( $command eq "generateMafs" ) { $result = MafSimulator::generateMafs( $opts ); }
elsif( $command eq "getSizeHisto" ) { $result = MafSimulator::getSizeHisto( $opts ); }
elsif( $command eq "mergeHistos" ) { $result = MafSimulator::mergeHistos( $opts ); }
else { print "\nError: Please specify a valid command.\n\n"; abort(); }


sub help {
    return <<EOF;


Usage:
  $0  --command getCoverage  --transcripts <csv_file>  --maf <maf_file>  [--geneList <hugo_gene1>,<hugo_gene2>,... | --geneListFile <gene_list>]  [--sites <outfile> [default: $sitesOutDefault]]  [--site-stats <outfile> [default: $statsOutDefault]]

  $0  --command randomize    --sites <infile>  --site-stats <infile>  [--distributions <outfile>]  [--geneList <hugo_gene1>,<hugo_gene2>,... |  --geneListFile <gene_list>]  [--num-shuffles <int> [default: $numShufflesDefault]]  [--random-unique]  [--offset <int> [default: $offsetDefault]]

  $0  --command generateMafs --source-maf <maf_file>  --distributions <infile>  [--out-maf-stem <filestem> [default: $mafStemOutDefault]]  [--geneList <hugo_gene1>,<hugo_gene2>,... |  --geneListFile <gene_list>]  [--offset <int> [default: $offsetDefault]]

  $0  --command getSizeHisto --clustersListFile <list_file>  --distributions <infile>  [--size-histo <outfile> [default: sizes.histo]]  [--clusterVarCounts <outfile>]

  $0  --command mergeHistos --histosListFile <list_file>  [--size-histo <outfile> [default: sizes.histo]]


EOF
}

sub abort {
    print help();
    exit 0;
}
