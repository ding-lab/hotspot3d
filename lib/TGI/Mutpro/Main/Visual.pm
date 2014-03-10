package TGI::Mutpro::Main::Visual;
#
#----------------------------------
# $Authors: Beifang Niu 
# $Date: 2014-01-14 14:34:50 -0500 (Tue Jan 14 14:34:50 CST 2014) $
# $Revision:  $
# $URL: $
# $Doc: $ proximity pairs visualization
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
    $this->{_PROXIMITY_FILE} = undef;
    $this->{_OUTPUT_DIR} = getcwd;
    $this->{_PDB_DIR} = getcwd;
    $this->{_PYMOL_DIR} = '/usr/bin/pymol';
    $this->{_MUT_COLOR} = 'red';
    $this->{_MUT_STYLE} = 'spheres';
    $this->{_BG_COLOR} = 'black';
    $this->{_NO_LABEL} = undef;
    $this->{_STAT} = undef;
    bless $this, $class;
    $this->process();
    return $this;
}

sub process {
    my $this = shift;
    my ( $help, $options );
    unless( @ARGV ) { die $this->help_text(); }
    $options = GetOptions (
        'proximity-file=s' => \$this->{_PROXIMITY_FILE},
        'output-dir=s' => \$this->{_OUTPUT_DIR},
        'pdb-dir=s' => \$this->{_PDB_DIR},
        'pymol-dir=s' =>\$this->{_PYMOL_DIR},
        'mut-color=s' => \$this->{_MUT_COLOR},
        'mut-style=s' => \$this->{_MUT_STYLE},
        'no-label' => \$this->{_NO_LABEL},
        'bg-color=s' => \$this->{_BG_COLOR},
        'help' => \$help,
    );
    if ( $help ) { print STDERR help_text(); exit 0; }
    unless( $options ) { die $this->help_text(); }
    unless( $this->{_PROXIMITY_FILE} ) { warn 'You must provide a proximity file ! ', "\n"; die $this->help_text(); }
    unless( -e $this->{_PROXIMITY_FILE} ) { warn ' proximity file is not exist  ! ', "\n"; die $this->help_text(); }
    map{ $this->{_STAT}{$_} = 0; } qw( num_muts pdb pairs );
    my $muts_to_show_ref = $this->closePairs( $this->{_PROXIMITY_FILE} );
    # download pdbs
    $this->downloadPDBs( $muts_to_show_ref );
    # mark mutations
    $this->markMutations( $muts_to_show_ref );
    return 1;
}

# reading close pairs 
sub closePairs {
    my ( $this, $closepairf ) = @_;
    my ( %mutsToShow );
    my $fh = new FileHandle;
    unless( $fh->open($closepairf) ) { die "Could not open proximity file\n" };
    while ( my $a = <$fh> ) {
        chomp($a);
        my @t = split /\t/, $a;
        my ($taac0, $tchain0, $pos0, $taac1, $tchain1, $pos1, $pdbs) = @t[4,5,6,13,14,15,19];
        my @pdbIds = map{ @_ = split /\s+/; $_[1]} split /\|/, $pdbs;
        my ($chain0) = $tchain0 =~ /\[(\D)\]/;
        my ($chain1) = $tchain1 =~ /\[(\D)\]/;
        next unless( $chain0 and $chain1 );
        my ($aac0) = $taac0 =~ /p\.(\S+)/; 
        my ($aac1) = $taac1 =~ /p\.(\S+)/;
        next unless( $aac0 and $aac1 );
        map{ $mutsToShow{$_}{$chain0}{$pos0}{$aac0} = 1; } @pdbIds;
        map{ $mutsToShow{$_}{$chain1}{$pos1}{$aac1} = 1; } @pdbIds;
    }
    $fh->close();
    return \%mutsToShow;
}
# download PDBs  
sub downloadPDBs {
    my ( $this, $mutstoshowf ) = @_;
    foreach my $pdbid ( keys %$mutstoshowf ) {
        my $tf = $this->{_PDB_DIR}."/$pdbid.pdb";
        next if ( -e $tf );
        print STDERR "Downloading $pdbid.pdb ...\n";
        my $fh = new FileHandle;
        unless( $fh->open(">$tf") ) { die "Could not creat PDB file\n" };
        my $pdbUrl = "http://www.rcsb.org/pdb/files/$pdbid.pdb";
        my $tcon = get($pdbUrl);
        print $fh $tcon;
        $fh->close();
    }
    return 1;
}
# mark mutations   
sub markMutations {
    my ( $this, $mutstoshowf ) = @_;
    foreach my $pdbid ( keys %$mutstoshowf ) {
        my $tf = $this->{_OUTPUT_DIR}."/$pdbid.pdb";
        warn "$tf is not exist! \n" unless( -e $tf );
        next unless( -e $tf );
        my $i = 1; my $PyMol = "/usr/bin/pymol";
        $PyMol = $this->{_PYMOL_DIR}; 
        my $pcmd = "$PyMol $tf -x -c -W 800 -H 600 -d \'bg_color $this->{_BG_COLOR}; set internal_gui=0; show_as cartoon; spectrum b;' ";
        foreach my $iterchain ( keys %{$mutstoshowf->{$pdbid}} ) {
            foreach my $iterloc ( keys %{$mutstoshowf->{$pdbid}->{$iterchain}} ) { my $aacs = "";
                foreach my $iteraac ( keys %{$mutstoshowf->{$pdbid}->{$iterchain}->{$iterloc}} ) { $aacs .= $iteraac.'|'; }
                chop( $aacs ) if ( $aacs ); my $one_mut = ""; 
                if ( $this->{_NO_LABEL} ) { $one_mut = " -d 'sele t$i, (resi $iterloc and chain $iterchain); color $this->{_MUT_COLOR}, selection=t$i; show $this->{_MUT_STYLE}, selection=t$i;'";
                } else { $one_mut = " -d 'sele t$i, (resi $iterloc and chain $iterchain); color $this->{_MUT_COLOR}, selection=t$i; show $this->{_MUT_STYLE}, selection=t$i; set label_font_id, 10; set label_position, (3,2,1); set label_size, 20; label t$i and name CA, \"$aacs\";'"; }
                $pcmd .= $one_mut;
                $i++;
            }
        }
        $pcmd .= " -d 'save $this->{_OUTPUT_DIR}/$pdbid.pse ' > /dev/null";
        system( $pcmd ); #print STDERR $pcmd; #print STDERR "\n";
    }
    return 1; 
}

sub help_text{
    my $this = shift;
        return <<HELP

Usage: 3dproximity visual [options]

--proximity-file        Input mutation proximity file ( mutations close each other/COSMIC/region of interest(ROI))
--pymol-dir             PyMoL program location
--output-dir		Output directory for PyMol scripts
--pdb-dir               PDB file directory

--bg-color              background color, default: black
--mut-color             mutation color, default: red
--mut-style             mutation style, default: spheres
--no-label              don't show amino acid change label

--help			this message

HELP

}

1;

