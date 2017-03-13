package TGI::Files::File;
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

use IO::File;
use FileHandle;
use List::MoreUtils qw( first_index );
use Scalar::Util qw( openhandle );
use TGI::Data::StringTemplate;

sub new {
	my $proto = shift;
    my $class = ref( $proto ) || $proto;
    my $this = {};
	$this->{'handle'} = new FileHandle;
	$this->{'file_name'} = undef;
	$this->{'header'} = {};
	$this->{'required'} = undef;#\@REQUIRED;
	$this->{'required_indices'} = undef;#\@REQUIREDINDICES;
	if ( @_ ) { $this->{'file_name'} = shift; }
    bless $this, $class;
    return $this;
}

sub mapColumns {
	my $this = shift;
	my ( $required , $exception );
	$this->open();
	seek( $this->{'handle'} , 0 , 0 );
	my $headline = $this->header( "\t" );
	if ( @_ ) { $required = shift; } else { die "No required columns provided\n"; }
	if ( @_ ) { $exception = shift; } else { $exception = ""; }
	my $mafcols;
	my $mafi = 0;
	#print $headline."\n";
	my %mafcols = map{ ( $_ , $mafi++ ) } split( /\t/ , $headline );
	foreach my $needColumn ( @{$required} ) {
		#print "need column= ".$needColumn."...";
		unless( exists $mafcols{$needColumn} ) { die $exception."\n"; }
		#print $mafcols{$needColumn}."\n";
		push @{$mafcols} , $mafcols{$needColumn};
	}
	return $mafcols;
}

sub isOpen {
	my $this = shift;
	#print "isOpen: "."is open= ";
	if ( defined openhandle( $this->{'handle'} ) ) {
		#print "isOpen: "."yes\n";
		return 1;
	}
	#print "isOpen: "."no\n";
	return 0;
}

sub open {
    my $this = shift;
	my $hasHeader = 0;
	if ( @_ ) {
		$hasHeader = shift;
	}
	#print "open: "."Has header= ".$hasHeader."\n";
	#print "open: "."file name= ";
    if ( defined $this->{'file_name'} ) {
		#print $this->{'file_name'}."\n";
		if ( not $this->isOpen() ) {
			#print "open: "."opening\n";
			$this->{'handle'}->open( $this->{'file_name'} , "r" );
		} else {
			seek( $this->{'handle'} , 0 , 0 );
		}
		#print "open: "."opened, now getting header if has header\n";
		if ( $hasHeader && (not $this->{'header'}) ) {
			#print "open: "."setting header\n";
			$this->setHeader();
		}
        return 1;
    }
	#print "open: "."cannot open\n";
    return 0;
}

sub close {
    my $this = shift;
	if ( $this->isOpen() ) {
		$this->{'handle'}->close();
	}
    return;
}

sub getlines {
    my $this = shift;
	$this->open();
	return $this->{'handle'}->getlines;
}

sub getline {
	my $this = shift;
	$this->open();
	return $this->{'handle'}->getline;
}

sub setHeader {
    my $this = shift;
	my $closed = 0;
	#print "setHeader: "."setting header\n";
    if ( $this->isOpen() ) {
		#print "setHeader: "."to first line\n";
		seek( $this->{'handle'} , 0 , 0 );
    } else {
		#print "setHeader: "."need to open file\n";
		$closed = 1;
		$this->open( 1 );
	}
    my $i = 0;
	#print "setHeader: "."now setting from first line\n";
	my $line = $this->getline();
	chomp( $line );
	#print "setHeader: ".$line."\n";
    foreach my $column ( split( /\t/ , $line ) ) {
		#print "setHeader: ".$column." => ".$i."\n";
        $this->{'header'}->{$column} = $i++;
    }
	if ( $closed ) {
		$this->close();
	}
    return;
}

sub getHeader {
	my $this = shift;
	#print "getHeader: "."getting header\n";
	#print "getHeader: "."is header set? ";
	if ( scalar keys %{$this->{'header'}} == 0 ) {
		#print "getHeader: "."no\n";
		$this->setHeader();
	}
	#print "getHeader: "."header is set\n";
	my @header;
	my %header;
	#print "getHeader: "."sorting header hash\n";
	foreach my $key ( keys %{$this->{'header'}} ) {
		my $index = $this->{'header'}->{$key};
		#print "getHeader: ".$key." => ".$index."\n";
		$header{$index} = $key;
	}
	foreach my $index ( sort { $a <=> $b } keys %header ) {
		#print "getHeader: ".$index." => ".$header{$index}."\n";
		push @header , $header{$index};
	}
	#print "getHeader: "."header is= ( ".join( ", " , @header )." )\n";
	#print "getHeader: "."header has= ( ".join( ", " , keys %{$this->{'header'}} )." )\n";

	return \@header;
}

sub header {
	my $this = shift;
	my $delim = ", ";
	#print "header: "."have header? ";
	if ( @_ ) { $delim = shift; }
	if ( scalar keys %{$this->{'header'}} > 0 ) {
		#print "yes\n";
		my $header = join( $delim , @{$this->getHeader()} );
		#print "header: ".$header."\n";
		return $header;
	} else {
		#print "header: "."no\n";
		$this->setHeader();
		#print "header: "."now header is set (".join( ", " , @{$this->getHeader()} )."\n";
		return $this->header();
	}
	return;
}

sub setColumnIndex {
    my $this = shift;
    my $columnHead = shift;
    my $index = shift;
    $this->{'header'}->{$columnHead} = $index;
    return;
}

sub getColumnIndex {
    my $this = shift;
    my $columnHead = shift;
    if ( exists $this->{'header'}->{$columnHead} ) {
        return $this->{'header'}->{$columnHead};
    }
    return;
}

sub getField {
	my $this = shift;
	my $column = shift;
	my $line = shift;
	#print "getField: ".$column."\n";
	#print "getField: ".join( "... " , @{$line} )."\n";
	if ( exists $this->{'header'}->{$column} ) {
		#print "getField: ".$this->getColumnIndex( $column )."\n";
        return $line->[$this->getColumnIndex( $column )];
    }
    return undef;
}

sub setCustomColumn {
    my $this = shift;
    if ( @_ ) {
        my $column = shift;
        $this->{$column} = undef;
        if ( @_ ) {
            $this->setHeaderIndex( $this->{$column} , shift );
        } else {
            $this->setHeaderIndex( $this->{$column} , -1 );
        }
    }
    return;
}

sub requireColumn {
    my $this = shift;
    if ( @_ ) {
        my $column = shift;
        push @{$this->{'required'}} , $column;
    } else {
        warn "TGI::Files::File::requireColumn warning: Tried to require a column, but no column name was given\n";
    }
    return $this->{'required'};
}

sub requireColumnIndex {
    my $this = shift;
    if ( @_ ) {
        my $index = shift;
        push @{$this->{'required_indices'}} , $this->{$index};
    } else {
        warn "TGI::Files::File::requireColumn warning: Tried to require a column, but no column name was given\n";
    }
    return $this->{'required'};
}

sub requireColumns {
    my $this = shift;
    if ( @_ ) {
        my $columns = shift;
        foreach my $column ( @{$columns} ) {
            my $found = first_index{ $column eq $_ } @{$this->getHeader()};
            if ( $found == -1 ) {
                warn "TGI::Files::File::requireColumns warning: a column needed to be required is not in the header (".$column.")\n";
            } else {
                $this->requireColumn( $column );
                $this->requireColumnIndex( $found );
            }
        }
    }
    return;
}

sub getRequired {
	my $this = shift;
	my %required;
	foreach my $required ( @{$this->{'required'}} ) {
		$required{$required} = $this->getColumnIndex( $required );
	}
	return \%required;
}

sub setColumnIndices {
    my $this = shift;
    my $exception = "";
	$exception = "TGI::Files::File::setColumnIndices error: required columns not present ";
	$exception .= "in this .maf (".$this->{'file_name'}."). Need the ";
	$exception .= "following columns:\n\t";
	$exception .= join( "\n\t" , @{$this->{'required'}} )."\n";
    $this->setHeader();
    my $header = $this->getHeader();
    my $mafcols = ();
    if ( not $this->{'required'} ) {
        $this->requireColumns( $this->getHeader() );
    }
	$mafcols = $this->mapColumns( $this->{'required'} , $exception );
	#print join( "-- " , @{$mafcols} )."\n"; 
    return $mafcols;
}

sub read {
    my $this = shift;
    my $keyTemplate = new TGI::Data::StringTemplate();
    my $valueTemplate = new TGI::Data::StringTemplate();
    my $haveKey = 0;
    my $haveValue = 0;
    if ( @_ ) {
        my $temp = shift;
        $keyTemplate->addToTemplate( $temp );
        $haveKey = 1;
    }
    if ( @_ ) {
        my $temp = shift;
        $valueTemplate->addToTemplate( $temp );
        $haveValue = 1;
    }
    print STDOUT "\nReading in ".$this->{'file_name'}."...\n";
    seek( $this->{'handle'} , 0 , 0 );
    my $indices = $this->setColumnIndices();
    my %output;
    my @output;
    #print "read: "."Have key?\n";
    if ( $haveKey ) {
        #print "read: "."yes\nHave value?\n";
        if ( $haveValue ) {
            #print "read: "."yes\n";
            foreach my $line ( $this->getlines() ) {
                chomp( $line );
                my @line = split( /\t/ , $line );
                next if ( exists $this->{'header'}->{$line[0]} );
                my $key = "";
                my $value = "";
                $key = $keyTemplate->constructFromColumns( \@line , $this->getRequired() );
                $value = $valueTemplate->constructFromColumns( \@line , $this->getRequired() );
                #print "read: ".$key." => ".$value."\n";
                $output{$key} = $value;
            }
        } else {
            #print "read: "."no\n";
            foreach my $line ( $this->getlines() ) {
                chomp( $line );
                my @line = split( /\t/ , $line );
                next if ( exists $this->{'header'}->{$line[0]} );
                my $key = "";
                $key = $keyTemplate->constructFromColumns( \@line , $this->getRequired() );
                $output{$key} += 1;
                #print "read: ".$key." => ".$output{$key}."\n";
            }
        }
        return \%output;
    } else {
        #print "read: "."no\n";
        my @output = $this->getlines();
        #print "read: ".join( "\t" , @output )."\n";
        return \@output;
    }
    return;
}


1;
