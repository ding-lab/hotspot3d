package TGI::Data::StringTemplate;
#
#----------------------------------
# $Authors: Adam Scott
# $Date: 2016-10-27
# $Revision: v0.0 $
# $URL: $
# $Doc: $ key/value template for perl hashes
# 
#----------------------------------
#
use strict;
use warnings;

use List::MoreUtils qw( first_index );

sub new {
	my $proto = shift;
    my $class = ref( $proto ) || $proto;
	my $this = {};
	$this->{'parts'} = ();
    if ( @_ ) { 
		my $part = shift;
		if ( ref( $part ) eq 'ARRAY' ) {
			push @{$this->{'parts'}} , @{$part};
		} else {
			push @{$this->{'parts'}} , $part;
		}
	}
    bless $this, $class;
    return $this;
}

sub addToTemplate {
	my $this = shift;
	if ( @_ ) {
		my $part = shift;
		push @{$this->{'parts'}} , @{$part};
	}
	return $this->{'parts'};
}

sub construct {
	my $this = shift;
	if ( @_ ) {
		return join( shift , @{$this->{'parts'}} );
	} else {
		return join( "" , @{$this->{'parts'}} );
	}
}

sub getTemplate {
	my $this = shift;
	return $this->{'parts'};
}

sub constructFromColumns {
	my $this = shift;
	my $template = "";
	if ( @_ ) {
		my $line = shift;
		my $required = shift;
		#print join( " | " , @{$line} )."\n";
		foreach my $part ( @{$this->getTemplate()} ) {
			#print "-----".$part;
			if ( $part =~ /^-?\d+$/ ) {
				#print "-----";
				#print $line->[$part]."\n";
				$template .= $line->[$part];
			} else {
				if ( exists $required->{$part} ) {
					#print "=====".$line->[$required->{$part}]."\n";
					$template .= $line->[$required->{$part}];
				} else {
					#print "+++++".$part."\n";
					$template .= $part;
				}
			}
		}
	}
	return $template;
}

1;
