use strict;
use warnings;

use Test::Most;
use TGI::Mutpro::Complicated;

TGI::Mutpro::Complicated::print_stuff();
ok(1, '1 is definitely okay');

done_testing;

