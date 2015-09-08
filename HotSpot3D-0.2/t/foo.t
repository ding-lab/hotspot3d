use strict;
use warnings;

use Test::Most;
use TGI::Mutpro::Preprocess::Complicated;

TGI::Mutpro::Preprocess::Complicated::print_stuff();
ok(1, '1 is definitely okay');

done_testing;

