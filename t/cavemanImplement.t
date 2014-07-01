##########LICENCE##########
#  Copyright (c) 2014 Genome Research Ltd.
#
#  Author: David Jones <cgpit@sanger.ac.uk>
#
#  This file is part of cavemanWrapper.
#
#  cavemanWrapper is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Affero General Public License as published by the Free
#  Software Foundation; either version 3 of the License, or (at your option) any
#  later version.
#
#  This program is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
#  FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
#  details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program. If not, see <http://www.gnu.org/licenses/>.
##########LICENCE##########

use strict;
use Test::More;
use File::Spec;
use Const::Fast qw(const);
use FindBin qw($Bin);

const my $MODULE => 'Sanger::CGP::Caveman::Implement';
const my $LINES_PER_FQ => 4;

my $test_data = "$Bin/../testData";

my $good_fq = File::Spec->catfile($test_data, '1_2.fq');
my $bad_fq = File::Spec->catfile($test_data, '1_2_bad.fq');

subtest 'Initialisation checks' => sub {
  use_ok($MODULE);
};

subtest 'Line count checks' => sub {
	is($LINES_PER_FQ,Sanger::CGP::Caveman::Implement::file_line_count($good_fq),'Line count correct');
	isnt($LINES_PER_FQ, Sanger::CGP::Caveman::Implement::file_line_count($bad_fq),'Line count correct (2)');
};

done_testing();
