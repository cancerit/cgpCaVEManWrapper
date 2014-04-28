use strict;
use Test::More;
use Test::Fatal;
use File::Spec;
use Try::Tiny qw(try catch finally);
use Const::Fast qw(const);

const my $MODULE => 'PCAP::Caveman';

my $test_data = "$Bin/../testData";

my $line_count_file = File::Spec->catfile($test_data, '1_2.fq');

subtest 'Initialisation checks' => sub {
  use_ok($MODULE);
};

subtest 'Line count checks' => sub {
	my $count = PCAP::Caveman::file_line_count($line_count_file);
	is(4,$count,'Line count correct');
}

done_testing();
