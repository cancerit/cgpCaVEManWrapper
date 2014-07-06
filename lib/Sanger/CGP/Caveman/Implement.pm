package Sanger::CGP::Caveman::Implement;

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
use warnings FATAL => 'all';
use File::Which qw(which);
use FindBin qw($Bin);
use autodie qw(:all);
use Const::Fast qw(const);

use Sanger::CGP::Caveman;
our $VERSION = Sanger::CGP::Caveman->VERSION;

use PCAP::Threaded;
use PCAP::Bam;

const my $CAVEMAN_SETUP => q{ setup -t %s -n %s -r %s -g %s -l %s -f %s -c %s -a %s};
const my $CAVEMAN_SPLIT => q{ split -i %d -f %s};
const my $CAVEMAN_MSTEP => q{ mstep -i %d -f %s};
const my $CAVEMAN_MERGE => q{ merge -c %s -p %s -f %s};
const my $CAVEMAN_ESTEP => q{ estep -i %d -e %s -j %s -k %f -g %s -o %s -v %s -w %s -f %s};
const my $MERGE_CAVEMAN_RESULTS => q{ mergeCavemanResults -o %s %s};

sub prepare {
  my $options = shift;
  $options->{'tumour_name'} = (PCAP::Bam::sample_name($options->{'tumbam'}))[0];
  $options->{'normal_name'} = (PCAP::Bam::sample_name($options->{'normbam'}))[0];
  return 1;
}

sub file_line_count {
	my $file = shift;
	my $contig_count = 0;
  {
    my $FH;
    open($FH, '<', $file) or die("Error trying to open $file: $!\n");
    	while(<$FH>){
    		my $line = $_;
    		$contig_count++;
    	}
    close($FH);
  }
  return $contig_count;
}

sub caveman_setup {
	# uncoverable subroutine
	my $options = shift;
	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);
	my $ref = $options->{'reference'};
	my $tumbam = $options->{'tumbam'};
	my $normbam = $options->{'normbam'};
	my $ignore = $options->{'ignore'};
	my $split_list_loc = $tmp."/splitList";
	my $results_loc = $tmp."/results";
	my $config = $options->{'cave_cfg'};
	my $alg_bean = $options->{'cave_alg'};


	my $command = _which('caveman') || die "Unable to find 'caveman' in path";

	$command .= sprintf($CAVEMAN_SETUP,
								$tumbam,
								$normbam,
								$ref,
								$ignore,
								$split_list_loc,
								$results_loc,
								$config,
								$alg_bean);

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);

	return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 'caveman_setup', 0);
}

sub caveman_split {
	# uncoverable subroutine
	my ($index,$options) = @_;

	return 1 if(exists $options->{'index'} && $index != $options->{'index'});
	my $tmp = $options->{'tmp'};
  my $config = $options->{'cave_cfg'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 'caveman_split', $index);

	my $command = _which('caveman') || die "Unable to find 'caveman' in path";
	$command .= sprintf($CAVEMAN_SPLIT,$index,$config);

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);
  	return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 'caveman_split', $index);
}

sub caveman_merge{
	# uncoverable subroutine
	my $options = shift;

	my $command = _which('caveman') || die "Unable to find 'caveman' in path";
	my $tmp = $options->{'tmp'};
	my $config = $options->{'cave_cfg'};
	my $prob_arr = $options->{'cave_parr'};
	my $cov_arr = $options->{'cave_carr'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 'caveman_merge', 0);

	$command .= sprintf($CAVEMAN_MERGE, $cov_arr, $prob_arr,$config);

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);

	return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 'caveman_merge', 0);
}

sub caveman_mstep{
	# uncoverable subroutine
	my ($index_in,$options) = @_;

	# first handle the easy bit, skip if limit not set
	return 1 if(!exists $options->{'limit'} && exists $options->{'index'} && $index_in != $options->{'index'});

	my @indicies = limited_xstep_indicies($options, $index_in);
  my $config = $options->{'cave_cfg'};
	my $tmp = $options->{'tmp'};
	for my $index(@indicies) {
    next if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 'caveman_mstep', $index);

    my $command = _which('caveman') || die "Unable to find 'caveman' in path";

    $command .= sprintf($CAVEMAN_MSTEP,
                    $index,
                    $config);

    PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);
    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 'caveman_mstep', $index);
  }
  return 1;
}

sub caveman_estep{
	# uncoverable subroutine
	my ($index_in,$options) = @_;

	# first handle the easy bit, skip if limit not set
	return 1 if(!exists $options->{'limit'} && exists $options->{'index'} && $index_in != $options->{'index'});

	my @indicies = limited_xstep_indicies($options, $index_in);
  my $config = $options->{'cave_cfg'};
	my $tmp = $options->{'tmp'};
	my $prob_arr = $options->{'cave_parr'};
	my $cov_arr = $options->{'cave_carr'};
	for my $index(@indicies) {
    next if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 'caveman_estep', $index);

    my $command = _which('caveman') || die "Unable to find 'caveman' in path";

    $command .= sprintf($CAVEMAN_ESTEP,
                    $index,
                    $options->{'normcn'},
                    $options->{'tumcn'},
                    $options->{'normcont'},
                    $cov_arr,
                    $prob_arr,
                    $options->{'species-assembly'},
                    $options->{'species'},
                    $config);

    PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);
    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 'caveman_estep', $index);
  }
  return 1;
}

sub caveman_merge_results {
	# uncoverable subroutine
	my $options = shift;
	my $tmp = $options->{'tmp'};
	my $out = $options->{'out_file'};
	my $target = $options->{'subvcf'};
	my $command = sprintf($MERGE_CAVEMAN_RESULTS,$out.".muts.vcf",$target);
	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0)
								unless (PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 'merge_muts', 0));
	PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 'merge_muts', 0);
	$target = $options->{'snpvcf'};
	$command = sprintf($MERGE_CAVEMAN_RESULTS,$out.".snps.vcf",$target);
	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0)
								unless (PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 'merge_snps', 0));
	PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 'merge_snps', 0);
	$target = $options->{'noanalysisbed'};
	$command = sprintf($MERGE_CAVEMAN_RESULTS,$out.".no_analysis.bed",$target);
	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0)
								unless (PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 'merge_no_analysis', 0));
	PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 'merge_no_analysis', 0);
	return 1;

}

sub limited_xstep_indicies {
  my ($options, $index_in) = @_;
  my @indicies;
  if(exists $options->{'limit'}) {
	  my $split_count = file_line_count($options->{'splitList'});
    # main script checks index is not greater than limit or < 1
	  my $base = $index_in;
	  while($base <= $split_count) {
	    push @indicies, $base;
	    $base += $options->{'limit'};
	  }
	}
	else {
	  push @indicies, $index_in;
	}
	return @indicies;
}

sub concat {
	# uncoverable subroutine
	my $options = shift;
	my $tmp = $options->{'tmp'};
	my $out = $options->{'out_file'};
	my $target = $options->{'target_files'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 'caveman_concat_split', 0);
	my $command = sprintf('cat %s > %s',$target,$out);
	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);

	return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 'caveman_concat_split', 0);

}

sub valid_index{
	# uncoverable subroutine
	my $options = shift;
	if($options->{'process'} =~ m/^split$/){
		return file_line_count($options->{'reference'});
	}elsif($options->{'process'} =~ m/^(mstep|estep)$/){
		return file_line_count($options->{'splitList'});
	}
	return 0;
}

sub _which {
  my $prog = shift;
  my $l_bin = $Bin;
  my $path = File::Spec->catfile($l_bin, $prog);
  $path = which($prog) unless(-e $path);
  return $path;
}

1;

__END__

=head1 NAME

Sanger::CGP::Caveman::Implement - Generate SNV calls from mapped bam tumour/normal pairs

=head2 Methods

=over 4

=item file_line_count

	Sanger::CGP::Caveman::Implement::file_line_count($file);

Counts the number of lines in a text file using open and count rather than 'wc -l' which is platform dependent.

=item caveman_setup

	Sanger::CGP::Caveman::Implement::caveman_setup($options);

Runs the caveman setup process for a pair of mapped bam files.

  options - Hashref, requires the following entries:
    -outdir     : working/output directory depending on application
    -tumbam     : Path to tumour bam file
    -normbam    : Path to normal bam file
    -reference  : Path to reference index fa.fai
    -ignore     : Path to ignore file (1 based bed file of regions not to analyse)

=item caveman_split

	Sanger::CGP::Caveman::Implement::caveman_split($index,$options);

Runs the caveman split process.

	options - Hashref, requires the following entries:
    -outdir     : working/output directory depending on application

=item concat

	Sanger::CGP::Caveman::Implement::concat($options);

Concatenates the files from caveman split to results

	options - Hashref, requires the following entries:
    -out_file      : working/output directory depending on application
    -target_files  : files to concatenate into output file.

=item caveman_mstep

	Sanger::CGP::Caveman::Implement::caveman_setup($index,$options);

Runs the caveman mstep

	options - Hashref, requires the following entries:
    -outdir     : working/output directory depending on application

=item caveman_estep

	Sanger::CGP::Caveman::Implement::caveman_estep($index,$options);

Runs the caveman estep

	options - Hashref, requires the following entries:
    -outdir     : working/output directory depending on application
    -tumcn      : Path to tumour copy number file
    -normcn     : Path to normal copy number file
    -normcont   : normal contamination value (float)

=item caveman_merge_results

	Sanger::CGP::Caveman::Implement::caveman_merge_results($options);

	options - Hashref, requires the following entries:
	  -out_file   : Preliminary name of output file. Will be appended with appropriate extension by script.
    -subvcf     : Pattern match of subs vcf output files passed to script
    -snpvcf     : Pattern match of snps vcf output files passed to script
    -noanalysisbed : Pattern match of no analysis bed output files passed to script

=back
