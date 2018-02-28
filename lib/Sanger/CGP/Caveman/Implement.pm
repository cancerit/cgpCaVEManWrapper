package Sanger::CGP::Caveman::Implement;

##########LICENCE##########
#  Copyright (c) 2014-2018 Genome Research Ltd.
#
#  Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
#
#  This file is part of cgpCaVEManWrapper.
#
#  cgpCaVEManWrapper is free software: you can redistribute it and/or modify it under
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
use Capture::Tiny qw(capture);
use List::Util qw(first);
use File::Basename;

use Sanger::CGP::Caveman;
our $VERSION = Sanger::CGP::Caveman->VERSION;

use PCAP::Threaded;
use PCAP::Bam;

const my $CAVEMAN_SETUP => q{ setup -t %s -n %s -r %s -g %s -l %s -f %s -c %s -a %s};
const my $CAVEMAN_SPLIT => q{ split -i %d -f %s -e %d};
const my $CAVEMAN_MSTEP => q{ mstep -i %d -f %s};
const my $CAVEMAN_MERGE => q{ merge -c %s -p %s -f %s};
const my $CAVEMAN_ESTEP => q{ estep -i %d -k %f -g %s -o %s -v %s -w %s -f %s -l %s -r %s};
const my $CAVEMAN_ESTEP_MUT_PRIOR_EXT => q{ -c %s};
const my $CAVEMAN_ESTEP_SNP_PRIOR_EXT => q{ -d %s};
const my $CAVEMAN_ESTEP_NPLATFORM_EXT => q{ -P %s};
const my $CAVEMAN_ESTEP_TPLATFORM_EXT => q{ -T %s};
const my $CAVEMAN_FLAG => q{ -i %s -o %s -s %s -m %s -n %s -b %s -g %s -umv %s -ref %s -t %s};
const my $MERGE_CAVEMAN_RESULTS => q{ mergeCavemanResults -s %s -o %s -f %s};
const my $CAVEMAN_VCF_IDS => q{ -i %s -o %s};
const my $CAVEMAN_MUT_PROB_CUTOFF => q{ -p %f};
const my $CAVEMAN_SNP_PROB_CUTOFF => q{ -q %f};
const my $CAVEMAN_DEBUG_MODE => q{ -s};
const my $CAVEMAN_VCF_SPLIT => q{ -i %s -o %s -s -l %d};
const my $CAVEMAN_VCF_FLAGGED_CONCAT => q{vcf-concat %s | vcf-sort > %s};
const my $FILE_COUNT => q{ls -1 %s | wc -l};

const my $FLAG_SCRIPT => q{cgpFlagCaVEMan.pl};
const my $IDS_SCRIPT => q{cgpAppendIdsToVcf.pl};
const my $VCF_SPLIT_SCRIPT => q{cgpVCFSplit.pl};

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

	if(exists($options->{'normcn'}) && defined($options->{'normcn'}))
	{
		$command .= " -j ".$options->{'normcn'};
	}
	if(exists($options->{'tumcn'}) && defined($options->{'tumcn'})){
		$command .= " -e ".$options->{'tumcn'};
	}
	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);

	return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
}

sub caveman_split {
	# uncoverable subroutine
	my ($index,$options) = @_;

	return 1 if(exists $options->{'index'} && $index != $options->{'index'});
	my $tmp = $options->{'tmp'};
  my $config = $options->{'cave_cfg'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);

	my $command = _which('caveman') || die "Unable to find 'caveman' in path";

	$command .= sprintf($CAVEMAN_SPLIT,
                      $options->{'valid_fai_idx'}->[$index-1], # only process the contigs we care about
                      $config,
                      $options->{'read-count'});

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);
  	return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
}

sub caveman_merge {
	# uncoverable subroutine
	my $options = shift;

	my $command = _which('caveman') || die "Unable to find 'caveman' in path";
	my $tmp = $options->{'tmp'};
	my $config = $options->{'cave_cfg'};
	my $prob_arr = $options->{'cave_parr'};
	my $cov_arr = $options->{'cave_carr'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

	$command .= sprintf($CAVEMAN_MERGE, $cov_arr, $prob_arr,$config);

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);

	return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
}

sub caveman_mstep {
	# uncoverable subroutine
	my ($index_in,$options) = @_;

	# first handle the easy bit, skip if limit not set
	return 1 if(exists $options->{'index'} && $index_in != $options->{'index'});

	my @indicies = limited_xstep_indicies($options, $index_in);
  my $config = $options->{'cave_cfg'};
	my $tmp = $options->{'tmp'};
	for my $index(@indicies) {
    next if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);

    my $command = _which('caveman') || die "Unable to find 'caveman' in path";

    $command .= sprintf($CAVEMAN_MSTEP,
                    $index,
                    $config);

    PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);
    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
  }
  return 1;
}

sub caveman_estep {
	# uncoverable subroutine
	my ($index_in,$options) = @_;

	# first handle the easy bit, skip if limit not set
	return 1 if(exists $options->{'index'} && $index_in != $options->{'index'});

	my @indicies = limited_xstep_indicies($options, $index_in);
  my $config = $options->{'cave_cfg'};
	my $tmp = $options->{'tmp'};
	my $prob_arr = $options->{'cave_parr'};
	my $cov_arr = $options->{'cave_carr'};
	my $normprot = $options->{'normprot'};
	my $tumprot = $options->{'tumprot'};

	for my $index(@indicies) {
    next if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);

    my $command = _which('caveman') || die "Unable to find 'caveman' in path";

    $command .= sprintf($CAVEMAN_ESTEP,
                    $index,
                    $options->{'normcont'},
                    $cov_arr,
                    $prob_arr,
                    $options->{'species-assembly'},
                    q{'}.$options->{'species'}.q{'},
                    $config,
                    $normprot,
                    $tumprot);

    if(exists($options->{'normdefcn'}) && defined($options->{'normdefcn'})){ #Add default normal cn
      $command .= ' -n '.$options->{'normdefcn'};
    }

    if(exists($options->{'tumdefcn'}) && defined($options->{'tumdefcn'})){ #Add default tumour cn
      $command .= ' -t '.$options->{'tumdefcn'};
    }

    if(exists($options->{'priorMut'}) && defined($options->{'priorMut'})){
      $command .= sprintf($CAVEMAN_ESTEP_MUT_PRIOR_EXT,$options->{'priorMut'});
    }

    if(exists($options->{'priorSnp'}) && defined($options->{'priorSnp'})){
      $command .= sprintf($CAVEMAN_ESTEP_SNP_PRIOR_EXT,$options->{'priorSnp'});
    }

    #Check for platform overrides.
    if(exists($options->{'nplat'}) && defined($options->{'nplat'})){
      $command .= sprintf($CAVEMAN_ESTEP_NPLATFORM_EXT,$options->{'nplat'});
    }

    if(exists($options->{'tplat'}) && defined($options->{'tplat'})){
      $command .= sprintf($CAVEMAN_ESTEP_TPLATFORM_EXT,$options->{'tplat'});
    }

    if(exists($options->{'mpc'}) && defined($options->{'mpc'})){
      $command .= sprintf($CAVEMAN_MUT_PROB_CUTOFF,$options->{'mpc'});
    }

    if(exists($options->{'spc'}) && defined($options->{'spc'})){
      $command .= sprintf($CAVEMAN_SNP_PROB_CUTOFF,$options->{'spc'});
    }

    if(exists($options->{'debug_cave'}) && defined($options->{'debug_cave'})){
      $command .= sprintf($CAVEMAN_DEBUG_MODE);
    }

    PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);
    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
  }
  return 1;
}

sub caveman_merge_results {
	# uncoverable subroutine
	my $options = shift;
	my $tmp = $options->{'tmp'};
	my $out = $options->{'out_file'};
	my $splitList = File::Spec->catfile($tmp, 'splitList');

	my $command = sprintf($MERGE_CAVEMAN_RESULTS,$splitList,$out.".muts.vcf",$options->{'subvcf'});
	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0)
								unless (PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 'merge_muts', 0));
	PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 'merge_muts', 0);

	$command = sprintf($MERGE_CAVEMAN_RESULTS,$splitList,$out.".snps.vcf",$options->{'snpvcf'});
	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0)
								unless (PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 'merge_snps', 0));
	PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 'merge_snps', 0);

  unless (PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 'merge_no_analysis', 0)) {
  	$command = sprintf($MERGE_CAVEMAN_RESULTS,$splitList,$out.".no_analysis.bed",$options->{'noanalysisbed'});
  	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);
    extend_no_analysis($options, $out.'.no_analysis.bed');
  }


	PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 'merge_no_analysis', 0);

	return 1;
}

sub caveman_add_vcf_ids {
	# uncoverable subroutine
	my ($options, $snps_or_muts) = @_;
	my $tmp = $options->{'tmp'};
	my $raw = $options->{'raw_file'};
	my $ids = $options->{'ids_file'};

	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $snps_or_muts);
	my $script = _which($IDS_SCRIPT) ||  die "Unable to find '$IDS_SCRIPT' in path";
	my $command = $^X.' '.$script;
	$command .= sprintf($CAVEMAN_VCF_IDS,
														$raw,
														$ids);

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $snps_or_muts);
	return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $snps_or_muts);
}

sub caveman_split_vcf {
  # uncoverable subroutine
	my $options = shift;
  my $tmp = $options->{'tmp'};
  my $infile = $options->{'for_split'};
  my $outstub = $options->{'split_out'};
  my $split_lines = $options->{'split_lines'};

  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'),0);

  my $script = _which($VCF_SPLIT_SCRIPT) ||  die "Unable to find '$VCF_SPLIT_SCRIPT' in path";
	my $command = $^X.' '.$script;
	$command .= sprintf($CAVEMAN_VCF_SPLIT,
														$infile,
														$outstub,
                            $split_lines);
  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);
	return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
}

sub count_files {
  # uncoverable subroutine
	my ($options,$match) = @_;
  my $tmp = $options->{'tmp'};
  my $command = sprintf($FILE_COUNT,$match);
  my ($stdout, $stderr, $exit) = capture {
    system($command);
  };
  die "ERROR: ($stderr) Encountered counting split files for flagging. Searching $match" unless($exit==0);
  chomp($stdout);
  $stdout =~ s/\s+//g;
  return $stdout;
}

sub caveman_flag {
  # uncoverable subroutine
  my ($index_in,$options) = @_;

  # first handle the easy bit, skip if limit not set
  return 1 if(exists $options->{'index'} && $index_in != $options->{'index'});

  my @indicies = limited_flag_indicies($options, $index_in);

	my $tmp = $options->{'tmp'};
	my $tumbam = $options->{'tumbam'};
	my $normbam = $options->{'normbam'};
	my $ref = $options->{'reference'};

  for my $index(@indicies) {
    next if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);
    my $script = _which($FLAG_SCRIPT) || die "Unable to find '$FLAG_SCRIPT' in path";
  	my $flag = $^X.' '.$script;
  	$flag .= sprintf($CAVEMAN_FLAG,
              $options->{'split_out'}.".$index",
              $options->{'flagged'}.".$index",
              q{'}.$options->{'species'}.q{'},
              $tumbam,
              $normbam,
              $options->{'flag-bed'},
              $options->{'germindel'},
              $options->{'unmatchedvcf'},
              $ref,
              $options->{'seqType'}
              );

    $flag .= ' -c '.$options->{'flagConfig'} if(defined $options->{'flagConfig'});
  	$flag .= ' -v '.$options->{'flagToVcfConfig'} if(defined $options->{'flagToVcfConfig'});
  	$flag .= ' -p '.$options->{'apid'} if(defined $options->{'apid'});
    $flag .= ' -ab '.$options->{'annot-bed'} if(defined $options->{'annot-bed'});

    PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $flag, $index);
    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);

  }

  return 1;
}

sub concat_flagged {
  # uncoverable subroutine
	my $options = shift;
	my $tmp = $options->{'tmp'};

  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'),0);
  my $command = sprintf($CAVEMAN_VCF_FLAGGED_CONCAT,
                  $options->{'flagged'}.".*",
                  $options->{'flagged'});

  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);
  return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
}

sub zip_flagged{
  # uncoverable subroutine
	my $options = shift;
	my $tmp = $options->{'tmp'};
  my $flagged = $options->{'flagged'};
  my $vcf_gz = $flagged.'.gz';
  my $bgzip = _which('bgzip');
  $bgzip .= sprintf ' -c %s > %s', $flagged, $vcf_gz;

  my $tabix = _which('tabix');
  $tabix .= sprintf ' -p vcf %s', $vcf_gz;

  my @commands = ($bgzip, $tabix);

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), \@commands, 0);
	return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
}

sub pre_cleanup_zip {
  # uncoverable subroutine
	my $options = shift;
	my $tmp = $options->{'tmp'};

	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

	my $vcf_muts_gz = $options->{'ids_muts_file'}.'.gz';
	my $bgzip_muts = _which('bgzip');
  $bgzip_muts .= sprintf ' -c %s > %s', $options->{'ids_muts_file'}, $vcf_muts_gz;

  my $tabix_muts = _which('tabix');
  $tabix_muts .= sprintf ' -p vcf %s', $vcf_muts_gz;

  my $vcf_snps_gz = $options->{'ids_snps_file'}.'.gz';
  my $bgzip_snps = _which('bgzip');
  $bgzip_snps .= sprintf ' -c %s > %s', $options->{'ids_snps_file'}, $vcf_snps_gz;

  my $tabix_snps = _which('tabix');
  $tabix_snps .= sprintf ' -p vcf %s', $vcf_snps_gz;
  my @commands = ($bgzip_muts, $tabix_muts, $bgzip_snps, $tabix_snps);

  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), \@commands, 0);
	return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
}

sub limited_flag_indicies {
  my ($options, $index_in) = @_;
  return limited_indices($options, $index_in, $options->{'vcf_split_count'});
}

sub limited_xstep_indicies {
  my ($options, $index_in) = @_;
	my $split_count = file_line_count($options->{'splitList'});
	return limited_indices($options, $index_in, $split_count);
}

sub limited_indices {
	my ($options, $index_in, $count) = @_;
  my @indicies;
  if(exists $options->{'limit'}) {
    # main script checks index is not greater than limit or < 1
	  my $base = $index_in;
	  while($base <= $count) {
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
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);
	my $command = sprintf('cat %s > %s',$target,$out);
	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);

	return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);

}

sub load_exclude {
  my $options = shift;
  my @exclude_patt;
  if(exists $options->{'exclude'}) {
    my @exclude = split /,/, $options->{'exclude'};
    for my $ex(@exclude) {
      $ex =~ s/%/.+/;
      push @exclude_patt, $ex;
    }
  }
  return @exclude_patt;
}

sub valid_seq_indexes {
  my $options = shift;

  my @exclude_patt = load_exclude($options);

  my @good_sq_idx;
  open my $FAI_IN, '<', $options->{'reference'};
  while(<$FAI_IN>) {
    my $sq = (split /\t/, $_)[0];
    # if doesn't match keep
    push @good_sq_idx, $. unless(first { $sq =~ m/^$_$/ } @exclude_patt);
  }
  close $FAI_IN;

  return \@good_sq_idx;
}

sub extend_no_analysis {
  my ($options, $no_analysis) = @_;
  my @exclude_patt = load_exclude($options);
  return if(@exclude_patt == 0);

  open my $na_fh, '>>', $no_analysis;
  open my $FAI_IN, '<', $options->{'reference'};
  while(<$FAI_IN>) {
    my ($sq, $len) = (split /\t/, $_)[0..1];
    # if matches then print
    printf $na_fh "%s\t0\t%d\n", $sq, $len if(first { $sq =~ m/^$_$/ } @exclude_patt);
  }
  close $FAI_IN;
  close $na_fh
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

=item caveman_add_vcf_ids

	Sanger::CGP::Caveman::implement::caveman_add_vcf_ids

Appends IDS to a vcf file

	options - Hashref, requires the following entries:
		-raw_file  : Path to the existing vcf file
		-ids_file  : Path to generate the vcf file with IDs


=item caveman_flag

	Sanger::CGP::Caveman::Implement::caveman_flag(,$options);

Runs flagging (post processing) over CaVEMan results

	options - Hashref, requires the following entries:
		-tumbam        : Path to tumour bam file
    -normbam       : Path to normal bam file
    -for_flagging  : Path to file to be flagged
	  -flagged       : Path to output flagged file
    -reference     : Path to reference index fa.fai
		-flag-bed      : Directory containing flagging related bed files
		-germindel     : path to germline indel file
		-unmatchedvcf  : Directory containing unmatched VCF files

=item caveman_merge_results

	Sanger::CGP::Caveman::Implement::caveman_merge_results($options);

	options - Hashref, requires the following entries:
	  -out_file   : Preliminary name of output file. Will be appended with appropriate extension by script.
    -subvcf     : Pattern match of subs vcf output files passed to script
    -snpvcf     : Pattern match of snps vcf output files passed to script
    -noanalysisbed : Pattern match of no analysis bed output files passed to script

=back
