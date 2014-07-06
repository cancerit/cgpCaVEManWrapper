#!/usr/bin/perl

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


BEGIN {
  use Cwd qw(abs_path);
  use File::Basename;
  push (@INC,dirname(abs_path($0)).'/../lib');
};

use strict;
use warnings FATAL => 'all';
use autodie qw(:all);

use File::Path qw(remove_tree make_path);
use File::Spec;
use Getopt::Long;
use Pod::Usage qw(pod2usage);
use List::Util qw(first);
use Const::Fast qw(const);

use PCAP::Cli;
use Sanger::CGP::Caveman::Implement;

const my @VALID_PROCESS => qw(setup split split_concat mstep merge estep merge_results);
const my $CAVEMAN_CONFIG => 'caveman.cfg.ini';
const my $CAVEMAN_ALG_BEAN => 'alg_bean';
const my $CAVEMAN_PROB_ARR => 'prob_arr';
const my $CAVEMAN_COV_ARR => 'cov_arr';

my %index_max = ( 'setup'  => 1,
									'split'  => -1,
									'mstep'  => -1,
									'merge'  => 1,
									'estep'  => -1,
									'merge_results'   => 1,
									'split_concat' => -1);

{
	my $options = setup();
	Sanger::CGP::Caveman::Implement::prepare($options);

	my $threads = PCAP::Threaded->new($options->{'threads'});
	&PCAP::Threaded::disable_out_err if(exists $options->{'index'});

  # register processes
	$threads->add_function('caveman_split', \&Sanger::CGP::Caveman::Implement::caveman_split);
	$threads->add_function('caveman_mstep', \&Sanger::CGP::Caveman::Implement::caveman_mstep);
  $threads->add_function('caveman_estep', \&Sanger::CGP::Caveman::Implement::caveman_estep);

	#Start processes in correct order, according to process (DEFAULT is caveman)

	#caveman process flow
	#Setup
	Sanger::CGP::Caveman::Implement::caveman_setup($options) if(!exists $options->{'process'} || $options->{'process'} eq 'setup');
	#Split
	#count the number of chromosomes/contigs in the fasta index

	if(!exists $options->{'process'} || $options->{'process'} eq 'split'){
		$options->{'out_file'} = $options->{'splitList'};
		my $contig_count = Sanger::CGP::Caveman::Implement::file_line_count($options->{'reference'});
		$threads->run($contig_count, 'caveman_split', $options);
	}

  if(!exists $options->{'process'} || $options->{'process'} eq 'split_concat'){
    $options->{'out_file'} = $options->{'splitList'};
    $options->{'target_files'} = $options->{'splitList'}.".*";
	  Sanger::CGP::Caveman::Implement::concat($options);
	}

	my $split_count = Sanger::CGP::Caveman::Implement::file_line_count($options->{'splitList'}) if(!exists $options->{'process'} || first { $options->{'process'} eq $_ } ('mstep', 'estep'));
	#Split & concatenate has succeeded in running, so now count the number of split files.
	if(!exists $options->{'process'} || $options->{'process'} eq 'mstep'){
		#Run the mstep with number of split jobs.
		$threads->run($split_count, 'caveman_mstep', $options);
	}

	#Run the merge step
	Sanger::CGP::Caveman::Implement::caveman_merge($options) if(!exists $options->{'process'} || $options->{'process'} eq 'merge');

	#Run the estep
	if(!exists $options->{'process'} || $options->{'process'} eq 'estep'){
		$threads->run($split_count, 'caveman_estep', $options);
	}

	#Now we have all the results... merge all the split results files into one for each type.
	if(!exists $options->{'process'} || $options->{'process'} eq 'merge_results'){
	  $options->{'out_file'} = File::Spec->catfile($options->{'outdir'},$options->{'tumour_name'}."_vs_".$options->{'normal_name'});
		Sanger::CGP::Caveman::Implement::caveman_merge_results($options);
	  #finally cleanup after ourselves by removing the temporary output folder, split files etc.
  	cleanup($options);
  }
}

sub cleanup{
	my $options = shift;
  #Make a directory for kept results

  #Move cov array
  #Move prob array
  #Move alg bean
  #move concatenated muts
  #move concatenated snps
  #move concatenated no analysis
  #move config file

  remove_tree (File::Spec->catdir($options->{'tmp'}, 'results'));
  remove_tree (File::Spec->catdir($options->{'tmp'}, 'logs'));
  remove_tree (File::Spec->catdir($options->{'tmp'}, 'progress'));
	return 0;
}


sub setup {
  my %opts;
  GetOptions(
  				'h|help' => \$opts{'h'},
					'm|man' => \$opts{'m'},
					'r|reference=s' => \$opts{'reference'},
					'o|outdir=s' => \$opts{'outdir'},
					'tb|tumour-bam=s' => \$opts{'tumbam'},
					'nb|normal-bam=s' => \$opts{'normbam'},
					'ig|ignore-file=s' => \$opts{'ignore'},
					'tc|tumour-cn=s' => \$opts{'tumcn'},
					'nc|normal-cn=s' => \$opts{'normcn'},
					't|threads=i' => \$opts{'threads'},
					'k|normal-contamination=f' => \$opts{'normcont'},
					's|species=s' => \$opts{'species'},
					'sa|species-assembly=s' => \$opts{'species-assembly'},
					'p|process=s' => \$opts{'process'},
					'i|index=i' => \$opts{'index'},
					'l|limit=i' => \$opts{'limit'},
  ) or pod2usage(2);

  pod2usage(-message => PCAP::license, -verbose => 2) if(defined $opts{'h'});
  pod2usage(-message => PCAP::license, -verbose => 1) if(defined $opts{'m'});

  # then check for no args:
  my $defined;
  for(keys %opts) { $defined++ if(defined $opts{$_}); }
  pod2usage(-msg  => "\nERROR: Options must be defined.\n", -verbose => 2,  -output => \*STDERR) unless($defined);

	pod2usage(-msg  => "\nERROR: Options must be defined.\n", -verbose => 2,  -output => \*STDERR) unless(defined($opts{'species'}) && defined($opts{'species-assembly'}));

  #check the reference is the fasta fai file.
  pod2usage(-msg  => "\nERROR: reference option (-r) does not appear to be a fasta index file.\n", -verbose => 2,  -output => \*STDERR) unless($opts{'reference'} =~ m/\.fai$/);

  #Check all files and dirs are readable and exist.
  PCAP::Cli::file_for_reading('reference',$opts{'reference'});
  PCAP::Cli::file_for_reading('tumour-bam',$opts{'tumbam'});
  PCAP::Cli::file_for_reading('normal-bam',$opts{'normbam'});
  #We should also check the bam indexes exist.
  my $tumidx = $opts{'tumbam'}.".bai";
  my $normidx = $opts{'normbam'}.".bai";
  PCAP::Cli::file_for_reading('tumour-bai',$tumidx);
  PCAP::Cli::file_for_reading('normal-bai',$normidx);
  PCAP::Cli::file_for_reading('ignore-file',$opts{'ignore'});
  PCAP::Cli::file_for_reading('tum-cn-file',$opts{'tumcn'});
  PCAP::Cli::file_for_reading('norm-cn-file',$opts{'normcn'});
  PCAP::Cli::out_dir_check('outdir', $opts{'outdir'});

  delete $opts{'process'} unless(defined $opts{'process'});
  delete $opts{'index'} unless(defined $opts{'index'});
  delete $opts{'limit'} unless(defined $opts{'limit'});

	if(exists $opts{'process'}) {
    PCAP::Cli::valid_process('process', $opts{'process'}, \@VALID_PROCESS);
    if(exists $opts{'index'}) {
      my $max = $index_max{$opts{'process'}};
      if($max==-1){
        if(exists $opts{'limit'}) {
          $max = $opts{'limit'};
        }
        else {
      	  $max = Sanger::CGP::Caveman::Implement::valid_index(\%opts);
      	}
      }

      die "ERROR: based on reference and exclude option index must be between 1 and $max\n" if($opts{'index'} < 1 || $opts{'index'} > $max);
      PCAP::Cli::opt_requires_opts('index', \%opts, ['process']);

      die "No max has been defined for this process type\n" if($max == 0);

      PCAP::Cli::valid_index_by_factor('index', $opts{'index'}, $max, 1);
    }
  }
  elsif(exists $opts{'index'}) {
    die "ERROR: -index cannot be defined without -process\n";
  }


	# now safe to apply defaults
	$opts{'threads'} = 1 unless(defined $opts{'threads'});

	$opts{'normcont'} = 0.1 unless(defined $opts{'normcont'});

	#Create the results directory in the output directory given.
	my $tmpdir = File::Spec->catdir($opts{'outdir'}, 'tmpCaveman');
	make_path($tmpdir) unless(-d $tmpdir);
	$opts{'tmp'} = $tmpdir;
	my $resultsdir = File::Spec->catdir($opts{'tmp'}, 'results');
	make_path($resultsdir) unless(-d $resultsdir);
	#directory to store progress reports
	my $progress = File::Spec->catdir($opts{'tmp'}, 'tmp');
   make_path($progress) unless(-d $progress);
	#Directory to store run logs.
	my $logs = File::Spec->catdir($opts{'outdir'}, 'tmp');
  make_path($logs) unless(-d $logs);
  $opts{'logs'} = $logs;
  my $config_file = File::Spec->catfile($opts{'tmp'},$CAVEMAN_CONFIG);
  $opts{'cave_cfg'} = $config_file;
  my $alg_bean = File::Spec->catfile($opts{'tmp'},$CAVEMAN_ALG_BEAN);
  $opts{'cave_alg'} = $alg_bean;
  my $prob_arr = File::Spec->catfile($opts{'tmp'},$CAVEMAN_PROB_ARR);
  $opts{'cave_parr'} = $prob_arr;
  my $cov_arr = File::Spec->catfile($opts{'tmp'},$CAVEMAN_COV_ARR);
  $opts{'cave_carr'} = $cov_arr;

  $opts{'splitList'} = File::Spec->catfile($opts{'tmp'},"splitList");
	#vcf concat subs & snps
  $opts{'subvcf'} = File::Spec->catfile($opts{'tmp'},"results/*/*.muts.vcf");
  $opts{'snpvcf'} = File::Spec->catfile($opts{'tmp'},"results/*/*.snps.vcf");
	#bed concat no_analysis
  $opts{'noanalysisbed'} = File::Spec->catfile($opts{'tmp'},"results/*/*.no_analysis.bed");


	return \%opts;
}

__END__

=head1 NAME

caveman.pl - Analyse aligned bam files for SNVs via CaVEMan using a single command.

=head1 SYNOPSIS

caveman.pl [options]

  	Required parameters:
    -outdir            -o   Folder to output result to.
    -reference         -r   Path to reference genome index file *.fai
	 -tumour-bam         -tb  Path to tumour bam file
    -normal-bam        -nb  Path to normal bam file
    -ignore-file       -ig  Path to ignored regions file
    -tumour-cn         -tc  Path to tumour copy number file
    -normal-cn         -nc  Path to normal copy number file
    -species           -s   Species name for (output in VCF)
    -species-assembly  -sa  Species assembly for (output in VCF)

   Optional parameters:
    -normal-contamination  -k   Normal contamination value (default 0.1)
    -threads               -t   Number of threads allowed on this machine (default 1)
    -limit                 -l   Limit the number of jobs required for m/estep (default undef)

   Targeted processing (further detail under OPTIONS):
    -process   -p   Only process this step then exit, optionally set -index
    -index     -i   Optionally restrict '-p' to single job

	Other:
    -help     -h   Brief help message.
    -man      -m   Full documentation.

=head1 OPTIONS

=over 8

=item B<-outdir>

Directory to write output to.  During processing a temp folder will be generated in this area,
should the process fail B<only delete this if> you are unable to resume the process.

Final output files are: muts.vcf.gz, snps.vcf.gz, no_analysis.bed.gz no_analysis.bed.gz.tbi

=item B<-reference>

Path to genome.fa.fai file and associated .fa file.

=item B<-tumour-bam>

Path to mapped, indexed, duplicate marked/removed tumour bam file.

=item B<-normal-bam>

Path to mapped, indexed, duplicate marked/removed normal bam file.

=item B<-ignore-file>

Path to ignore file. 1-based first coordinate bed style format of regions for caveman not to analyse.

=item B<-tumour-cn>

Path to tumour copy number file (1-based first coordinate bed style format). All analysed bases must have a CN assigned.

=item B<-normal-cn>

Path to normal copy number file (1-based first coordinate bed style format). All analysed bases must have a CN assigned.

=item B<-species>

Species name for (output in VCF) e.g HUMAN

=item B<-species-assembly>

Species assembly for (output in VCF) e.g. 37

=item B<-process>

Used to restrict to a single process. Valid processes are [setup|split|split_concat|mstep|merge|estep|merge_results]. Use in conjunction with -index to restrict to a single job in a process.

=item B<-index>

Use in conjunction with -process to restrict to a single job index in a process.

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<caveman.pl> will attempt to run all caveman steps automatically including collation of output files.

=cut

