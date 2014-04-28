#!/usr/bin/perl

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
use Const::Fast qw(const);


use PCAP::Cli;
use Sanger::CGP::Caveman::Implement;
use PCAP::Bedtools;
use PCAP::Vcftools;
use PCAP::Tabix;

{
	my $options = setup();
	Sanger::CGP::Caveman::Implement::prepare($options);
	
	my $threads = PCAP::Threaded->new($options->{'threads'});

  # register processes
	$threads->add_function('caveman_setup', \&Sanger::CGP::Caveman::Implement::caveman_setup);
	$threads->add_function('caveman_split', \&Sanger::CGP::Caveman::Implement::caveman_split);
	$threads->add_function('caveman_split_cat', \&Sanger::CGP::Caveman::Implement::concat);
	$threads->add_function('caveman_merge', \&Sanger::CGP::Caveman::Implement::caveman_merge);
	$threads->add_function('caveman_mstep', \&Sanger::CGP::Caveman::Implement::caveman_mstep);
  $threads->add_function('caveman_estep', \&Sanger::CGP::Caveman::Implement::caveman_estep);
  $threads->add_function('caveman_merge_results', \&Sanger::CGP::Caveman::Implement::caveman_merge_results);
	  
	#Start processes in correct order, according to process (DEFAULT is caveman)
	
	#caveman process flow	
	#Setup 
	$threads->run(1, 'caveman_setup', $options);
	#Split
	#count the number of chromosomes/contigs in the fasta index
	my $contig_count = Sanger::CGP::Caveman::Implement::file_line_count($options->{'reference'});
	$threads->run($contig_count, 'caveman_split', $options);
	$options->{'out_file'} = $options->{'splitList'};
	$options->{'target_files'} = File::Spec->catfile($options->{'outdir'},$options->{'splitList'}.".*");
	$threads->run(1, 'caveman_split_cat', $options);

	#Split & concatenate has succeeded in running, so now count the number of split files.
	my $split_file = File::Spec->catfile($options->{'outdir'},"splitList");
	my $split_count = Sanger::CGP::Caveman::Implement::file_line_count($split_file);	

	#Run the mstep with number of split jobs.
	$threads->run($split_count, 'caveman_mstep', $options);

	#Run the merge step
	$threads->run(1,'caveman_merge',$options);

	#Run the estep
	$threads->run($split_count, 'caveman_estep', $options);
	
	#Now we have all the results... merge all the split results files into one for each type.
	$options->{'out_file'} = File::Spec->catfile($options->{'outdir'},"."$options->{'tumour_name'}."_vs_".$options->{'normal_name'});
	$threads->run(1, 'caveman_merge_results', $options);
	
	#finally cleanup after ourselves by removing the temporary output folder, split files etc.
  &cleanup($options);  
}

sub cleanup{
	my $options = shift;
	#Get the splitFiles.* and remove them
	
	unlink($options->{'subvcf'});
	unlink($options->{'snpvcf'});
	unlink($options->{'noanalysisbed'});
  remove_tree File::Spec->catdir($options->{'outdir'}, 'results');
  remove_tree File::Spec->catdir($options->{'outdir'}, 'logs');
  remove_tree File::Spec->catdir($options->{'outdir'}, 'progress');
	return 0;
}


sub setup {
  my %opts;
  GetOptions( 	'h|help' => \$opts{'h'},
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
  ) or pod2usage(2);

  pod2usage(-message => PCAP::license, -verbose => 2) if(defined $opts{'h'});
  pod2usage(-message => PCAP::license, -verbose => 1) if(defined $opts{'m'});

  # then check for no args:
  my $defined;
  for(keys %opts) { $defined++ if(defined $opts{$_}); }
  pod2usage(-msg  => "\nERROR: Options must be defined.\n", -verbose => 2,  -output => \*STDERR) unless($defined);

  #check the reference is the fasta fai file.
  pod2usage(-msg  => "\nERROR: reference option (-r) does not appear to be a fasta index file.\n", -verbose => 2,  -output => \*STDERR) unless($opts{'reference'} =~ m/\.fai$/);

  #Check all files and dirs are readable and exist.  
  PCAP::Cli::file_for_reading('reference',$opts{'reference'});
  PCAP::Cli::file_for_reading('tumour-bam',$opts{'tumbam'});
  PCAP::Cli::file_for_reading('normal-bam',$opts{'normbam'});
  #We should also check the bam indexes exist.
  my $tumidx = File::Spec->catfile($opts{'tumbam'},".bai");
  my $normidx = File::Spec->catfile($opts{'normbam'},".bai");
  PCAP::Cli::file_for_reading('tumour-bai',$tumidx);
  PCAP::Cli::file_for_reading('normal-bai',$normidx);  
  PCAP::Cli::file_for_reading('ignore-file',$opts{'ignore'});
  PCAP::Cli::file_for_reading('tum-cn-file',$opts{'tumcn'});
  PCAP::Cli::file_for_reading('norm-cn-file',$opts{'normcn'});
  PCAP::Cli::out_dir_check('outdir', $opts{'outdir'});

	# now safe to apply defaults
	$opts{'threads'} = 1 unless(defined $opts{'threads'});

	delete $opts{'normcont'} unless(defined $opts{'normcont'});
  
	#Create the results directory in the output directory given.
	my $tmpdir = File::Spec->catdir($opts{'outdir'}, 'results');
	make_path($tmpdir) unless(-d $tmpdir); 
	#directory to store progress reports 
	my $progress = File::Spec->catdir($tmpdir, 'progress');
   make_path($progress) unless(-d $progress);
	#Directory to store run logs.
	my $logs = File::Spec->catdir($opts{'outdir'}, 'logs');
   make_path($logs) unless(-d $logs);
   
  $opts{'splitList'} = File::Spec->catfile($opts{'outdir'},"splitList");
	#vcf concat subs & snps
  $opts{'subvcf'} = File::Spec->catfile($opts{'outdir'},"*.muts.vcf");
  $opts{'snpvcf'} = File::Spec->catfile($opts{'outdir'},"*.snps.vcf");
	#bed concat no_analysis
  $opts{'noanalysisbed'} = File::Spec->catfile($opts{'outdir'},"*.no_analysis.bed");	

	return \%opts; 
}

__END__

=head1 NAME

caveman.pl - Analyse aligned bam files for SNVs via CaVEMan using a single command.

=head1 SYNOPSIS

caveman.pl [options]

  	Required parameters:
    -outdir       -o   Folder to output result to.
    -reference    -r   Path to reference genome index file *.fai
	 -tumour-bam   -tb  Path to tumour bam file
    -normal-bam   -nb  Path to normal bam file
    -ignore-file  -ig  Path to ignored regions file
    -tumour-cn    -tc  Path to tumour copy number file
    -normal-cn    -nc  Path to normal copy number file
    
   Optional parameters:
    -normal-contamination  -k   Normal contamination value (default 0.1)

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

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<caveman.pl> will attempt to run all caveman steps automatically including collation of output files.

=cut

