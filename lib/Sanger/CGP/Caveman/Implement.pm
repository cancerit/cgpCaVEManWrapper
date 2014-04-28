package Sanger::CGP::Caveman::Implement;

use strict;
use warnings FATAL => 'all';
use autodie qw(:all);
use Const::Fast qw(const);

use Sanger::CGP::Caveman;
our $VERSION = Sanger::CGP::Caveman->VERSION;

use PCAP::Threaded;
use PCAP::Bam;

const my $CAVEMAN_SETUP => q{ setup -t %s -n %s -r %s -g %s};
const my $CAVEMAN_SPLIT => q{ split -i %d};
const my $CAVEMAN_MSTEP => q{ mstep -i %d};
const my $CAVEMAN_MERGE => q{ merge};
const my $CAVEMAN_ESTEP => q{ estep -i %d -e %s -j %s -k %f};
const my $MERGE_CAVEMAN_RESULTS => q{ caveman_merge_results.pl -o %s %s};
const my $CONCAT_CAVEMAN_SPLIT => q{ mergeCavemanResults -o %s %s};

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
	my $tmp = $options->{'outdir'};
	my $ref = $options->{'reference'};
	my $tumbam = $options->{'tumbam'};
	my $normbam = $options->{'normbam'};
	my $ignore = $options->{'ignore'};
	
	my $command = which('caveman') || die "Unable to find 'caveman' in path";
	
	$command .= sprintf($CAVEMAN_SETUP, 
								$options->{'tumbam'},
								$options->{'normbam'},
								$options->{'reference'},
								$options->{'ignore'});
								
	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);	
	
	return 1;
}

sub caveman_split {
	# uncoverable subroutine
	my ($index,$options) = @_;
	
	return 1 if(exists $options->{'index'} && $index != $options->{'index'});
	my $tmp = $options->{'outdir'};
	
	return if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);
	
	my $command = which('caveman') || die "Unable to find 'caveman' in path";
	$command .= sprintf($CAVEMAN_SPLIT,
									$index);
									
	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);
  	return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 'caveman_split', $index);	
}

sub caveman_merge{
	# uncoverable subroutine
	my $options = shift;
	
	my $command = which('caveman') || die "Unable to find 'caveman' in path";
	my $tmp = $options->{'outdir'};
	$command .= sprintf($CAVEMAN_MERGE);	
	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);	
	
	return 1;
}

sub caveman_mstep{
	# uncoverable subroutine
	my ($index,$options) = @_;
	return 1 if(exists $options->{'index'} && $index != $options->{'index'});
	my $tmp = $options->{'outdir'};
	return if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);
	
	my $command = which('caveman') || die "Unable to find 'caveman' in path";
	
	$command .= sprintf($CAVEMAN_MSTEP,
									$index);
	
	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);
  	return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 'caveman_mstep', $index);		
}

sub caveman_estep{
	# uncoverable subroutine
	my ($index,$options) = @_;
	
	return 1 if(exists $options->{'index'} && $index != $options->{'index'});
	my $tmp = $options->{'outdir'};
	return if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);
	
	my $command = which('caveman') || die "Unable to find 'caveman' in path";

	$command .= sprintf($CAVEMAN_ESTEP,
									$index,
									$options->{'normcn'},
									$options->{'tumcn'},
									$options->{'normcont'});
									
	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);
  	return PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 'caveman_estep', $index);		
}

sub caveman_merge_results {
	# uncoverable subroutine
	my $options = shift;
	my $out = $options->{'out_file'};
	my $target = $options->{'subvcf'}.' '.$options->{'snpvcf'}.' '.$options->{'noanalysisbed'};
	my $command = sprintf($MERGE_CAVEMAN_RESULTS,$out,$target);
	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);	
	
	return 1;
	
}

sub concat {
	# uncoverable subroutine
	my $options = shift;
	my $out = $options->{'out_file'};
	my $target = $options->{'target_files'};
	my $command = sprintf('cat %s > %s',$target,$out);
	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);	
	
	return 1;
	
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