#!/usr/bin/perl
# Copyright (c) 2014-2022
#
# Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>
#
# This file is part of cgpCaVEManWrapper.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# 1. The usage of a range of years within a copyright statement contained within
# this distribution should be interpreted as being equivalent to a list of years
# including the first and last year specified and all consecutive years between
# them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
# 2009, 2011-2012’ should be interpreted as being identical to a statement that
# reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
# statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
# identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
# 2009, 2010, 2011, 2012’.
#

BEGIN {
  use Cwd qw(abs_path);
  use File::Basename;
  push (@INC,dirname(abs_path($0)).'/../lib');
};

use strict;
use warnings FATAL => 'all';
use autodie qw(:all);

use File::Which qw(which);
use Getopt::Long;
use Try::Tiny;
use Carp;
use Pod::Usage qw(pod2usage);

{
  my $options = setup();
  merge_vcf($options->{'out'}.'.snps', $options->{'snp_vcf'});
  merge_vcf($options->{'out'}.'.subs', $options->{'sub_vcf'});
  merge_bed($options->{'out'}.'.no_analysis', $options->{'na_bed'});
}

sub merge_vcf {
  my ($path_prefix, $vcf_files) = @_;
  my $new_vcf = $path_prefix.'.vcf';
  system(qq{grep '^#' $vcf_files->[0] > $new_vcf});
  system(qq{cat @{$vcf_files} | grep -v '^#' | sort -k1,1 -k2,2n >> $new_vcf});

  my $vcf_gz = $new_vcf.'.gz';
  my $command = which('bgzip');
  $command .= sprintf ' -c %s > %s', $new_vcf, $vcf_gz;
  system($command);

  $command = which('tabix');
  $command .= sprintf ' -p vcf %s', $vcf_gz;
  system($command);
  return 1;
}

sub merge_bed {
	my ($path_prefix, $bed_files) = @_;
  my $new_bed = $path_prefix.'.bed';
  my $new_tmp_bed = $path_prefix.'.tmp.bed';
  system(qq{grep '^#' $bed_files->[0] > $new_tmp_bed});
  system(qq{cat @{$bed_files} | grep -v '^#' | sort -k1,1 -k2,2n >> $new_tmp_bed});

	#Now merge the bed file with bedtools merge so we have a much smaller output file.
	my $command = which('bedtools');
	$command .= sprintf ' merge -i %s > %s', $new_tmp_bed, $new_bed;
	system($command);

  unlink($new_tmp_bed);
  return 1;
}


sub setup {
  my %opts;
  $opts{'cmd'} = join " ", $0, @ARGV;
  GetOptions( 'h|help' => \$opts{'h'},
              'm|man' => \$opts{'m'},
              'o|out=s' => \$opts{'out'},
  ) or pod2usage(2);

  my $version = Sanger::CGP::Caveman->VERSION;

  if(defined $opts{'v'}){
    print "Version: $version\n";
    exit;
  }

  pod2usage(-verbose => 1) if(defined $opts{'h'});
  pod2usage(-verbose => 2) if(defined $opts{'m'});

  # the files from the command line pattern match
  my @bad_files;
  for my $file(@ARGV) {
    if($file =~ m/\.snps\.vcf$/) {
      push @{$opts{'snp_vcf'}}, $file;
    }elsif($file =~ m/\.muts\.vcf$/){
    	push @{$opts{'snp_vcf'}}, $file;
    }
    elsif($file =~ m/\.no_analysis\.bed$/) {
      push @{$opts{'na_bed'}}, $file;
    }
    else {
      push @bad_files, $file;
    }
  }

  die "ERROR: The following unexpected files were presented: \n\t".(join "\n\t", @bad_files)."\n"
    if(scalar @bad_files);

  return \%opts;
}

__END__

=head1 NAME

caveman_merge_results.pl - Merges provided VCF and BED files into final result files.

=head1 SYNOPSIS

caveman_merge_results.pl [options] files...

  Input files must follow the expected naming convention of:
    *.snps.vcf        - the unmerged SNP VCF files
    *.muts.vcf        - the unmerged mut VCF files
    *.no_analysis.bed - the unmerged no analysis bed files

  This matches the standard output of pindel_2_combined_vcf.pl.

  Required parameters:
    -out       -o   Output stub for final files e.g.
                    somepath/sample_vs_sample, results in:

                      VCF+index
                        somepath/sample_vs_sample.snps.vcf.gz
                        somepath/sample_vs_sample.snps.vcf.gz.tbi
                        
                      VCF+index
                        somepath/sample_vs_sample.muts.vcf.gz
                        somepath/sample_vs_sample.muts.vcf.gz.tbi

                      BED for no analysis sections
                        somepath/sample_vs_sample.no_analysis.bed

  Other:
    -help      -h   Brief help message.
    -man       -m   Full documentation.
    -version   -v   Prints the version number.

  Example:
   caveman_merge_results.pl -o someloc/tum_vs_norm in/*.snps.vcf in/*.muts.vcf in/*no_analysis.bed

