#!/usr/bin/perl

use ExtUtils::MakeMaker;

WriteMakefile(
  NAME          => 'cgpCaVEManWrapper',
  VERSION_FROM  => 'lib/Sanger/CGP/Caveman.pm',
  EXE_FILES     => [qw( bin/caveman.pl)],
  PREREQ_PM     => {
                     'Const::Fast' => 0.014,
                     'Try::Tiny' => 0.19,
                     'File::Which' => 0.05,
                     'List::Util' => 1.38,
                     'Test::Fatal' => 0.013,
                     'Devel::Cover' => 1.09,
                     'Pod::Coverage' => 0.23,
                   }
);
