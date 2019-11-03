#!/bin/bash

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

CAVEMAN_CORE="https://github.com/cancerit/CaVEMan/archive/1.13.16.tar.gz"

get_distro () {
  if hash curl 2>/dev/null; then
    curl -sS -o $1.tar.gz -L $2
  else
    wget -nv -O $1.tar.gz $2
  fi
  mkdir -p $1
  tar --strip-components 1 -C $1 -zxf $1.tar.gz
}

if [ ! -z ${CAVE_C_REMOTE_TAR+x} ] ; then
  CAVEMAN_CORE=$CAVE_C_REMOTE_TAR
fi

if [ "$#" -ne "1" ] ; then
  echo "Please provide an installation path  such as /opt/ICGC"
  exit 0
fi

CPU=`cat /proc/cpuinfo | egrep "^processor" | wc -l`
echo "Max compilation CPUs set to $CPU"

INST_PATH=$1

# get current directory
INIT_DIR=`pwd`

# cleanup inst_path
mkdir -p $INST_PATH/bin
cd $INST_PATH
INST_PATH=`pwd`
cd $INIT_DIR

# make sure that build is self contained
unset PERL5LIB
ARCHNAME=`perl -e 'use Config; print $Config{archname};'`
PERLROOT=$INST_PATH/lib/perl5
export PERL5LIB="$PERLROOT"

#create a location to build dependencies
SETUP_DIR=$INIT_DIR/install_tmp
mkdir -p $SETUP_DIR


# log information about this system
echo '============== System information ===='
set -x
lsb_release -a
uname -a
sw_vers
system_profiler
grep MemTotal /proc/meminfo
set +x
echo

CHK=`perl -le 'eval "require $ARGV[0]" and print $ARGV[0]->VERSION' PCAP`
if [[ "x$CHK" == "x" ]] ; then
  echo "PREREQUISITE: Please install PCAP-core before proceeding: https://github.com/ICGC-TCGA-PanCancer/PCAP-core/releases"
  exit 1;
fi


if [ ! -e $PERLROOT/Sanger/CGP/CavemanPostProcessor.pm ]; then
  echo "PREREQUISITE: Please install cgpCaVEManPostProcessing before proceeding: https://github.com/cancerit/cgpCaVEManPostProcessing/releases"
  exit 1;
fi


echo -n "Building CaVEMan ..."
if [ -e $SETUP_DIR/caveman.success ]; then
  echo -n " previously installed ...";
else
  cd $SETUP_DIR
  set -xe
  if [ ! -e caveman ]; then
    get_distro "caveman" $CAVEMAN_CORE
  fi
  cd caveman
  ./setup.sh $INST_PATH
  cp scripts/mergeCavemanResults $INST_PATH/bin/.
  chmod a+x $INST_PATH/bin/mergeCavemanResults
  touch $SETUP_DIR/caveman.success
fi

#add bin path for PCAP install tests
export PATH="$INST_PATH/bin:$PATH"

cd $INIT_DIR

echo -n "Installing Perl prerequisites ..."
if ! ( perl -MExtUtils::MakeMaker -e 1 >/dev/null 2>&1); then
  echo
  echo "WARNING: Your Perl installation does not seem to include a complete set of core modules.  Attempting to cope with this, but if installation fails please make sure that at least ExtUtils::MakeMaker is installed.  For most users, the best way to do this is to use your system's package manager: apt, yum, fink, homebrew, or similar."
fi
set -x
perl $INIT_DIR/bin/cpanm -v --mirror http://cpan.metacpan.org --notest -l $INST_PATH/ --installdeps . < /dev/null
set +x

echo -n "Installing cgpCaVEManWrapper ..."
cd $INIT_DIR
perl Makefile.PL INSTALL_BASE=$INST_PATH
make
make test
make install

# cleanup all junk
rm -rf $SETUP_DIR

echo
echo
echo "Please add the following to beginning of path:"
echo "  $INST_PATH/bin"
echo "Please add the following to beginning of PERL5LIB:"
echo "  $PERLROOT"
echo
