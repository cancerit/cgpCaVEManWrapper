cgpCaVEManWrapper
=================

cgpCaVEManWrapper provides a simplified usage implementation for the complete Cancer Genome Project processing flow of the algorithm CaVEMan.

For details of the underlying algorithm please see the [CaVEMan](http://cancerit.github.io/CaVEMan/) site.

For details of the filtering process please see the [cgpCaVEManPostProcessing](http://cancerit.github.io/cgpCaVEManPostProcessing/) site.

| Master | Dev |
|---|---|
| [![Build Status](https://travis-ci.org/cancerit/cgpCaVEManWrapper.svg?branch=master)](https://travis-ci.org/cancerit/cgpCaVEManWrapper)  | [![Build Status](https://travis-ci.org/cancerit/cgpCaVEManWrapper.svg?branch=dev)](https://travis-ci.org/cancerit/cgpCaVEManWrapper) |

---

### Dependencies/Install
Please install the following first:

* [PCAP-core](http://github.com/ICGC-TCGA-PanCancer/PCAP-core/releases)
* [cgpCaVEManPostProcessing](http://github.com/cancerit/cgpCaVEManPostProcessing/releases)

Please see these for any child dependencies.

Once complete please run:

./setup.sh /some/install/location

This will automatically get the appropriate version of the core [CaVEMan](http://cancerit.github.io/CaVEMan/) algorithm.

---

## Creating a release
#### Preparation
* Commit/push all relevant changes.
* Pull a clean version of the repo and use this for the following steps.

#### Cutting the release
1. Update `perl/lib/Sanger/CGP/Caveman.pm` to the correct version.
2. Run `./prerelease.sh`
3. Check all tests and coverage reports are acceptable.
4. Commit the updated docs tree and updated module/version.
5. Push commits.
6. Use the GitHub tools to draft a release.

LICENCE
=======
Copyright (c) 2014-2017 Genome Research Ltd.

Author: Cancer Genome Project <cgpit@sanger.ac.uk>

This file is part of cgpCaVEManWrapper.

cgpCaVEManWrapper is free software: you can redistribute it and/or modify it under
the terms of the GNU Affero General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option) any
later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

1. The usage of a range of years within a copyright statement contained within
this distribution should be interpreted as being equivalent to a list of years
including the first and last year specified and all consecutive years between
them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
2009, 2011-2012’ should be interpreted as being identical to a statement that
reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
2009, 2010, 2011, 2012’."
