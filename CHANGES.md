# CHANGES

## 1.14.1

* Use CaVEMan core [1.13.16](https://github.com/cancerit/CaVEMan/releases/tag/1.13.16)

## 1.14.0

* Created Docker file to build a containder for the wrapper. 
* Added quay.io badge for the image to the README.

## 1.13.14

* Adds `-noclean` option to allow resumption following `-noflag` execution.

## 1.13.13

* Build with CaVEMan [1.13.15](https://github.com/cancerit/CaVEMan/releases/tag/1.13.15)

## 1.13.12

* Minor errors in CaVEMan 1.13.13 force update.

## 1.13.11

* Add AMPLICON and TARGETED to accepted protocol list. Change RNA-seq, but also keep RNA.
* https://www.ebi.ac.uk/ena/submit/reads-library-strategy
* In line with Update of CaVEMan:
* Build with CaVEMan [1.13.13](https://github.com/cancerit/CaVEMan/releases/tag/1.13.13)

## 1.13.10

* Fix command line parsing of short option names, np/NP and tp/TP only ever set uppercase versions.

## 1.13.9

* Build with CaVEMan [1.13.12](https://github.com/cancerit/CaVEMan/releases/tag/1.13.12)

## 1.13.8

* Added libdb-dev and libgd-dev to travis.yml

## 1.13.7

* Build with CaVEMan [1.13.11](https://github.com/cancerit/CaVEMan/releases/tag/1.13.11)

## 1.13.6

* Build with CaVEMan [1.13.10](https://github.com/cancerit/CaVEMan/releases/tag/1.13.10)
* **See CaVEMan changes below**

### Behaviour change

**Where the proper pair filter flag is used, this code now checks that the paired-end orientation is also used.**
**This will mean that mate-pair orientation (F/F or R/R) will be rejected**
**If you wish to use mate-pair data, please use previous version**

* Where a proper pair filter is used, now check for the correct paired-end orientation of F/R.
* If this is not met the read is ignored.

## 1.13.5

* Build with CaVEMan [1.13.9](https://github.com/cancerit/CaVEMan/releases/tag/1.13.9)

## 1.13.4

* Build with CaVEMan [1.13.8](https://github.com/cancerit/CaVEMan/releases/tag/1.13.8)

## 1.13.3

* Build with CaVEMan [1.13.2](https://github.com/cancerit/CaVEMan/releases/tag/1.13.2)

## 1.13.1

* Build with CaVEMan [1.13.1](https://github.com/cancerit/CaVEMan/releases/tag/1.13.1)

## 1.13.0

* Overlapping read support
* Support for csi indexing of bam/cram files
* Scatter gather implememnted for flagging
* Updates to README, license dates and contact information
* Uses CaVEMan core overlapping reads enabled [1.13.0](https://github.com/cancerit/CaVEMan/releases/tag/1.13.0)

## 1.12.0

* Update to use CaVEMan 1.12.0
* Deal with CaVEMan 1.12.0 intermediate output now being gzipped
