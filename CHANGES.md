# CHANGES

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
