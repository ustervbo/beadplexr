## Resubmission
This is a resubmission. In this version I have extended the 'description' field of DESCRIPTION to provide more information. The 'description' was changed from

* Functions and wrappers for automatic analysis of fcs files resulting from a CBA, LEGENDplex, or MACSPlex experiment.

to

* Reproducible and automated analysis of multiplex bead assays such as CBA, LEGENDplex, and MACSPlex. The package provides functions for streamlined reading of fcs  files, and identification of bead clusters and analyte expression. It eases the  calculation of a standard curve and the subsequent calculation of the analyte concentration.

## Test environments

### Local
* Ubuntu 15.10, R 3.3.1

### Rbub (05.Mar.2018)

* [Pass] Windows, R release: http://builder.r-hub.io/status/beadplexr_0.0.0.9000.tar.gz-8a4d51955c9c43d5b3de780a167f68f6
* [Pass] Windows, R devel: http://builder.r-hub.io/status/beadplexr_0.0.0.9000.tar.gz-0d7139a7001c4b52a8eb61c48ebde03a
* [Pass] Debian, R release: http://builder.r-hub.io/status/beadplexr_0.0.0.9000.tar.gz-24afd509009f4a1589502e48ff4f43d6
* [Pass] Debian, R devel: http://builder.r-hub.io/status/beadplexr_0.0.0.9000.tar.gz-1b717b454b944e879c9439cda8fdd102
* [Fail] Ubuntu, R relase: http://builder.r-hub.io/status/beadplexr_0.0.0.9000.tar.gz-0b397672d0384a94ad70881b3e9b8f6f
* [Pass] Ubuntu, R devel: http://builder.r-hub.io/status/beadplexr_0.0.0.9000.tar.gz-872fdf827a9144f59fb7fd5beca47a1d
* [Pass] Fedora, R devel: http://builder.r-hub.io/status/beadplexr_0.0.0.9000.tar.gz-52d38adc50864e958b3984060b953a5e
* [Pass] MacOS, R release: http://builder.r-hub.io/status/beadplexr_0.0.0.9000.tar.gz-644f35591ef64c4ab938ec0ef0cfde33

## R CMD check results
With the exception of a single platform, there were no ERRORs, WARNINGs, or NOTEs.

The problematic platform was Ubuntu 16.04 LTR, R release. Here two problems raised an error:

* 'fpc' package
  Error in dyn.load(file, DLLpath = DLLpath, ...) :
  unable to load shared object '/home/docker/R/kernlab/libs/kernlab.so':
  libRlapack.so: cannot open shared object file: No such file or directory

* 'flowCore' package
  Error in dyn.load(file, DLLpath = DLLpath, ...) :
  unable to load shared object '/home/docker/R/rrcov/libs/rrcov.so':
  libRlapack.so: cannot open shared object file: No such file or directory

However, the problem seems to be with `libRlapack.so`, rather than with `beadplexr`
