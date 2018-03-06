# Rhub Check
library(rhub)

platform_to_check <- c("windows-x86_64-release",
                       "windows-x86_64-devel",
                       "debian-gcc-release",
                       "debian-gcc-devel",
                       "ubuntu-gcc-release",
                       "ubuntu-gcc-devel",
                       "fedora-clang-devel",
                       "macos-elcapitan-release")

# Ubuntu fails with:
# flowCore error: 'unable to load shared object '/home/docker/R/rrcov/libs/rrcov.so''
# fpc: 'unable to load shared object '/home/docker/R/kernlab/libs/kernlab.so''
#
# Consequently, check_for_cran() fails.
#
# check_on_ubuntu() # http://builder.r-hub.io/status/beadplexr_0.0.0.9000.tar.gz-163d7675f0084c7aa116add851a241e8

check(platform = platform_to_check)

