## Overview

BreakDancerMax has the following dependencies:

* git
* cmake 2.8+ ([cmake.org](http://cmake.org))
* zlib
* ncurses
* A C++ compiler

The associated helper program `bam2cfg.pl` also requires the following Perl modules:
* Statistics::Descriptive
* GD::Graph

## Build Instructions
### Obtaining Build dependencies

* For APT-based systems (Debian, Ubuntu), install the following packages:

```
sudo apt-get install build-essential git-core cmake zlib1g-dev libncurses-dev
```

* For RPM-based systems (Fedora, CentOS, RHEL), install the following packages instead:

```
sudo yum groupinstall "Development tools" 
sudo yum install zlib-devel ncurses-devel cmake
```

Note that for some RPM based systems (like RHEL), you will need to install cmake 2.8 or greater yourself as the packaged version is much older.

### Obtain the BreakDancer sourcecode

Clone the git repository

```
git clone https://github.com/genome/breakdancer.git
```

This creates a breakdancer directory in the current working directory. We refer
to this directory as `$BD_ROOT` from now on.

### Build BreakDancer

BreakDancer does not support in-source builds. So create a subdirectory, enter it, build, and run tests:

```
mkdir $BD_ROOT/build
cd $BD_ROOT/build
cmake ..
make deps
make
make test
```

The binary `breakdancer-max` can then be found under `$BD_ROOT/build/bin`. If you have administrative rights, then run `sudo make install` to install the tool for all users under `/usr/bin`.

### Obtaining Perl modules
The procedure of downloading CPAN modules is described on the [CPAN Web site]
(http://www.cpan.org/modules/INSTALL.html).
