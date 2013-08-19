BreakDancer uses CMake which is a cross-platform build tool. Basically it will
generate a Makefile so you can use `make`. The requirements are samtools,
gcc, gmake, cmake 2.8+. Beginning with version 1.4.0, BreakDancer requires 
samtools 0.1.19 or higher. Older versions require samtools 0.1.6.

Here's the steps to build:

    # --recrusive option is important so that it gets the submodules too
    $ git clone --recursive https://github.com/genome/breakdancer.git
    ...
    Resolving deltas: 100% (38/38), done.
    Submodule path 'build-common': checked out '...'

    $ cd breakdancer
    $ mkdir build
    $ cd build

    $ wget "http://downloads.sourceforge.net/project/samtools/samtools/0.1.19/samtools-0.1.19.tar.bz2"
    $ tar -jxvf samtools-0.1.19.tar.bz2
    $ cd samtools-0.1.19
    $ make
    $ export SAMTOOLS_ROOT=$(pwd)
    $ cd ..

    $ cmake .. -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=/usr/local
    ...
    -- Build files have been written to: .../breakdancer/build

    $ make
    ...
    Linking CXX executable ../../../../bin/breakdancer-max
    [100%] Built target breakdancer-max

    $ sudo make install
