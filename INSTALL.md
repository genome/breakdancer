BreakDancer uses CMake which is a cross-platform build tool. Basically it will
generate a Makefile so you can use `make`. The requirements are the zlib,
development library, gcc, gmake, cmake 2.8+. Beginning with version 1.4.4, 
BreakDancer includes samtools as part of the build process.

Here are the steps to build:

    # --recursive option is important so that it gets the submodules too
    $ git clone --recursive https://github.com/genome/breakdancer.git
    ...
    Resolving deltas: 100% (38/38), done.
    Submodule path 'build-common': checked out '...'

    $ cd breakdancer
    $ mkdir build
    $ cd build

    $ cmake .. -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=/usr/local
    ...
    -- Build files have been written to: .../breakdancer/build

    $ make
    ...
    Linking CXX executable ../../../../bin/breakdancer-max
    [100%] Built target breakdancer-max

    $ sudo make install
