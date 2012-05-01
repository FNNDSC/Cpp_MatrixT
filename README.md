Cpp_MatrixT
==========

A C++ Matrix library.

### Dependencies ###

The cppmatrixt library depends on the following libraries:

 * FFTW3
 * GSL
 * cmake 2.6

### Ubuntu Linux ###

Install dependent packages:

    $ sudo apt-get install fftw-dev
    $ sudo apt-get install libgsl0-dev
    $ sudo apt-get install cmake

To build, in root directory:

    cppmatrix2/$ mkdir build
    cppmatrix2/$ cd build
    cppmatrix2/build$ cmake ../
    cppmatrix2/build$ make

### Mac OS X ###

Install dependent packages:

    $ sudo port install fftw
    $ sudo port install gsl
    $ sudo port install cmake

To build, in root directory:

    cppmatrix2/$ mkdir build
    cppmatrix2/$ cd build
    cppmatrix2/build$ cmake ../
    cppmatrix2/build$ make
    
    