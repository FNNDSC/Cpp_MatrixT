
AM_CPPFLAGS=-DGSL_USE -I/opt/local/include

bin_PROGRAMS = cppmatrixt TSVol TSfftf TSfft TSsimpleMaths TSmanipulate speedComp TSconstructors
cppmatrixt_SOURCES = cppmatrixt.cpp

# set the include path found by configure
INCLUDES= $(all_includes)

# the library search path.
cppmatrixt_LDFLAGS = $(all_libraries)

libCMatrix_la_SOURCES =  cmatrix.cpp cmatrix.h C_mkl.h math_misc.cpp math_misc.h
lib_LTLIBRARIES =  libCMatrix.la










speedComp_LDADD = $(top_builddir)/src/libCMatrix.la
TSfft_LDADD = $(top_builddir)/src/libCMatrix.la
TSsimpleMaths_LDADD = $(top_builddir)/src/libCMatrix.la
TSmanipulate_LDADD = $(top_builddir)/src/libCMatrix.la
TSconstructors_LDADD = $(top_builddir)/src/libCMatrix.la
cppmatrixt_LDADD = $(top_builddir)/src/libCMatrix.la

TSfftf_LDADD = $(top_builddir)/src/libCMatrix.la
TSVol_SOURCES =  TSvol.cpp
TSVol_LDADD = $(top_builddir)/src/libCMatrix.la
speedComp_SOURCES =  speedComp.cpp
TSconstructors_SOURCES =  TSconstructors.cpp
TSmanipulate_SOURCES =  TSmanipulate.cpp
TSsimpleMaths_SOURCES =  TSsimpleMaths.cpp
TSfft_SOURCES =  TSfft.cpp
TSfftf_SOURCES =  TSfftf.cpp

libCMatrix_la_LIBADD = -L/opt/local/lib -lgsl -lgslcblas -lpthread -l$(LIB_FFTW_A)

