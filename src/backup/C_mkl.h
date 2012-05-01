/***************************************************************************
                          C_mkl.h  -  description
                             -------------------
    begin                : Wed 22 April 2003
    copyright            : (C) 2003 by Rudolph Pienaar
    email                : rudolph@nmr.mgh.harvard.edu
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/*
//
// NAME
//
//	C_mkl.h
//
// DESCRIPTION
//
//      'C_mkl.h' header for a thin C-based wrapper around core MKL routines.
//
//      Since the MKL functions are C-based, and strictly C++ compatible
//      headers are *not* supplied in the MKL distribution, relevant functions
//      are explicitly declared here. By wrapping them around an explicit
//      extern "C" structure, we prevent C++ name-mangling.
//
// HISTORY
// 22 April 2003
// o Initial design and coding.
//
// 26 August 2003
// o Added zfft2dc(...)
*/

#ifdef __cplusplus
        extern  "C" {
#endif

#ifndef __C_MKL_H__
#define __C_MKL_H__

	  //#include <mkl/mkl_fft.h>
//#include </opt/intel/mkl/include/mkl_fftc_ln.h>

void zfft1dc(double* r, double* i, int n, int isign, double* wsave);
void zfft2dc(double* r, double* i, int m, int n, int isign);

void cfft1dc(float* r, float* i, int n, int isign, float* wsave);
void cfft2dc(float* r, float* i, int m, int n, int isign);

#endif //__C_MKL_H__

#ifdef __cplusplus
}
#endif

