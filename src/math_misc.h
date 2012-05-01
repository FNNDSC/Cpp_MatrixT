/***************************************************************************
                          math_misc.h  -  description
                             -------------------
    begin                : Thu Feb 10 2000
    copyright            : (C) 2000 by Rudolph Pienaar
    email                : pienaar@bme.ri.ccf.org
                           rudolph@nmr.mgh.harvard.edu
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
//
// NAME
//
//	math_misc.h	$Id: math_misc.h,v 1.1.1.1 2004/06/16 15:53:07 rudolph Exp $
//
// DESCRIPTION
//
//	`math_misc.h' contains some miscellaneous mathematical-related
//	routines.
//
//	These routines are not associated with any particular class
//	or object implementation.	
//
// HISTORY
//
// 09-02-1999
// o quantize routine
//
// 02-24-2000
// o radix conversion routines
//
// 1 November 2000
// o Added bit randomisation
//
// 2 November 2000
// o Added number range randomisation
//
// 14 November 2000
// o Added <type>Compare functions
//
// 21 November 2000
// o Added simple (scalar) normalisation function
//
// 02 December 2000
// o Added simple max/min functions
//
// 10 July 2001
// o Added sqr()
//
// 07 September 2003
// o Added isEven() / isOdd()
//
// 25 February 2004
// o Added ((i)fft)shiftIndex(...)
//
// DEVELOPMENT NOTES
//  Might want to consider templetising some of these functions
//

#ifndef __MATH_MISC_H__
#define __MATH_MISC_H__

#include <cstdlib>
#include <string>
using namespace std;


float   quantize(   float   f_num,
                    float   f_grain,
                    float   f_lowerBound	= 0.0,
                    float   f_upperBound 	= 0.0);
float   normalise(  float	f_num,
                    float	f_maxVal,
                    float	f_lowerBound	= 0.0,
                    float	f_upperBound	= 1.0);

float   fmax(   float   af_a,
                float   af_b);
float   fmin(   float   af_a,
                float   af_b);

double  sqr(    double  af);

int     doubleCompareAscending( const void* apv_A,
                                const void* apv_B);
int     floatCompareAscending(  const void* apv_A,
                                const void* apv_B);
int     intCompareAscending(    const void* apv_A,
                                const void* apv_B);
int     doubleCompareDescending(const void* apv_A,
                                const void* apv_B);
int     floatCompareDescending( const void* apv_A,
                                const void* apv_B);
int     intCompareDescending(   const void* apv_A,
                                const void* apv_B);

bool    bit_random();
int     num_random( int  a_range);

bool    isEven(     int a_num);
bool    isOdd(      int a_num);

string  a2binstr(       int a_number);
void    powersOf2(
                        int     a_number,
                        int&    a_ones,
                        int&    a_highestPower
                    );


int     shiftIndex(
                        int     a_dimensionLength,
			int     a_index,
			int     a_pivotOffset
			);

int 	ifftshiftIndex(
                        int     a_dimensionLength,
			int     a_index
			);

int 	fftshiftIndex(
                        int     a_dimensionLength,
			int     a_index
			);


#endif // __MATH_MISC_H__

