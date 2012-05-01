/***************************************************************************
                          math_misc.cpp  -  description
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
//	math_misc.cpp	$Id: math_misc.cpp,v 1.1.1.1 2004/06/16 15:53:07 rudolph Exp $
//
// DESCRIPTION
//
//	`math_misc.cpp' contains some miscellaneous mathematical-related
//	routines.
//
//	These routines are not associated with any particular class
//	or object implementation.	
//
// HISTORY
//
// o See the header file, `math_misc.h' for history information
//

#include <iostream>
#include <string>
#include <sstream>
using namespace std;

#include <math.h>
#include "math_misc.h"

void
math_misc_error(
    string          astr_proc,
    string          astr_msg,
    int             code
) {
    //
    // ARGS
    // astr_proc    in      method being evaluated
    // astr_msg     in      error message
    // code         in      exit code
    //
    // DESC
    // Simple error function.
    //
    // HISTORY
    // 22 November 2000
    // o Initial version simply terminates process. Subsequent
    //   versions should throw exceptions.
    //

    cerr << "\nFatal error encountered in <math_misc> execution.\n";
    cerr << "\tCurrent function: "	<< astr_proc	<< endl;
    cerr << "\t" << astr_msg << "\n";
    cerr << "Aborting with code " << code << "\n\n";
    exit(code);
}

float
quantize(
        float   f_num,
	float   f_grain,
	float   f_lowerBound    /*= 0.0         */,
	float   f_upperBound    /*= 0.0         */
) {
    //
    // ARGS
    // f_num		in	number to be quantized
    // f_grain		in	grain factor (expressed as < 1)
    // f_lowerBound	in	lower bound
    // f_upperBound	in	upper bound
    // f_quantized	return	quantized number
    //
    // DESC
    // Quantizes a number with given grain and passed
    // boundaries.
    //
    // If lower and upper bounds are equal to zero (default)
    // then no hardlimiting is implemented; otherwise quantized
    // numbers such that
    //
    // (quantized < f_lowerBound) ? quantized = f_lowerBound
    // (quantized > f_upperBound) ? quantized = f_upperBound
    //

    float f_invGrain = 1 / f_grain;
    float f_quantized = f_num;

    f_quantized *= f_invGrain;
    f_quantized = (int) f_quantized;
    f_quantized /= f_invGrain;

    if(f_lowerBound && f_upperBound) {
	if(f_quantized < f_lowerBound)
	    f_quantized = f_lowerBound;
	if(f_quantized > f_upperBound)
	    f_quantized = f_upperBound;
    }
    return f_quantized;
}

int
doubleCompareAscending(
        const void*     apv_A,
        const void*     apv_B
) {
    //
    // ARGS
    // apv_A            in              number
    // apv_B            in              number
    //
    // DESC
    // Compares the numbers pointed to by the
    // void pointers.
    //
    // Returns:
    //          -1      "A"  < "B"
    //          0       "A" == "B"
    //          1       "A"  > "B"
    //
    // HISTORY
    // 14 November 2000
    // o Initial design and coding
    //

    const double        *a      = (const double*) apv_A;
    const double        *b      = (const double*) apv_B;
    double              diff    = *a - *b;

    return ((diff >= 0.0) ? ((diff > 0.0) ? +1 : 0) : -1);
}

int
doubleCompareDescending(
        const void*     apv_A,
        const void*     apv_B
) {
    //
    // ARGS
    // apv_A            in              number
    // apv_B            in              number
    //
    // DESC
    // Compares the numbers pointed to by the
    // void pointers.
    //
    // Returns:
    //          +1      "A"  < "B"
    //          0       "A" == "B"
    //          -1      "A"  > "B"
    //
    // HISTORY
    // 14 November 2000
    // o Initial design and coding
    //

    const double        *a      = (const double*) apv_A;
    const double        *b      = (const double*) apv_B;
    double              diff    = *a - *b;

    return ((diff >= 0.0) ? ((diff > 0.0) ? -1 : 0) : +1);
}

int
floatCompareAscending(
        const void*     apv_A,
        const void*     apv_B
) {
    //
    // ARGS
    // apv_A            in              number
    // apv_B            in              number
    //
    // DESC
    // Compares the numbers pointed to by the
    // void pointers.
    //
    // Returns:
    //          -1      "A"  < "B"
    //          0       "A" == "B"
    //          1       "A"  > "B"
    //
    // HISTORY
    // 14 November 2000
    // o Initial design and coding
    //

    const float         *a      = (const float*) apv_A;
    const float         *b      = (const float*) apv_B;
    float               diff    = *a - *b;

    return ((diff >= 0.0) ? ((diff > 0.0) ? +1 : 0) : -1);
}

int
floatCompareDescending(
        const void*     apv_A,
        const void*     apv_B
) {
    //
    // ARGS
    // apv_A            in              number
    // apv_B            in              number
    //
    // DESC
    // Compares the numbers pointed to by the
    // void pointers.
    //
    // Returns:
    //          +1      "A"  < "B"
    //          0       "A" == "B"
    //          -1      "A"  > "B"
    //
    // HISTORY
    // 14 November 2000
    // o Initial design and coding
    //

    const float         *a      = (const float*) apv_A;
    const float         *b      = (const float*) apv_B;
    float               diff    = *a - *b;

    return ((diff >= 0.0) ? ((diff > 0.0) ? -1 : 0) : +1);
}

int
intCompareAscending(
        const void*     apv_A,
        const void*     apv_B
) {
    //
    // ARGS
    // apv_A            in              number
    // apv_B            in              number
    //
    // DESC
    // Compares the numbers pointed to by the
    // void pointers.
    //
    // Returns:
    //          -1      "A"  < "B"
    //          0       "A" == "B"
    //          1       "A"  > "B"
    //
    // HISTORY
    // 14 November 2000
    // o Initial design and coding
    //

    const int           *a      = (const int*) apv_A;
    const int           *b      = (const int*) apv_B;
    int                 diff    = *a - *b;

    return ((diff >= 0.0) ? ((diff > 0.0) ? +1 : 0) : -1);
}

int
intCompareDescending(
        const void*     apv_A,
        const void*     apv_B
) {
    //
    // ARGS
    // apv_A            in              number
    // apv_B            in              number
    //
    // DESC
    // Compares the numbers pointed to by the
    // void pointers.
    //
    // Returns:
    //          +1      "A"  < "B"
    //          0       "A" == "B"
    //          -1      "A"  > "B"
    //
    // HISTORY
    // 14 November 2000
    // o Initial design and coding
    //

    const int           *a      = (const int*) apv_A;
    const int           *b      = (const int*) apv_B;
    int                 diff    = *a - *b;

    return ((diff >= 0.0) ? ((diff > 0.0) ? -1 : 0) : +1);
}

bool
bit_random() {
    //
    // DESC
    // Returns a random bit.
    //
    // HISTORY
    // 1 November 2000
    // o Initial design and coding.
    //

    // Choose a random number between 1 and 10 using high order
    // bits. See (Linux) man:rand
    int randomNum       = 1+(int) (10.0*rand()/(RAND_MAX+1.0));

    return (randomNum <=5 ? false : true);
}

int
num_random(
        int     a_range
) {
    //
    // ARGS
    // a_range          in              the upper limit(+1) of
    //                                          the inclusive
    //                                          range.
    // DESC
    // Returns a random number between 0 and a_range-1 (inclusive).
    //
    // HISTORY
    // 2 November 2000
    // o Initial design and coding, using high order bits
    //

    // Choose a random number between  0 and a_range using high order
    // bits. See (Linux) man:rand
    int randomNum       = (int) ((float)a_range*rand()/(RAND_MAX+1.0));
    return randomNum;
}

float
normalise(
        float	f_number,
        float	f_maxVal,
        float	f_lowerBound	        /*= 0.0         */,
        float   f_upperBound	        /*= 1.0         */
) {

    //
    // ARGS
    // f_number         in          number to be normalised
    // f_maxVal         in          the maximum value that f_number
    //                                  might be
    // f_lowerBound     in          ower bound of normalisation
    // f_upperBound     in          upper bound of normalisation
    //
    // DESC
    // Given a number and its maximum value, normalises it between the
    // given range.
    //
    // HISTORY
    // 22 November 2000
    // o Initial design and coding
    //

    string  str_proc    = "normalise";

    float   f_range     =   f_upperBound    - f_lowerBound;
    if(f_range <= 0)
        math_misc_error(
            str_proc,
            "Error with determining range for normalisation",
            1);

    float   f_normPure  = f_number / f_maxVal;
    float   f_norm      = (f_normPure * f_range) + f_lowerBound;

    return(f_norm);
}

float
fmax(
        float   af_a,
        float   af_b
) {
    //
    // ARGS
    // af_a     in      number
    // af_b     in      number
    //
    // POSTCONDITION
    // return max( af_a, af_b)
    //

    return af_a > af_b ? af_a : af_b;
}

float
fmin(
        float   af_a,
        float   af_b
) {
    //
    // ARGS
    // af_a     in      number
    // af_b     in      number
    //
    // POSTCONDITION
    // return min( af_a, af_b)
    //

    return af_a < af_b ? af_a : af_b;
}

double
sqr(
    double	af
) {
    //
    // ARGS
    //  af      in      value to square
    //
    // DESC
    //  Simply returns the square of its argument
    //
    // HISTORY
    // 10 July 2001
    // o Initial design and coding
    //

    return(af*af);
}

bool
isEven(
    int a_num
) {
    //
    // ARGS
    //  a_num       in      number to check
    //
    // DESC
    //  Checks if the passed number is even.
    //
    // POSTCONDITIONS
    //  o If 'even', return true, else return false.
    //
    // HISTORY
    // 07 September 2003
    //  o Intial design and coding.
    //

    if(a_num%2)
        return false;
    else
        return true;

};

bool
isOdd(
    int a_num
) {
    //
    // ARGS
    //  a_num       in      number to check
    //
    // DESC
    //  Checks if the passed number is odd.
    //
    // POSTCONDITIONS
    //  o If 'odd', return true, else return false.
    //
    // HISTORY
    // 07 September 2003
    //  o Intial design and coding.
    //

    if(a_num%2)
        return true;
    else
        return false;

};

string
a2binstr(
    int a_number
) {
    //
    // ARGS
    //  a_number        in      number to return as binary "string"
    //
    // DESC
    //  Converts its argument to a binary "string" format.
    //
    // NOTE
    //  The stringstream buffer seems to lose the very first character
    //  written to it. In order to correct for this strange behaviour,
    //  the input number is first left shifted (multiplied by 2) before
    //  the main algorithm processes it.
    //
    // HISTORY
    // 10 September 2003
    //  o Initial design and coding.
    //

    stringstream    sout("");
    int             remainder;

    if(!a_number)
        return "0";

    a_number <<= 1;

    while(a_number > 0) {
        remainder   = a_number%2;
        a_number  >>= 1;
        sout    << remainder;
    }

    string          str_inverted    = sout.str();
    string          str_ret         = str_inverted;

    for(int i=0; i<str_inverted.length(); i++)
        str_ret[i]  = str_inverted[str_inverted.length()-i];

    return str_ret;
}

void
powersOf2(
    int     a_number,
    int&    a_ones,
    int&    a_highestPower
) {
    //
    // ARGS
    //  a_number        in      number to analyse
    //  a_ones          in/out  number of "ones" in binary decomposition
    //  a_highestPower  in/out  highest power of 2 in number
    //
    // DESC
    //  Analyses a given number and returns the number of "ones" in the
    //  binary decomposition, as well as the highest power of two in the
    //  number.
    //
    // HISTORY
    // 10 September 2003
    //  o Initial design and coding.
    //

    a_ones          = 0;
    a_highestPower  = 0;

    while(a_number > 0) {
        if(a_number%2) a_ones++;
        a_number  >>= 1;
        a_highestPower++;
    }
    a_highestPower--;
}

int 
shiftIndex(
    int	                a_dimensionLength,
    int			a_index,
    int			a_pivotOffset
) {
    //
    // ARGS
    //	a_dimensionLength 
    //                  in		length (span) of dimension along which
    //                                          to shift
    //	a_index		in		index along dimension that is to be
    //						shifted
    //  a_pivotOffset    in              offset to pivot point
    //                                          about indices are shifted.
    //
    // DESC
    //	This method calculates the shift in an index position, and is
    //	often used when indices need to be re-calculated for fft/ifft
    //	operations
    //
    // PRECONDITIONS
    //  o In the interests of speed, no checking is done on whether
    //	  or not a_index falls within [0, 1, ..., a_dimensionLength-1].
    //	  If you send ridiculous values, you'll get ridiculous results!	  
    //
    // POSTCONDITIONS
    //  o A shifted index is returned. 
    //    For fftshift,  a_pivotOffset is +1
    //	  Fot ifftshift, a_pivotOffset is -1
    //
    // HISTORY
    // 25 February 2004
    //  o Initial design and coding, based around the inPlace space
    //	  algorithm.
    //
    
    
	
    int indexDelta;			// the value that an index "jumps"
    int	indexNext;			// "next" value after jump
	
    if(isOdd(a_dimensionLength))
	indexDelta	= (a_dimensionLength - a_pivotOffset)/2;
    else	
	indexDelta	= a_dimensionLength/2;
	
	    
    indexNext = a_index + indexDelta;
    
    // Check for boundary violation in indices
    if(indexNext >= a_dimensionLength)
	indexNext -= a_dimensionLength;
    return indexNext;
}

int 	
ifftshiftIndex(
    int     a_dimensionLength,
    int     a_index
) { 
    return shiftIndex(
	a_dimensionLength,
	a_index,	
	-1);
}

int
fftshiftIndex(
    int     a_dimensionLength,
    int     a_index
) { 
    return shiftIndex(
	a_dimensionLength,
	a_index,	
	+1);
}

