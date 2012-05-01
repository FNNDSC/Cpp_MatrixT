/***************************************************************************
                          TSsimpleMaths.cpp  -  description
                             -------------------
    begin                : Fri 18 April 2003
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
//
// NAME
//
//	TSsimpleMaths - CMatrix simple maths test suite
//
// SYNOPSOIS
//
//	TSsimpleMaths
//
// DESCRIPTION
//
//	'TSsimpleMaths' is a simple test suite designed to test the various
//	CMatrix operator overload methods.
//
// HISTORY
// 18 April 2003
// o Initial design and coding.
//	
 
using namespace std;

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <stdlib.h>
#include <complex>
#include <cmath>

#include "cmatrix.h"

int main(int argc, char *argv[])
{
    GSL_complex		gz_a(10, 10);
    GSL_complex		gz_b(2, 2);

    cout << -gz_a;
    
    // Use one of the following blocks of code - make sure the
    //	second block is commented out!
        
    CMatrix<GSL_complex>		M_A(10, 10, gz_a);
    CMatrix<GSL_complex>		M_B(10, 10, gz_b);
    CMatrix<GSL_complex>		M_C(10, 10);
    
//    CMatrix<double>		M_A(10, 10, 10);
//    CMatrix<double>		M_B(10, 10, 2);
//    CMatrix<double>		M_C(10, 10);
    
    M_A.print("M_A");
    M_B.print("M_B");
    M_C.print("M_C");

    cout << "Unary negation of M_A" 	<< endl;
    (-M_A).print("-M_A");    
    
    M_C		= M_A * M_B;		// * (matrix)
    M_C.print("M_C = M_A * M_B");
    
    M_C 	= M_A * (GSL_complex)(4);	// * (scalar)
    M_C.print("M_C = M_A * 4");
    
    M_A *= M_B;				// *= (matrix)
    M_A.print("M_A *= M_B");
    M_A /= M_B;				// */ (matrix)
    M_A.print("M_A /= M_B");
        
    M_A *= (GSL_complex)30;		// *= (scalar)
    M_A.print("M_A *= 30");
    M_A /= (GSL_complex)30;				// /= (scalar)
    M_A.print("M_A /= 30");

    M_C		= M_A / M_B;		// / (matrix)
    M_C.print("M_C = M_A / M_B");
    
    M_C 	= M_A / (GSL_complex)3;		// / (scalar)
    M_C.print("M_C = M_A / 3");

              
    M_C		= M_A + M_B;		// + (matrix)
    M_C.print("M_C = M_A + M_B");
    
    M_C 	= M_A + (GSL_complex)3;		// + (scalar)
    M_C.print("M_C = M_A + 3");
    
    M_A += M_B;				// += (matrix)
    M_A.print("M_A += M_B");
    M_A -= M_B;				// -= (matrix)
    M_A.print("M_A -= M_B");

    M_C		= M_A - M_B;		// - (matrix)
    M_C.print("M_C = M_A - M_B");
    
    M_C 	= M_A - (GSL_complex)3;		// - (scalar)
    M_C.print("M_C = M_A - 3");

    GSL_complex                         vz(10, 10);
    CMatrix<GSL_complex>		M_Ci(10, 10, vz);
    
    M_Ci.print("M_Ci");
    M_Ci = M_C.inverse();
    M_Ci.print("M_Ci");
    (M_Ci * M_C).print("M_Ci * M_C");
                      
    return EXIT_SUCCESS;
}
