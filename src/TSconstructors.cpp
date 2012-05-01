/***************************************************************************
                          TSconstructors.cpp  -  description
                             -------------------
    begin                : Tue 1 April 2003
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
//	TSConstructors - CMatrix constructor test suite
//
// SYNOPSOIS
//
//	TSconstructors
//
// DESCRIPTION
//
//	'TSconstructors' is a simple test suite designed to test the various
//	CMatrix constructor methods.
//
// HISTORY
// 01 April 2003
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

    double      pf_dataSet1[] = {
        -1.221,		1.434,		4.336,		6.778,
        -5.432,		4.254,		8.543,		9.346
    };
    double	pf_dataSetRow[] = {
    	1,	 	2,		3,		4	
    };
    double	pf_dataSetCol[] = {
	-1,
	-2    
    };

    // Constructor examples
    CMatrix<double>	M_I10("I", 10);
    		// Create a 10x10 identity matrix
    CMatrix<double>     M_A(10, 10, -1);
                // Create a matrix, M_A, of size 10x10
                //      with values set to -1 (default
                //      is 0.0)
    CMatrix<double>     M_B(M_A);
                // Create a matrix, M_B, initialised
                //      from M_A
    M_B.randomize(0, 100);
    CMatrix<double>	M_Aa(10, 10, 0.0);
    CMatrix<double>     M_C(2, 4, pf_dataSet1);
                // Create a matrix, M_C, of size 2x4
                //      with values read from float array
                //      pf_dataSet1
    CMatrix<double>	M_D(4, 2, pf_dataSet1);
                // Create a matrix, M_D, of size 4x2
                //      with values read from float array
                //      pf_dataSet1
    CMatrix<double>*    pM_E;
    pM_E        = new CMatrix<double>(10,5,0.0);
                // A pointer-based example.
    pM_E->randomize(0, 100);
    CMatrix<double>	M_c(1, 4, pf_dataSetRow);
    CMatrix<double>	M_d(2, 1, pf_dataSetCol);
    CMatrix<double>	M_ct(4, 1);
    
    CMatrix<double>	V_100(1, 100, 0.0);
    for(int i=0; i<V_100.cols_get(); i++)
    	V_100(0, i) = i;
    CMatrix<double>	M_Vr(10, 10, V_100, e_row);
    CMatrix<double>	M_Vc(10, 10, V_100, e_column);
    
    CMatrix<double>	V_r100(M_Vr, e_row, e_rowVector);
    CMatrix<double>	V_c100(M_Vr, e_column, e_columnVector);
    
    // Print examples
    M_I10.print("Matrix M_I10");
    cout << endl;
    
    M_Vr.print("Row dominant creation from vector");
    M_Vc.print("Column dominant creation from vector");
    V_r100.print("V_r100"); cout << endl;
    V_c100.print("V_c100"); cout << endl;
    
    M_A.print("Matrix A");
    cout << M_A.sprint("Matrix A (using sprint)")				<< endl;

    M_Aa.randomly_fill(1, 0.1);
    M_Aa.print("M_Aa randomly filled");
    cout									<< endl;

    M_B.print("Matrix B");
    M_C.print("Matrix C");
    M_D.print("Matrix D");
    pM_E->print("Matrix E");
    cout << endl;
    
    cout << "Inserting row into M_C..." 					<< endl;
    M_C.row_insert(M_c, 1);
    M_C.print("New Matrix C");
    cout << endl;
    
    cout << "Inserting col into M_D..." 					<< endl;
    M_ct = !M_c;
    M_D.col_insert(M_ct, 2);
    M_D.print("New Matrix D");
    cout << endl;
    
    cout << "Saving Matrix B to file \'4x5random.mat\'"				<< endl;
    M_B.fprint("4x5random.mat", "Matrix B from \'simpleConstructors.cpp\'");
    cout << endl;
    cout << "Readig Matrix Bfile from file \'4x5random.mat\'"			<< endl;
    CMatrix<double>	M_Bfile("4x5random.mat");
    M_Bfile.print("M_Bfile");
    M_Bfile.hardlimit(20, 80, 20, 80);
    M_Bfile.print("M_Bfile hardlimited");
    M_Bfile.quantize(0.1);
    M_Bfile.print("M_Bfile quantized to 0.1");
    M_Bfile.delta_add(100, 21, 79);
    M_Bfile.print("M_Bfile.delta_add");
    
    V_c100.quantize_linearly(-50, 50);
    V_c100.print("V_c100.quantize_linearly");
    
    V_c100.quantize_linearlyShifted(-50, 50);
    V_c100.print("V_c100.quantize_linearlyShifted");

    cout << "Diagonalising M_c..."						<< endl;
    M_c.diagonalise();
    M_c.print("M_c");
    cout << endl;

    CMatrix<double> M_Cr(1,1);
    CMatrix<double> M_Cc(1,1);

    cout << "Replicating M_C row-wise twice..." 				<< endl;
    M_Cr.copy(M_C);
    M_Cr.replicate(2, e_row);
    M_Cr.print("Replicated M_C row-wise");
    cout << endl;

    cout << "Replicating M_C column-wise twice..." 				<< endl;
    M_Cc.copy(M_C);
    M_Cc.replicate(2, e_column);
    M_Cc.print("Replicated M_C column-wise");
    cout << endl;

    GSL_complex			gz_A(100, 10);
    CMatrix<GSL_complex>	Mc_gsl0(5,5, -10);
    CMatrix<GSL_complex>	Mc_gsl1(5,5, gz_A);
    CMatrix<GSL_complex>	Mc_gsl2(5,5);
    Mc_gsl0.print("Mc_gsl0");
    Mc_gsl1.print("Mc_gsl1");
    Mc_gsl2.print("Mc_gsl2");

    Mc_gsl2 = Mc_gsl0 * Mc_gsl1;
    
    Mc_gsl2.print("Mc_gsl2");
//    Mc_gsl2.fprint("complex.mat");
    
    CMatrix<GSL_complex>	Mc_z0("complex.mat");
    Mc_z0.print("Mc_z0");
            
    delete pM_E;
    return EXIT_SUCCESS;
}
