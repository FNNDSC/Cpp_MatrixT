/***************************************************************************
                          TSmanipulate.cpp  -  description
                             -------------------
    begin                : Wed April 2 EST 2003
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
//	TSmanipulate - demonstrate some CMatrix class usage
//
// SYNOPSOIS
//
//	TSmanipulate
//
// DESCRIPTION
//
//	'TSmanipulate' shows how to use the CMatrix class in some
//	common manipulation operations.
//
//	It's primary purpose is to test the manipulate methods and check
//	for memory leaks. This is accomplished by running the program and
//	checking system usage concurrently with a monitor program such
//	as 'top'.
//
// HISTORY
// 02 April 2003
// o Initial design and coding.
//	
 
using namespace std;

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <stdlib.h>

#include "cmatrix.h"

int main(int argc, char *argv[])
{

    // Constructor examples
    CMatrix<double>     M_A(10, 10, 0.0);
                // Create a matrix, M_A, of size 10x10
                //      with values set to 0.0
    
    int 		i 	= 0;		// arbitrary variable
    const int		LOOP	= 100000;	// loop limit
    
    for(int row=0; row<M_A.rows_get(); row++)
    	for(int col=0; col<M_A.cols_get(); col++)
	    M_A(row, col) = i++;

    M_A.print("M_A");
    
    // Matrix remove and replace
    for(i=0; i<LOOP; i++) {
    	CMatrix<double>*	pM_Asub	= new CMatrix<double>(1,1);
    	M_A.matrix_remove(pM_Asub, 0, 0, 5, 5);
	M_A.matrix_replace(4, 4, *pM_Asub);
	delete pM_Asub;
    }
    M_A.print("M_A after matrix remove / replace");
    
    // Row remove / insert / replace
    for(i=0; i<LOOP; i++) {
    	CMatrix<double>*	pV_row	= new CMatrix<double>(1,1);
	CMatrix<double>*	pM_Aa	= new CMatrix<double>(1,1);
	pM_Aa->copy(M_A);
	pM_Aa->row_remove(pV_row, 0);
	pM_Aa->row_insert(*pV_row, 0);
	pM_Aa->row_replace(*pV_row, 9);
	if(i==LOOP-1)
	   M_A.copy(*pM_Aa);
	delete	pV_row;
	delete 	pM_Aa;
    }
    
    M_A.print("M_A after row remove / insert / replace");
    
    // Col remove / insert / replace
    for(i=0; i<LOOP; i++) {
    	CMatrix<double>*	pV_col	= new CMatrix<double>(1,1);
	CMatrix<double>*	pM_Aa	= new CMatrix<double>(1,1);
	pM_Aa->copy(M_A);
	pM_Aa->col_remove(pV_col, 0);
	pM_Aa->col_insert(*pV_col, 0);
	pM_Aa->col_replace(*pV_col, 9);
	if(i==LOOP-1)
	   M_A.copy(*pM_Aa);
	delete	pV_col;
	delete 	pM_Aa;
    }
    
    M_A.print("M_A after col remove / insert / replace");
    
    // diagonal replace
    for(i=0; i<LOOP; i++) {
    	CMatrix<double>*	pV_col	= new CMatrix<double>(1,1);
	CMatrix<double>*	pM_Aa	= new CMatrix<double>(1,1);
	pM_Aa->copy(M_A);
	pM_Aa->col_remove(pV_col, 0);
	(*pV_col) *= 0.0;
	(*pV_col) += -1;
	pM_Aa->diagonal_replace(*pV_col);
	pM_Aa->diagonal_replace(*pV_col, false);
	if(i==LOOP-1)
	   M_A.copy(*pM_Aa);
	delete	pV_col;
	delete 	pM_Aa;
    }
    
    M_A.print("M_A after diagonal_replace (major and minor)");
    
    // diagonal replace
    CMatrix<double>*	pM_Aa	= new CMatrix<double>(1,1);
    for(i=0; i<LOOP; i++) {
    	CMatrix<double>*	pV_row	= new CMatrix<double>(1,1);
	M_A.row_remove(pV_row, 0);
	
	pV_row->diagonalise();	
	pM_Aa->copy(*pV_row);
	delete	pV_row;
    }
    
    pM_Aa->print("pM_Aa diagonalised");
    delete pM_Aa;
    
    for(i=0; i<LOOP; i++) {
	CMatrix<double>*	pM_Aar	= new CMatrix<double>(1,1);
	CMatrix<double>*	pM_Aac	= new CMatrix<double>(1,1);
    	M_A.matrix_remove(pM_Aar, 0, 0, 2, 2);
    	M_A.matrix_remove(pM_Aac, 0, 0, 2, 2);
	pM_Aar->replicate(2, e_row);
	pM_Aac->replicate(2, e_column);
	if(i==LOOP-1) {
	    pM_Aar->print("pM_Aar replicated row-wise");	
	    pM_Aac->print("pM_Aac replicated col-wise");	
	}
	delete pM_Aar;
	delete pM_Aac;
    }

    for(i=0; i<LOOP; i++) {
	CMatrix<double>*	pM_Asorted	= new CMatrix<double>(1,1);
	CMatrix<double>*	pM_AsortedD	= new CMatrix<double>(1,1);
	pM_Asorted->copy(M_A.sort());
	pM_AsortedD->copy(M_A.sort(e_descending));
	if(i==LOOP-1) {
	    pM_Asorted->print("pM_Aar sorted");	
	    pM_AsortedD->print("pM_Aar sorted descending");	
	}
	delete pM_Asorted;
    }
                
    return EXIT_SUCCESS;
}
