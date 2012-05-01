/***************************************************************************
                          speedComp.cpp  -  description
                             -------------------
    begin                : Tue Sep 03 18:31:02 EST 2002
    copyright            : (C) 2002 by Rudolph Pienaar
    email                : pienaar@bme.ri.ccf.org
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
//      speedComp
//
// SYNOPSIS
//      speedComp <numberOfTrials> [<defaultSize=100>]
//
// DESC
//      `speedComp' shows how to dynamically create two matrices, multiply
//      them together, and store the result in a third.
//
//      While demonstrating a straightforward C-style approach to achieve
//      this, and comparing it with a C++-style CMatrix based solution,
//      the real purpose of this program is to compare the speed differences
//      between the two approaches.
//
// NOTE
//      When analysing time taken for multiplications, explore the effect
//      of different optimisation flags. For g++, setting CFLAGS=-O3 seemed
//      to result in the most speed.
//
// HISTORY
// 03 September 2002
// o Initial design and coding.
//
// 04 March 2003
// o Added <defaultSize=100> argument
// o Added gsl support
//
// 05 March 2003
// o Changed matrix cell elements to a causal function of indices
//	(this is to simplify comparison across different math libraries)
//

using namespace std;

#ifdef HAVE_CONFIG
#include <config.h>
#endif

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <sys/times.h>

//
// gsl includes
//
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics.h>

#include "cmatrix.h"

int main(int argc, char *argv[])
{

   float                f_totaltime = 0.0;
   struct tms           st_start, st_stop;

//////
////// C-style method
//////
    int			DEFSIZE	= 100;
    if(argc==3)
    	DEFSIZE			= atoi(argv[argc-1]);
    cout << "Default size " << DEFSIZE << endl;
    const int           ROWS    = DEFSIZE;  	// Number of rows in matrix
    const int           COLS    = DEFSIZE;	// Number of cols in matrix
    int                 TRIALS;         	// Number of trials
    int                 i,j,k;          	// Misc for loop counting variable
    int                 row	= 0;
    int			col	= 0;    	// Misc matrix variables
    char                pch_text[1024]; 	// Misc 1k char buffer
    double              f_sum;          	// Needed for matrix mult
//    int                 numRand = 0;

    if(argc==1 || argc>3) {
        printf("Synopsis:\n\t./speedComp <trials> [<default size=100>]\n");
        exit(1);
    }

    TRIALS = atoi(argv[1]);

    // Create two blocks of memory containing random double data
    //  of size 100 rows by 100 cols
    printf("Creating memory block A with (%d x %d) values...\n", DEFSIZE, DEFSIZE);
    double**    ppf_blockA      = (double**) calloc(ROWS, sizeof(double*));
    if(ppf_blockA==NULL) {
        errno   = ENOMEM;
        perror("Cannot allocate memory for entire block A");
        exit(1);
    }
    for(i=0; i<ROWS; i++) {
        ppf_blockA[i]           = (double*) calloc(COLS, sizeof(double));
        if(ppf_blockA[i]==NULL) {
            sprintf(pch_text,
                        "Cannot allocate memory for block A, row %d", i);
            errno   = ENOMEM;
            perror(pch_text);
            exit(1);
        }
        for(j=0; j<COLS; j++) {
            ppf_blockA[i][j]    = i+j;
        }
    }

    printf("Creating memory block B with (%d x %d) values...\n", DEFSIZE, DEFSIZE);
    double**    ppf_blockB      = (double**) calloc(ROWS, sizeof(double*));
    if(ppf_blockB==NULL) {
        errno   = ENOMEM;
        perror("Cannot allocate memory for entire block B");
        exit(1);
    }
    for(i=0; i<ROWS; i++) {
        ppf_blockB[i]           = (double*) calloc(COLS, sizeof(double));
        if(ppf_blockB[i]==NULL) {
            sprintf(pch_text,
                        "Cannot allocate memory for block B, row %d", i);
            errno   = ENOMEM;
            perror(pch_text);
            exit(1);
        }
        for(j=0; j<COLS; j++) {
            ppf_blockB[i][j]    = i-j;
        }
    }

    // Create a third block to contain results
    printf("Creating memory block C...\n");
    double**    ppf_blockC      = (double**) calloc(ROWS, sizeof(double*));
    if(ppf_blockC==NULL) {
        errno   = ENOMEM;
        perror("Cannot allocate memory for entire block C");
    }
    for(i=0; i<ROWS; i++) {
        ppf_blockC[i]           = (double*) calloc(COLS, sizeof(double));
        if(ppf_blockC[i]==NULL) {
            sprintf(pch_text,
                        "Cannot allocate memory for block C, row %d", i);
            errno   = ENOMEM;
            perror(pch_text);
        }
    }

    // Now perform the actual matrix multiplication
    printf("\nImplementing %d trials of block A multiplied by block B...\n", TRIALS);
    times(&st_start);
    for (k=0; k<TRIALS; k++)
        for (row=0; row < ROWS; row++)
            for (col=0; col < COLS; col++) {
                f_sum = 0.0;
                for (int i=0; i<COLS; i++)
                    f_sum += ppf_blockA[row][i] * ppf_blockB[i][col];
                ppf_blockC[row][col] = f_sum;
            }
    times(&st_stop);
    f_totaltime = difftime(st_stop.tms_utime, st_start.tms_utime) / 100;
    printf("Total time for C-style multiplication: %f seconds.\n", f_totaltime);

    free(ppf_blockA);
    free(ppf_blockB);
    free(ppf_blockC);

//////
////// C++-style method using CMatrix
//////
    
    // Create some CMatrix matrices from the already-declared structures
    CMatrix<double>     M_A(ROWS, COLS);
    CMatrix<double>     M_B(ROWS, COLS);
    CMatrix<double>     M_C(ROWS, COLS);
    
    for(i=0; i<ROWS; i++) {
    	for(j=0; j<COLS; j++) {
            M_A.val(i, j) = i+j+2;
            M_B.val(i, j) = i-j;
	}
    }

    cout << "\nImplementing " << TRIALS << " trials of M_C = M_A * M_B..." << endl;
    times(&st_start);
    for(k=0; k<TRIALS; k++)
        M_C = M_A * M_B;
    times(&st_stop);
    f_totaltime = difftime(st_stop.tms_utime, st_start.tms_utime) / 100;
    cout << "Total time for CMatrix multiplication: ";
    cout <<  f_totaltime << " seconds." << endl << endl;

    bool	b_sample = true;
    cout << "Max A =\t\t\t\t" << M_A.max() << endl;
    cout << "Max B =\t\t\t\t" << M_B.max() << endl << endl;
    cout << "Mean of M_A:\t\t\t"                << M_A.mean()   << endl;
    cout << "Standard deviation of M_A:\t"      << M_A.std(b_sample)    << endl << endl;
    cout << "Mean of M_B:\t\t\t"                << M_B.mean()   << endl;
    cout << "Standard deviation of M_B:\t"      << M_B.std(b_sample)    << endl << endl;
    cout << "Mean of M_C:\t\t\t"                << M_C.mean()   << endl;
    cout << "Standard deviation of M_C:\t"      << M_C.std(b_sample)    << endl << endl;

//////
////// C++-style method using complex CMatrix
//////
    
    // Create some CMatrix matrices from the already-declared structures
    CMatrix<GSL_complex>     cM_A(ROWS, COLS);
    CMatrix<GSL_complex>     cM_B(ROWS, COLS);
    CMatrix<GSL_complex>     cM_C(ROWS, COLS);
    
    GSL_complex		z0;
    GSL_complex		z1;
    
    for(i=0; i<ROWS; i++) {
    	for(j=0; j<COLS; j++) {
	    GSL_SET_COMPLEX(&z0, i+j+2, i-j);
	    GSL_SET_COMPLEX(&z1, i-j, i+j+2);
            cM_A.val(i, j) = z0;
            cM_B.val(i, j) = z1;
	}
    }

    cout << "\nImplementing " << TRIALS << " trials of (complex) cM_C = cM_A * cM_B..." << endl;
    times(&st_start);
    for(k=0; k<TRIALS; k++)
        cM_C = cM_A * cM_B;
    times(&st_stop);
    f_totaltime = difftime(st_stop.tms_utime, st_start.tms_utime) / 100;
    cout << "Total time for CMatrix multiplication: ";
    cout <<  f_totaltime << " seconds." << endl << endl;
    
//
// gsl method
//
    gsl_matrix* gsl_A = gsl_matrix_alloc(ROWS, COLS);
    gsl_matrix* gsl_B = gsl_matrix_alloc(ROWS, COLS);
    gsl_matrix* gsl_C = gsl_matrix_alloc(ROWS, COLS);

    for(i=0; i<ROWS; i++) {
    	for(j=0; j<COLS; j++) {
            gsl_matrix_set(gsl_A,i,j,i+j+2);
            gsl_matrix_set(gsl_B,i,j,i-j);
	}
    }

    cout << "\nImplementing " << TRIALS << " trials of gsl_C = gsl_A * gsl_B..." << endl;
    times(&st_start);
    for(k=0; k<TRIALS; k++)
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                  1.0, gsl_A, gsl_B,
                  0.0, gsl_C);
    
    times(&st_stop);
    f_totaltime = difftime(st_stop.tms_utime, st_start.tms_utime) / 100;
    cout << "Total time for gsl multiplication: ";
    cout <<  f_totaltime << " seconds." << endl << endl;
    	    
    cout << "Max gsl_A =\t\t\t";
    cout << gsl_stats_max(gsl_A->data, 1, gsl_A->size1*gsl_A->size2) << endl;
    cout << "Max gsl_B =\t\t\t";
    cout << gsl_stats_max(gsl_B->data, 1, gsl_B->size1*gsl_B->size2) << endl << endl;

    cout << "Mean of gsl_A:\t\t\t";
    cout << gsl_stats_mean(gsl_A->data, 1, gsl_A->size1*gsl_A->size2) << endl;
    cout << "Standard deviation of gsl_A:\t";
    cout << gsl_stats_sd(gsl_A->data, 1, gsl_A->size1*gsl_A->size2) << endl << endl;

    cout << "Mean of gsl_B:\t\t\t";
    cout << gsl_stats_mean(gsl_B->data, 1, gsl_B->size1*gsl_B->size2) << endl;
    cout << "Standard deviation of gsl_B:\t";
    cout << gsl_stats_sd(gsl_B->data, 1, gsl_B->size1*gsl_B->size2) << endl << endl;

    cout << "Mean of gsl_C:\t\t\t";
    cout << gsl_stats_mean(gsl_C->data, 1, gsl_C->size1*gsl_C->size2) << endl;
    cout << "Standard deviation of gsl_C:\t";
    cout << gsl_stats_sd(gsl_C->data, 1, gsl_C->size1*gsl_C->size2) << endl << endl;
    
    return EXIT_SUCCESS;
}
