/***************************************************************************
                          TSvol.cpp  -  description
                             -------------------
    begin                : Wed 04 February 2004
    copyright            : (C) 2004 by Rudolph Pienaar
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
//        TSvol   - simple volume structure test suite
//
// SYNOPSOIS
//
//      See synopsis_show()
//
// DESCRIPTION
//
//      `TSvol' is a simple test suite designed to test the various
//      "compound" (high dimensional) CVol* classes.
//
// HISTORY
// 04 February 2004
//  o Initial design and coding.
//

using namespace std;

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <stdlib.h>
#include <complex>
#include <cmath>

#include <sys/times.h>
#include "cmatrix.h"
#include "math_misc.h"

#include <getopt.h>

// Some global variables
string		Gstr_version	="$Id";
string		Gstr_processName;	// Name of this process

static struct option const longopts[] =
{
  {"cols", 		required_argument, 	NULL, 'c'},
  {"rows", 		required_argument, 	NULL, 'r'},
  {"slices", 		required_argument, 	NULL, 's'},
  {"trials", 		required_argument, 	NULL, 't'},
  {"dim4",              required_argument,      NULL, '4'},
  {"dim5",              required_argument,      NULL, '5'},
  {NULL, 0, NULL, 0}
};

void
synopsis_show() {
    //
    // DESC
    //	Show the usage of process
    //
    // HISTORY
    // 27 August 2003
    // o Initial design and coding.
    //

    cout << "Usage: " << Gstr_processName	<< " [OPTIONS]" << endl;
    cout << "\
A straightforward platform on which to test high dimension volume       \n\
performance of the CMatrix library.                                     \n\
                                                                        \n\
OPTIONS:								\n\
                                                                        \n\
REQUIRED OPTIONS							\n\
                                                                        \n\
        There are no *required* options.                                \n\
                                                                        \n\
        If no options are passed, internal defaults will be used.       \n\
                                                                        \n\
\"OPTIONAL\" OPTIONS							\n\
									\n\
        --cols          <int:cols>                                      \n\
                                                                        \n\
                        This specifies the number of elements in        \n\
                        a single 1D vector, or the number of columns    \n\
                        for a 2D matrix.                                \n\
                                                                        \n\
        --rows          <int:rows>                                      \n\
                                                                        \n\
                        Specifies the number of rows in a 2D matrix.    \n\
                                                                        \n\
        --slices        <int:slices>                                    \n\
                                                                        \n\
                        Numer of slices in a 3D volume. Meaningful      \n\
                        only in combination with 'cols' and 'rows'      \n\
                        specified.                                      \n\
                                                                        \n\
        --dim4          <int:dim4>                                      \n\
                                                                        \n\
                        Size of the 4th volume dimension.               \n\
                                                                        \n\
        --dim5          <int:dim5>                                      \n\
                                                                        \n\
                        Size of the 5th volume dimension.               \n\
                                                                        \n\
";
}

CVol<GSL_complex_float>
dataSet3D_create(
        CVol<GSL_complex_float>*        pVol,
        int                             rows,
        int                             columns,
        int                             slices,
        int                             start           = 1,
        int                             inc             = 1
) {
    //
    // ARGS
    //  columns                 in              number of columns
    //  rows                    in              number of rows
    //  slices                  in              number of slices
    //  start                   in              start value of first slice
    //                                                  upper-left quad
    //  inc                     in              increment between adjacent
    //                                                  cell values
    //
    // DESC
    //  Create a simple 3D volume.
    //
    // PRECONDITIONS
    //  o None.
    //
    // POSTCONDITIONS
    //  o A single volume conforming to the above description is created.
    //
    // HISTORY
    // 05 February 2004
    //  o Initial design and coding.
    //

    int                                 i, j, k, l, m;
    GSL_complex_float                   z_tmp;
    GSL_SET_COMPLEX(&z_tmp, start, start);

    m = start;
    l = -start;

    for(k=0; k<slices; k++) {
        for(i=0; i<rows; i++)
	    for(j=0; j<columns; j++) {
                GSL_SET_COMPLEX(&z_tmp, m, l);
                m += inc;
                l -= inc;
		pVol->val(i, j, k)	= z_tmp;
	    }
    }
    return *pVol;
}

CVol4D<GSL_complex_float>
dataSet4D_create(
        CVol4D<GSL_complex_float>*      pVol,
        int                             rows,
        int                             columns,
        int                             slices,
        int                             dim4,
        int                             start           = 1,
        int                             inc             = 1
) {
    //
    // ARGS
    //  columns                 in              number of columns
    //  rows                    in              number of rows
    //  slices                  in              number of slices
    //  dim4                    in              extent of 4th dimension
    //  start                   in              start value of first slice
    //                                                  upper-left quad
    //  inc                     in              increment between adjacent
    //                                                  cell values
    //
    // DESC
    //  Create a simple 3D volume.
    //
    // PRECONDITIONS
    //  o None.
    //
    // POSTCONDITIONS
    //  o A single volume conforming to the above description is created.
    //
    // HISTORY
    // 05 February 2004
    //  o Initial design and coding.
    //

    int                                 i, j, k, l, m;
    int                                 L, M;
    GSL_complex_float                   z_tmp;
    GSL_SET_COMPLEX(&z_tmp, start, start);


    for(l=0; l<dim4; l++) {
        M =  start * (int)pow(10., l);
        L = -start * (int)pow(10., l);
        for(k=0; k<slices; k++) {
            for(i=0; i<rows; i++)
	        for(j=0; j<columns; j++) {
                    GSL_SET_COMPLEX(&z_tmp, M, L);
                    M += inc;
                    L -= inc;
		    pVol->val(i, j, k, l)	= z_tmp;
	        }
        }
    }
    return *pVol;
}

CVol5D<GSL_complex_float>
dataSet5D_create(
        CVol5D<GSL_complex_float>*      pVol,
        int                             rows,
        int                             columns,
        int                             slices,
        int                             dim4,
        int                             dim5,
        int                             start           = 1,
        int                             inc             = 1
) {
    //
    // ARGS
    //  columns                 in              number of columns
    //  rows                    in              number of rows
    //  slices                  in              number of slices
    //  dim4                    in              extent of 4th dimension
    //  start                   in              start value of first slice
    //                                                  upper-left quad
    //  inc                     in              increment between adjacent
    //                                                  cell values
    //
    // DESC
    //  Create a simple 3D volume.
    //
    // PRECONDITIONS
    //  o None.
    //
    // POSTCONDITIONS
    //  o A single volume conforming to the above description is created.
    //
    // HISTORY
    // 05 February 2004
    //  o Initial design and coding.
    //

    int                                 i, j, k, l, m, n;
    int                                 L, M;
    GSL_complex_float                   z_tmp;
    GSL_SET_COMPLEX(&z_tmp, start, start);


    for(m=0; m<dim5; m++)
        for(l=0; l<dim4; l++) {
            M =  start * (int)pow(10., l);
            L = -start * (int)pow(10., l);
            for(k=0; k<slices; k++) {
                for(i=0; i<rows; i++)
	            for(j=0; j<columns; j++) {
                        GSL_SET_COMPLEX(&z_tmp, M, L);
                        M += inc;
                        L -= inc;
		        pVol->val(i, j, k, l, m)	= z_tmp;
	            }
            }
        }
    return *pVol;
}


int main(int argc, char *argv[]) {

    // Set default values for command-line specified
    //  variables.
    int                 rows            = 4;
    int                 columns         = 4;
    int                 slices          = 3;
    int                 dim4            = 0;
    int                 dim5            = 0;
    int                 trials          = 1;
    bool                b_1D            = false;
    bool                b_2D            = false;
    bool                b_3D            = false;
    bool                b_4D            = false;
    bool                b_5D            = false;

    float               f_totaltime = 0.0;
    struct tms          st_start, st_stop;
    int                 i=0, j=0, k=0;

    Gstr_processName	= argv[0];

    while(1) {
        int opt;
        int optionIndex = 0;
        opt = getopt_long(argc, argv, "", longopts, &optionIndex);
        if( opt == -1)
            break;

        switch(opt) {
            case 'c':
                columns	        = atoi(optarg);
            break;
            case 'r':
                rows	        = atoi(optarg);
            break;
            case 's':
                slices	        = atoi(optarg);
            break;
            case 't':
                trials	        = atoi(optarg);
            break;
            case '4':
                dim4            = atoi(optarg);
            break;
            case '5':
                dim5            = atoi(optarg);
            break;
            case '?':
                synopsis_show();
                exit(1);
            break;
            default:
                cout << "?? getopt returned character code " << opt << endl;
        }
    }

    if(columns) {
        b_1D    = true;
        b_2D    = false;
        b_3D    = false;
        b_4D    = false;
        b_5D    = false;
    }

    if(columns && rows) {
        b_1D    = false;
        b_2D    = true;
        b_3D    = false;
        b_4D    = false;
        b_5D    = false;
    }

    if(b_2D && slices) {
        b_1D    = false;
        b_2D    = false;
        b_3D    = true;
        b_4D    = false;
        b_5D    = false;
    }

    if(b_3D && dim4) {
        b_1D    = false;
        b_2D    = false;
        b_3D    = false;
        b_4D    = true;
        b_5D    = false;
    }

    if(b_4D && dim5) {
        b_1D    = false;
        b_2D    = false;
        b_3D    = false;
        b_4D    = false;
        b_5D    = true;
    }

    cout << "Rows               = " << rows             << endl;
    cout << "Columns            = " << columns          << endl;
    cout << "Slices             = " << slices           << endl;
    cout << "Dim4               = " << dim4             << endl;
    cout << "Dim5               = " << dim5             << endl;
    cout << "trials             = " << trials           << endl;
    cout << "b_1D               = " << b_1D             << endl;
    cout << "b_2D               = " << b_2D             << endl;
    cout << "b_3D               = " << b_3D             << endl;
    cout << "b_4D               = " << b_4D             << endl;
    cout << "b_5D               = " << b_5D             << endl;

    CVol<GSL_complex_float>*            pVl3D_input     = new CVol<GSL_complex_float>(1, 1, 1);
    CVol<GSL_complex_float>*            pVl3D_output    = new CVol<GSL_complex_float>(1, 1, 1);

    CVol4D<GSL_complex_float>*		pVl4D_input	= new CVol4D<GSL_complex_float>(1, 1, 1, 1);
    CVol4D<GSL_complex_float>*		pVl4D_output	= new CVol4D<GSL_complex_float>(1, 1, 1, 1);
    CVol4D<GSL_complex_float>*          pVl4D_test      = new CVol4D<GSL_complex_float>(1, 1, 1, 1);

    CVol5D<GSL_complex_float>*		pVl5D_input	= new CVol5D<GSL_complex_float>(1, 1, 1, 1, 1);
    CVol5D<GSL_complex_float>*		pVl5D_output	= new CVol5D<GSL_complex_float>(1, 1, 1, 1, 1);
    CVol5D<GSL_complex_float>*          pVl5D_test      = new CVol5D<GSL_complex_float>(1, 1, 1, 1, 1);

    if(b_1D) {
        cerr << "This program should only be used on 4 (and higher) dimensional volumes.";
        exit(1);
    }

    if(b_2D) {
        cerr << "This program should only be used on 4 (and higher) dimensional volumes.";
        exit(1);
    }

    if(b_3D) {
	delete 			pVl3D_input;
	delete			pVl3D_output;
	pVl3D_input		= new CVol<GSL_complex_float>(rows, columns, slices);
	pVl3D_output		= new CVol<GSL_complex_float>(rows, columns, slices);
        dataSet3D_create(pVl3D_input, rows, columns, slices);
    }

    if(b_4D) {
	delete 			pVl4D_input;
	delete			pVl4D_output;
        delete                  pVl4D_test;
	pVl4D_input		= new CVol4D<GSL_complex_float>(rows, columns, slices, dim4);
	pVl4D_output		= new CVol4D<GSL_complex_float>(rows, columns, slices, dim4);
        dataSet4D_create(pVl4D_input, rows, columns, slices, dim4);

        pVl4D_test              = new CVol4D<GSL_complex_float>(*pVl4D_input);
        pVl4D_test              = pVl4D_output;

    }

    if(b_5D) {
	delete 			pVl5D_input;
	delete			pVl5D_output;
        delete                  pVl5D_test;
	pVl5D_input		= new CVol5D<GSL_complex_float>(rows, columns, slices, dim4, dim5);
	pVl5D_output		= new CVol5D<GSL_complex_float>(rows, columns, slices, dim4, dim5);
        dataSet5D_create(pVl5D_input, rows, columns, slices, dim4, dim5);

        pVl5D_test              = new CVol5D<GSL_complex_float>(*pVl5D_input);
        pVl5D_test              = pVl5D_output;

    }


//    pVl_input->print("Input volume");

    cout << "\nImplementing " << trials << " loop trial(s)."    << endl;
    times(&st_start);

    for(k=0; k<trials; k++) {
	if(b_3D) {
           pVl3D_output->copy(*pVl3D_input);
           pVl3D_output->ifftshift(true);
	}
	if(b_4D) {
           pVl4D_output->copy(*pVl4D_input);
           pVl4D_output->vol3D(0).fftshift(true);
	}
	if(b_5D) {
           pVl5D_output->copy(*pVl5D_input);
           //pVl5D_output->vol3D(0, 1) = pVl5D_input->vol3D(0, 1).fftshift(false);
	   pVl5D_output->vol3D(0, 1).fftshift(true);
	}

    }
    //pVl_output->saveBinary("test.bin");
    times(&st_stop);
    f_totaltime = difftime(st_stop.tms_utime, st_start.tms_utime) / 100;
    cout << "Total (CPU) time for fft: ";
    cout <<  f_totaltime << " seconds, which is " << f_totaltime/trials << " per trial." << endl << endl;

    pVl3D_input->print("Input volume after operation");
    pVl3D_output->print("Output volume after operation");
    
    CMatrix<int>	M_dim3(1, 3);
    CMatrix<int>	M_index(1, 3);
    CMatrix<int>	M_shift(1, 3);
    M_dim3(0)	= rows;
    M_dim3(1)	= columns;
    M_dim3(2)	= slices;
    
    for(int i=0; i<rows; i++)
	for(int j=0; j<columns; j++)
	    for(int k=0; k<slices; k++) 
	    {
   
		M_index(0)	= i;
		M_index(1)	= j;
		M_index(2)	= k;
		M_shift = M_dim3.fftshiftIndex(M_index);
		
		M_index.print("\nindex");
		M_shift.print("shift");
	    }	

    delete	pVl3D_input,    pVl4D_input,    pVl5D_input;
    delete	pVl3D_output,   pVl4D_output,   pVl5D_output;

    cout << endl;
    return EXIT_SUCCESS;
}
