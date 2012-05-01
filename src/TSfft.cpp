/***************************************************************************
                          TSfft.cpp  -  description
                             -------------------
    begin                : Mon 22 April 2003
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
//        TSfft   - CMatrix simple FFT test suite
//
// SYNOPSOIS
//
//      See synopsis_show()
//
// DESCRIPTION
//
//      'TSfft' is a simple test suite designed to test the various
//      CMatrix complex FFT functions.
//
//      Originally, the program started out to test 1D FFTs, but has
//      since been expanded to accommodate both 2D and 3D FFTs. These
//      are specified by command line switches.
//
// HISTORY
// 18 April 2003
// o Initial design and coding.
//
// 26 August 2003
// o Added a cumbersome wrap for 2D.
//
// 01 October 2003
// o Removed "in place" parameter.
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
string		Gstr_version	="$Id: TSfft.cpp,v 1.1.1.1 2004/06/16 15:53:07 rudolph Exp $";
string		Gstr_processName;	// Name of this process

static struct option const longopts[] =
{
  {"cols", 		required_argument, 	NULL, 'c'},
  {"rows", 		required_argument, 	NULL, 'r'},
  {"slices", 		required_argument, 	NULL, 's'},
  {"trials", 		required_argument, 	NULL, 't'},
  {"MKL",		no_argument,		NULL, 'm'},
  {"filesDump",	        no_argument, 		NULL, 'f'},
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
Performs some straightforward 1D, 2D, and 3D FFTs. Internal data is     \n\
created by this process in a manner analogous to FFT routines           \n\
bundled with the Intel Math Kernel Libraries (MKL). This allows         \n\
for easy comparison between output values for this process and          \n\
MKL examples.                                                           \n\
                                                                        \n\
Note however that 3D FFTs have no directly bundled MKL example          \n\
In order to test, a custom MKL analogue should be used.                 \n\
                                                                        \n\
OPTIONS:								\n\
                                                                        \n\
REQUIRED OPTIONS							\n\
                                                                        \n\
        There are no *required* options                                 \n\
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
                        The presence of this parameter implies a 2D FFT \n\
                        and its absence implies a 1D FFT.               \n\
                                                                        \n\
                                                                        \n\
        --slices        <int:slices>                                    \n\
                                                                        \n\
                        Numer of slices in a 3D volume. Meaningful      \n\
                        only in combination with 'cols' and 'rows'      \n\
                        specified.                                      \n\
                                                                        \n\
                NOTE: In order to use the MKL routines, the above       \n\
                dimensions must be lengths of powers of 2.              \n\
                                                                        \n\
	--trials	<int:trials>					\n\
                                                                        \n\
			Number of trials to perform. The core FFT       \n\
                        routine is called *trials* number of times.     \n\
                        This allows for more accurate timing, as well   \n\
                        as monitoring dynamic memory usage of the       \n\
                        process, and comparing speed for \"in place\"   \n\
                        vs. non-\"in place\" function calls.            \n\
                                                                        \n\
        The following options have no arguments, but function as        \n\
        boolean type flag toggles.                                      \n\
                                                                        \n\
	--MKL                                                           \n\
                                                                        \n\
                        If specified, wrap around MKL functions         \n\
                        instead of the default GSL routines.            \n\
                                                                        \n\
	--filesDump                                                     \n\
                                                                        \n\
			If specified, dump the input and output         \n\
                        data streams to files (pzM_input.mat and        \n\
                        pzM_output.mat).                                \n\
                                                                        \n\
	--inPlace		                                        \n\
                        If specified, perform an inplace FFT. The       \n\
                        contents of the original input stream are lost. \n\
                        This could be useful for speedup purposes.      \n\
	                                                                \n\
				EXAMPLES				\n\
									\n\
	TSfft --cols=1024 --t--inPlace	rials=1                         \n\
                                                                        \n\
                Here, a 1D FFT is performed on a 1024 length vector.    \n\
                The MKL is *not* used.                                  \n\
                                                                        \n\
        TSfft --cols=2048 --rows=1024 --trials=10 --MKL --filesDump     \n\
                                                                        \n\
                Perform 10 2D FFT trials.                               \n\
                Wrap around the MKL and dump the input and output       \n\
                data to file.                                           \n\
                                                                        \n\
        TSfft --cols=2048 --rows=1024 --slices=8 --trials=10 --MKL      \n\
              --filesDump                                               \n\
                                                                        \n\
                Perform 10 3D FFT trials.                               \n\
                Wrap around the MKL and dump the input and output       \n\
                data to numbered files.                                 \n\
                                                                        \n\
";
}

    // Create a sample data 1D and 2D sets.
    //  These are based on the GetInputArray functions provided in the MKL
    //  example FFT suite, and allows for direct comparison with MKL
    //  example programs on the same basic data set.

void
dataSet1D_create(
        CMatrix<GSL_complex>*&  pzM_input,
        int                     columns) {
    //
    // ARGS
    //  pzM_input               in/out          holder for input vector
    //  columns                 in              number of columns
    //
    // DESC
    //  Create a simple vector that will be used as input to 1D FFT routine.
    //
    // PRECONDITIONS
    //  o pzM_input must be non-NULL.
    //
    // POSTCONDITIONS
    //  o If pzM_input is not compatible with the columns argc, the
    //    memory is deleted and a new compatible vector is created.
    //
    // HISTORY
    // 02 september 2003
    //  o Initial design and coding.
    //

    double                      v_step  = -M_PI;
    double                      v_step0;
    int                         i;
    GSL_complex                 z_tmp;

    if(!pzM_input->compatible(1, columns)) {
        delete pzM_input;
        pzM_input       = new CMatrix<GSL_complex>(1, columns);
    }

    v_step0 = (2*M_PI)/(double)columns;
    for(i=0; i<columns; i++) {
        GSL_SET_COMPLEX(&z_tmp,
                (sin(v_step)*sqrt(3.0))/2.0,
                sin(v_step)/sqrt(3.0)
        );
        v_step += v_step0;
        pzM_input->val(0, i)   = z_tmp;
    }
}

void
dataSet2D_create(
        CMatrix<GSL_complex>*&  pzM_input,
        int                     columns,
        int                     rows,
        bool                    b_ordered  = false) {
    //
    // ARGS
    //  pzM_input               in/out          holder for input matrix
    //  columns                 in              number of columns
    //  rows                    in              number of rows
    //  b_ordered               in/opt          if true, create data set as
    //                                                  a linearly increasing
    //                                                  element array.
    //
    // DESC
    //  Create a simple matrix that will be used as input to 2D FFT routine.
    //
    // PRECONDITIONS
    //  o pzM_input must be non-NULL.
    //
    // POSTCONDITIONS
    //  o If pzM_input is not compatible with the columns argc, the
    //    memory is deleted and a new compatible vector is created.
    //
    // HISTORY
    // 02 september 2003
    //  o Initial design and coding.
    //
    // 22 January 2003
    //  o Added b_linear
    //

    double                      v_step  = -M_PI;
    double                      v_step0;
    int                         row, col;
    GSL_complex                 z_tmp;

    if(!pzM_input->compatible(rows, columns)) {
        delete pzM_input;
        pzM_input       = new CMatrix<GSL_complex>(rows, columns);
    }

    if(!b_ordered) {
    v_step0 = (2*M_PI)/(double)rows;
        for(row=0; row<rows; row++) {
            GSL_SET_COMPLEX(&z_tmp,
                (sin(v_step)*sqrt(3.0))/2.0,
                sin(v_step)/sqrt(3.0)
            );
            for(col=0; col<columns; col++) {
                pzM_input->val(row, col)        = z_tmp;
            }
            v_step += v_step0;
        }
    } else {
        int value = 0;
        for(row=0; row<rows; row++)
            for(col=0; col<columns; col++) {
                GSL_SET_COMPLEX(&z_tmp,
                        ++value,
                        -value
                );
                pzM_input->val(row, col)        = z_tmp;
            }
    }
}

CVol<GSL_complex>
dataSetD_3Dcreate(
        int                     rows,
        int                     columns,
        int                     slices,
        bool                    b_ordered   = false) {
    //
    // ARGS
    //  columns                 in              number of columns
    //  rows                    in              number of rows
    //  slices                  in              number of slices
    //  b_ordered               in              if true, set each element
    //                                              to an increasing count
    //                                              value.
    //
    // DESC
    //  Create a simple volume that will be used as input to 3D FFT routine.
    //  The volume is composed of a series of matrix. In each matrix, a
    //  single column (at index col/2) is set to contain '1'. All other
    //  values are zero.
    //
    //  In the special case when the optional b_ordered is true, the volume
    //  will have each cell set to an increasing count value.
    //
    // PRECONDITIONS
    //  o None.
    //
    // POSTCONDITIONS
    //  o A single volume conforming to the above description is created.
    //
    // HISTORY
    // 02 september 2003
    //  o Initial design and coding.
    //

    int                         i, j, k, l, m;
    GSL_complex                 z_tmp;
    GSL_SET_COMPLEX(&z_tmp, 0, 0);
    CVol<GSL_complex>		Vol(rows, columns, slices);

    m = 1;
    l = 1;

    for(k=0; k<slices; k++) {
        GSL_SET_COMPLEX(&z_tmp, 0, 0);
	if(!b_ordered) {
	    for(i=0; i<rows; i++) {
		GSL_SET_COMPLEX(&z_tmp, 1, 0);
		Vol(i, columns/2, k)	= z_tmp; 
	    }
	} else {
	    for(i=0; i<rows; i++)
		for(j=0; j<columns; j++) {
                    GSL_SET_COMPLEX(&z_tmp, m++, l--);
		    Vol(i, j, k)	= z_tmp;
		}
	}
    }
    return Vol;
}

int main(int argc, char *argv[]) {

    // Set default values for command-line specified
    //  variables.
    int                 columns         = 512;
    int                 rows            = 0;
    int                 slices          = 0;
    int                 trials          = 1;
    bool                b_filesDump     = false;
    bool                b_MKLwrap       = false;
    bool                b_1D            = false;
    bool                b_2D            = false;
    bool                b_3D            = false;

    float                f_totaltime = 0.0;
    struct tms           st_start, st_stop;
    int                 k;

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
            case 'm':
                b_MKLwrap       = true;
            break;
            case 'f':
                b_filesDump     = true;
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
    }

    if(columns && rows) {
        b_1D    = false;
        b_2D    = true;
        b_3D    = false;
    }

    if(columns && rows && slices) {
        b_1D    = false;
        b_2D    = false;
        b_3D    = true;
    }

    cout << "Columns            = " << columns          << endl;
    cout << "Rows               = " << rows             << endl;
    cout << "Slices             = " << slices           << endl;
    cout << "trials             = " << trials           << endl;
    cout << "b_filesDump        = " << b_filesDump      << endl;
    cout << "b_MKLwrap          = " << b_MKLwrap        << endl;
    cout << "b_1D               = " << b_1D             << endl;
    cout << "b_2D               = " << b_2D             << endl;
    cout << "b_3D               = " << b_3D             << endl;

    CMatrix<GSL_complex>*       pzM_input       = new CMatrix<GSL_complex>(1, 1);
    CMatrix<GSL_complex>*       pzM_output      = new CMatrix<GSL_complex>(1, 1);
    
    CVol<GSL_complex>*		pVl_input	= new CVol<GSL_complex>(1, 1, 1);
    CVol<GSL_complex>*		pVl_output	= new CVol<GSL_complex>(1, 1, 1);

    CMatrix<GSL_complex>**      ppzM_input      = NULL;
    CMatrix<GSL_complex>**      ppzM_output     = NULL;

    if(b_1D) {
        delete                  pzM_input;
        pzM_input               = new CMatrix<GSL_complex>(1, columns);
        dataSet1D_create(pzM_input, columns);
        pzM_output->copy(*pzM_input);
    }

    if(b_2D) {
        delete                  pzM_input;
        pzM_input               = new CMatrix<GSL_complex>(rows, columns);
        dataSet2D_create(pzM_input, columns, rows, true);
        pzM_output->copy(*pzM_input);
    }

    if(b_3D) {
	delete 			pVl_input;
	delete			pVl_output;
	pVl_input		= new CVol<GSL_complex>(rows, columns, slices);
	pVl_output		= new CVol<GSL_complex>(rows, columns, slices);
	*pVl_input              = dataSetD_3Dcreate(rows, columns, slices, true);
    }

    pVl_input->print("Input volume");

    cout << "\nImplementing " << trials << " trial(s) of ";
    cout <<
        (b_2D ?  "pzM_input->fft2D()..."  :
        b_3D ?  "fft3D(ppzM_input)"      :
                "pzM_input->fft1D()...") << endl;
    times(&st_start);
    delete	pVl_output;
    pVl_output	= new CVol<GSL_complex>(rows, columns, slices);
    for(k=0; k<trials; k++) {
	//*pVl_output = pVl_input->zeroPad(e_row, 2);

	if(b_1D) {
	    pzM_output->copy(*pzM_input);
	    pzM_output->fft1D(e_forward, b_MKLwrap);
	}
	if(b_2D) {
	    pzM_output->copy(*pzM_input);
	    //pzM_output->fft2D(e_forward, b_MKLwrap);
            pzM_output->fftshift(true);
	}
	if(b_3D) {
	    pVl_output->copy(*pVl_input);
	    pVl_output->fft3D(e_inverse, b_MKLwrap);
	}
    }
    //pVl_output->saveBinary("test.bin");
    times(&st_stop);
    f_totaltime = difftime(st_stop.tms_utime, st_start.tms_utime) / 100;
    cout << "Total (CPU) time for fft: ";
    cout <<  f_totaltime << " seconds, which is " << f_totaltime/trials << " per trial." << endl << endl;

    pzM_input->print("Input matrix after operation");
    pzM_output->print("Output matrix after operation");

    pVl_input->print("Input volume after operation");
    pVl_output->print("Output volume after operation");

    if(b_filesDump) {
        if(!b_3D) {
            pzM_input->fprint("TSfft_files/pzM_input.mat");
            pzM_output->fprint("TSfft_files/pzM_output.mat");
        } else {
            for(k=0; k<slices; k++) {
                char    pch_num[16];
                sprintf(pch_num, "%d", k);
                string  str_inputName   = "TSfft_files/ppzM_input";
                string  str_outputName  = "TSfft_files/ppzM_output";
                ppzM_input[k]->fprint(str_inputName + pch_num + ".mat");
                ppzM_output[k]->fprint(str_outputName + pch_num + ".mat");
            }
        }
    }

    delete      pzM_input;
    delete      pzM_output;
    
    delete	pVl_input;
    delete	pVl_output;

    cout << endl;
    return EXIT_SUCCESS;
}
