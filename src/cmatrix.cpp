/***************************************************************************
                          cmatrix.cpp  -  description
                             -------------------
    begin                : Thu Feb 10 2000
    copyright            : (C) 2000 by Rudolph Pienaar
    email                : pienaar@bme.ri.ccf.org
                           rudolph@nmr.mgh.harvard.edu
 ***************************************************************************/
//
//                                NOTE
// Portions of this library are based on source code orginally written by
// Bruce Eckel. These portions are copyright either Osborne/McGraw-Hill or
// copyright Bruce Eckel. See http://www.mindview.net
//
// Additional contributors include:
// Antonio        (alrl1@alu.um.es)
//
//
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
//      cmatrix.cpp     $Id: cmatrix.cpp,v 1.5 2004/07/28 16:01:47 rudolph Exp $
//
// DESCRIPTION
//
//      `cmatrix.cpp' defines the class structures used in a matrix/
//      vector handling application.
//
// HISTORY
//
//  (See cmatrix.h for history information)
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <new>
#include <complex>
using namespace std;

#define TINY 1e-20

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <math_misc.h>
#include <time.h>

#ifdef __APPLE__
#include <dfftw3.h>
#else
#include <fftw3.h>
#endif

#include "cmatrix.h"
#include "C_mkl.h"

#ifdef GSL_USE
std::ostream& operator<< (std::ostream& o, const GSL_complex& rval) {

    return o << "( " << rval.dat[0] << "; " << rval.dat[1] << "i )";
}

std::ostream& operator<< (std::ostream& o, const GSL_complex_float& rval) {

    return o << "( " << rval.dat[0] << "; " << rval.dat[1] << "i )";
}


#endif

const type_info&        Ginfo_intMatrix                 = typeid(CMatrix<int>);
const type_info&        Ginfo_floatMatrix               = typeid(CMatrix<float>);
const type_info&        Ginfo_doubleMatrix              = typeid(CMatrix<double>);
const type_info&        Ginfo_complexFloatMatrix        = typeid(CMatrix<GSL_complex_float>);
const type_info&        Ginfo_complexDoubleMatrix       = typeid(CMatrix<GSL_complex>);

template<typename _CMDATA>
int CMatrix<_CMDATA>::createdCounter        = 0;

template<typename _CMDATA>
void
CMatrix<_CMDATA>::_error (
        char *apch_proc,
        char *apch_msg,
        int   code) {
    //
    // ARGS
    //  apch_proc               in              (optional) function name
    //                                                  in which error has
    //                                                  occured
    //  apch_msg                in              (optional) error message
    //                                                  text
    //  code                    in              (optional) error code
    //
    // DESC
    //  Main error reporting method for the CMatrix class.
    //
    // TODO
    //        o Wrap exception throwing around #ifdefs?
    //
    // HISTORY
    // 12 March 2003
    //        o Added id, rows, and col output
    //

    cerr << "\nFatal error encountered.\n";
    cerr << "\tCMatrix object id: " << ps_mat->id << endl;
    cerr << "\tRows: " << ps_mat->rows << ", cols: " << ps_mat->cols << endl;
    cerr << "\tCurrent function: " << "CMatrix::" << apch_proc << "\n";
    cerr << "\t" << apch_msg << "\n";
    cerr << "Throwing an exception to (this) with code " << code << endl;
    throw (this);
}

//////---------->
////// Non Class functions
//////---------->

void
CM_error(
        char*   apch_procname,
        char*   apch_errorMessage,
        int     errorCode) {
    cerr << "\nFatal error encountered.\n";
    cerr << "\tCMatrix object non class function\n";
    cerr << "\tCurrent function: "                 << apch_procname        << endl;
    cerr << "\t"                                 << apch_errorMessage    << endl;
    cerr << "Exiting to system with code "         << errorCode            << endl;
    exit(errorCode);
}

template<typename _CMDATA>
CMatrix<_CMDATA>
linearSystem_solve(
        CMatrix<_CMDATA>&        aM_A,
        CMatrix<_CMDATA>&        aV_b
) {
    //
    // ARGS
    //  aM_A            in              the 'A' matrix
    //  aV_b            in              the 'b' vector
    //
    // DESC
    //  Determines the vector x that solves the expression:
    //          A*x = b
    //
    // PRECONDITIONS
    //  o A must be a square matrix
    //  o b must be a column vector
    //
    // POSTCONDTIONS
    //  o creates a vector, x, and passes this vector out.
    //
    // HISTORY
    //  14 February 2002
    //  o Incorporated patch from Antonio.
    //
    // 03 April 2003
    //        o Templatization.
    //
    
#ifndef CMATRIX_NORANGE
    char* pch_proc = "linearSystem_solve";
    if(!aM_A.matrix_type() != e_matrix)
        CM_error(pch_proc, "\'A\' needs to be a matrix.");
    if(!aV_b.is_column_vector())
        CM_error(pch_proc, "\'b\' needs to be a column vector");
#endif

    CMatrix<_CMDATA>     V_linearSolution(aM_A.cols_get(), 1, 0.0);

    V_linearSolution = aM_A.inverse() * aV_b;
    return V_linearSolution;
}

template<typename _CMDATA>
CMatrix<_CMDATA>
leastSquaresError_find(
        CMatrix<_CMDATA>&        aM_A,
        CMatrix<_CMDATA>&        aV_b
) {
    //
    // ARGS
    //  aM_A            in              the 'A' matrix
    //  aV_b            in              the 'b' vector
    //
    // DESC
    //  Determines the vector x that minimsed the expression:
    //          (A*x-b)^2
    //
    // PRECONDITIONS
    //  o A must be a square matrix
    //  o b must be a column vector
    //
    // POSTCONDTIONS
    //  o creates a vector, x, and passes this vector out.
    //
    // HISTORY
    //  14 February 2002
    //  o Incorporated patch from Antonio.
    //
    // 03 April 2003
    //        o Templatization.
    //

#ifndef CMATRIX_NORANGE
    char* pch_proc = "leastSquaresError_find";
    if(aM_A.matrix_type() != e_matrix)
        CM_error(pch_proc, "\'A\' needs to be a matrix.");
    if(!aV_b.is_column_vector())
        CM_error(pch_proc, "\'b\' needs to be a column vector");
#endif

    CMatrix<_CMDATA>        V_leastSquareError(aM_A.cols_get(), 1, 0.0);

    // this is matrix for linear system coeficients (of the variables).
    CMatrix<_CMDATA>        M_lsc(aM_A.cols_get(), aM_A.cols_get(),0.0);
    // and this will be linear system's objectives matrix.
    CMatrix<_CMDATA>        V_lsv(aM_A.cols_get(), 1, 0.0);

    // we compute objectives of the equations.
    for (int k=0; k<aM_A.rows_get(); k++)
        for (int i=0; i<aM_A.cols_get(); i++)
            V_lsv.val(i,0) += aM_A.val(k,i)*aV_b.val(k,0);

    // we just process triangular matrix of coefs of the equations
    // because ...
    for (int k=0; k<aM_A.rows_get(); k++)
        for (int i=0; i<M_lsc.rows_get(); i++)
            for (int j=0; j<=i; j++)
                M_lsc.val(i,j) += aM_A.val(k,j) * aM_A.val(k,i);

    // ... elements of the equations are symmetrical, so we save about half
    // of the cpu time in computing the linear equations.
    for (int i=0; i<M_lsc.rows_get(); i++)
        for (int j=i+1; j<M_lsc.cols_get(); j++)
            M_lsc.val(i,j) =        M_lsc.val(j,i);

    V_leastSquareError = linearSystem_solve(M_lsc, V_lsv);
    return V_leastSquareError;
}


//////---------->
////// Constructor, destructor, and miscellaneous initialization
////// routines
//////---------->

template<typename _CMDATA>
void
CMatrix <_CMDATA>::coreMatrix_construct(
        int                arows,
        int                acols
) {
    // ARGS
    //        arows                in                        number of rows in matrix
    //        acols                in                        number of cols in matrixtype
    //
    // DESCRIPTION
    //        This method is a consolidation of the main initialization code
    //        that is shared across all the *matrix* constructors.
    //
    // PRECONDITIONS
    //        o Should only be called from a constructor method.
    //
    // POSTCONDITIONS
    //        o A gsl-aware ps_mat structure is created.
    //
    // HISTORY
    // 19 March 2003
    //        o Initial design and coding.
    //

    // Use RTTI to determine datatype
    const type_info&        Tinfo = typeid(*this);                        // RTTI data

    // create the structure
    try {
        ps_mat = new matstruct;
        ps_mat->data    = (_CMDATA**) 0x0;
    } catch (bad_alloc xa) {
        _error("base constructor",
                "could not allocate memory for matstruct (heap exhausted?)");
    }
    ps_mat->rows = arows;
    ps_mat->cols = acols;
    // allocate memory for the structure

#ifndef GSL_USE
    try {
        ps_mat->data = new _CMDATA *[ps_mat->rows];
    } catch (bad_alloc xa) {
        _error("base constructor",
                "could not allocate memory for matrix rows (heap exhausted?)");
    }
    for (int y = 0; y < arows; y++) {
            try {
        ps_mat->data[y] = new _CMDATA[ps_mat->cols];
        } catch (bad_alloc xa) {
        _error("base constructor",
                "could not allocate memory for matrix cols (heap exhausted?)");
        }
    }
#else
    //ps_mat->ps_gsl = (gslwrap*)malloc(sizeof(gslwrap));
    if( Tinfo!=Ginfo_intMatrix                  &&
        Tinfo!=Ginfo_floatMatrix                &&
        Tinfo!=Ginfo_doubleMatrix               &&
        Tinfo!=Ginfo_complexDoubleMatrix        &&
        Tinfo!=Ginfo_complexFloatMatrix
        )
        _error("base constructor", "GSL interface not defined for requested type");

    ps_mat->ps_gsl                              = new gslwrap;
    ps_mat->ps_gsl->matrix_int                  = 0x0;
    ps_mat->ps_gsl->matrix_float                = 0x0;
    ps_mat->ps_gsl->matrix_double               = 0x0;
    ps_mat->ps_gsl->matrix_complex              = 0x0;
    ps_mat->ps_gsl->matrix_complex_float        = 0x0;

    if(Tinfo==Ginfo_intMatrix) {
        ps_mat->ps_gsl->matrix_int              = gsl_matrix_int_alloc(arows, acols);
        if(!ps_mat->ps_gsl->matrix_int)
            _error("base constructor", "gsl_matrix_int_alloc error");
        ps_mat->data = (_CMDATA**) ps_mat->ps_gsl->matrix_int->data;
    }
    if(Tinfo==Ginfo_floatMatrix) {
        ps_mat->ps_gsl->matrix_float            = gsl_matrix_float_alloc(arows, acols);
        if(!ps_mat->ps_gsl->matrix_float)
            _error("base constructor", "gsl_matrix_float_alloc error");
        ps_mat->data =  (_CMDATA**)ps_mat->ps_gsl->matrix_float->data;
    }
    if(Tinfo==Ginfo_doubleMatrix) {
        ps_mat->ps_gsl->matrix_double           = gsl_matrix_alloc(arows, acols);
        if(!ps_mat->ps_gsl->matrix_double)
            _error("base constructor", "gsl_matrix_alloc error");
        ps_mat->data = (_CMDATA**) ps_mat->ps_gsl->matrix_double->data;
    }
#endif

    // set first reference/id to this data
    ps_mat->refcnt = 1;
    ps_mat->id           = CMatrix::createdCounter++;
}

template <typename _CMDATA>
CMatrix<_CMDATA>::CMatrix(
        int             mrows,
        int             mcols,
        _CMDATA         initvalue) {
    //
    // ARGS
    //  mrows                   in              number of rows for matrix
    //  mcols                   in              number of cols for matrix
    //  initvalue               in              initial value for matrix
    //                                                  elements
    //
    // DESC
    //  Matrix constructor
    //
    // PRECONDITIONS
    //  o Make sure that mrows and mcols are not negative numbers!
    //
    // HISTORY
    // 18 March 2003
    //        o Began gsl integration
    //

    coreMatrix_construct(mrows, mcols);
     
    for (int i = 0; i < mrows; i++)
        for (int j = 0; j < mcols; j++)
            _mval (i, j) = initvalue;
}

template <typename _CMDATA>
CMatrix<_CMDATA>::CMatrix(
        int             mrows,
        int             mcols,
        const _CMDATA*  p_initvalues) {
    //
    // ARGS
    //  mrows                   in              number of rows for matrix
    //  mcols                   in              number of cols for matrix
    //  p_initvalues            in              pointer to block of memory
    //                                                  containing initial
    //                                                  values
    //
    // DESC
    //  Matrix constructor
    //
    // PRECONDITIONS
    //  o Make sure that mrows and mcols are not negative numbers!
    //
    // HISTORY
    // 19 March 2003
    //        o Integration with gsl/coreMatrix_construct
    //

    coreMatrix_construct(mrows, mcols);
    
    int        c = 0;
    for (int i = 0; i < mrows; i++)
        for (int j = 0; j < mcols; j++)
            _mval (i, j) = p_initvalues[c++];
}

template<typename _CMDATA>
CMatrix<_CMDATA>::CMatrix(
        int                             mrows,
        int                             mcols,
        CMatrix<_CMDATA>&                aV,
        e_DOMINANCE                     e_adominance) {
    //
    // ARGS
    //  mrows           in              rows of matrix
    //  mcols           in              cols of matrix
    //  aV              in              Vector whose data will be used to
    //                                          construct matrix
    //  e_adominance    in              Will vector be built from row-first
    //                                          or col-first from this matrix
    //
    // DESC
    //  Builds a matrix from the passed vector.
    //
    //  The `dominance' variable defines whether the vector is
    //  row or column dominant, i.e. the values from the vector
    //  are assumed organized as [row1, row2, row3, ... , rowN]
    //  or as [col1, col2, col3, ... , colN]
    //
    // POSTCONDITIONS
    //  o A *matrix* is built.
    //
    // HISTORY
    // 31 March 2003
    //        o Templatization.
    //

    char*        pch_name         = "CMatrix (construct from vector)";
    char*        pch_errorSize         = "Base vector cannot build target matrix";

    int                row, col, y, i = 0;

    if (aV.cols_get() * aV.rows_get() != mrows * mcols)
        _error (pch_name, pch_errorSize);

    // create the structure
    coreMatrix_construct(mrows, mcols);
    
    switch (e_adominance) {
    case e_row:
        for (row = 0; row < mrows; row++)
            for (col = 0; col < mcols; col++) {
                switch (aV.matrix_type()) {
                case e_rowVector:
                    _mval (row, col) = aV.val (0, i++);
                    break;
                case e_columnVector:
                    _mval (row, col) = aV.val (i++, 0);
                    break;
                case e_matrix:
                case e_square:
                    break;
                }
            }
        break;
    case e_column:
        for (col = 0; col < mcols; col++)
            for (row = 0; row < mrows; row++) {
                switch (aV.matrix_type ()) {
                case e_rowVector:
                    _mval (row, col) = aV.val (0, i++);
                    break;
                case e_columnVector:
                    _mval (row, col) = aV.val (i++, 0);
                    break;
                case e_matrix:
                case e_square:
                    break;
                }
            }
    }
}

template<typename _CMDATA>
CMatrix<_CMDATA>::CMatrix(
        CMatrix<_CMDATA>&       aM,
        e_DOMINANCE             e_adominance,
        e_MATRIXTYPE            e_atype) {
    //
    // ARGS
    //  aM              in              Matrix whose data will be used to
    //                                          construct vector
    //  e_adominance    in              Will vector be built from row-first
    //                                          or col-first from this matrix
    //  e_atype         in              type of resultant vector (row or column)
    //
    // DESC
    //  Builds a vector from the passed matrix.
    //
    //  The `dominance' variable defines whether the vector is
    //  row or column dominant, i.e. the values from the vector
    //  are assumed organized as [row1, row2, row3, ... , rowN]
    //  or as [col1, col2, col3, ... , colN]
    //
    // POSTCONDITIONS
    //  o A *vector* is built.
    //
    // HISTORY
    // 01 April 2003
    //        o Templatization.
    //

    char*        pch_name         = "CMatrix (construct from matrix)";
    char*        pch_error         = "Resultant vector must be either row or column";
    int                row, col, y, i = 0;
    int                rows, cols;

    // create the structure
    switch (e_atype) {
    case e_rowVector:
        rows = 1;
        cols = aM.rows_get () * aM.cols_get ();
        break;
    case e_columnVector:
        cols = 1;
        rows = aM.rows_get () * aM.cols_get ();
        break;
    case e_matrix:
    case e_square:
        _error (pch_name, pch_error);
    }
    coreMatrix_construct(rows, cols);

    switch (e_adominance) {
    case e_row:
        for (row = 0; row < aM.rows_get (); row++)
            for (col = 0; col < aM.cols_get (); col++) {
                switch (e_atype) {
                case e_columnVector:
                    val (i++, 0) = aM.val (row, col);
                    break;
                case e_rowVector:
                    val (0, i++) = aM.val (row, col);
                    break;
                case e_matrix:
                case e_square:
                    break;
                }
            }
        break;
    case e_column:
        for (col = 0; col < aM.cols_get (); col++)
            for (row = 0; row < aM.rows_get (); row++) {
                switch (e_atype) {
                case e_columnVector:
                    val (i++, 0) = aM.val (row, col);
                    break;
                case e_rowVector:
                    val (0, i++) = aM.val (row, col);
                    break;
                case e_matrix:
                case e_square:
                    break;
                }
            }
        break;
    }
}

template<typename _CMDATA>
CMatrix<_CMDATA>::CMatrix (
        char            *flag,
        int             dimension)
{
    //
    // AEGS
    //  *flag           in              Currently only supports 'I'
    //                                          as flag
    //  dimension       in              size of identity matrix to
    //                                          create
    //
    // DESC
    //  Constructor.
    //
    //  Creates an identity matrix, of size (dimension X dimension).
    //
    // POSTCONDITIONS
    //  o An identity matrix is returned
    //
    // HISTORY
    // 31 March 2003
    //        o Templatization.
    //
    
    if (flag[0] != 'I')
        _error ("To create an identity matrix, use",
                 "CMatrix(\"I\", dimension)");
    
    coreMatrix_construct(dimension, dimension);
        
    for (int i = 0; i < dimension; i++)
        for (int j = 0; j < dimension; j++)
            _mval (i, j) = (i == j ? (_CMDATA)1 : (_CMDATA)0);
}

// error message used when failing to read a "non standard" matrix file
static char
    nonstandard[] =
    "\nis a 'non-standard' MATRIX file. A 'standard' file must\n"
    "start with the dimensions of the matrix, i.e.:\n"
    "\t rows=rr columns=cc\nor abbreviated:\n\tr=rr c=cc\n"
    "Note: rows appear before columns, and chars are lowercase.\n"
    "Comments follow '#' signs to the end of the line, data follows :::::::\n";

template<typename _CMDATA>
CMatrix<_CMDATA>::CMatrix (
        char            *matfile)
{
    //
    // ARGS
    //  matfile         in              name of file containing the
    //                                          matrix to be parsed
    //
    // DESC
    //  Constructor.
    //
    //  Creates a matrix from the specified filename.
    //
    // PRECONDITIONS
    //  o The information in *matfile must be in the format that
    //    this method expects!
    //
    // POSTCONDITIONS
    //  o Matrix is created.
    //
    // HISTORY
    // 31 March 2003
    //        o Templatization.
    //
    
    const int        BUFSIZE = 120;
    FILE*        infile;
    char            buf[BUFSIZE], *cp, *cp2;
    int                rfound = 0, cfound = 0, colonsfound = 0;
    int                rows = 0, cols = 0;

    if ((infile = fopen (matfile, "r")) == 0)
        _error ("File open error! Cannot open file:-", matfile);
    while (fgets (buf, BUFSIZE, infile)) {        // parse file initialization header
        // for each header line remove comments
        if ((cp = strpbrk (buf, "#")) != NULL)        // if comment found,
            *cp = '\0';                // terminate string at comment
        if ((cp = strpbrk (buf, "r")) != NULL)
            if ((cp2 = strpbrk (cp, "=")) != NULL)
                if ((cp = strpbrk (cp2, "0123456789")) != NULL) {
                    rows = atoi (cp);
                    rfound++;
                }
        if ((cp = strpbrk (buf, "c")) != NULL)
            if ((cp2 = strpbrk (cp, "=")) != NULL)
                if ((cp = strpbrk (cp2, "0123456789")) != NULL) {
                    cols = atoi (cp);
                    cfound++;
                }
        if (strstr (buf, ":::::::") != NULL) {
            colonsfound++;
            break;                // out of while loop
        }
    }
    if (!rfound || !cfound || !colonsfound) {
        fprintf (stderr, "%s %s", matfile, nonstandard);
        exit (1);
    }

    coreMatrix_construct(rows, cols);
            
    for (int row = 0; row < ps_mat->rows; row++)
        for (int col = 0; col < ps_mat->cols; col++) {
            char                info[20];
            fscanf (infile, "%s", info);
            _mval (row, col) = (_CMDATA) atof (info);
            if (ferror (infile))
                _error ("Input file error! Problem accessing file:-",
                         matfile);
        }
    fclose (infile);
}

template <typename _CMDATA>
CMatrix<_CMDATA>::CMatrix(const CMatrix<_CMDATA> & z)
{
    //
    // DESC
    //  Copy constructor - simple pointer based.
    //
    // POSTCONDITIONS
    //  o NB!! NB!! NB!!
    //    The copy constructor merely increases the reference count of the
    //    "source" data, and directs the target to point to the source.
    //    This is *not* a deepcopy!. Although allowing for fast copies between
    //    matrices, it can potentially suffer from problems relating to scope 
    //    local variable variable problems. If used in function arguments
    //    or as returns out of functions, then it is not really a problem.
    //
    // HISTORY
    // 27 January 2004
    //  o refcnt_inc()
    //

//    z.ps_mat->refcnt++;                // adding another reference
    z.refcnt_inc();                     // adding another reference
    ps_mat = z.ps_mat;                  // point to the new matstruct
}

template <typename _CMDATA>
void
CMatrix<_CMDATA>::coreMatrix_destruct() {
    //
    // DESC
    //        Explicitly destroys the core data stuctures of a gsl-aware *matrix*
    //        object.
    //
    // HISTORY
    // 20 March 2003
    //        o Initial design and coding.
    //
    // 28 January 2004
    //  o Corrected non-GSL destruction using delete [] properly.
    //

    // Use RTTI to determine datatype
    const type_info&        Tinfo = typeid(*this);                        // RTTI data

#ifdef GSL_USE
    if(Tinfo==Ginfo_intMatrix)
        gsl_matrix_int_free(ps_mat->ps_gsl->matrix_int);
    if(Tinfo==Ginfo_floatMatrix)
        gsl_matrix_float_free(ps_mat->ps_gsl->matrix_float);
    if(Tinfo==Ginfo_doubleMatrix)
        gsl_matrix_free(ps_mat->ps_gsl->matrix_double);

    delete ps_mat->ps_gsl;
#else
    for (int y = 0; y < ps_mat->rows; y++)
        delete [] ps_mat->data[y];
    delete [] ps_mat->data;
#endif
    delete ps_mat;

}

template <typename _CMDATA>
CMatrix<_CMDATA>::~CMatrix () {
    //
    // DESC
    //  Destructor.
    //
    // POSTCONDITIONS
    //  Decrements reference count. If counter reaches zero, frees allocated
    //  memory.
    //
    // HISTORY
    // 27 January 2004
    //  o Added refcnt_dec()
    //

    //if (--ps_mat->refcnt == 0)
    if(refcnt_dec()==0)
             coreMatrix_destruct();
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::reinitialise (
        const _CMDATA   *apf_data,
        int             a_elements)
{
    //
    // ARGS
    //  apf_data        in              array of new values
    //  a_elements      in              (optional) number of elements in
    //                                          array. Ignored
    //                                          if zero.
    //
    // DESC
    //  Overwrites the internal data of the matrix with values
    //  copied from the passed array.
    //
    //  If a_elements is non-zero and positive, it is assumed to
    //  indicate the number of elements in the passed array; otherwise
    //  it is assumed that the number of elements is the same as size1D().
    //
    //  It is a very good idea (and much safer) to pass a valid element
    //  value in a_elements. Basically, you are simply pointing to an
    //  unstructured block of memory and "sucking" it into this matrix - the
    //  a_elements allows at least a basic error check in case there is
    //  a mismatch between the size of this matrix and the block of
    //  memory.
    //
    // PRECONDITIONS
    //  o Make sure that the number of elements passed in array are equal
    //    to the number of elements in matrix itself.
    //
    // POSTCONDITIONS
    //  o Matrix reinitialised.
    //
    // HISTORY
    //  17 January 2001
    //  o Initial design, coding, and testing
    //
    // 31 March 2003
    //         o Templatization.
    //

    int                elements = size1D ();
    if (a_elements && (a_elements != elements))
        _error ("reinitialise", "number of elements mismatch");

    int                dataptr = 0;
    for (int i = 0; i < rows_get (); i++)
        for (int j = 0; j < cols_get (); j++)
            val (i, j) = apf_data[dataptr++];
    return *this;
}

template<typename _CMDATA>
int
CMatrix<_CMDATA>::positionVectorise(
        CMatrix<_CMDATA>&        aM_dim ) {
    //
    // ARGS
    //  aM_dim           in              vector of dimensions and sizes
    //
    // DESC
    //  Given a vector describing size and extent of a space, this
    //  method takes its current *data as an index into that space
    //  and returns an offset into a 1D concatenation of M_dim.
    //
    //  Used really for a N-Dimensional (where N is the size of aM_dim)
    //  to 1 dimensional offset. This is conceptually similar to doing
    //  base radix transforms on numbers, where aM_dim is a number in
    //  a particular radix, and the returned value is the base ten
    //  equivalent of this number.
    //
    // PRECONDITIONS
    //  o *this matrix's val() _must_ be within the space defined by
    //    M_dim.
    //
    // POSTCONDITIONS
    //  o returns an integer offset into a 1D "vectorised" version of
    //    the N dimensional space.
    //
    // HISTORY
    //  19 January 2001
    //  o Initial design and coding.
    //
    // 31 March 2003
    //         o Templatization.
    //

    char*        pch_name                 = "positionVectorise";
    int                i, j, partialsum         = 0;
    int                 partialprod                 = 1;

    if (size1D () != aM_dim.cols_get ())
        _error (pch_name, "dimension / internal size mismatch");

    for (i = aM_dim.cols_get () - 1; i >= 0; i--) {
        partialprod = (int) val (0, i);
        for (j = i - 1; j >= 0; j--) {
            partialprod *= (int) aM_dim.val (0, j);
        }
        partialsum += partialprod;
    }
    return partialsum;
}

template<typename _CMDATA>
void
CMatrix<_CMDATA>::randomize (
        _CMDATA                lower,
        _CMDATA         upper,
        bool                 b_asInt)
{
    //
    // ARGS
    //        lower                in                lower bound of range
    //        upper                in                upper bound of range
    //         b_asInt                in                (optional) integer or _CMDATA values
    //
    // DESC
    //        This method fills the contents of *this matrix with
    //        random values between lower and upper (inclusive).
    //
    //        If b_asInt is true, the random values are truncated to
    //        form integers.
    //
    // PRECONDITIONS
    //        o Matrix must exist
    //        o upper must be bigger than or equal to lower
    //
    // POSTCONDITIONS
    //        o Matrix randomly filled with (asInt==0 ? integer : _CMDATA)
    //          values between lower and upper.
    //  o If lower == upper, i.e. range = 0, all matrix values
    //          are set to lower.
    //
    // HISTORY
    //  13 September 2001
    //  o Added code for range == 0.
    //
    //  12 February 2002
    //  o Removed random seed setting ability
    //
    //         31 March 2003
    //         o Templatization.
    //

    char*        pch_proc = "randomize";
    _CMDATA        f_range = upper - lower;
    _CMDATA        f_normalized, f_num;
    int                num;
    int                row, col;

    for (row = 0; row < ps_mat->rows; row++) {
        for (col = 0; col < ps_mat->cols; col++) {
            if (f_range > 0) {
                    if (!b_asInt) {
                    num = rand ();
                     f_normalized = (_CMDATA) num / RAND_MAX;
                    if ((f_normalized < 0) || (f_normalized > 1))
                            _error (pch_proc, "Out of range random int normalized");
                    f_num = (f_normalized * f_range + lower);
                    val (row, col) = (_CMDATA)f_num;
                    } else
                    val (row, col) = (int)        ((_CMDATA) rand () /
                                                ((_CMDATA) (RAND_MAX) / f_range)) +
                                                    lower;
            } else
                    val(row, col) = lower;
        }
    }
}

template<typename _CMDATA>
void
CMatrix<_CMDATA>::randomly_fill (
        _CMDATA         af_fillValue,
        _CMDATA         af_density) {
    //
    // ARGS
    //  af_fillValue    in      value randomly filled into matrix
    //  af_density      in      density of fill_value
    //                          (0.. 1)
    //
    // DESC
    //  Fills a matrix with fill_value at given density.
    //
    // NOTE
    //  This routine is not guaranteed to fill a matrix with exactly
    //  the required density (although I suppose in the limit as the
    //  matrix size tends to infinity it would). Simply stated, each
    //  element of the matrix is cycled through, and its particular
    //  value is changed to fill_value with probability of `density'.
    //
    //  Some rudimentary error checking is done on bounds of
    //  density. I should actually properly implement exception
    //  handling... :-P
    //
    // HISTORY
    //  6 September 2000
    //  o Added random seed parameters
    //
    //  12 February 2002
    //  o Removed random seeding ability
    //
    //         31 March 2003
    //         o Templatization.
    //

    _CMDATA        af_randValue;
    int                row, col;

    if (af_density < 0)
        af_density = 0;
    if (af_density > 1)
        af_density = 1;

    for (row = 0; row < ps_mat->rows; row++)
        for (col = 0; col < ps_mat->cols; col++) {
            af_randValue = (_CMDATA) rand () / (_CMDATA) (RAND_MAX);
            if (af_randValue >= 1 - af_density)
                val (row, col) = af_fillValue;
        }
}

template<typename _CMDATA>
int
CMatrix<_CMDATA>::hardlimit (
        _CMDATA         low,
        _CMDATA         high,
        _CMDATA         tolow,
        _CMDATA         tohigh) {
    //
    // ARGS
    //  low             in              the lower boundary
    //  high            in              the upper boundary
    //  tolow           in              values less than `low' are set to `tolow'
    //  tohigh          in              values greater than `high' are set to `tohigh'
    //
    // DESC
    //  This method implements a hardlimiting filter. A target matrix
    //  is scanned for entries that are out of bounds
    //  (entry < low || entry > high) and limits them to the
    //  passed bounds.
    //
    // POSTCONDITIONS
    //  o return the number of `out of bounds' cells in the method name.
    //
    // NOTE
    //  o The comparison for lower is <= low, and for higher is > high
    //
    // HISTORY
    //  05 October 2000
    //  o Added equality in lower comparison
    //
    //         31 March 2003
    //         o Templatization.
    //

    int                row, col;
    int                outOfBounds         = 0;

    for (row = 0; row < ps_mat->rows; row++)
        for (col = 0; col < ps_mat->cols; col++) {
            if (val (row, col) <= low) {
                outOfBounds++;
                val (row, col) = tolow;
            }

            if (val (row, col) > high) {
                outOfBounds++;
                val (row, col) = tohigh;
            }
        }
    return outOfBounds;
}

template<typename _CMDATA>
void
CMatrix<_CMDATA>::quantize(
        _CMDATA                 f_grain) {
    //
    // ARGS
    //  f_grain                 in      grain factor (0.1, 0.01, 0.001) etc.
    //
    // DESC
    //  Quantize matrix contents. Each element is quantized to grain size
    //  f_grain.
    //
    // HISTORY
    //         31 March 2003
    //         o Templatization.
    //

    int                row, col;
    _CMDATA        f_num, f_invGrain;

    f_invGrain = 1 / f_grain;
    for (row = 0; row < ps_mat->rows; row++)
        for (col = 0; col < ps_mat->cols; col++) {
            f_num = val (row, col);
            f_num *= f_invGrain;
            f_num = (int) f_num;
            f_num /= f_invGrain;
            val (row, col) = f_num;
        }
}

template<typename _CMDATA>
bool 
CMatrix<_CMDATA>::quantizeIndex(
        CMatrix<_CMDATA>&                M_lookup,
        CMatrix<_CMDATA>*&               pM_index
) {
    //
    // ARGS
    //  M_lookup                in      lookup table of
    //                                  quantized values
    //  pM_index                out     index containing
    //                                  the "quantization"
    //                                  of the base vector
    //
    // DESC
    //  Assuming a continuous variable *row* vector, and given
    //  a lookup table, this method returns an index into the
    //  lookup table.
    //
    //  The lookup table must have the same number of rows as
    //  the underlying vector, and each column of the lookup
    //  table defines the range of quantized values for each
    //  dimension component of this underlying vector. Hmmm...
    //  if I try real hard, maybe I can make it more obtuse than
    //  even this.
    //
    // PRECONDITIONS
    //  o *this must be a *row* vector containing continuous values
    //  o M_lookup must have same number of rows as *this
    //
    // POSTCONDITIONS
    //  o pM_index defines the "position" in M_lookup that is "closest"
    //   to the continuous value of *this.
    //
    // RETURN
    //  b_outOfBounds           out     error value of method:
    //                                  0:      all values were "cleanly"
    //                                          quantised within M_lookup
    //                                  1:      one (or more) values in *this
    //                                          were outOfBounds and
    //                                          hardlimited to the closest
    //                                          boundary
    //
    // HISTORY
    //  08 Aug 2000
    //  o Initial design and coding
    //
    //  04 April 2001
    //  o Changed the meaning of the return value from b_error to
    //    b_outOfBounds
    //
    //         31 March 2003
    //         o Templatization.
    //

    static CMatrix<_CMDATA>*        pM_colIndex         = new CMatrix<_CMDATA>(1, 2, (_CMDATA) 0.0);
    static CMatrix<_CMDATA>*        pV_col                 = new CMatrix<_CMDATA>(1, 1, (_CMDATA) 0.0);
    bool   b_outOfBounds                         = false;

    if (!is_vectorRow ())
        _error ("quantizedIndex", "Base must be a *row* vector", 1);
    if (cols_get () != M_lookup.cols_get ())
        _error ("quantizedIndex", "column amount mismatch", 1);

    if (!pM_index->compatible (*this)) {
        delete                    pM_index;
        pM_index         = new CMatrix<_CMDATA> (1, this->cols_get (), (_CMDATA) 0.0);
    }

    for (int col = 0; col < cols_get (); col++) {
        M_lookup.col_remove (pV_col, col);
        b_outOfBounds |= pV_col->find_quantized (val (0, col), pM_colIndex);
        pM_index->val (0, col) = pM_colIndex->val (0, 0);
    }
    return b_outOfBounds;
}

template<typename _CMDATA>
void
CMatrix<_CMDATA>::quantize_linearly(
        _CMDATA                 af_start,
        _CMDATA                 af_stop)
{
    //
    // ARGS
    //  af_start                in      start value
    //  af_stop                 in      stop value
    //
    // DESC
    //  This method accepts a vector and quantizes the values
    //  it contains between af_start and af_stop. The quantum
    //  step size is a linear function of the number of elements
    //  in the vector.
    //
    // PRECONDITIONS
    //        o column vector
    //        o stop > start
    //
    // POSTCONDITIONS
    //        o The left-most value of the vector  = af_start
    //        o The right-most value of the vector = af_stop
    //
    // HISTORY
    //  04 Aug 2000
    //  o Initial design and coding. For the sake of getting
    //    a working version out of the door, the initial version
    //    of this method assumes an underlying column vector.
    //
    //         31 March 2003
    //         o Templatization.
    //

    if (!is_vectorColumn())
        _error ("quantize_linearly", "Needs a column vector", 1);

    if (af_stop < af_start)
        _error ("quantize_linearly", "Invalid range", 1);

    _CMDATA        f_range = af_stop - af_start;
    _CMDATA        f_quantum = f_range / (rows_get ());

    val (0, 0) = af_start;
    for (int row = 1; row < rows_get (); row++)
        val (row, 0) = val (row - 1, 0) + f_quantum;
}

template<typename _CMDATA>
void
CMatrix<_CMDATA>::quantize_linearlyShifted(
        _CMDATA                        af_start,
        _CMDATA                        af_stop
) {
    //
    // ARGS
    //         af_start                in      start value
    //         af_stop                        in      stop value
    //
    // DESC
    //         This method accepts a vector and quantizes the values
    //         it contains between af_start and af_stop. The quantum
    //         step size is a linear function of the number of elements
    //         in the vector.
    //
    //  The difference between this method and quantize_linearly()
    //        is that the step between val(0) and val(1) is one-half a
    //        quantum, such that the rightmost lookup value is one-half
    //        a quantum less than af_stop.
    //
    //        This means that using quantiseIndex() with a vector returned
    //        from this method classifies values close to the right edge
    //        correctly.
    //
    // PRECONDITIONS
    //        o column vector
    //        o stop > start
    //
    // POSTCONDITIONS
    //        The left-most value of the vector  = af_start
    //        The right-most value of the vector = af_stop - 0.5 * quantum
    //        At least three quantisation levels
    //
    // HISTORY
    // 04 Aug 2000
    // o Initial design and coding. For the sake of getting
    //   a working version out of the door, the initial version
    //   of this method assumes an underlying column vector.
    //
    // 12 September 2001
    //        Adaptation from quantize_linearly()
    //
    //         31 March 2003
    //         o Templatization.
    //
    
    if (!is_vectorColumn())
        _error ("quantize_linearly", "Needs a column vector", 1);

    if (af_stop < af_start)
        _error ("quantize_linearly", "Invalid range", 1);

    _CMDATA        f_range = af_stop - af_start;
    _CMDATA        f_quantum = f_range / (rows_get () - 1);

    val(0, 0) = af_start;
    val(1, 0) = af_start + f_quantum / 2;
    for (int row = 2; row < rows_get (); row++)
        val(row, 0) = val (row - 1, 0) + f_quantum;
}

template<typename _CMDATA>
void
CMatrix<_CMDATA>::delta_add(
        _CMDATA         af_delta,
        _CMDATA         af_lower,
        _CMDATA         af_upper)
{
    //
    // ARGS
    //  af_delta        in      delta value to add to appropriate
    //                           cells
    //  af_lower        in      lower bound of cells that receive delta
    //  af_upper        in      upper bound of cells that receive delta
    //
    // DESC
    //  This routine merely adds the `af_delta' value to cells in the
    //  matrix that satisfy af_lower <= cell <= af_upper.
    //
    // HISTORY
    //  05 June 2000
    //  o Initial design and coding.
    //
    //         31 March 2003
    //         o Templatization.
    //

    for (int row = 0; row < rows_get (); row++)
        for (int col = 0; col < cols_get (); col++) {
            if ((val (row, col) >= af_lower) && (val (row, col) <= af_upper))
                val (row, col) += af_delta;
        }
}

//////---------->
////// Access and set internal values
//////---------->

template<typename _CMDATA>
_CMDATA &        
CMatrix<_CMDATA>::_mval(  
        int        row,
        int        col)  const {
    //
    // ARGS
    //        row                in                row to access
    //        col                in                col to access
    //
    // DESC
    //         Element selection: can be used to read or write.
    //
    // PRECONDITIONS
    //         This routine does *not* use range checking! It is an
    //        interal access routine that is used by class methods
    //        and is not available for "public" use.
    //
    // HISTORY
    // 18 March 2003
    //        o Moved from header file into main source body as part 
    //          of gsl integration.
    //

    const type_info&        Tinfo = typeid(*this);                        // RTTI data

#ifndef GSL_USE
    return(ps_mat->data[row][col]);
#else
    if(Tinfo==Ginfo_intMatrix)
        return( (_CMDATA&)*gsl_matrix_int_ptr(ps_mat->ps_gsl->matrix_int, row, col));
    if(Tinfo==Ginfo_floatMatrix)
        return( (_CMDATA&)*gsl_matrix_float_ptr(ps_mat->ps_gsl->matrix_float, row, col));
    if(Tinfo==Ginfo_doubleMatrix)
        return( (_CMDATA&)*gsl_matrix_ptr(ps_mat->ps_gsl->matrix_double, row, col));
#endif
}

template<typename _CMDATA>
_CMDATA &
CMatrix<_CMDATA>::val(
        int         row,
        int         col)
{
    //
    // ARGS
    //        row                in                row to access
    //        col                in                col to access
    //
    // DESC
    //         Element selection: can be used to read or write.
    //
    // PRECONDITIONS
    //        o Range checking can be expensive, particularly in tightly
    //          nested iterative loops. By defining CMATRIX_NORANGE range
    //          checking is disabled. This will result in tremendous performance
    //          boost, but at the risk of a single out of bound access
    //          killing the entire process.
    //
    // HISTORY
    //  3 Decemeber 2001
    //        o Added compiler directives to turn off range checking.
    //
    // 30 January 2004
    //  o Added GSL support for int/float
    //

#ifndef CMATRIX_NORANGE
    char *
        pch_name = "val";
    char *
        pch_errormsg = "Addressing error. Index out of range";
    if (!(row >= 0 && row < ps_mat->rows && col >= 0 && col < ps_mat->cols))
        _error (pch_name, pch_errormsg);
    else
#endif
    return (_mval (row, col));
}

template<typename _CMDATA>
_CMDATA &
CMatrix<_CMDATA>::operator()(
        int         row,
        int         col)
{
    //
    // ARGS
    //        row                in                row to access
    //        col                in                col to access
    //
    // DESC
    //         Element selection: can be used to read or write.
    //
    // PRECONDITIONS
    //        o Range checking can be expensive, particularly in tightly
    //          nested iterative loops. By defining CMATRIX_NORANGE range
    //          checking is disabled. This will result in tremendous performance
    //          boost, but at the risk of a single out of bound access
    //          killing the entire process.
    //
    // HISTORY
    //  3 Decemeber 2001
    //        o Added compiler directives to turn off range checking.
    //

#ifndef CMATRIX_NORANGE
    char *
        pch_name = "val";
    char *
        pch_errormsg = "Addressing error. Index out of range";
    if (!(row >= 0 && row < ps_mat->rows && col >= 0 && col < ps_mat->cols))
        _error (pch_name, pch_errormsg);
    else
#endif
    return (_mval (row, col));
}

template<typename _CMDATA>
_CMDATA &
CMatrix<_CMDATA>::val(
        int         element)
{
    //
    // ARGS
    //  element            in                element to access
    //
    // DESC
    //  Element selection: can be used to read or write. This method
    //	is for use with *vectors*.
    //
    // PRECONDITIONS
    //	o Only use with *vectors*, not matrices. 
    //  o Range checking can be expensive, particularly in tightly
    //    nested iterative loops. By defining CMATRIX_NORANGE range
    //    checking is disabled. This will result in tremendous performance
    //    boost, but at the risk of a single out of bound access
    //    killing the entire process.
    //
    // HISTORY
    //  3 Decemeber 2001
    //  o Added compiler directives to turn off range checking.
    //
    // 24 September 2003
    //	o Designed for *vectors*.
    //

#ifndef CMATRIX_NORANGE
    char* pch_name = "val";
    char* pch_errormsg = "Addressing error. Index out of range";
    if ((!is_vectorRow() && !is_vectorColumn()))
        _error (pch_name, pch_errormsg);
#endif
    if(is_vectorRow())
	return (val (0, element));
    if(is_vectorColumn())
	return (val(element, 0));
}

template<typename _CMDATA>
_CMDATA &
CMatrix<_CMDATA>::operator()(
        int         element)
{
    //
    // ARGS
    //  element            in                element to access
    //
    // DESC
    //  Element selection: can be used to read or write. This method
    //	is for use with *vectors*.
    //
    // PRECONDITIONS
    //	o Only use with *vectors*, not matrices.
    //  o Range checking can be expensive, particularly in tightly
    //    nested iterative loops. By defining CMATRIX_NORANGE range
    //    checking is disabled. This will result in tremendous performance
    //    boost, but at the risk of a single out of bound access
    //    killing the entire process.
    //
    // HISTORY
    //  3 Decemeber 2001
    //  o Added compiler directives to turn off range checking.
    //
    // 24 September 2003
    //	o Designed for *vectors*.
    //

#ifndef CMATRIX_NORANGE
    char* pch_name = "val";
    char* pch_errormsg = "Addressing error. Index out of range";
    if ((!is_vectorRow() && !is_vectorColumn()))
        _error (pch_name, pch_errormsg);
#endif
    if(is_vectorRow())
	return (val (0, element));
    if(is_vectorColumn())
	return (val(element, 0));
}

//////---------->
////// Output routines
//////---------->

template<typename _CMDATA>
void
CMatrix<_CMDATA>::print (
        char*                   apch_msg        /*= "ans"       */,
        int                     a_precision     /*= 6           */,
        int                     a_width         /*= 12          */,
        ios::fmtflags           a_userFormat    /*= 0           */,
        int                     a_marginOffset  /*= 0           */,
        char                    ach_left        /*= (char) 0    */,
        int                     a_tabOffset     /*= 0           */,
        bool                    ab_fancy        /*= 1           */) {
    //
    // ARGS
    //  apch_msg                in      message that prepends matrix dump
    //  a_precision             in      the precision of numerical output
    //  a_width                 in      the width of the numerical field
    //  a_userFormat            in      if non-zero, set stream format to arg
    //  a_marginOffset          in      the amount of tab spaces prepending ch_left
    //  ach_left                in      also somewhat primitive. If non-zero, will
    //                                          print (char) left_char as left most character
    //                                          of each new line. Useful for when the matrix is
    //                                          part of a larger output dump
    //  a_tabOffset             in      the amount of tab spaces prepending dump
    //  ab_fancy                in      still rather primitive. If false, will print
    //                                          *only* the numbers, nothing else
    //
    // DESC
    //        Print a matrix in a variety of ways. The a_userFormat parameter allows the
    //        user to send extra formatting flags to the underlying stream. Consult a
    //        C++ manual for flags.
    //
    // NOTE
    //         The `marginOffset' and `tabOffset' are used in conjunction
    //         with ch_left as follows:
    //
    //                 [<--- marginOffset --->]<ch_left>[<--- tabOffset --->]
    //
    // HISTORY
    //  18 February 2002
    //        o Changed all c-style printf's to c++-style cout
    //

    ios::fmtflags
            streamFlags        =  cout.flags();        // get current stream flags
    char        pch_tab[64]        = "";
    char        pch_margin[64]        = "";
    char        pch_indent[64]        = "";
    int                i;

    for (i = 0; i < a_marginOffset; i++)
        strcat (pch_margin, "\t");

    for (i = 0; i < a_tabOffset; i++)
        strcat (pch_tab, "\t");

    sprintf (pch_indent, "%s%c%s", pch_margin, ach_left, pch_tab);

    if (ab_fancy) {
        if (ps_mat->rows > 1)
            cout << pch_tab << apch_msg << " = [" << endl << pch_indent;
        else
            cout << pch_tab << apch_msg << " = [ ";
    }
    for (int row = 0; row < ps_mat->rows; row++) {
        for (int col = 0; col < ps_mat->cols; col++) {
//            cout.width(a_width);
            cout.precision(a_precision);
            cout.setf(ios::left);
            cout.setf(ios::skipws);
            if(a_userFormat)
                    cout.flags(a_userFormat);
            cout << _mval(row, col) << "\t";
            //if(ab_fancy) cout << "| ";
        }
        if (ps_mat->rows > 1)
            cout << endl << pch_indent;
    }
    if (ab_fancy)
        if (ps_mat->rows > 1)
            cout << "]" << endl << pch_indent;
        else
            cout << "] ";
    cout.flags(streamFlags);                // Restore stream flags
}

template<typename _CMDATA>
string
CMatrix<_CMDATA>::sprint (
        char*                   apch_msg        /*= "ans"       */,
        int                     a_precision     /*= 6           */,
        int                     a_width         /*= 12          */,
        ios::fmtflags           a_userFormat    /*= 0           */,
        int                     a_marginOffset  /*= 0           */,
        char                    ach_left        /*= (char) 0    */,
        int                     a_tabOffset     /*= 0           */,
        bool                    ab_fancy        /*= 1           */)
{
    //
    // ARGS
    //  apch_msg                in      message that prepends matrix dump
    //  a_precision             in      the precision of numerical output
    //  a_width                 in      the width of the numerical field
    //  a_userFormat            in      if non-zero, set stream format to arg
    //  a_marginOffset          in      the amount of tab spaces prepending ch_left
    //  ach_left                in      also somewhat primitive. If non-zero, will
    //                                          print (char) left_char as left most character
    //                                          of each new line. Useful for when the matrix is
    //                                          part of a larger output dump
    //  a_tabOffset             in      the amount of tab spaces prepending dump
    //  ab_fancy                in      still rather primitive. If false, will print
    //                                          *only* the numbers, nothing else
    // DESC
    //        Print a matrix in a variety of ways. The a_userFormat parameter allows the
    //        user to send extra formatting flags to the underlying stream. Consult a
    //        C++ manual for flags.
    //
    // NOTE
    //  The `marginOffset' and `tabOffset' are used in conjunction
    //  with ch_left as follows:
    //
    //  [<--- marginOffset --->]<ch_left>[<--- tabOffset --->]
    //
    // HISTORY
    // 18 February 2002
    //        o Changed all c-style printf's to c++-style cout
    //
    // 01 April 2003
    //        o Templatization.
    //

    stringstream        sstream("");

    ios::fmtflags            streamFlags        = cout.flags();        // get current stream flags
    char                pch_tab[64]        = "";
    char                pch_margin[64]        = "";
    char                pch_indent[64]        = "";
    int                        i;

    for (i = 0; i < a_marginOffset; i++)
        strcat (pch_margin, "\t");

    for (i = 0; i < a_tabOffset; i++)
        strcat (pch_tab, "\t");

    sprintf (pch_indent, "%s%c%s", pch_margin, ach_left, pch_tab);

    if (ab_fancy) {
        if (ps_mat->rows > 1)
            sstream << pch_tab << apch_msg << " = [" << endl << pch_indent;
        else
            sstream << pch_tab << apch_msg << " = [ ";
    }
    for (int row = 0; row < ps_mat->rows; row++) {
        for (int col = 0; col < ps_mat->cols; col++) {
//            sstream.width(a_width);
            sstream.precision(a_precision);
            sstream.setf(ios::internal);
            sstream.setf(ios::skipws);
            if(a_userFormat)
                    sstream.flags(a_userFormat);
            sstream << _mval(row, col) << "\t";
        }
        if (ps_mat->rows > 1)
            sstream << endl << pch_indent;
    }
    if (ab_fancy)
        if (ps_mat->rows > 1)
            sstream << "]" << endl << pch_indent;
        else
            sstream << "] ";
    sstream.flags(streamFlags);                // Restore stream flags
    return(sstream.str());
}

template<typename _CMDATA>
void
CMatrix<_CMDATA>::fprint(
        string                  astr_filename,
        char*                   apch_msg        /*= ""          */,
        int                     a_precision     /*= 6           */,
        int                     a_width         /*= 12          */,
        ios::fmtflags           a_userFormat    /*= 0           */
) {
    //
    // ARGS
    //  astr_filename           in              filename to write matrix to
    //  apch_msg                in              message that prepends matrix dump
    //  a_precision             in              the precision of numerical output
    //  a_width                 in              the width of the numerical field
    //  a_userFormat            in              if non-zero, set stream format to arg
    //
    // DESC
    //        Print a matrix to file in a variety of ways. The a_userFormat parameter allows
    //        the user to send extra formatting flags to the underlying stream. Consult a
    //        C++ manual for flags.
    //
    // HISTORY
    // 18 February 2002
    //        o Changed all c-style printf's to c++-style cout
    //
    // 01 April 2003
    //        o Templatization.
    //

    char*        pch_proc        = "fprint";

    ofstream        sout(astr_filename.c_str());
    if(!sout)
            _error(pch_proc, "Could not create output file");

    time_t         time_now        = time(NULL);
    string        str_time        = ctime(&time_now);
    // strip trailing \n
    str_time.erase(str_time.length()-1);
    string        str_hostname(getenv("HOSTNAME"));
    string        str_user(getenv("USER"));

    sout << "#"                                        << endl;
    sout << "# Standard CMatrix save file."        << endl;
    sout << "#        Created by\t"         << str_user         << endl;
    sout << "#        Date stamp:\t"         << str_time        << endl;
    sout << "#        Machine name:\t"<< str_hostname        << endl;
    if(strlen(apch_msg))
    sout << "#        Matrix name:\t"        << apch_msg        << endl;
    sout << "#"                                        << endl;
    sout << "rows="         << rows_get();
    sout << " columns=" << cols_get()                 << endl;
    sout << ":::::::"                                << endl;

    string        str_data        = sprint("",
                                            a_precision,
                                            a_width,
                                            a_userFormat,
                                            0, 0, 0,
                                            0);

    sout << str_data;
    sout.close();
}

//////---------->
////// Matrix information routines
//////---------->

template<typename _CMDATA>
bool
CMatrix<_CMDATA>::is_outOfBounds(
        CMatrix<_CMDATA>&                aM_boundaryUpper,
        CMatrix<_CMDATA>&                aM_boundaryLower,
        bool&                                ab_violate
) {
    //
    // ARGS
    //        aM_boundaryUpper        in                the upper boundary
    //        aM_boudnaryLower        in                the lower boundary
    //        ab_violate              out                if a boundary has been
    //                                                        violated, this
    //                                                        values describes
    //                                                        (true) upper
    //                                                        (false) lower
    // DESC
    //        This method checks if any elements in *this are larger than
    //        corresponding values in aM_boundary. If any are larger, true
    //        is returned.
    //
    // PRECONDITIONS
    //        o An element-by-element comparison is made - thus all matrices
    //          (this, boundaryLower, and boundaryUpper) must all have identical
    //          sizes.
    //          - boundaryUpper violation:
    //                 element in *this is >= corresponding element in boundaryUpper
    //          - boundaryLower violation:
    //                element in *this is <= corresponding element in boundaryLower
    //
    // POSTCONDITIONS
    //        o if any values in *this are "out of bounds" with respect to
    //          the passed upper and lower boundary matrices, true is returned
    //          and the ab_violate defines which boundary has been violated.
    //
    // HISTORY
    //  04 October 2001
    //  o Initial design and coding.
    //
    // 09 October 2001
    //  o Ooops! Discovered an idiotic bug! Revamped method to cater for
    //          upper and lower boundary violations.
    //

    if(!compatible(aM_boundaryLower))
            _error("is_outOfBounds", "boundaryLower must be same size as *this", -1);
    if(!compatible(aM_boundaryUpper))
            _error("is_outOfBounds", "boundaryUpper must be same size as *this", -1);
    for(int row=0; row<rows_get(); row++)
            for(int col=0; col<cols_get(); col++)
                if(        val(row, col) >= aM_boundaryUpper.val(row, col) ||
                        val(row, col) <= aM_boundaryLower.val(row, col)) {
                        if(val(row, col) >= aM_boundaryUpper.val(row, col))
                           ab_violate = true;
                        else
                           ab_violate = false;
                        return true;
                }
    return false;
}

template<typename _CMDATA>
bool 
CMatrix<_CMDATA>::is_vector() const
{
    //
    // Checks if object is a vector
    //
    return (is_vectorRow() || is_vectorColumn());
}

template<typename _CMDATA>
bool
CMatrix<_CMDATA>::is_vectorColumn() const
{
    //
    // Checks if object is column vector
    //
    return ((ps_mat->cols == 1) && (ps_mat->rows > 1));
}

template<typename _CMDATA>
bool
CMatrix<_CMDATA>::is_vectorRow() const
{
    //
    // Checks if object is row vector
    //
    return ((ps_mat->rows == 1) && (ps_mat->cols > 1));
}

template<typename _CMDATA>
bool
CMatrix<_CMDATA>::is_square() const
{
    //
    // Checks if matrix is square
    //
    return (ps_mat->cols == ps_mat->rows);
}

template<typename _CMDATA>
int
CMatrix<_CMDATA>::size1D ()	const
{
    //
    // DESC
    //  Returns the 1D size of a matrix, i.e. the product of the rows
    //  and columns
    //
    // HISTORY
    //  14 November 2000
    //  o Initial design and coding.
    //

    return (rows_get () * cols_get ());
}

template<typename _CMDATA>
CMatrix<int>
CMatrix<_CMDATA>::size()	const
{
    //
    // DESC
    //  Returns a matrix containing the size of *this.
    //
    // HISTORY
    //  23 September 2000
    //  o Initial design and coding
    //
    //  14 November 2000
    //  o Possible dangling / lost pointer!
    //
    // 31 March 2003
    //         o Templatization.
    //

    CMatrix<int>        iM_size(1, 2);

    // rows
    iM_size.val (0, 0) = rows_get ();
    //cols
    iM_size.val (0, 1) = cols_get ();

    return iM_size;
}

template<typename _CMDATA>
e_MATRIXTYPE
CMatrix<_CMDATA>::matrix_type () {
    //
    // Returns the matrix type
    // of its argument
    //

    if (is_vectorColumn())
        return e_columnVector;
    if (is_vectorRow())
        return e_rowVector;
    if (is_square())
        return e_square;
    return e_matrix;
}

template<typename _CMDATA>
bool
CMatrix<_CMDATA>::compatible (
        int                     dim1,
        int                     dim2) {
    //
    // ARGS
    //  dim1, dim2              in              dimensions under consideration
    //
    // DESC
    //  Checks if base is compatible with passed arguments
    //
    // 26 March 2003
    //         o Templatization.
    //

    if (ps_mat->rows != dim1 || ps_mat->cols != dim2)
        return false;
    return true;
}

template<typename _CMDATA>
bool
CMatrix<_CMDATA>::compatible (
        const CMatrix&          M)
{
    //
    // ARGS
    //
    //  M                       in      Matrix under consideration
    //
    // DESC
    //  Checks if base is compatible with passed arguments
    //
    // 26 March 2003
    //         o Templatization.
    //

    int                dim1, dim2;

    dim1         = M.rows_get();
    dim2         = M.cols_get();
    if (ps_mat->rows != dim1 || ps_mat->cols != dim2)
        return false;
    return true;
}

//////---------->
////// Matrix structural manipulation routines
//////---------->

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::copy (
        const CMatrix&          source)
{
    // ARGS
    //  source                  in              matrix whose contents to copy
    //                                                  into *this
    //
    // DESC
    //  Explicit deepcopy.
    //
    // POSTCONDITIONS.
    //  o The data in source is copied over into *this. Space for this data
    //    is created as in necessary.
    //
    // HISTORY
    // 20 March 2003
    //        o coreMatrix_[construct/destruct] integration.
    //
    // 26 March 2003
    //         o Templatization.
    //
    
    int                row, col, y;
        
    if (!compatible (source)) {
            coreMatrix_destruct();
        coreMatrix_construct(source.rows_get(), source.cols_get());
    }
    for (row = 0; row < ps_mat->rows; row++)
        for (col = 0; col < ps_mat->cols; col++)
            _mval (row, col) = source._mval (row, col);
    return *this;
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::matrix_remove(
        CMatrix<_CMDATA>*&      apM,
        int                     a_trow,
        int                     a_tcol,
        int                     a_rows,
        int                     a_cols) {
    //
    // ARGS
    //  apM             in/out          structure to contain removed
    //                                          sub-matrix
    //  a_trow, a_tcol  in              top left coordinate in base matrix
    //                                          of submatrix
    //  a_rows, a_cols  in              relative to trow and tcol, number
    //                                          of rows and cols to extract
    //                                          into submatrix
    //
    // DESC
    //  Removes a submatrix from a base matrix
    //
    // PRECONDITIONS
    //  o 0 < a_trow < rows_get()
    //  o 0 < a_tcol < cols_get()
    //  o a_trow + a_rows <= rows_get() + 1
    //  o a_tcol + a_cols <= cols_get() + 1
    //
    // POSTCONDITIONS
    //  o Returns submatrix (in name) and pointer to submatrix (in argument list).
    //
    // HISTORY
    //  03 December 2000
    //  o Changed calling parameters to explicitly allow for
    //    holding memory structure.
    //
    //         31 March 2003
    //         o Templatization.
    //

    char*        pch_name = "remove_matrix";
    int                row, col;

    if (a_trow + a_rows > ps_mat->rows + 1 ||
        a_tcol + a_cols > ps_mat->cols + 1)
            _error (pch_name, "Specified target dimensions invalid!");

    if (!apM->compatible (a_rows, a_cols)) {
        delete        apM;
        apM         = new CMatrix<_CMDATA>(a_rows, a_cols, (_CMDATA) 0.0);
    }

    for (row = 0; row < a_rows; row++)
        for (col = 0; col < a_cols; col++) {
            apM->val (row, col) = val (row + a_trow, col + a_tcol);
        }
    return *apM;
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::matrix_replace(
        int                     trow,
        int                     tcol,
        CMatrix<_CMDATA>&        aM_replacement)
{
    //
    // ARGS
    //  trow, tcol      in      top left coordinate in base matrix
    //                                  where replacement will be added.
    //  aM_replacement  in      replacement matrix
    //
    // DESC
    //  Overwrites part of a matrix with another matrix.
    //
    // 01 April 2003
    //         o Templatization.
    //

    char*        pch_name         = "matrix_replace";
    char*        pch_errorSize         = "Replacement cannot fit into base.";
    char*        pch_errorCoord         = "Specified insert point is out of range.";
    int                row, col;

    if ((trow < 0) || (tcol < 0) ||
        (trow > ps_mat->rows) || (tcol > ps_mat->cols))
            _error (pch_name, pch_errorCoord);
    if ((trow + aM_replacement.rows_get () > ps_mat->rows) ||
        (tcol + aM_replacement.cols_get () > ps_mat->cols))
        _error (pch_name, pch_errorSize);
    for (row = 0; row < aM_replacement.rows_get (); row++) {
        for (col = 0; col < aM_replacement.cols_get (); col++) {
            val (row + trow, col + tcol) = aM_replacement.val (row, col);
        }
    }
    return *this;
}

template<typename _CMDATA>
CMatrix<_CMDATA> 
CMatrix<_CMDATA>::row_remove(
        CMatrix<_CMDATA>*&        apV,
        int                     a_rownum) {
    //
    // ARGS
    //  apV             in/out          vector to contain removed
    //                                          row
    //  a_rownum        in              target row to be removed
    //
    // DESC
    //  Removes a specified row.
    //
    // PRECONDITIONS
    //  o 0 < rownum < rows_get()
    //
    // POSTCONDITIONS
    //  o returns row and pointer to row.
    //
    // HISTORY
    //  03 December 2000
    //  o Changed calling parameters to explicitly allow for
    //    holding memory structure.
    //
    // 01 April 2003
    //         o Templatization.
    //

    int                col;
    if (!apV->compatible (1, ps_mat->cols)) {
        delete        apV;
        apV         = new CMatrix<_CMDATA> (1, ps_mat->cols, (_CMDATA) 0.0);
    }

    for (col = 0; col < ps_mat->cols; col++)
        apV->val (0, col) = val (a_rownum, col);
    return *apV;
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::row_insert(
        CMatrix<_CMDATA>&       v,
        int                     trow) {
    //
    // ARGS
    //  v               in              row vector to be inserted
    //  trow            in              row number for insertion
    //
    // DESC
    //  Inserts a row into *this matrix
    //
    // POSTCONDITIONS
    // o The original (this) matrix needs to be destroyed
    //   and enlarged to accommodate the new row
    //
    // HISTORY
    // 01 April 2003
    //         o Templatization.
    //

    char*                pch_name         = "row_insert";
    char*                pch_errorNotRow = "Vector must be a *row* vector.";
    char*                pch_errorSize         = "Vector is incorrect size.";
    char*                pch_errortrow         = "Specified target row is out of range.";
    int                        y, row, col;
    int                        mrows, mcols;
    CMatrix<_CMDATA>*   pM                 = 
            new CMatrix<_CMDATA>(ps_mat->rows, ps_mat->cols, (_CMDATA) 0.0);

    // Make a deepcopy of myself
    pM->copy (*this);

    // Track the current matrix size...
    mrows = ps_mat->rows;
    mcols = ps_mat->cols;

    // Then destroy and create a new matrix in its core structures
    
    coreMatrix_destruct();
    coreMatrix_construct(mrows+1, mcols);
    
    if (v.cols_get () < 2)
        _error (pch_name, pch_errorNotRow);
    if (v.cols_get () != ps_mat->cols)
        _error (pch_name, pch_errorSize);
    if (trow > ps_mat->rows - 1)
        _error (pch_name, pch_errortrow);
    for (row = 0; row < ps_mat->rows; row++)
        for (col = 0; col < ps_mat->cols; col++) {
            if (row < trow)
                val (row, col) = pM->val (row, col);
            if (row == trow)
                val (row, col) = v.val (0, col);
            if (row > trow)
                val (row, col) = pM->val (row - 1, col);
        }
    
    delete pM;
    return *this;
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::row_replace(
        CMatrix<_CMDATA>&       v,
        int                     trow) {
    //
    // DESC
    //  Replaces a row in *this matrix
    //
    // ARGS
    //  v               in      row vector to be inserted
    //  trow            in      row number for insertion
    //
    // POSTCONDITIONS
    //  o *this matrix is returned - internally the target row has
    //    been replaced.
    //
    // 01 April 2003
    //         o Templatization.
    //

    char*        pch_name         = "row_replace";
    char*        pch_errorNotRow = "Vector must be a *row* vector.";
    char*        pch_errorSize         = "Vector is incorrect size.";
    char*        pch_errortrow         = "Specified target row is out of range.";
    int                col;

    if (v.cols_get () < 2)
        _error (pch_name, pch_errorNotRow);
    if (v.cols_get () != ps_mat->cols)
        _error (pch_name, pch_errorSize);
    if (trow > ps_mat->rows - 1)
        _error (pch_name, pch_errortrow);
    for (col = 0; col < ps_mat->cols; col++)
        val (trow, col) = v.val (0, col);
    return *this;
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::col_remove(
        CMatrix<_CMDATA>*&      apV,
        int                     a_colnum) {
    //
    // ARGS
    //  apV             in/out          vector to contain removed
    //                                          column
    //  a_colnum        in              target column to be removed
    //
    // DESC
    //  Removes a specified row.
    //
    // PRECONDITIONS
    //  o 0 < colnum < cols_get()
    //
    // POSTCONDITIONS
    //  o returns column and pointer to column.
    //
    // HISTORY
    // 03 December 2000
    //  o Changed calling parameters to explicitly allow for
    //    holding memory structure.
    //
    // 01 April 2003
    //         o Templatization.
    //

    int                row;
    if (!apV->compatible (ps_mat->rows, 1)) {
        delete        apV;
        apV         = new CMatrix<_CMDATA>(ps_mat->rows, 1, (_CMDATA) 0.0);
    }

    for (row = 0; row < ps_mat->rows; row++)
        apV->val (row, 0) = val (row, a_colnum);
    return *apV;
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::col_insert(
        CMatrix<_CMDATA>&       aV,
        int                     tcol) {
    //
    // ARGS
    //  aV              in              column vector to be inserted
    //  tcol            in              column number for insertion
    //
    // DESC
    //  Inserts a column into a matrix
    //
    // POSTCONDITIONS
    //  o The original (this) matrix needs to be destroyed
    //     and enlarged to accommodate the new column
    // 
    // HISTORY
    // 20 March 2003
    //        o coreMatrix_[construct/destruct] additions
    //
    // 01 April 2003
    //         o Templatization.
    //

    char*                pch_name         = "col_insert";
    char*                pch_errorNotCol = "Vector must be a *column* vector.";
    char*                pch_errorSize         = "Vector is incorrect size.";
    char*                pch_errortcol        = "Specified target column is out of range.";
    int                        y, row, col;
    int                        mrows, mcols;
    CMatrix<_CMDATA>    M (ps_mat->rows, ps_mat->cols, (_CMDATA) 0.0);

    // Make a deepcopy of myself
    M.copy (*this);

    // Track the current matrix size...
    mrows = ps_mat->rows;
    mcols = ps_mat->cols;
    
    // Then destroy and create a new matrix in its core structures
    coreMatrix_destruct();
    coreMatrix_construct(mrows, mcols+1);

    if (aV.rows_get () < 2)
        _error (pch_name, pch_errorNotCol);
    if (aV.rows_get () != ps_mat->rows)
        _error (pch_name, pch_errorSize);
    if (tcol > ps_mat->cols - 1)
        _error (pch_name, pch_errortcol);
    for (row = 0; row < ps_mat->rows; row++)
        for (col = 0; col < ps_mat->cols; col++) {
            if (col < tcol)
                val (row, col) = M.val (row, col);
            if (col == tcol)
                val (row, col) = aV.val (row, 0);
            if (col > tcol)
                val (row, col) = M.val (row, col - 1);
        }
    return *this;
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::col_replace(
        CMatrix<_CMDATA>&       aV,
        int                     tcol) {
    //
    // ARGS
    //  v               in      column vector to be inserted
    //  colnum          in      column number for replacement
    //
    // DESC
    //  Replaces a column in a matrix
    //
    // POSTCONDITIONS
    //  o returns *this matrix containing the "new" column
    //
    // HISTORY
    // 01 April 2003
    //         o Templatization.
    //

    char*        pch_name = "col_replace";
    char*        pch_errorNotCol = "Vector must be a *column* vector.";
    char*        pch_errorSize = "Vector is incorrect size.";
    char*        pch_errortcol = "Specified target column is out of range.";
    int                row;

    if (aV.rows_get () < 2)
        _error (pch_name, pch_errorNotCol);
    if (aV.rows_get () != ps_mat->rows)
        _error (pch_name, pch_errorSize);
    if (tcol > ps_mat->cols - 1)
        _error (pch_name, pch_errortcol);
    for (row = 0; row < ps_mat->rows; row++)
        val (row, tcol) = aV.val (row, 0);
    return *this;
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::diagonal_replace(
        CMatrix<_CMDATA>&       v,
        bool                    ab_dominant        /* = true          */)
{
    //
    // ARGS
    //  v               in      diagonal vector (row or column)
    //  dominant        in      flag describing which diagonal
    //                                  to replace
    //
    // DESC
    //  Replaces a diagonal.
    //
    // Defaults to the dominant diagonal. If `dominant' is false,
    // `anti-diagonal' is replaced.
    //
    // POSTCONDITIONS
    //  o *this matrix has its diagonal replaced and returned.
    //
    // 01 April 2003
    //         o Templatization.
    //

    int                row;
    char*        pch_name         = "diagonal_replace";
    char*        pch_errorRow         = "Vector must be a *column* vector";
    char*        pch_errorVct         = "Diagonal is not a vector";
    char*        pch_errorFit         = "Vector does not fit into matrix";
    char*        pch_errorSqu         = "Matrix must be square to have diagonal replaced";

    if (!v.is_vector ())
        _error (pch_name, pch_errorVct);

    if (v.is_vectorRow ())
        _error (pch_name, pch_errorRow);

    if (v.rows_get () != ps_mat->rows)
        _error (pch_name, pch_errorFit);

    if (!is_square ())
        _error (pch_name, pch_errorSqu);

    if (ab_dominant) {
        for (row = 0; row < rows_get (); row++)
            val (row, row) = v.val (row, 0);
    } else {
        for (row = 0; row < rows_get (); row++)
            val (row, rows_get () - row - 1) = v.val (row, 0);
    }
    return *this;
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::diagonalise(
        _CMDATA         af_offdiag      /*= 0.0         */)
{
    //
    // ARGS
    //  af_offdiag      in              value of off-diagonal elements
    //
    // DESC
    //  Creates a NxN (N=dimension of *this) matrix with *this as the
    //  diagonal (results are stored in pM_target).
    //  Off diagonal elements are set at af_offdiag.
    //
    //  Of course, *this must be a vector!
    //
    // HISTORY
    // 23 September 2000
    //  o Initial design and coding
    //
    // 24 September 2000
    //  o Changed calling interface. Method now creates diagonal
    //    matrix and returns this matrix.
    //
    // 21 October 2000
    //  o Fundamental error in original design. There is no way
    //    to delete the diagonalised matrix that is created by this
    //    method. As as result, method is changed to completely
    //    destroy original matrix, and recreate a new one.
    //
    // 20 March 2003
    //        o coreMatrix_[construct/destruct] integration.
    //
    // 01 April 2003
    //         o Templatization.
    //

    char*        pch_name         = "diagonalise";
    char*        pch_errorNotRow = "passed vector must be a row vector";

    if (!this->is_vectorRow ())
        _error (pch_name, pch_errorNotRow);

    int                        y;
    int                        mrows, mcols;
    CMatrix<_CMDATA>        M (ps_mat->rows, ps_mat->cols, (_CMDATA) 0.0);

    // Make a deepcopy of myself
    M.copy (*this);

    // Track the current matrix size...
    mrows = ps_mat->cols;
    mcols = ps_mat->cols;
    
    // Then destroy and create a new matrix in its core structures
    coreMatrix_destruct();
    coreMatrix_construct(mrows, mcols);
    
    for (int row = 0; row < cols_get (); row++)
        for (int col = 0; col < cols_get (); col++) {
            if (row == col)
                val (row, col) = M.val (0, col);
            else
                val (row, col) = af_offdiag;
        }
    return *this;
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::replicate(
        int             times,
        e_DOMINANCE     dir) {
    //
    // ARGS
    //  times           in      number of times the matrix is replicated
    //  dir             in      direction of replication
    //
    // DESC
    //  Replicates a matrix
    //
    //  This routine replaces the base matrix with a much larger
    //  copy of itself, made up of (times) instances of the original
    //  base. The replications occurs in either a column- or row-dominant
    //  direction.
    //
    //  In other words, if we have an [N x M] matrix that we wish
    //  to replicate p times in a row-dominant fashion, this routine
    //  replaces the base [N x M] matrix with a new matrix:
    //
    //          [ [ N x M] [N x M] ... [N x M] ]
    //
    //  (The column dominant is simply the transpose of the above)
    //
    // PRECONDITIONS
    //        o *this is destroyed and rebuilt!
    //
    // HISTORY
    // 01 April 2003
    //         o Templatization.
    //


    int                y, i;
    int                mrows, mcols;
    CMatrix        M (ps_mat->rows, ps_mat->cols, (_CMDATA) 0.0);

    // Make a deepcopy of myself
    M.copy (*this);

    // Track the current matrix size...
    mrows = ps_mat->rows;
    mcols = ps_mat->cols;

    // Then destroy and create a new matrix in its core structures
    coreMatrix_destruct();
    switch (dir) {
    case e_row:
        coreMatrix_construct(mrows, mcols*times);
        break;
    case e_column:
            coreMatrix_construct(mrows*times, mcols);
        break;
    }

    for (i = 0; i < times; i++) {
        switch (dir) {
        case e_row:
            matrix_replace (0, i * mcols, M);
            break;
        case e_column:
            matrix_replace (i * mrows, 0, M);
            break;
        }
    }
    return *this;
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::sort(
        e_SORTTYPE      e_sorttype      /*= e_ascending         */,
        e_DOMINANCE     e_dominance     /*= e_row               */) {
    //
    // ARGS
    //         e_sorttype        in      direction in which to sort
    //         e_dominance        in      for matrices, row dominant
    //                                  or column dominant
    //                                  sorting
    // DESC
    //         Sorts a matrix (using qsort()).
    //
    //         The e_sorttype should be self-explanatory.
    //         Internally, the matrix is converted to a vector,
    //         and qsort is applied to this vector. The e_dominance
    //         defines how this vector is transformed back into
    //         the original matrix.
    //
    //         This method works most efficiently on row vectors. Ths
    //         bulk of its bloat is to accommodate column vectors and
    //         matrices.
    //
    // HISTORY
    // 14 November 2000
    //  o Initial design and coding.
    //  o testing ok - no leaks over 1e6 iterations
    //
    // 17 February 2002
    //  o Changed behaviour to preserve original matrix contents -
    //          a sorted copy is returned out.
    //
    // 01 April 2003
    //         o Templatization.
    //

    CMatrix<_CMDATA>*        pV_sorted = 0x0;        // Internal row vector that is
                                                //        sorted by quicksort
        
    CMatrix<_CMDATA>                                // The sorted version of *this
            M_sorted(rows_get(), cols_get(), (_CMDATA) 0.0);

    // If *this is a vector, vectorise into a row vector
    if (!is_vector ()) {
        pV_sorted         = new CMatrix<_CMDATA>(*this, e_row, e_rowVector);
    }
    // If *this is a column vector, transpose to a row vector
    if (is_vectorColumn ()) {
        pV_sorted         = new CMatrix<_CMDATA>(1, size1D (), (_CMDATA) 0.0);
        *pV_sorted         = !(*this);                // ! implies transpose
    }

    if (is_vectorRow ()) {
        pV_sorted         = new CMatrix<_CMDATA>(1, size1D (), (_CMDATA) 0.0);
        *pV_sorted         = *this;
    }

    _CMDATA*        pf_data = new _CMDATA[size1D ()];
    for (int i = 0; i < size1D (); i++)
        pf_data[i] = pV_sorted->val (0, i);

    switch (e_sorttype) {
    case e_ascending:
#if _CMDATA == double
        qsort (pf_data, size1D (), sizeof (_CMDATA), doubleCompareAscending);
#elif _CMDATA == float
        qsort (pf_data, size1D (), sizeof (_CMDATA), floatCompareAscending);
#elif _CMDATA == int
        qsort (pf_data, size1D (), sizeof (_CMDATA), intCompareAscending);
#else
        qsort (pf_data, size1D (), sizeof (_CMDATA), doubleCompareAscending);
#endif
        break;
    case e_descending:
#if _CMDATA == double
        qsort (pf_data, size1D (), sizeof (_CMDATA), doubleCompareDescending);
#elif _CMDATA == float
        qsort (pf_data, size1D (), sizeof (_CMDATA), floatCompareDescending);
#elif _CMDATA == int
        qsort (pf_data, size1D (), sizeof (_CMDATA), intCompareDescending);
#else
        qsort (pf_data, size1D (), sizeof (_CMDATA), doubleCompareDescending);
#endif
        break;
    }

    for (int i = 0; i < size1D (); i++)
        pV_sorted->val (0, i) = pf_data[i];

    if (!is_vector ()) {
        int            row, col;
        int            i = 0;
        switch (e_dominance) {
        case e_row:
            for (row = 0; row < rows_get (); row++)
                for (col = 0; col < cols_get (); col++)
                    M_sorted._mval (row, col) = pV_sorted->val (0, i++);
            break;
        case e_column:
            for (col = 0; col < cols_get (); col++)
                for (row = 0; row < rows_get (); row++)
                    M_sorted._mval (row, col) = pV_sorted->val (0, i++);
            break;
        }
    }

    if (is_vectorColumn ()) {
        for (int row = 0; row < size1D (); row++)
            M_sorted._mval (row, 0) = pV_sorted->val (0, row);
    }

    delete[]        pf_data;
    delete          pV_sorted;
    return          M_sorted;
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::fft1D(
    e_FFTDIR        e_direction     /*= e_forward       */,
    bool            b_MKLwrap       /*= false		*/
) {
    
    char* pch_proc = "CMatrix<_CMDATA>::fft1D";
    _error(pch_proc, "1D FFT not yet defined for this data type", 2);
    
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::fft2D(
    e_FFTDIR        e_direction     /*= e_forward       */,
    bool            b_MKLwrap       /*= false		*/
) {
    char* pch_proc = "CMatrix<_CMDATA>::fft2D";
    _error(pch_proc, "2D FFT not yet defined for this data type", 2);
}    

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::nextPowerOf2() {
    //
    // DESC
    //	This method assumes that *this matrix is a row vector describing
    //	dimension sizes. For each value (i.e. dimension), if not a pure
    //	power of two, it determines the next highest power of two.
    //
    // PRECONDITIONS
    //	o Must be a row vector.
    //	o Although templatised, this method really only makes sense for
    //	  <int> matrices.
    //
    // POSTCONDITIONS
    //	o This is a destructive operation!
    //
    // HISTORY
    // 25 March 2004
    //	o Initial design and coding.
    //
    // 18 June 2004
    //	o Painfully embarrassing bug discovered. The dimensions
    //	  should loop over *cols* (not *rows* as the buggy version
    //	  did). This is particularly annoying - more so since this bug was 
    //	  masked by a test vector that had incorrect data ordering -
    //	  the test vector *happened* to work in its initial ordering.
    //
        
    int     dimension           = -1;
    int     ones                = -1;
    int     highestPower        = -1;     // for powerOf2() analysis
    
    if(!is_vectorRow())
    	_error("nextPowerOf2", "*this must be a row vector", 1);
        
    for(dimension=0; dimension<cols_get(); dimension++) {
        powersOf2((int)val(dimension), ones, highestPower);
        if(ones>1) 
	    val(dimension) = 1 << (highestPower+1);
    }
    return *this;
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::zeroPad(
    e_DOMINANCE         dir,
    int                 size) {
    //
    // ARGS
    //  dir             in          the dimension along which to zero pad
    //  size            in          size of padding along a dimension
    //
    // DESC
    //  This method "zeropads" a matrix with 2 submatrices of padding
    //  length 'size' along the direction 'dir'.
    //
    //  Zero padding is necessary before calling fft-type methods on
    //  non-power of two lengthed dimensions.
    //
    // PRECONDITIONS
    //  o None.
    //
    // POSTCONDITIONS
    //  o Non destructive
    //
    // HISTORY
    // 08 September 2003
    //  o Initial design and coding.
    //

    _CMDATA                 z_zero = (_CMDATA) (0);
    int                     rows            = rows_get();
    int                     cols            = cols_get();
    int                     i, j;
    int                     trow, tcol;

    CMatrix<_CMDATA>        M(rows, cols, (_CMDATA) 0.0);
    CMatrix<_CMDATA>*       pM_zeroPadded;

    // Make a deepcopy of myself
    M.copy (*this);

    switch (dir) {
    case e_column:
        pM_zeroPadded   = new CMatrix<_CMDATA>(rows, cols+2*size, z_zero);
        trow    = 0;
        tcol    = size;
        break;
    case e_row:
        pM_zeroPadded   = new CMatrix<_CMDATA>(rows+2*size, cols, z_zero);
        trow    = size;
        tcol    = 0;
        break;
    }

    // and insert the original data in the center
    pM_zeroPadded->matrix_replace(trow, tcol, M);
    CMatrix<_CMDATA>    M_padded(*pM_zeroPadded);
    delete  pM_zeroPadded;
    return M_padded;
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::ifftshiftNewMem() {
    //
    // DESC
    //  This method implements an ifftshift operation, analogous to the MatLAB
    //  function of the same name.
    //
    //  Note that although "fft" appears in the name, no actual fft is performed;
    //  this method merely "reorganises" the internal data contents, shifting the
    //  center of k-space.
    //
    //  In the case of a matrix, the quadrants are labelled as follows:
    //
    //      [B] [C]
    //      [D] [A]
    //
    //  and are "shifted" to
    //
    //      [A] [D]
    //      [C] [B]
    //
    //  where [A] is swapped with [B], and [C] is swapped with [D].
    //  For vectors, this is [B] [A] and [[B] [A]]' (transpose) with
    //  [A] and [B] swapped.
    //
    // POSTCONDITIONS
    //  o This is a non-destructive operation. An internal matrix is
    //    created locally and returned (using the copy constructor).
    //    Make certain to capture the returned matrix!
    //
    // HISTORY
    // 06 September 2003
    //  o Initial design and coding.
    //

    CMatrix<_CMDATA> M_dispatch(*this);

    M_dispatch  = shiftNewMem(-1);

    return(M_dispatch);
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::fftshiftNewMem() {
    //
    // DESC
    //  This method implements an fftshift operation, analogous to the MatLAB
    //  function of the same name.
    //
    //  Note that although "fft" appears in the name, no actual fft is performed;
    //  this method merely "reorganises" the internal data contents, shifting the
    //  center of k-space.
    //
    //  In the case of a matrix, the quadrants are labelled as follows:
    //
    //      [B] [C]
    //      [D] [A]
    //
    //  and are "shifted" to
    //
    //      [A] [D]
    //      [C] [B]
    //
    //  where [A] is swapped with [B], and [C] is swapped with [D].
    //  For vectors, this is [B] [A] and [[B] [A]]' (transpose) with
    //  [A] and [B] swapped.
    //
    // POSTCONDITIONS
    //  o This is a non-destructive operation. An internal matrix is
    //    created locally and returned (using the copy constructor).
    //    Make certain to capture the returned matrix!
    //
    // HISTORY
    // 06 September 2003
    //  o Initial design and coding.
    //

    CMatrix<_CMDATA> M_dispatch(*this);

    M_dispatch  = shiftNewMem(+1);

    return(M_dispatch);
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::shiftNewMem(
    int     a_pivotOffset
) {
    //
    // ARGS
    //  a_pivotOffset       in          offset to pivot point
    //                                      about which matrix
    //                                      elements are swapped.
    //
    // DESC
    //  This method performs the actual innards of both the
    //  fftshift and iffshift operation. Conceptually, both
    //  methods are identical, differing only in the pivot
    //  point offset.
    //
    // PRECONDITIONS
    //  o None.
    //
    //
    // POSTCONDITIONS
    //  o This is a non-destructive operation. An internal matrix is
    //    created locally and returned (using the copy constructor).
    //    Make certain to capture the returned matrix!
    //
    // HISTORY
    // 08 September 2003
    //  o Initial design and coding.
    //

    CMatrix<_CMDATA>*   pM_shifted = new CMatrix<_CMDATA>(1, 1);
    pM_shifted->copy(*this);

    if(is_vectorRow()) {
        // shift for row vectors
        int     Atcol, Acols;
        int     Btcol, Bcols;
        int     Pcol;     // Insertion pivot
        int     cols            = cols_get();

        if(isOdd(cols)) {
            Atcol   = (cols+a_pivotOffset)/2;
            Pcol    = Atcol-a_pivotOffset;
        } else {
            Atcol   = cols/2;
            Pcol    = Atcol;
        }
        Acols   = cols-Atcol;
        Bcols   = Atcol;
        Btcol   = 0;

        CMatrix<_CMDATA>* pM_quadA  = new CMatrix<_CMDATA>(1, Acols);
        CMatrix<_CMDATA>* pM_quadB  = new CMatrix<_CMDATA>(1, Bcols);

        matrix_remove(pM_quadA, 0, Atcol, 1, Acols);
        matrix_remove(pM_quadB, 0, Btcol, 1, Bcols);

        pM_shifted->matrix_replace(0,   0,      *pM_quadA);
        pM_shifted->matrix_replace(0,   Pcol,   *pM_quadB);

        delete pM_quadA;
        delete pM_quadB;

    }

    if(is_vectorColumn()) {
        // shift for col vectors
        int     Atrow, Arows;
        int     Btrow, Brows;
        int     Prow;     // Insertion pivot
        int     rows            = rows_get();

        if(isOdd(rows)) {
            Atrow   = (rows+a_pivotOffset)/2;
            Prow    = Atrow-a_pivotOffset;
        } else {
            Atrow   = rows/2;
            Prow    = Atrow;
        }
        Arows   = rows-Atrow;
        Brows   = Atrow;
        Btrow   = 0;

        CMatrix<_CMDATA>* pM_quadA  = new CMatrix<_CMDATA>(Arows, 1);
        CMatrix<_CMDATA>* pM_quadB  = new CMatrix<_CMDATA>(Brows, 1);

        matrix_remove(pM_quadA, Atrow, 0, Arows, 1);
        matrix_remove(pM_quadB, Btrow, 0, Brows, 1);

        pM_shifted->matrix_replace(0,       0,      *pM_quadA);
        pM_shifted->matrix_replace(Prow,    0,      *pM_quadB);

        delete pM_quadA;
        delete pM_quadB;

    }

    if
    (!is_vector()) {
        // shift for matrices
        //  We need to determine the top left coords and lengths of
        //  each quadrant. This is also dependant on whether or not
        //  the dimension length is even or odd.
        int     Atrow, Atcol, Arows, Acols;
        int     Btrow, Btcol, Brows, Bcols;
        int     Ctrow, Ctcol, Crows, Ccols;
        int     Dtrow, Dtcol, Drows, Dcols;

        int     Prow, Pcol;     // Insertion pivot

        int     rows    = rows_get();
        int     cols    = cols_get();

        Btrow   = 0;
        Btcol   = 0;
        if(isOdd(cols_get())) {
            Atcol   = (cols+a_pivotOffset)/2;
            Pcol    = Atcol-a_pivotOffset;
        } else {
            Atcol   = cols/2;
            Pcol    = Atcol;
        }
        Acols   =  cols-Atcol;
        Bcols   = Atcol;
        Ctcol   = Atcol;
        Ccols   = Acols;
        Dtcol   = 0;
        Dcols   = Bcols;

        if(isOdd(rows_get())) {
            Atrow   = (rows+a_pivotOffset)/2;
            Prow    = Atrow -a_pivotOffset;
        } else {
            Atrow   = rows/2;
            Prow    = Atrow;
        }
        Arows   = rows-Atrow;
        Brows   = Atrow;
        Ctrow   = 0;
        Crows   = Brows;
        Dtrow   = Atrow;
        Drows   = Arows;

        CMatrix<_CMDATA>* pM_quadA  = new CMatrix<_CMDATA>(Arows, Acols);
        CMatrix<_CMDATA>* pM_quadB  = new CMatrix<_CMDATA>(Brows, Bcols);
        CMatrix<_CMDATA>* pM_quadC  = new CMatrix<_CMDATA>(Crows, Ccols);
        CMatrix<_CMDATA>* pM_quadD  = new CMatrix<_CMDATA>(Drows, Dcols);

        matrix_remove(pM_quadA, Atrow, Atcol, Arows, Acols);
        matrix_remove(pM_quadB, Btrow, Btcol, Brows, Bcols);
        matrix_remove(pM_quadC, Ctrow, Ctcol, Crows, Ccols);
        matrix_remove(pM_quadD, Dtrow, Dtcol, Drows, Dcols);

        pM_shifted->matrix_replace(0,       0,      *pM_quadA);
        pM_shifted->matrix_replace(Prow,    Pcol,   *pM_quadB);
        pM_shifted->matrix_replace(0,       Pcol,   *pM_quadD);
        pM_shifted->matrix_replace(Prow,    0,      *pM_quadC);

        delete pM_quadA;
        delete pM_quadB;
        delete pM_quadC;
        delete pM_quadD;

    }
    CMatrix M_shifted(*pM_shifted);
    delete pM_shifted;
    return(M_shifted);
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::fftshift(
        bool    ab_inPlace
) {
    //
    // ARGS
    //  ab_inPlace              in              if true, branch to the
    //                                                  inPlace() method
    //                                          else, branch to the
    //                                                  newMem() method
    //
    // DESCRIPTION
    //  shift*() dispatching layer. Depending on the value of ab_inPlace,
    //  either an inPlace or new memory operation is performed.
    //
    // HISTORY
    //  22 January 2003
    //  o Initial design and coding.
    //

    if(ab_inPlace) {
        shiftInPlace(+1);
        return *this;
    } else {
        CMatrix<_CMDATA> M_dispatch(*this);
        M_dispatch  = shiftNewMem(+1);
        return(M_dispatch);
    }
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::ifftshift(
        bool    ab_inPlace
) {
    //
    // ARGS
    //  ab_inPlace              in              if true, branch to the
    //                                                  inPlace() method
    //                                          else, branch to the
    //                                                  newMem() method
    //
    // DESCRIPTION
    //  shift*() dispatching layer. Depending on the value of ab_inPlace,
    //  either an inPlace or new memory operation is performed.
    //
    // HISTORY
    //  22 January 2003
    //  o Initial design and coding.
    //

    if(ab_inPlace) {
        shiftInPlace(-1);
        return *this;
    } else {
        CMatrix<_CMDATA> M_dispatch(*this);
        M_dispatch  = shiftNewMem(-1);
        return(M_dispatch);
    }
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::fftshiftInPlace()
{
    shiftInPlace(+1);
    return *this;
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::ifftshiftInPlace()
{
    shiftInPlace(-1);
    return *this;
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::ifftshiftIndex( CMatrix<_CMDATA>&   aM_index)
{ 
    CMatrix<_CMDATA>	M_copy(*this);
    M_copy = shiftIndex(aM_index, -1);
    return(M_copy);
}

template<typename _CMDATA>
CMatrix<_CMDATA>	
CMatrix<_CMDATA>::fftshiftIndex(  CMatrix<_CMDATA>&   aM_index)
{ 
    CMatrix<_CMDATA>	M_copy(*this);
    M_copy = shiftIndex(aM_index, +1); 
    return(M_copy);
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::shiftIndex(
    CMatrix<_CMDATA>&	aM_index,
    int			a_pivotOffset
) {
    //
    // ARGS
    //	aM_index	in		index of position to shift
    //  a_pivotOffset   in              offset to pivot point
    //                                      about which matrix
    //                                      elements are swapped.
    //
    // DESC
    //	This method calculates the shift in an index position,
    //	aM_index. The *this pointer denotes a matrix containing the
    //	dimensions of the space (either 2D or 3D).
    //
    // PRECONDITIONS
    //  o *this is a row vector (1x2 or 1x3) denoting the space inwhich 
    //	  aM_index is to be shifted.
    //
    //
    // POSTCONDITIONS
    //  o A vector is returned denoting the shift of aM_index in 
    //	  *this space.
    //
    // HISTORY
    // 23 February 2004
    //  o Initial design and coding, based around the inPlace space
    //	  algorithm.
    //
    
    char*	pch_name = "shiftIndex(...)";
    if(!compatible(aM_index))
	_error(	pch_name,
		"aM_index must have same dimensions as *this",
		-1);
    
    CMatrix<_CMDATA>    M_indexDelta(1,2);      // the "jump" shift
                                                //	note that is 
                                                // 	2D! A biphase
                                                //	calculation is
                                                //	needed for 3D:
                                                //	first *within* a slice
                                                //      and then *along* the
                                                //	slices
    CMatrix<_CMDATA>	M_indexCurrent(1, 2);
    CMatrix<_CMDATA>	M_indexNext(1,2);
    CMatrix<_CMDATA>	M_shifted(1, 1);
    
    M_shifted.copy(*this);

    long int            rows;
    long int            cols;
    long int		slices		= cols_get()==3 ? (long int)val(2) : 0;
    int			loop		= 1;
    
    if(slices) 		loop            = 2;	// if 3D, loop once for slice
                                                //	plane, and once for
                                                //	slices
    
    for(int i=0; i<loop; i++) {
	if(!i) {
	    rows			= (long int) val(0);
	    cols			= (long int) val(1);
	    M_indexCurrent(0)		= aM_index(0);
	    M_indexCurrent(1)		= aM_index(1);
	} else {
	    rows			= 1;
	    cols			= slices;
	    M_indexCurrent(0)		= 0;
	    M_indexCurrent(1)		= aM_index(2);
	}
	
	if(isOdd(rows))     M_indexDelta(0)   =  (rows-a_pivotOffset)/2;
	    else            M_indexDelta(0)   =  rows/2;
	if(isOdd(cols))     M_indexDelta(1)   =  (cols-a_pivotOffset)/2;
            else            M_indexDelta(1)   =  cols/2;
	    
	M_indexNext = M_indexCurrent + M_indexDelta;
    
	// Check for boundary violation in indices
	if(M_indexNext(0)>=rows)        M_indexNext(0)-=rows;   // rows
	if(M_indexNext(0)<0)            M_indexNext(0)+=rows;
	if(M_indexNext(1)>=cols)        M_indexNext(1)-=cols;   // cols
	if(M_indexNext(1)<0)            M_indexNext(1)+=cols;
	
	if(!i) {
	    M_shifted.val(0)	= M_indexNext.val(0);
	    M_shifted.val(1)	= M_indexNext.val(1);
	} else
	    M_shifted.val(2)	= M_indexNext.val(1);
 
    }

    return M_shifted;
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::shiftInPlace(
    int                 a_pivotOffset
) {
    //
    // ARGS
    //  a_pivotOffset       in          offset to pivot point
    //                                      about which matrix
    //                                      elements are swapped.
    //
    // DESC
    //  This method performs the actual innards of both the
    //  fftshift and iffshift operation. Conceptually, both
    //  methods are identical, differing only in the pivot
    //  point offset.
    //
    // PRECONDITIONS
    //  o None.
    //
    //
    // POSTCONDITIONS
    //  o This is a destructive (in place) operation. The original matrix
    //    has its internal contents shifted, and is not preserved.
    //
    // HISTORY
    // 20 January 2004
    //  o Initial design and coding.
    //    Usage of the non-destructive shift quickly becomes prohibitive with
    //    when using large volumes. As a matter of fact, the apparent simplicity
    //    of a non-destructive operation has not only a memory penalty, but a
    //    time penalty as well: the time required to create a new block of
    //    memory can be quite large (especially if the operating system starts
    //    to swap).
    //

    //  For odd lengthed dimensions, the "swap" is asymmetrical. An element is
    //  not simply swapped with its "shifted" target, but this target
    //  displaces a different element.
    //
    //  This displaced element in turn "jumps" forward, displacing
    //  another (different) element, which displaces another element
    //  etc. This chain continues until the jumps
    //  terminate again at the original start point.
    //
    CMatrix<int>                M_indexDelta(1,2);      // the "jump" shift

    long int            rows            = rows_get();
    long int            cols            = cols_get();

    // The following keep track of the indices in the matrix index space that
    //  are being shifted. M_startPos begins at (0, 0) and moves through the
    //  indices as governed by M_currentDelta. The shiftCounter counts how
    //  many cells have been shifted. Once this shiftCounter equals the number
    //  of matrix elements, the shift operation terminates.
    CMatrix<int>        M_indexStart(1, 2);             // start position of
                                                        //      current shift;
    CMatrix<int>        M_indexCurrent(1, 2);           // current shift index
    CMatrix<int>        M_indexNext(1, 2);              // next shift index
    CMatrix<int>        M_indexStartDelta(1, 2);        // delta for start index
    long int            shiftCounter    = 0;            // running counter of
                                                        //      shifted cells
    long int            totalElements   = 0;
    long int            juggleSeqLength = 0;
    _CMDATA             cellCurrentValue        = (_CMDATA) 0;
    _CMDATA             cellNextValue           = (_CMDATA) 0;
    _CMDATA             cellJuggleValue         = (_CMDATA) 0;

    // The indexStartDelta is used to update the "head" of a new shift vector in
    //  the matrix index space. The matrix elements are updated *orthogonally*
    //  to the indexDelta direction
    M_indexStartDelta(0)        = isEven(rows);
    M_indexStartDelta(1)        = isEven(cols);
    // If matrix dimensions are odd/odd,
    if(!M_indexStartDelta(0) && !M_indexStartDelta(1)) {
        M_indexStartDelta(0) = -1;      // "up" one row
        M_indexStartDelta(1) =  1;      // "over" one column
    }

    if(isOdd(rows))     M_indexDelta(0)   =  (rows-a_pivotOffset)/2;
        else            M_indexDelta(0)   =  rows/2;
    if(isOdd(cols))     M_indexDelta(1)   =  (cols-a_pivotOffset)/2;
        else            M_indexDelta(1)   =  cols/2;

    // start at (0, 0) in the index space
    M_indexStart(0)         = 0;            M_indexStart(1)   = 0;
    M_indexCurrent(0)       = 0;            M_indexCurrent(1) = 0;

    totalElements   = rows*cols;
    while(++shiftCounter <= totalElements) {
        M_indexNext = M_indexCurrent + M_indexDelta;

        // Check for boundary violation in indices
        if(M_indexNext(0)>=rows)        M_indexNext(0)-=rows;   // rows
        if(M_indexNext(0)<0)            M_indexNext(0)+=rows;
        if(M_indexNext(1)>=cols)        M_indexNext(1)-=cols;   // cols
        if(M_indexNext(1)<0)            M_indexNext(1)+=cols;

        if(!juggleSeqLength++)
            cellJuggleValue    = val(M_indexCurrent(0), M_indexCurrent(1));

        cellNextValue       = val(M_indexNext(0), M_indexNext(1));
        val(M_indexNext(0), M_indexNext(1))     = cellJuggleValue;

        M_indexCurrent      = M_indexNext;
        cellJuggleValue     = cellNextValue;

        if(M_indexNext.equal(M_indexStart)) {
            M_indexStart += M_indexStartDelta;
            // Check for boundary violation in indices
            if(M_indexStart(0)>=rows)   M_indexStart(0)-=rows;  // rows
            if(M_indexStart(0)<0)       M_indexStart(0)+=rows;
            if(M_indexStart(1)>=cols)   M_indexStart(1)-=cols;  // cols
            if(M_indexStart(1)<0)       M_indexStart(1)+=cols;

            if(M_indexStartDelta(0)==1 && M_indexStartDelta(1)==1) {
                // In this case the matrix dimensions are even/even
                //  and the shift is completely symmetrical: i.e.
                //  a straight swap of elements. In this case, we simply
                //  advance the M_indexStart linearly across the matrix
                int row     = 0;
                int col     = 0;
                                                // Since there at least two
                row         = (shiftCounter/2)/cols;    // jumps required to
                col         = (shiftCounter/2)%cols;    // return to start point
                                                // we divide shiftCounter by 2
                M_indexStart(0)     = row;
                M_indexStart(1)     = col;
            }
            M_indexCurrent = M_indexStart;
            juggleSeqLength = 0;
        }

    }
    return *this;
}

//////---------->
////// Search/replace/fill routines
//////---------->

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::findAll(
        _CMDATA         what,
        int&            occurences) {
    //
    // ARGS
    //  what            in              target number to search for
    //  occurences      in              number of occurences
    //
    // DESC
    //  searches for instances of `what'
    //
    // POSTCONDITIONS
    //  o returns a matrix of (row, col) of instances in a ->col matrix <-
    //
    // 02 April 2003
    //         o Templatization.
    //

    int                row, col, count;

    occurences = 0;
    for (row = 0; row < ps_mat->rows; row++)
        for (col = 0; col < ps_mat->cols; col++)
            if (val (row, col) == what)
                occurences++;
    if (!occurences) {
        CMatrix<_CMDATA>        occur(1, 1, (_CMDATA) 0.0);
        return (occur);
    }
    if (occurences) {
        count = 0;
        CMatrix        occur (occurences, 2, (_CMDATA) 0.0);
        for (row = 0; row < ps_mat->rows; row++)
            for (col = 0; col < ps_mat->cols; col++)
                if (val (row, col) == what) {
                    occur.val (count, 0) = row;
                    occur.val (count++, 1) = col;
                }
        return (occur);
    }
    // The following code is never executed - merely added so that
    // compiling with -Wall doesn't produce warnings!!
    CMatrix<_CMDATA>        occur(1, 1, (_CMDATA) 0.0);
    return (occur);
}

template<typename _CMDATA>
int
CMatrix<_CMDATA>::findAndReplace (
        _CMDATA         target,
        _CMDATA         source)
{
    //
    // ARGS
    //  target          in              _CMDATA value that will be replaced
    //  source          in              _CMDATA value that will be entered into replace cells
    //
    // DESC
    //  Find and replace all values that are equal to target with source.
    //
    // POSTCONDITIONS
    //  o return the number of instances replaced
    //
    // 01 April 2003
    //         o Templatization.
    //

    int                row, col;
    int                num                 = -1;

    for (row = 0; row < ps_mat->rows; row++)
        for (col = 0; col < ps_mat->cols; col++) {
            if (val (row, col) == target) {
                num++;
                val (row, col) = source;
            }
        }
    return num;
}

template<typename _CMDATA>
bool
CMatrix<_CMDATA>::find_quantized(
        _CMDATA                 af_what,
        CMatrix<_CMDATA>*       pM_index) {
    //
    // ARGS
    //  af_what         in              search value
    //  pM_index        out             index in matrix
    //
    // DESC
    //  Searches for `af_what' and returns its position in the
    //  matrix. Note that a *quantized* search is performed,
    //  which means that the matrix position "closest" to the search term
    //  is returned. The concept of "closest" is defined within
    //  the context of the matrix itself, i.e. the span between two
    //  successive elements of the matrix.
    //
    // PRECONDITIONS
    // o Currently the input matrix, *this, must be a column vector
    //
    // POSTCONDITIONS
    // o If the value, `af_what', falls completely out of range, it is
    //   "hardlimited" as it were to one of the end indices.
    // o Search is over the range [af_start, af_stop] - not the square
    //   brackets, not parenthesis. End points are included. This may
    //   have an impact on deciding if af_what falls in the right-most or
    //   next-to-rightmost position.
    // o Intervals over indices are as follows:
    //          [..] (..] (..] ... (..]
    //            0    1    2        n-1
    //
    // RETURN
    // o b_outOfBounds    return          0 - value found inside range
    //                                    1 - value found outside range
    //
    // HISTORY
    // 06 Aug 2000
    //  o Initial design and coding
    //  o NOTE: column vector is assumed in these initial versions
    //
    // 28 Nov 2000
    //  o Added <= condition on af_what <= val(row...)
    //    if equality is not tested for, error will be returned even
    //    when af_what is in valid range.
    //  o Expanded logic such that quantised lookup includes end points
    //    Intervals over indices are as follows:
    //          [..] (..] (..] ... (..]
    //            0    1    2        n-1
    //
    // 04 April 2001
    //  o Added "hardlimiting" of out-of-range values to endpoints of
    //    spectrum
    //  o Return value changed to b_outOfBounds.
    //
    // 001 April 2003
    //        o Templatization.
    //

    if (!is_vectorColumn())
        _error ("find_quantized", "Input vector must be a column vector", 1);

    if (!pM_index->compatible (1, 2)) {
        delete                    pM_index;
        pM_index         = new CMatrix<_CMDATA>(1, 2, (_CMDATA) 0.0);
    }
    // First, check for extrema
    // These are values that fall outside the defined range, either to
    // the left or the right. Such values are "hardlimited" to the ends
    // of the quantized range
    if (af_what < val (0, 0)) {
        pM_index->val (0, 0) = 0;
        pM_index->val (0, 1) = 0;
        return true;
    }
    if (af_what > val (rows_get () - 1, 0)) {
        pM_index->val (0, 0) = rows_get () - 1;
        pM_index->val (0, 1) = 0;
        return true;
    }
    // At this point, the af_what value must lie somewhere within
    // the quantised range. Search through this range...
    for (int row = 0; row < rows_get () - 1; row++) {
        bool            b_found = false;
        if (row == 0) {
            // left most quantised index
            if ((af_what >= val (row, 0)) && (af_what <= val (row + 1, 0)))
                b_found = true;
        } else {
            // remaining indices
            if ((af_what > val (row, 0)) && (af_what <= val (row + 1, 0)))
                b_found = true;
        }
        if (b_found) {
            pM_index->val (0, 0) = row;
            pM_index->val (0, 1) = 0;
            return false;
        }
    }
    // If processing ever gets to this point, something bizarre has
    // gone wrong
    _error ("find_quantized", "Completely out-of-state condition reached",
             1);
    // This line is never reached - only added to force clean compile...
    return false;
}

template<typename _CMDATA>
void
CMatrix<_CMDATA>::perimeter_set(
        _CMDATA         af_value,
        int             offset          /*= 0           */) {
    //
    // ARGS
    //  af_value        in              the value all "perimeter" cells
    //                                          are set to
    //  offset          in              "offset" from edge of matrix
    //
    // DESC
    //  Fills a "square" perimeter with af_value.
    //
    // NOTE
    //  offset unused
    //
    // PRECONDITIONS
    //  o Assumes that perimeter can indeed be set
    //
    // POSTCONDITIONS
    //  o Perimeter is set to af_value
    //
    // HISTORY
    // 02 April 2001
    //  o Initial design and coding.
    //
    // 01 April 2003
    //         o Templatization.
    //

    if (rows_get () - 2 * offset < 1)
        _error ("perimeter_set", "row offset invalid", -1);
    if (cols_get () - 2 * offset < 1)
        _error ("perimeter_set", "col offset invalid", -1);
    int                rows = rows_get ();
    int                cols = cols_get ();

    for (int row = 0; row < rows_get () - 2 * offset; row++) {
        val (row + offset, offset) = af_value;
        val (row + offset, cols - offset - 1) = af_value;
    }
    for (int col = 0; col < cols_get () - 2 * offset; col++) {
        val (offset, col + offset) = af_value;
        val (rows - offset - 1, col + offset) = af_value;
    }
}

//////---------->
////// Mathematical considerations
//////---------->

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::functionApply (
        double (*func)(double)
) {
    //
    // ARGS
    //  (*func)(double)         in      pointer to a function that
    //                                          has a single double
    //                                          argument and returns
    //                                          a double
    //
    // DESC
    //  This method applies the function (*func) to each element
    //  of the matrix.
    //
    // PRECONDITIONS
    //  o Make sure that (*func) is a valid function
    //
    // POSTCONDITIONS
    //  o returns *this
    //
    // HISTORY
    //  09 july 2001
    //  o Initial design and coding.
    //
    // 26 March 2003
    //         o Templatization.
    //

    for (int row = 0; row < rows_get (); row++)
        for (int col = 0; col < cols_get (); col++)
            _mval (row, col) = (_CMDATA) (*func) (_mval (row, col));
    return (*this);
}

//////---------->
////// Simple statistics        
//////---------->

template<typename _CMDATA>
_CMDATA
CMatrix<_CMDATA>::mean ()
{
    //
    // DESC
    //  Determines the mean of an entire matrix.
    //
    // PRECONDITIONS
    //  o None
    //
    // POSTCONDITIONS
    //  o A single _CMDATA describing the mean is returned
    //
    // HISTORY
    //  06 July 2001
    //  o Expanded from old definition
    //
    // 26 March 2003
    //         o Templatization.
    //

    _CMDATA                sum = 0;
    int                        r = ps_mat->rows;
    int                        c = ps_mat->cols;
    int                        col = 0;
    int                        row = 0;

    for (row = 0; row < r; row++)
        for (col = 0; col < c; col++)
            sum += _mval (row, col);
    return sum / (row * col);
}

template<typename _CMDATA>
_CMDATA
CMatrix<_CMDATA>::mean (
        e_DOMINANCE     ae_type,
        int             a_index)
{
    //
    // ARGS
    //  ae_type         in              type of mean to determine:
    //                                          row wise or
    //                                          col wise
    //  a_index         in              row of column for which the mean is
    //                                  calculated
    //
    // DESC
    //  Determines the mean of a particular column.row.
    //
    // PRECONDITIONS
    //  o *this MUST be of type e_matrix.
    //
    // POSTCONDITIONS
    //  o the mean of the index'd row/column is returned.
    //
    // HISTORY
    // 06 July 2001
    //         o Initial design and coding.
    //
    // 02 April 2003
    //        o Templatization.
    //
    
    char*        pch_proc         = "mean(e_DOMINANCE ae_type, int a_index)";
    _CMDATA        f_sum                 =  (_CMDATA) 0.0;
    _CMDATA        f_mean                 =  (_CMDATA) 0.0;
    int                rows = rows_get ();
    int                cols = cols_get ();

    // some initial error checking
    if (matrix_type () != e_matrix && matrix_type() != e_square)
        _error (pch_proc, "Matrix must not be a vector!");

    switch (ae_type) {
    case e_row:
        for (int col = 0; col < cols; col++)
            f_sum += _mval (a_index, col);
        f_mean = f_sum / cols;
        break;
    case e_column:
        for (int row = 0; row < rows; row++)
            f_sum += _mval (row, a_index);
        f_mean = f_sum / rows;
        break;
    }
    return f_mean;
}

template<typename _CMDATA>
CMatrix<_CMDATA> 
CMatrix<_CMDATA>::mean (
        e_DOMINANCE     ae_type)
{
    //
    // ARGS
    //  ae_type         in              governs how the mean is calculated.
    //
    // DESC
    //  Determines the mean of matrix elements. The ae_type parameter
    //  governs whether the row- or column-wise elements
    //  are used in determining the mean.
    //
    //  If ae_type == e_row,    then the row-wise mean is determined;
    //  if ae_type == e_column, then the col-wise mean is determined.
    //
    //  In both cases a vector is returned whose row/col element is the mean
    //  of the corresponding row/col vector.
    //
    // PRECONDITIONS
    //  o None
    //
    // POSTCONDITIONS
    //  o a vector is returned.
    //  o the return vector is created locally and passed out.
    //
    // HISTORY
    //  06 July 2001
    //  o Expanded from old definition
    //
    // 12 February 2002
    //  o Fixed minor memory leak
    //
    // 02 April 2003
    //        o Templatization.
    //

    char*        pch_proc         = "mean(e_DOMINANCE ae_type)";
    _CMDATA        f_sum                 = (_CMDATA)0.0;
    int                rows                 = rows_get ();
    int                cols                 = cols_get ();
    CMatrix        V_mean (1, 1, (_CMDATA) 0.0);
    CMatrix*        pV = 0x0;

    // some initial error checking
    if (matrix_type () != e_matrix && matrix_type() != e_square)
        _error (pch_proc, "Matrix must not be a vector!");

    switch (ae_type) {
    case e_row:
        pV = new CMatrix<_CMDATA>(rows_get (), 1, (_CMDATA) 0.0);
        for (int row = 0; row < rows; row++) {
            f_sum = (_CMDATA) 0.0;
            for (int col = 0; col < cols; col++)
                f_sum += _mval (row, col);
            pV->_mval (row, 0) = f_sum / cols;
        }
        break;
    case e_column:
        pV = new CMatrix<_CMDATA>(1, cols_get (), (_CMDATA) 0.0);
        for (int col = 0; col < cols; col++) {
            f_sum = (_CMDATA) 0.0;
            for (int row = 0; row < rows; row++)
                f_sum += _mval (row, col);
            pV->_mval (0, col) = f_sum / rows;
        }
        break;
    }
    V_mean.copy (*pV);
    delete pV;
    return V_mean;
}

template<typename _CMDATA>
_CMDATA
CMatrix<_CMDATA>::std (
        bool                    ab_sample       /*= false               */)
{
    //
    // ARGS
    //  ab_sample               in              if true, use the "sample"
    //                                                  standard deviation
    //
    // DESC
    //  Determines the standard deviation of an entire matrix.
    //
    // PRECONDITIONS
    //  o None
    //
    // POSTCONDITIONS
    //  o A single _CMDATA describing the standard deviation is returned
    //
    // HISTORY
    //  10 July 2001
    //  o Expanded from mean definition
    //
    // 26 March 2003
    //         o Templatization.
    //

    _CMDATA                f_std         = (_CMDATA) 0.0;
    _CMDATA                f_mean         = mean();
    int                        N = ps_mat->rows * ps_mat->cols;
    CMatrix<_CMDATA>        M_std(ps_mat->rows, ps_mat->cols, (_CMDATA) 0.0);
    M_std                 = (*this - f_mean);
    M_std.functionApply(sqr);
    _CMDATA                f_sum = M_std.sum ();

    if (ab_sample)
        N--;
    f_std = (_CMDATA) sqrt((double) f_sum / (_CMDATA) N);

    return f_std;
}

template<typename _CMDATA>
_CMDATA
CMatrix<_CMDATA>::std (
        e_DOMINANCE     ae_type,
        int             a_index,
        bool            ab_sample       /*= false               */)
{
    //
    // ARGS
    //  ae_type         in              type of std to determine:
    //                                          row wise or
    //                                          col wise
    //  a_index         in              row of column for which the std is
    //                                          calculated
    //  ab_sample       in              if true, determine the "sample"
    //                                          standard deviation
    //
    // DESC
    //  Determines the standard deviation of a particular column or row.
    //
    // PRECONDITIONS
    //  o *this MUST be of type e_matrix.
    //
    // POSTCONDITIONS
    //  o the std of the index'd row/column is returned.
    //
    // HISTORY
    // 09 July 2001
    //  o Initial design and coding.
    //
    // 02 April 2003
    //         o Templatization.
    //
    
    char*        pch_proc = "std(e_DOMINANCE ae_type, int a_index, bool ab_sample)";
    int                N;
    N = (ae_type == e_column) ? rows_get () : cols_get ();
    if (ab_sample)
        N--;

    // some initial error checking
    if (matrix_type () != e_matrix && matrix_type() != e_square)
        _error (pch_proc, "Matrix must not be a vector!");

    CMatrix*        pV = new CMatrix<_CMDATA>(1, 1, (_CMDATA) 0.0);
    if (ae_type == e_column)
        col_remove (pV, a_index);
    else
        row_remove (pV, a_index);

    _CMDATA        f_std =
        (_CMDATA)sqrt ((double)(*pV - pV->mean ()).functionApply (sqr).sum () / (_CMDATA) N);

    delete        pV;
    return f_std;
}

template<typename _CMDATA>
CMatrix<_CMDATA> 
CMatrix<_CMDATA>::std (
        e_DOMINANCE     ae_type,
        bool            ab_sample       /*= false               */)
{
    //
    // ARGS
    //  ae_type         in              governs how the std dev is calculated.
    //  ab_sample       in              if true, use the "sample" standard
    //                                          deviation
    //
    // DESC
    //  Determines the std dev of matrix elements. The optional ae_type
    //  defaults to e_matrix, implying that all elements in the matrix
    //  are used in determining the mean.
    //
    //  If ae_type == e_row,    then the row-wise std  is determined;
    //  if ae_type == e_column, then the col-wise std  is determined.
    //
    //  In both cases a vector is returned whose row/col element is the std
    //  of the corresponding row/col vector.
    //
    // PRECONDITIONS
    //  o None
    //
    // POSTCONDITIONS
    //  o a vector is returned.
    //  o the return vector is created locally and passed out.
    //
    // HISTORY
    // 09 July 2001
    //  o Expanded from mean definition
    //
    // 13 February 2002
    //  o Noticed that original design really only assumed "column"
    //    dominance. Finished coding for row dominance as well.
    //
    // 02 April 2003
    //        o Templatization.
    //

    char*        pch_proc = "std(e_DOMINANCE ae_type, bool ab_sample)";
    int                N;

    // some initial error checking
    if (matrix_type () != e_matrix && matrix_type() != e_square)
        _error (pch_proc, "Matrix must not be a vector!");

    N = (ae_type == e_column) ? rows_get () : cols_get ();

    CMatrix<_CMDATA>        M_mean(1, 1, (_CMDATA) 0.0);
    M_mean.copy (mean (ae_type));
    M_mean.replicate (N, ae_type);
    CMatrix<_CMDATA>        M_sumSqrd(M_mean.rows_get (), M_mean.cols_get (), (_CMDATA) 0.0);
    M_sumSqrd = *this - M_mean;
    M_sumSqrd.functionApply (sqr);
    CMatrix<_CMDATA>        V_std(1, 1, (_CMDATA) 0.0);
    CMatrix<_CMDATA>        *pV;

    if(ae_type == e_row)
        pV = new CMatrix(rows_get(), 1, (_CMDATA) 0.0);
    else
        pV = new CMatrix(1, cols_get(), (_CMDATA) 0.0);

    V_std.copy(*pV);
    V_std = M_sumSqrd.sum (ae_type);

    if (ab_sample)
        N--;
    V_std.scale (1 / (_CMDATA) N);
    V_std.functionApply (sqrt);
    delete pV;
    return (V_std);
}

template<typename _CMDATA>
_CMDATA
CMatrix<_CMDATA>::sum (
        bool            ab_abs  /*= false       */)
{
    //
    // ARGS
    //  ab_abs          in              if true, uses the absolute value of
    //                                          each matrix element
    //
    // DESC
    //  Determines the sum of an entire matrix.
    //
    // PRECONDITIONS
    //  o None
    //
    // POSTCONDITIONS
    //  o A single _CMDATA describing the sum is returned
    //
    // HISTORY
    //  10 July 2001
    //  o Expanded from old definition
    //
    // 26 March 2003
    //         o Templatization.
    //

    _CMDATA                sum         = 0;
    int                        r         = ps_mat->rows;
    int                        c         = ps_mat->cols;
    int                        col;
    int                        row;
    for (row = 0; row < r; row++)
        for (col = 0; col < c; col++) {
            if (ab_abs)
                sum += (_CMDATA) fabs ((double)_mval (row, col));
            else
                sum += _mval (row, col);
        }
    return sum;
}

template<typename _CMDATA>
_CMDATA
CMatrix<_CMDATA>::sum (
        e_DOMINANCE     ae_type,
        int             a_index,
        bool            ab_abs          /*= false       */)
{
    //
    // ARGS
    //  ae_type         in              type of sum to determine:
    //                                          row wise or
    //                                          col wise
    //  a_index         in              row of column for which the sum is
    //                                          calculated
    //  ab_abs          in              if true, uses the absolute value of
    //                                          each matrix element
    //
    // DESC
    //  Determines the sum of a particular column/row.
    //
    // PRECONDITIONS
    //  o *this MUST be of type e_matrix.
    //
    // POSTCONDITIONS
    //  o the mean of the index'd row/column is returned.
    //
    // HISTORY
    // 10 July 2001
    //         o Initial design and coding.
    //
    // 02 April 2003
    //        o Templatization.
    //

    char*        pch_proc         = "sum(e_DOMINANCE ae_type, int a_index, bool ab_abs)";
    _CMDATA        f_sum                 = (_CMDATA)0.0;
    int                rows = rows_get ();
    int                cols = cols_get ();

    // some initial error checking
    if (matrix_type () != e_matrix && matrix_type() != e_square)
        _error (pch_proc, "Matrix must not be a vector!");

    switch (ae_type) {
    case e_row:
        for (int col = 0; col < cols; col++) {
            if (ab_abs)
                f_sum += (_CMDATA)fabs ((double)_mval (a_index, col));
            else
                f_sum += _mval (a_index, col);
        }
        break;
    case e_column:
        for (int row = 0; row < rows; row++) {
            if (ab_abs)
                f_sum += (_CMDATA)fabs ((double)_mval (row, a_index));
            else
                f_sum += _mval (row, a_index);
        }
        break;
    }
    return f_sum;
}

template<typename _CMDATA>
CMatrix<_CMDATA> 
CMatrix<_CMDATA>::sum (
        e_DOMINANCE     ae_type,
        bool            ab_abs          /*= false               */)
{
    //
    // ARGS
    //  ae_type         in              governs how the sum is calculated.
    //  ab_abs          in              if true, uses the absolute value of
    //                                          each matrix element
    //
    // DESC
    //  Determines the sum of matrix elements. The ae_type argument specifies
    //  whether to determine the row- or column-wise sum.
    //
    //  If ae_type == e_row,    then the row-wise sum is determined;
    //  if ae_type == e_column, then the col-wise sum is determined.
    //
    //  In both cases a vector is returned whose row/col element is the sum
    //  of the corresponding row/col vector.
    //
    // PRECONDITIONS
    //  o None
    //
    // POSTCONDITIONS
    //  o a vector is returned.
    //  o the return vector is created locally and passed out.
    //
    // HISTORY
    // 10 July 2001
    //         o Expanded from old definition
    //
    // 02 April 2003
    //        o Templatization.
    //

    char*                        pch_proc         = "sum(e_DOMINANCE ae_type, bool ab_abs)";
    _CMDATA                        f_sum                 = (_CMDATA)0.0;
    int                                rows = rows_get ();
    int                                cols = cols_get ();
    CMatrix<_CMDATA>            V_sum (1, 1, (_CMDATA) 0.0);
    CMatrix<_CMDATA>*                pV = 0x0;

    // some initial error checking
    if (matrix_type () != e_matrix && matrix_type() != e_square)
        _error (pch_proc, "Matrix must not be a vector!");

    switch (ae_type) {
    case e_row:
        pV = new CMatrix<_CMDATA>(rows_get (), 1, (_CMDATA) 0.0);
        for (int row = 0; row < rows; row++) {
            f_sum = (_CMDATA)0.0;
            for (int col = 0; col < cols; col++) {
                if (ab_abs)
                    f_sum += (_CMDATA)fabs ((double)_mval (row, col));
                else
                    f_sum += _mval (row, col);
            }
            pV->_mval (row, 0) = f_sum;
        }
        break;
    case e_column:
        pV = new CMatrix<_CMDATA>(1, cols_get (), (_CMDATA) 0.0);
        for (int col = 0; col < cols; col++) {
            f_sum = (_CMDATA)0.0;
            for (int row = 0; row < rows; row++) {
                if (ab_abs)
                    f_sum += (_CMDATA)fabs ((double)_mval (row, col));
                else
                    f_sum += _mval (row, col);
            }
            pV->_mval (0, col) = f_sum;
        }
        break;
    }
    V_sum.copy (*pV);
    delete        pV;
    return V_sum;
}

template<typename _CMDATA>
_CMDATA
CMatrix<_CMDATA>::min ()
{
    //
    // DESC
    //  Finds the value of the minimum element of *this.
    //
    // 02 April 2003
    //        o Templatization.
    //

    if (rows_get () <= 0 || cols_get () <= 0)
        _error ("Matrix Error!",
                 "One of the dimensions has a negative length.");
    _CMDATA        min = _mval (0, 0);
    int                r = rows_get ();
    int                c = cols_get ();
    for (int row = 0; row < r; row++)
        for (int col = 0; col < c; col++)
            if (_mval (row, col) < min)
                min = _mval (row, col);
    return (min);
}

template<typename _CMDATA>
_CMDATA
CMatrix<_CMDATA>::min(
        e_DOMINANCE     ae_type,
        int             a_index)
{
    //
    // ARGS
    //  ae_type         in              type of min to determine:
    //                                          row wise or
    //                                          col wise
    //  a_index         in              row of column for which the min is
    //                                  calculated
    //
    // DESC
    //  Determines the min of a particular column or row.
    //
    // PRECONDITIONS
    //  o *this MUST be of type e_matrix.
    //
    // POSTCONDITIONS
    //  o the min of the index'd row/column is returned.
    //
    // HISTORY
    // 12 February 2002
    //         o Initial design and coding.
    //
    // 02 April 2003
    //        o Templatization.
    //

    char*        pch_proc = "min(e_DOMINANCE ae_type, int a_index)";
    _CMDATA        f_min = _mval(0, 0);
    int                rows = rows_get ();
    int                cols = cols_get ();

    // some initial error checking
    if (matrix_type () != e_matrix && matrix_type() != e_square)
        _error (pch_proc, "Matrix must not be a vector!");

    switch (ae_type) {
    case e_row:
        f_min = _mval(a_index, 0);
        for (int col = 0; col < cols; col++)
            if(f_min > _mval(a_index, col))
                f_min = _mval(a_index, col);
        break;
    case e_column:
        f_min = _mval(0, a_index);
        for (int row = 0; row < rows; row++)
            if(f_min > _mval(row, a_index))
                f_min = _mval(row, a_index);
        break;
    }
    return f_min;
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::min (
        e_DOMINANCE     ae_type)
{
    //
    // ARGS
    //  ae_type         in              governs how the min is calculated.
    //
    // DESC
    //  Determines the min of matrix elements. The ae_type parameter
    //  governs whether the row- or column-wise elements
    //  are used in determining the min.
    //
    //  If ae_type == e_row,    then the row-wise min is determined;
    //  if ae_type == e_column, then the col-wise min is determined.
    //
    //  In both cases a vector is returned whose row/col element is the min
    //  of the corresponding row/col vector.
    //
    // PRECONDITIONS
    //  o None
    //  o This is a "simple" min, abs values are not considered.
    //
    // POSTCONDITIONS
    //  o a vector is returned.
    //  o the return vector is created locally and passed out.
    //
    // HISTORY
    // 12 February 2002
    //  o Expanded from old definition
    //
    // 02 April 2003
    //        o Templatization.
    //

    char*                        pch_proc = "min(e_DOMINANCE     ae_type)";
    _CMDATA                        f_min;
    int                                rows = rows_get ();
    int                                cols = cols_get ();
    CMatrix<_CMDATA>            V_min(1, 1, (_CMDATA)0.0);
    CMatrix<_CMDATA>*                pV = 0x0;

    // some initial error checking
    if (matrix_type () != e_matrix && matrix_type() != e_square)
        _error (pch_proc, "Matrix must not be a vector!");

    switch (ae_type) {
    case e_row:
        pV = new CMatrix<_CMDATA>(rows_get (), 1, (_CMDATA) 0.0);
        for (int row = 0; row < rows; row++) {
            f_min = _mval(row, 0);
            for (int col = 0; col < cols; col++)
                if(f_min > _mval(row, col))
                    f_min = _mval(row, col);
            pV->_mval (row, 0) = f_min;
        }
        break;
    case e_column:
        pV = new CMatrix<_CMDATA>(1, cols_get (), (_CMDATA) 0.0);
        for (int col = 0; col < cols; col++) {
            f_min = _mval(0, col);
            for (int row = 0; row < rows; row++)
                if(f_min > _mval(row, col))
                    f_min = _mval(row, col);
            pV->_mval (0, col) = f_min;
        }
        break;
    }
    V_min.copy (*pV);
    delete pV;
    return V_min;
}

template<typename _CMDATA>
_CMDATA
CMatrix<_CMDATA>::min_filter (
        CMatrix<_CMDATA>&                M_mask)
{
    //
    // ARGS
    //  M_mask          in              vector of values that will
    //                                          *not* be considered
    //                                          in determining the min
    // DESC
    //  Finds the minimum value in a matrix (the values in the M_mask
    //  matrix are not considered in determining the minimum).
    //
    // HISTORY
    //  1 Aug 2000
    //  o Initial design and coding
    //
    // 02 April 2003
    //        o Templatization.
    //
    
    _CMDATA        min = max();
    _CMDATA        val = min;
    int                notInMask = 1;
    for (int row = 0; row < rows_get (); row++)
        for (int col = 0; col < cols_get (); col++) {
            val = _mval (row, col);
            notInMask = 1;
            for (int i = 0; i < M_mask.cols_get (); i++)
                if (M_mask.val (0, i) == val)
                    notInMask = 0;
            if (notInMask && val < min)
                min = val;
        }
    return min;
}

template<typename _CMDATA>
_CMDATA
CMatrix<_CMDATA>::min_epsilon (
        _CMDATA                 af_center,
        _CMDATA                 af_range)
{
    //
    // ARGS
    //  af_center               in              center of values for range
    //                                                  calculation.
    //  af_range                in              range from center.
    //
    // DESC
    //  Finds the minium in an `epsilon' area in a matrix, i.e.
    //  all elements with absolute values greater than af_range -
    //  af_center are not considered in finding the min.
    //
    // NOTE
    //  af_range *must* be positive!
    //
    // HISTORY
    // 19 June 2000
    //  o Initial design and coding.
    //
    // 02 April 2003
    //        o Templatization.
    //

    _CMDATA        f_min = _mval (0, 0);

    if (af_range < 0)
        _error ("Matrix error!", "`af_range' must be positive.");
    for (int row = 0; row < rows_get (); row++)
        for (int col = 0; col < cols_get (); col++) {
            if (fabs ((double)_mval (row, col)) <= (af_range - af_center))
                if (_mval (row, col) < f_min)
                    f_min = _mval (row, col);
        }
    return (f_min);
}

template<typename _CMDATA>
_CMDATA
CMatrix<_CMDATA>::min_find (
        CMatrix<_CMDATA>&                where)
{
    //
    // ARGS
    //  where   in              Matrix describing the locations
    //                                  containing the min value
    // DESC
    //  Finds the location of the minimum value in a matrix.
    //
    // POSTCONDITION
    //  o Returns the minimum value and their locations in where.
    //
    // HISTORY
    // 02 April 2003
    //        o Templatization.
    //
    
    int                occurences;
    _CMDATA        f_min = min ();
    where.copy (findAll (f_min, occurences));
    return f_min;
}

template<typename _CMDATA>
_CMDATA
CMatrix<_CMDATA>::minabs ()
{
    //
    // DESC
    //  Finds the absolute minimum value in a matrix.
    //
    // POSTCONDITIONS
    //  o Returns the value (type _CMDATA) of the absolute minimum.
    //
    // HISTORY
    // 02 April 2003
    //        o Templatization.
    //
    
    if (rows_get () <= 0 || cols_get () <= 0)
        _error ("Matrix Error!",
                 "One of the dimensions has a negative length.");
    _CMDATA        min         = (_CMDATA)fabs ((double)_mval (0, 0));
    int                r         = rows_get ();
    int                c         = cols_get ();
    for (int row = 0; row < r; row++)
        for (int col = 0; col < c; col++)
            if (fabs ((double)_mval (row, col)) < min)
                min = _mval (row, col);
    return (min);
}

template<typename _CMDATA>
_CMDATA
CMatrix<_CMDATA>::max ()
{
    //
    // DESC
    //  Finds the absolute maximum value in a matrix.
    //
    // POSTCONDITIONS
    //  o Returns the value (type _CMDATA) of the absolute maximum.
    //
    // 26 March 2003
    //         o Templatization.
    //

    if (ps_mat->rows<=0 || ps_mat->cols<= 0)
        _error ("Matrix Error!",
                 "One of the dimensions has a negative length.");
    _CMDATA                max         = _mval (0, 0);
    int                        r         = ps_mat->rows;
    int                        c         = ps_mat->cols;
    for (int row = 0; row < r; row++)
        for (int col = 0; col < c; col++)
            if (_mval (row, col) > max)
                max = _mval (row, col);
    return (max);
}

template<typename _CMDATA>
_CMDATA
CMatrix<_CMDATA>::max(
        e_DOMINANCE     ae_type,
        int             a_index)
{
    //
    // ARGS
    //  ae_type         in              type of max to determine:
    //                                          row wise or
    //                                          col wise
    //  a_index         in              row of column for which the max is
    //                                  calculated
    //
    // DESC
    //  Determines the max of a particular column or row.
    //
    // PRECONDITIONS
    //  o *this MUST be of type e_matrix.
    //
    // POSTCONDITIONS
    //  o the max of the index'd row/column is returned.
    //
    // HISTORY
    // 12 February 2002
    //         o Initial design and coding.
    //
    // 02 April 2003
    //        o Templatization.
    //

    char*        pch_proc = "max";
    _CMDATA        f_max = _mval(0, 0);
    int                rows = rows_get ();
    int                cols = cols_get ();

    // some initial error checking
    if (matrix_type () != e_matrix && matrix_type() != e_square )
        _error (pch_proc, "Matrix must not be a vector!");

    switch (ae_type) {
    case e_row:
        f_max = _mval(a_index, 0);
        for (int col = 0; col < cols; col++)
            if(f_max < _mval(a_index, col))
                f_max = _mval(a_index, col);
        break;
    case e_column:
        f_max = _mval(0, a_index);
        for (int row = 0; row < rows; row++)
            if(f_max < _mval(row, a_index))
                f_max = _mval(row, a_index);
        break;
    }
    return f_max;
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::max (
        e_DOMINANCE     ae_type)
{
    //
    // ARGS
    //  ae_type         in              governs how the max is calculated.
    //
    // DESC
    //  Determines the max of matrix elements. The ae_type parameter
    //  governs whether the row- or column-wise elements
    //  are used in determining the max.
    //
    //  If ae_type == e_row,    then the row-wise max is determined;
    //  if ae_type == e_column, then the col-wise max is determined.
    //
    //  In both cases a vector is returned whose row/col element is the max
    //  of the corresponding row/col vector.
    //
    // PRECONDITIONS
    //  o None
    //  o This is a "simple" max, abs values are not considered.
    //
    // POSTCONDITIONS
    //  o a vector is returned.
    //  o the return vector is created locally and passed out.
    //
    // HISTORY
    // 12 February 2002
    //  o Expanded from old definition
    //
    // 02 April 2003
    //        o Templatization.
    //

    char*                pch_proc = "max(e_DOMINANCE ae_type)";
    _CMDATA                f_max;
    int                        rows = rows_get ();
    int                        cols = cols_get ();
    CMatrix<_CMDATA>    V_max (1, 1, (_CMDATA) 0.0);
    CMatrix<_CMDATA>*        pV = 0x0;

    // some initial error checking
    if (matrix_type () != e_matrix && matrix_type() != e_square)
        _error (pch_proc, "Matrix must not be a vector!");

    switch (ae_type) {
    case e_row:
        pV = new CMatrix<_CMDATA>(rows_get (), 1, (_CMDATA) 0.0);
        for (int row = 0; row < rows; row++) {
            f_max = _mval(row, 0);
            for (int col = 0; col < cols; col++)
                if(f_max < _mval(row, col))
                    f_max = _mval(row, col);
            pV->_mval (row, 0) = f_max;
        }
        break;
    case e_column:
        pV = new CMatrix<_CMDATA>(1, cols_get (), (_CMDATA) 0.0);
        for (int col = 0; col < cols; col++) {
            f_max = _mval(0, col);
            for (int row = 0; row < rows; row++)
                if(f_max < _mval(row, col))
                    f_max = _mval(row, col);
                pV->_mval (0, col) = f_max;
        }
        break;
    }
    V_max.copy (*pV);
    delete pV;
    return V_max;
}

template<typename _CMDATA>
_CMDATA
CMatrix<_CMDATA>::max_filter (
        CMatrix<_CMDATA>&                M_mask)
{
    //
    // ARGS
    //  M_mask          in              vector of values that will
    //                                          *not* be considered
    //                                          in determining the max
    // DESC
    //  Finds the maximum value in a matrix (the values in the M_mask
    //  matrix are not considered in determining the minimum).
    //
    // HISTORY
    // 1 Aug 2000
    //  o Initial design and coding
    //
    // 02 April 2003
    //        o Templatization.
    //

    _CMDATA        max = min ();
    _CMDATA        val = max;
    int                notInMask = 1;
    for (int row = 0; row < rows_get (); row++)
        for (int col = 0; col < cols_get (); col++) {
            val = _mval (row, col);
            notInMask = 1;
            for (int i = 0; i < M_mask.cols_get (); i++)
                if (M_mask.val (0, i) == val)
                    notInMask = 0;
            if (notInMask && val > max)
                max = val;
        }
    return max;
}

template<typename _CMDATA>
_CMDATA
CMatrix<_CMDATA>::max_epsilon (
        _CMDATA         af_center,
        _CMDATA         af_range)
{
    //
    // ARGS
    //  af_center               in              center of values for range
    //                                                  calculation.
    //  af_range                in              range from center.
    //
    // DESC
    //  Finds the maxium in an `epsilon' area in a matrix, i.e.
    //  all elements with absolute values greater than af_range -
    //  af_center are not considered in finding the max.
    //
    // PRECONDITIONS
    //  o af_range *must* be positive!
    //
    // HISTORY
    // 19 June 2000
    //         o Initial design and coding.
    //
    // 02 April 2003
    //        o Templatization.
    //

    _CMDATA        f_max = _mval (0, 0);

    if (af_range < 0)
        _error ("Matrix error!", "`af_range' must be positive.");
    for (int row = 0; row < rows_get (); row++)
        for (int col = 0; col < cols_get (); col++) {
            if (fabs ((double)_mval (row, col)) <= (af_range - af_center))
                if (_mval (row, col) > f_max)
                    f_max = _mval (row, col);
        }
    return (f_max);
}

template<typename _CMDATA>
_CMDATA
CMatrix<_CMDATA>::max_find (
        CMatrix<_CMDATA>&        where)
{
    //
    // ARGS
    //  where           in/out          Matrix describing the locations
    //                                          containing the max value
    // DESC
    //  Finds the location of the maximum value in a matrix.
    //
    // POSTCONDITIONS
    //  o The locations of the maximum values are returned in [row col]
    //    format in `where'.
    //  o The maximum value is returned in method name.
    //
    // HISTORY
    // 02 April 2003
    //        o Templatization.
    //

    int                occurences;
    _CMDATA        f_max = max ();
    where.copy (findAll (f_max, occurences));
    return f_max;
}

template<typename _CMDATA>
_CMDATA
CMatrix<_CMDATA>::maxabs ()
{
    //
    // DESC
    //  Finds the absolute maximum value in a matrix.
    //
    // POSTCONDITIONS
    //  o The maximum value is returned in method name.
    //
    // HISTORY
    // 02 April 2003
    //  o Templatization.
    //

    if (rows_get () <= 0 || cols_get () <= 0)
        _error ("Matrix Error!",
                 "One of the dimensions has a negative length.");
    _CMDATA        max = (_CMDATA)fabs ((double)_mval (0, 0));
    int                r = rows_get ();
    int                c = cols_get ();
    for (int row = 0; row < r; row++)
        for (int col = 0; col < c; col++)
            if (fabs ((double)_mval (row, col)) > max)
                max = _mval (row, col);
    return (max);
}

//////---------->
////// Miscellaneous Maths
//////---------->

template<typename _CMDATA>
_CMDATA
CMatrix<_CMDATA>::dot (
        const CMatrix<_CMDATA>&                rval)
{
    //
    // ARGS
    //  rval            in              right hand argument
    //
    // DESC
    //  Determines the dot product, i.e.
    //          A . B = (a0*bo + a1*b1 + ... an*bn)
    //
    // HISTORY
    // 10 November 2000
    //  o Initial design and coding.
    //
    // 03 April 2003
    //        o Templatization.
    //

#ifndef CMATRIX_NORANGE    
    // check for compatibility
    if (!compatible (rval))
        _error ("dot product error", "matrices are not the same size");
#endif

    _CMDATA        f_sum = (_CMDATA) 0.0;

    for (int row = 0; row < rows_get (); row++)
        for (int col = 0; col < cols_get (); col++)
            f_sum += _mval (row, col) * rval._mval (row, col);
    return f_sum;
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::cross (
        const CMatrix<_CMDATA>&                rval)
{
    //
    // ARGS
    //  rval            in              right hand argument
    //
    // DESC
    //  Determines the cross product, i.e. A X B
    //
    // PRECONDITIONS
    // o Inputs must be three element vectors.
    //
    // POSTCONDITIONS
    // o New *row* matrix is created and passed out.
    //
    // HISTORY
    // 21 June 2004
    // o Initial design and coding.
    //
    // 23 June 2004
    // o Accommodates row / col vectors. Return vector dominance is the
    //	 same as operand dominance.
    //

#ifndef CMATRIX_NORANGE    
    // check for compatibility
    if (!compatible (rval))
        _error ("cross product error", "vectors are not the same size");
    if (!rval.is_vector())
    	_error ("cross product error", "operands must be vectors");
    if (rval.size1D() != 3)
    	_error ("cross product error", "operands must be 3 elemet vectors");
#endif

    CMatrix<_CMDATA>*	pV_C;
    if(is_vectorColumn()) {
    	pV_C	= new CMatrix<_CMDATA>(3, 1);
	pV_C->val(0)=  _mval(1,0)*rval._mval(2,0)-rval._mval(1,0)*_mval(2,0);
	pV_C->val(1)=-(_mval(0,0)*rval._mval(2,0)-rval._mval(0,0)*_mval(2,0));
	pV_C->val(2)=  _mval(0,0)*rval._mval(1,0)-rval._mval(0,0)*_mval(1,0);
    } else {
    	pV_C	= new CMatrix<_CMDATA>(1, 3);
	pV_C->val(0)=  _mval(0,1)*rval._mval(0,2)-rval._mval(0,1)*_mval(0,2);
	pV_C->val(1)=-(_mval(0,0)*rval._mval(0,2)-rval._mval(0,0)*_mval(0,2));
	pV_C->val(2)=  _mval(0,0)*rval._mval(0,1)-rval._mval(0,0)*_mval(0,1);
    }
    CMatrix<_CMDATA>	V_C(*pV_C);
    delete		pV_C;
    return V_C;
}


template<typename _CMDATA>
_CMDATA
CMatrix<_CMDATA>::sqrdDistance (
        CMatrix<_CMDATA>&                M_to)
{
    //
    // ARGS
    //  M_to            in              matrix to find distance to
    //
    // DESC
    //  Finds the squared "distance" between *this matrix and
    //  M_to.
    //
    // HISTORY
    //  27 September 2000
    //  o Initial design and coding
    //
    // 03 April 2003
    //        o Templatization.
    //

    char*        pch_name         = "sqrdDistance";
    _CMDATA        f_distance         = (_CMDATA) 0.0;

    if (!compatible (M_to))
        _error (pch_name, "Matrices must have identical sizes");

    for (int row = 0; row < rows_get (); row++)
        for (int col = 0; col < cols_get (); col++)
            f_distance += (val (row, col) - M_to.val (row, col)) *
                (val (row, col) - M_to.val (row, col));
    return f_distance;
}

template<typename _CMDATA>
_CMDATA
CMatrix<_CMDATA>::distance (
        CMatrix<_CMDATA>&                M_to)
{
    //
    // ARGS
    //  M_to            in              matrix to find distance to
    //
    // DESC
    //  Finds the "distance" between *this matrix and M_to.
    //
    // HISTORY
    //  27 September 2000
    //  o Initial design and coding
    //
    // 03 April 2003
    //        o Templatization.
    //
    
    return ((_CMDATA)sqrt ((double)sqrdDistance (M_to)));
}

template<typename _CMDATA>
_CMDATA
CMatrix<_CMDATA>::determinant ()
{
    //
    // DESC
    //  Finds the determinant of *this matrix.
    //
    // POSTCONDITIONS
    //  o Minor error checking is done.
    //
    // HISTORY
    // 03 April 2003
    //        o Templatization.
    //
    
    CMatrix<_CMDATA>            indx (ps_mat->cols);
    CMatrix<_CMDATA>            B (ps_mat->cols);
    CMatrix<_CMDATA>                decomp(ps_mat->rows, ps_mat->cols);
    int                                d = 0, i;
    _CMDATA                        determinant;

    if (rows_get () != cols_get ())
        _error ("Determinant Error!",
                 "Matrix must be square for determinant.");
    decomp = lu_decompose (indx, d);
    determinant = d;
    for (i = 0; i < cols_get (); i++)
        determinant *= decomp._mval (i, i);
    return (determinant);
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::inverse ()
{
    //
    // DESC
    //  Finds the inverse of *this matrix.
    //
    // POSTCONDITIONS
    //  o A new matrix containing the result is created and returned.
    //
    // 03 April 2003
    //        o Templatization.
    //

    CMatrix<_CMDATA>                Y ("I", ps_mat->rows);
    CMatrix<_CMDATA>            trans (ps_mat->cols, ps_mat->rows, (_CMDATA) 0.0);
    CMatrix<_CMDATA>                B (ps_mat->cols);
    CMatrix<_CMDATA>            indx (ps_mat->cols);
    CMatrix<_CMDATA>                decomp(ps_mat->rows, ps_mat->cols, (_CMDATA) 0.0);
    int                                d = 0, row, col;

    if (rows_get () != cols_get ())
        _error ("Inverse Error!", "Matrix must be square for inverse.");
    decomp = lu_decompose (indx, d);
    for (col = 0; col < cols_get (); col++) {
        B._column_copy (Y, col, 0);
        decomp.lu_back_subst (indx, B);
        Y._column_copy (B, 0, col);
    }
    for (row = 0; row < ps_mat->rows; row++)
        for (col = 0; col < ps_mat->cols; col++)
            trans._mval (col, row) = Y.val (row, col);
    return (trans);
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::normalise (
        e_NORMALISATIONTYPE     ae_type          /* = e_MAGNITUDE       */,
	bool                    ab_inPlace       /* = false             */)
{
    //
    // ARGS
    //  ae_type                 in              the type of normalisation to
    //                                                  calculate
    //	ab_inPlace		in		memory management: if false,
    //							create a new matrix,
    //							leaving original intact.
    //							If true, apply operation
    //							to current matrix.
    //
    // DESC
    //  Normalises the contents of a matrix. There are several different
    //  types of normalisation (as specified by ae_type).
    //
    // PRECONDITIONS
    //  None. Works on all matrix/vector types
    //
    // POSTCONDITIONS
    //	Original matrix contents are affected by ab_inPlace. If false,
    //	this method creates a new matrix which is the normalise() of *this,
    //	else, *this is normalised.
    //
    // HISTORY
    // 24 September 2000
    // o Initial design and coding
    //
    // 07 November 2000
    // o Expanded method to allow for normalisation relative to
    //          - magnitude of matrix
    //          - maximum element in matrix
    //          - sum of all elements in matrix
    //
    // 12 February 2002
    // o Added code for MEANSTD e_type normalisation
    //
    // 03 April 2003
    // o Templatization.
    //
    // 05 March 2004
    ///	o Added ab_inPlace - prepped for extending functionality to 
    //	  volumes (of GSL_complex type).
    //	o Extended memory management to cover ab_inPlace
    //
    
    char*         pch_name         = "e_NORMALISATIONTYPE ae_type";
    
    _CMDATA     f_normalisation = (_CMDATA)0.0;
    _CMDATA     f_mean          = (_CMDATA)0.0; // The mean and std variables
    _CMDATA     f_std           = (_CMDATA)0.0; //      are really only used
                                                //      for e_MEANSTD type
        
    switch (ae_type) {
    case e_MAGNITUDE:
        f_normalisation = mag ();
        break;
    case e_MAXIMUM:
        f_normalisation = max ();
        break;
    case e_SUM:
        f_normalisation = innerSum ();
        break;
    case e_MEANSTD:
        f_normalisation = std();                // Just for error checking
        f_mean          = mean();
        f_std           = std();
        break;
    default:
        f_normalisation = mag ();
        break;
    }

    if (!f_normalisation)
        _error (pch_name, "Cannot normalise - Matrix seems singular.");
    
    CMatrix<_CMDATA>*	pM_normalised;
        
    if(!ab_inPlace)
	pM_normalised	= new CMatrix<_CMDATA>(	rows_get(), cols_get(),
						(_CMDATA)0.0);

    for (int row = 0; row < rows_get (); row++)
	for (int col = 0; col < cols_get (); col++) {
	    if(!ab_inPlace)
		(ae_type != e_MEANSTD) ?
		    pM_normalised->_mval (row, col) = _mval(row, col)/f_normalisation 
						    :
		    pM_normalised->_mval (row, col) = (_mval(row, col) - f_mean)/f_std;
	    else
		(ae_type != e_MEANSTD) ?
		    _mval (row, col) = _mval(row, col)/f_normalisation 
						    :
		    _mval (row, col) = (_mval(row, col) - f_mean)/f_std;
	}
    if(!ab_inPlace) {
	CMatrix<_CMDATA>	M_return(*pM_normalised);
	delete			pM_normalised;
	return M_return;
    } else
	return *this;
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::normalise(
        e_NORMALISATIONTYPE     ae_type,
        e_DOMINANCE             ae_dominance
) {
    //
    // ARGS
    //  ae_type         in              type of normalisation
    //  ae_dominance    in              normalise in either a row or column
    //                                          wise manner
    //
    // DESC
    //  This method returns a normalised version of *this matrix. The
    //  normalisation is either row or column dominant, which means that
    //  each row or column "vector" of the matrix is normalised relative
    //  to a row-wise or column-wise normalisation vector.
    //
    //  ae_dominance == e_column:
    //          A row vector containing the normalised values of
    //          each column is calculated, and each row is normalised
    //          to this row.
    //  ae_dominance == e_row:
    //          A column vector containing the normalised values of
    //          each row is calculated, and each column is normalised
    //          to this column.
    //
    // PRECONDITIONS
    //  o *this must be a matrix
    //
    // POSTCONDITIONS
    //  o a normalised matrix is returned. The original matrix is
    //    not affected.
    //  o note that this is subtly different from the mag(), mean(), etc.
    //    series of methods in that these methods return a *vector* while
    //    this method returns a matrix.
    //
    // HISTORY
    // 13 February 2002
    //  o Initial design and coding
    //
    // 03 April 2003
    //        o Templatization.
    //

    CMatrix<_CMDATA>                M_normalised(rows_get(), cols_get());
    CMatrix<_CMDATA>                 V(1, 1);
    CMatrix<_CMDATA>*                pV;

    if(ae_dominance == e_column)                // Depending on dominance,
        pV = new CMatrix(1, cols_get());        //      create either a row
    else                                        //      or a column vector
        pV = new CMatrix(rows_get(), 1);        //      that will contain
    V.copy(*pV);                                //      the normalised data

    switch(ae_type) {
        case e_MAGNITUDE:
            V = mag(ae_dominance);
            break;
        case e_MAXIMUM:
            V = max(ae_dominance);
            break;
        case e_SUM:
            V = sum(ae_dominance);
            break;
        case e_MEANSTD:
            V = mean(ae_dominance);
            break;
    }

    // Now replicate this vector to create a normalisation matrix
    //  - this allows us to express the normalisation operations
    //  in matrix maths rather than needing to work in an
    //  element by element manner.

    if(ae_dominance == e_row)
        V.replicate(cols_get(), e_row);
    else
        V.replicate(rows_get(), e_column);

    // For all normalisation cases other than e_MEANSTD, the
    //  desired output is merely the element by element division
    //  of *this with the created normalised matrix.
    if(ae_type != e_MEANSTD)
        M_normalised = (*this /= V);
    else {
        // For the case of e_MEANSTD, we need to subtract the mean
        //      from each element and perform an element by element
        //      division with the std. The mean matrix is contained in V
        //
        // We also need to create a std matrix analogous to V above

        CMatrix<_CMDATA>         M_std(1, 1, (_CMDATA) 0.0);
        CMatrix<_CMDATA>*        pM_std;

        if(ae_dominance == e_row) {
            pM_std = new CMatrix(1, cols_get(), (_CMDATA) 0.0);
            pM_std->copy(std(ae_dominance));
            pM_std->replicate(cols_get(), e_row);
        } else {
            pM_std = new CMatrix(rows_get(), 1, (_CMDATA) 0.0);
            pM_std->copy(std(ae_dominance));
            pM_std->replicate(rows_get(), e_column);
        }
        M_std.copy(*pM_std);
        M_normalised = ((*this - V) /= M_std);
        delete pM_std;
    }

    delete pV;
    return M_normalised;
}

template<typename _CMDATA>
void
CMatrix<_CMDATA>::flip (
        _CMDATA         af_prob,
        _CMDATA         af_false        /*= 0.0         */,
        _CMDATA         af_true         /* = 1.0        */)
{
    //
    // ARGS
    //  af_prob         in              bit-flip probability
    //  af_false        in              default for false
    //  af_true         in              default for true
    //
    // DESC
    //  Flips the contents of a matrix with af_prob probability.
    //  Naturally, it only really makes sense on <bool> type matrix
    //  data.
    //
    //  The af_false and af_true parameters were added since the definition
    //  of true and false for a given matrix might not necessarily be (1, 0).
    //  In particular, Hopfield neural networks hardlimit values to (-1, 1),
    //  in which case false should be -1 and not zero.
    //
    // HISTORY
    //  09 October 2000
    //  o Initial design and coding.
    //
    //  10 October 2000
    //  o Added af_false and af_true
    //
    // 03 April 2003
    //        o Templatization.
    //
    
    _CMDATA                f_rand;

    for (int row = 0; row < rows_get (); row++)
        for (int col = 0; col < cols_get (); col++) {
            f_rand = (_CMDATA) rand () / (_CMDATA) RAND_MAX;
            if (f_rand < af_prob) {
                if (val (row, col) == af_false)
                    val (row, col) = af_true;
                else
                    val (row, col) = af_false;
            }
        }
}

template<typename _CMDATA>
_CMDATA
CMatrix<_CMDATA>::mag ()
{
    //
    // DESC(double)
    // Find the magnitude of the matrix
    //
    // HISTORY
    // 24 Septemner 2000
    // o Initial design and coding
    //
    // 03 April 2003
    //        o Templatization.
    //
    
    _CMDATA                f_abs = (_CMDATA) 0.0;
    for (int row = 0; row < rows_get (); row++)
        for (int col = 0; col < cols_get (); col++)
            f_abs += _mval (row, col) * val (row, col);
    return (_CMDATA)sqrt ((double)f_abs);
}

template<typename _CMDATA>
_CMDATA
CMatrix<_CMDATA>::mag(
        e_DOMINANCE     ae_type,
        int             a_index)
{
    //
    // ARGS
    //  ae_type         in              type of mag to determine:
    //                                          row wise or
    //                                          col wise
    //  a_index         in              row of column for which the mag is
    //                                          calculated
    //
    // DESC
    //  Determines the mag of a particular column or row.
    //
    // PRECONDITIONS
    //  o *this MUST be of type e_matrix.
    //
    // POSTCONDITIONS
    //  o the mag of the index'd row/column is returned.
    //
    // HISTORY
    // 12 February 2002
    // o Initial design and coding.
    //
    // 03 April 2003
    //        o Templatization.
    //

    char*        pch_proc         = "max(e_DOMINANCE ae_type, int a_index)";
    _CMDATA        f_mag                 = (_CMDATA) 0.0;
    int                rows                 = rows_get ();
    int                cols                 = cols_get ();

    // some initial error checking
    if (matrix_type () != e_matrix && matrix_type() != e_square )
        _error (pch_proc, "Matrix must not be a vector!");

    switch (ae_type) {
    case e_column:
        for (int row = 0; row < rows; row++)
            f_mag += _mval(row, a_index) * _mval(row, a_index);
        break;
    case e_row:
        for (int col = 0; col < cols; col++)
            f_mag += _mval(a_index, col) * _mval(a_index, col);
        break;
    }
    return (_CMDATA)sqrt((double)f_mag);
}

template<typename _CMDATA>
CMatrix<_CMDATA> 
CMatrix<_CMDATA>::mag (
        e_DOMINANCE     ae_type)
{
    //
    // ARGS
    //  ae_type         in              governs how the mag is calculated.
    //
    // DESC
    //  Determines the mag of matrix elements. The ae_type parameter
    //  governs whether the row- or column-wise elements
    //  are used in determining the mag.
    //
    //  If ae_type == e_row,    then the row-wise mag is determined
    //          i.e. a column vector is returned, each element of which has
    //          the magnitude of the corresponding row.
    //  if ae_type == e_column, then the col-wise mag is determined
    //          i.e. a row vector is returned, each element of which has
    //          the magnitude of the corresponding column.
    //
    //  In both cases a vector is returned whose row/col element is the mag
    //  of the corresponding row/col vector.
    //
    // PRECONDITIONS
    //  o None
    //
    // POSTCONDITIONS
    //  o a vector is returned.
    //  o the return vector is created locally and passed out.
    //
    // HISTORY
    // 12 February 2002
    //  o Expanded from old definition
    //
    // 03 April 2003
    //        o Templatization.
    //
    
    char*                pch_proc         = "max(e_DOMINANCE ae_type)";
    _CMDATA                f_mag                 = (_CMDATA) 0.0;
    int                        rows                 = rows_get ();
    int                        cols                 = cols_get ();
    CMatrix<_CMDATA>    V_max(1, 1, (_CMDATA) 0.0);
    CMatrix<_CMDATA>*        pV = 0x0;

    // some initial error checking
    if (matrix_type () != e_matrix && matrix_type() != e_square)
        _error (pch_proc, "Matrix must not be a vector!");

    switch (ae_type) {
    case e_row:
        pV = new CMatrix (rows_get (), 1);
        for (int row = 0; row < rows; row++) {
            f_mag = (_CMDATA) 0.0;
            for (int col = 0; col < cols; col++)
                f_mag += _mval(row, col) * _mval(row, col);
            pV->_mval (row, 0) = (_CMDATA)sqrt((double)f_mag);
        }
        break;
    case e_column:
        pV = new CMatrix (1, cols_get ());
        for (int col = 0; col < cols; col++) {
            f_mag = (_CMDATA)0.0;
            for (int row = 0; row < rows; row++)
                f_mag += _mval(row, col) * _mval(row, col);
            pV->_mval (0, col) = (_CMDATA)sqrt((double)f_mag);
        }
        break;
    }
    V_max.copy (*pV);
    delete pV;
    return V_max;
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::abs (
    bool		ab_inPlace	/* = false */
){
    //
    // ARGS
    //	ab_inPlace		in		memory management: if false,
    //							create a new matrix,
    //							leaving original intact.
    //							If true, apply operation
    //							to current matrix.
    //
    // DESC
    //  Finds the absolute value of a matrix.
    //
    // POSTCONDITIONS
    //  o Memory management is handled by ab_inPlace. Either a new
    //	  matrix is returned, or *this is changed.
    //
    // HISTORY
    // 03 April 2003
    //   o Templatization.
    //
    // 04 March 2004
    //	o Added ab_inPlace memory management.
    //

    if (rows_get () <= 0 || cols_get () <= 0)
        _error ("abs()",
                 "One of the dimensions has a negative length.");
    
    CMatrix<_CMDATA>*	pM_positive;
    
    if(!ab_inPlace)
	pM_positive	= new CMatrix<_CMDATA>(	rows_get(), cols_get());
    
    for (int row = 0; row < rows_get(); row++)
        for (int col = 0;  col < cols_get(); col++)
	    if(!ab_inPlace)
		pM_positive->_mval (row, col) = (_CMDATA)fabs((double)_mval (row, col));
            else
		_mval (row, col) = (_CMDATA)fabs((double)_mval (row, col));
    if(!ab_inPlace) {
	CMatrix M_return(*pM_positive);
	delete pM_positive;
	return M_return;
    } else 
    return *this;
}

template<typename _CMDATA>
_CMDATA
CMatrix<_CMDATA>::variance ()
{
    //
    // DESC
    //  Compute the statistical variance of the matrix.
    //
    // POSTCONDITIONS
    //  o A single value of type _CMDATA is returned.
    //
    // HISTORY
    // 03 April 2003
    //        o Templatization.
    //
    
    _CMDATA        s_squared         = (_CMDATA)0;
    _CMDATA        mn                 = mean ();
    int                r                 = rows_get ();
    int                c                 = cols_get ();
    int                col                 = 0;
    int                row                 = 0;
    
    for (row = 0; row < r; row++)
        for (col = 0; col < c; col++) {
            _CMDATA temp = _mval (row, col) - mn;
            temp *= temp;
            s_squared += temp;
        }
    s_squared /= (row * col) - 1;
    return s_squared;
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::lu_decompose (
        CMatrix<_CMDATA>&       indx,
        int&                    d)
{
    //
    // ARGS
    //  indx            in/out          output vector that records the row
    //                                          permutation effected by the
    //                                          partial pivoting
    //  d               out             output as +-1 depending on whether the
    //                                          number of row interchanges was
    //                                          even or odd respectively.
    //
    // DESC
    //  Determine the Lower-Upper decomposition of a matrix.
    //
    //  Used in conjunction with lu_back_subst to solve linear equations or
    //  invert a matrix.
    //
    // POSTCONDITIONS
    //  o Creates and returns a matrix.
    //
    // HISTORY
    // 03 April 2003
    //        o Templatization.
    //   

    if (rows_get () != cols_get ())
        _error ("lu_decompose",
                 "Matrix must be square for L-U decomposition.");
    d                         = 1;                        // parity check
    int         row         = 0;
    int         col         = 0;
    int         k           = 0;
    int         col_max        = 0;

    _CMDATA                dum;
    _CMDATA                sum;
    _CMDATA                max;
    CMatrix<_CMDATA>    lu_decomp (ps_mat->rows, ps_mat->cols);
    CMatrix<_CMDATA>        scale_vector(ps_mat->rows);
    _deepcopy (*this, lu_decomp);
    scale_vector         = lu_decomp._scale ();
    int                r = rows_get ();
    int                c = cols_get ();

    for (row = 0; row < r; row++) {
        if (row > 0) {
            for (col = 0; col <= row - 1; col++) {
                sum = lu_decomp._mval (row, col);
                if (col > 0) {
                    for (k = 0; k <= col - 1; k++)
                        sum -=
                            lu_decomp._mval (row, k) * lu_decomp._mval (k, col);
                    lu_decomp._mval (row, col) = sum;
                }
            }
        }
        max = 0;
        for (col = row; col <= c - 1; col++) {
            sum = lu_decomp._mval (row, col);
            if (row > 0) {
                for (k = 0; k <= row - 1; k++)
                    sum -=
                        lu_decomp._mval (k, col) * lu_decomp._mval (row, k);
                lu_decomp._mval (row, col) = sum;
            }
            dum = (_CMDATA) (scale_vector._mval (col, 0) * fabs ((double)sum));
            if (dum >= max) {
                col_max = col;
                max = dum;
            }
        }
        if (row != col_max) {
            lu_decomp._column_switch (col_max, row);
            d *= -1;
            dum = scale_vector._mval (col_max, 0);
            scale_vector._mval (col_max, 0) = scale_vector._mval (row, 0);
            scale_vector._mval (row, 0) = dum;
        }
        indx._mval (row, 0) = col_max;
        if (row != r - 1) {
            if (lu_decomp._mval (row, row) == 0)
                lu_decomp._mval (row, row) = (_CMDATA)TINY;
            dum = 1 / lu_decomp._mval (row, row);
            for (col = row + 1; col <= c - 1; col++)
                lu_decomp._mval (row, col) *= dum;
        }
    }
    if (lu_decomp._mval (r - 1, c - 1) == 0)
        lu_decomp._mval (r - 1, c - 1) = (_CMDATA)TINY;
    return lu_decomp;
}

template<typename _CMDATA>
void
CMatrix<_CMDATA>::lu_back_subst (
        CMatrix<_CMDATA>&        indx,
        CMatrix<_CMDATA>&        b)
{
    //
    // ARGS
    //  indx            in/out          input as the permutation vector
    //                                          returned from lu_decompose().         
    //  b               in/out          input as the right-hand side vector
    //                                          b and returns as the solution
    //                                          vector X.
    // DESC
    //  Solves the set of N linear equations AX=B. Here "this" is the LU-decomposition
    //  of the matrix A, determined by the routine lu_decompose(). Indx is input as the
    //  permutation vector returned by lu_decompose(). B is input as the right-hand
    //  side vector B, and returns with the solution vector X. This routine takes into
    //  account the possiblity that B will begin with many zero elements, so it is
    //  efficient for use in matrix inversion.
    //
    //
    // POSTCONDITIONS
    //  o Output is returned in b.
    //
    // HISOTORY
    // 03 April 2003
    //        o Templatization.
    //
    
    if (rows_get () != cols_get ())
        _error ("Decompose Error!",
                 "Matrix must be square for lu_back_subst()");
    if (rows_get () != b.rows_get ())
        _error ("Decompose Error!",
                 "Wrong size B vector passed to lu_back_subst()");
    if (rows_get () != indx.rows_get ())
        _error ("Decompose Error!",
                 "Wrong size indx vector passed to lu_back_subst()");
    
    int                row, col, Il;
    int                ii = 0;
    _CMDATA        sum;
    
    for (col = 0; col < cols_get (); col++) {
        Il = (int) indx._mval (col, 0);
        sum = b._mval (Il, 0);
        b._mval (Il, 0) = b._mval (col, 0);
        if (ii >= 0)
            for (row = ii; row <= col - 1; row++)
                sum -= _mval (row, col) * b._mval (row, 0);
        else if (sum != 0)
            ii = col;
        b._mval (col, 0) = sum;
    }
    for (col = cols_get () - 1; col >= 0; col--) {
        sum = b._mval (col, 0);
        if (col < cols_get () - 1)
            for (row = col + 1; row <= rows_get () - 1; row++)
                sum -= _mval (row, col) * b._mval (row, 0);
        b._mval (col, 0) = sum / _mval (col, col);
    }
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::scale (
        const _CMDATA           rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements scalar multiplication, i.e. each matrix element is
    //  is multiplied by rval.
    //
    // POSTCONDITIONS
    //  o the contents of *this are overwritten.
    //
    // HISTORY
    // 20 March 2003
    // o Made GSL and CBLAS aware
    //
    // 02 April 2003
    // o Templatization.
    //
    // 30 January 2004
    //  o GSL float / int
    //

    const type_info&        Tinfo = typeid(*this);                        // RTTI data

#ifndef GSL_USE
    for (int row = 0; row < ps_mat->rows; row++)
        for (int col = 0; col < ps_mat->cols; col++)
            _mval (row, col) = _mval (row, col) * rval;
#else
    if(Tinfo==Ginfo_intMatrix)
        gsl_matrix_int_scale(   ps_mat->ps_gsl->matrix_int, rval);
    if(Tinfo==Ginfo_floatMatrix)
        gsl_matrix_float_scale( ps_mat->ps_gsl->matrix_float, rval);
    if(Tinfo==Ginfo_doubleMatrix)
        gsl_matrix_scale(       ps_mat->ps_gsl->matrix_double, rval);
#endif
    return (*this);
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::scale (
        const CMatrix<_CMDATA>&        aM)
{
    //
    // ARGS
    //  aM              in      matrix which will scale *this
    //
    // DESC
    //  This method "scales" *this matrix by aM - effectively
    //  performing an element-for-element multiplication on
    //  matrix cells
    //
    // POSTCONDITIONS
    //        o Matrix contents are overwritten.
    //
    // HISTORY
    //  27 September 2000
    //  o Initial design and coding
    //
    // 20 March 2003
    //  o Made GSL and CBLAS aware
    //
    // 02 April 2003
    //  o Templatization.
    //
    // 30 January 2004
    //  o GSL float / int
    //

    const type_info&        Tinfo = typeid(*this);                        // RTTI data

    char*        pch_name = "scale(const CMatrix<_CMDATA>& aM)";

    if (!compatible (aM))
        _error (pch_name, "Matrices must have identical sizes");

#ifndef GSL_USE
    for (int row = 0; row < rows_get (); row++)
        for (int col = 0; col < cols_get (); col++)
            _mval(row, col) = _mval(row, col) *  aM._mval(row, col);
#else
    if(Tinfo==Ginfo_intMatrix)
        gsl_matrix_int_mul_elements (   ps_mat->ps_gsl->matrix_int,
                                        aM.ps_mat->ps_gsl->matrix_int);
    if(Tinfo==Ginfo_floatMatrix)
        gsl_matrix_float_mul_elements ( ps_mat->ps_gsl->matrix_float,
                                        aM.ps_mat->ps_gsl->matrix_float);
    if(Tinfo==Ginfo_doubleMatrix)
        gsl_matrix_mul_elements (       ps_mat->ps_gsl->matrix_double,
                                        aM.ps_mat->ps_gsl->matrix_double);
#endif
    return *this;
}

template<typename _CMDATA>
_CMDATA
CMatrix<_CMDATA>::innerProd ()
{
    //
    // DESC
    //  Finds the inner product (product of all cells)
    //
    // POSTCONDITIONS
    //  o Returns a single value of type _CMDATA
    //

    _CMDATA        prod =         (_CMDATA)1.0;
    int                row, col;

    for (row = 0; row < rows_get (); row++)
        for (col = 0; col < cols_get (); col++)
            prod *= val (row, col);
    return prod;
}

template<typename _CMDATA>
_CMDATA
CMatrix<_CMDATA>::innerSum ()
{
    //
    // DESC
    //  Finds the inner summation (sum of all cells).
    //
    // POSTCONDITIONS
    //  o Returns a single value of type _CMDATA
    //
    // HISTORY
    //  07 November 2000
    //  o Creation
    //
    // 03 April 2003
    //        o Templatization.
    //

    _CMDATA        sum =         (_CMDATA) 0.0;
    int                row, col;

    for (row = 0; row < rows_get (); row++)
        for (col = 0; col < cols_get (); col++)
            sum += val (row, col);
    return sum;
}

template<typename _CMDATA>
_CMDATA
CMatrix<_CMDATA>::rms (
        CMatrix&         rval)
{
    //
    // ARGS
    //  rval            in                      right hand value
    //
    // DESC
    //  Determines the Root Mean Square (rms) value between this and rval.
    //
    //  Note that this rms is not a `true' rms, but is normalized to the values
    //  over both matrices.
    //
    //            ___________________________________________
    //           / 1          2          2          2
    //          / --- [(a - A)  + (b - B)  + (c - C)  + ... ]
    //        \/   N
    // rms = -----------------------------------------------
    //             1
    //            ---- ( a + B + b + B + .... )
    //             2N
    //
    //
    // HISTORY
    // 03 April 2003
    //        o Templatization.
    //

    int                row, col, N;
    _CMDATA        f_rms                 =  (_CMDATA) 0.0;

    N = ps_mat->rows * ps_mat->cols;
    for (row = 0; row < ps_mat->rows; row++)
        for (col = 0; col < ps_mat->cols; col++)
            f_rms += (val (row, col) - rval.val (row, col)) *
                (val (row, col) - rval.val (row, col));
    f_rms = (_CMDATA) (2 * N * sqrt ((double) f_rms / N) / (rval.sum (1) + this->sum (1)));
    return f_rms;
}

template<typename _CMDATA>
int
CMatrix<_CMDATA>::b10_convertTo (
        int             radix)
{
    //
    // ARGS
    //  radix           in              radix of number
    //
    // DESC
    //  number_radix -> number_10
    //
    //  This routine converts from a number in a given base radix to a
    //  number in base 10.
    //
    // POSTCONDITIONS
    //  o The number in base 10 is returned.
    //
    // HISTORY
    // 03 April 2003
    //        o Templatization.
    //    

    int                i;
    int                num10 = 0;

    char*        pch_name         = "b10_convertTo";
    char*        pch_errorRow         = "The passed number must be in row vector form.";
    char*        pch_errorRad         = "Radix must be greater than 2";

    if (radix < 2)
        CM_error (pch_name, pch_errorRad);
    if (matrix_type () != e_rowVector)
        CM_error (pch_name, pch_errorRow);
    for (i = 0; i < cols_get (); i++)
        num10 += (int) val (0, cols_get () - i - 1) * (int) pow ((float)radix, i);
    return num10;
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::b10_convertFrom (
        int             num,
        int             radix,
        int             forcelength     /*= 0                   */)
{
    //
    // ARGS
    //  num             in              number to convert from base 10
    //  radix           in              radix of number
    //  forcelength     in              if true, defines the size of vector
    //                                          containing the conversion.
    //                                          High order elements will be
    //                                          zero if padding is necessary.
    //
    // DESC
    //  number_10 -> number_radix
    //
    // POSTCONDITIONS
    //  o A vector that represents the number in base radix is created and
    //    returned.
    //
    // HISTORY
    //  06 October 2000
    //  o Added forcelength data
    //
    // 20 March 2003
    //        o coreMatrix_[construct/destruct] integration.
    //
    // 03 April 2003
    //        o Templatization.
    //
    // 27 January 2004
    //  o Added refcnt_dec()
    //

    char*        pch_name         = "b10_convertFrom";
    char*        pch_errorRad         = "Radix must be greater than 2";
    int                i                 = 0;
    int                k                 = 0;
    int                j, numm;


    if (radix < 2)
        CM_error (pch_name, pch_errorRad);

    // Cycle up in powers of radix until the largest exponent is found.
    while (pow ((float)radix, i++) <= num);
    i--;

    if (forcelength && forcelength < i)
        _error (pch_name, "forcelength is too small for returned matrix");

    if (num == 0 || i < 1) {
        //if (--ps_mat->refcnt == 0)
        if (refcnt_dec() == 0)
            coreMatrix_destruct();
        // Now create a new matrix ([1 X i])
        // allocate memory for the structure
        int length = forcelength && forcelength > i ? forcelength : 1;
        coreMatrix_construct(1, length);
        // first reference to this data
        ps_mat->refcnt = 1;
        scale((_CMDATA)0.0);
        return *this;
    }

    int length = forcelength && forcelength > i ? forcelength : i;

    // Check if base matrix is compatible with converted result
    if (!compatible (1, length)) {
        //if (--ps_mat->refcnt == 0)
        if (refcnt_dec() == 0)
            coreMatrix_destruct();
        // Now create a new matrix ([1 X i])
        // allocate memory for the structure
        coreMatrix_construct(1, length);
        // first reference to this data
        ps_mat->refcnt = 1;
    }

    scale((_CMDATA)0.0);
    // and implement the conversion
    numm = num;
    if (forcelength)
        k = forcelength - i;
    for (j = i - 1; j >= 0; j--) {
        val (0, k++) = (int) (numm / pow ((float)radix, j));
        numm = numm % (int) pow ((float)radix, j);
    }
    return *this;
}

template<typename _CMDATA>
bool
CMatrix<_CMDATA>::equal (
        CMatrix<_CMDATA>&                aM_B)
{
    //
    // ARGS
    //  aM_B            in              matrix that self is compared with
    //
    // DESC
    // Am I numerically equal to B?
    //
    // POSTCONDITIONS
    // o returns:
    //          0                not equal
    //          1                equal
    //
    // HISTORY
    // 03 April 2003
    //        o Templatization.
    //    

    int                row, col;

    if (!compatible (aM_B))
        return 0;

    for (row = 0; row < rows_get (); row++)
        for (col = 0; col < cols_get (); col++)
            if (val (row, col) != aM_B.val (row, col))
                return false;
    return true;
}

template<typename _CMDATA>
_CMDATA
CMatrix<_CMDATA>::distanceTo (
        CMatrix<_CMDATA>&                aM_B)
{
    //
    // ARGS
    //  aM_B            in              vector to which distance
    //                                          is to be determined
    // DESC
    //  Determines the Euclidean distance between (this) and aM_B
    //
    // HISTORY
    // 03 April 2003
    //        o Templatization.
    //
    

    char*        pch_name         = "distanceTo";
    char*        pch_noRow         = "can only operate on row vectors.";
    char*        pch_diffSize         = "vectors are different sizes";

    _CMDATA        f_distance         = (_CMDATA) 0;
    _CMDATA        f_square         = (_CMDATA) 0;

    if (!is_vectorRow() || !aM_B.is_vectorRow())
        _error (pch_name, pch_noRow);
    if (cols_get () != aM_B.cols_get ())
        _error (pch_name, pch_diffSize);
    for (int col = 0; col < cols_get (); col++)
        f_square += (val (0, col) - aM_B.val (0, col)) *
            (val (0, col) - aM_B.val (0, col));
    f_distance = (_CMDATA) sqrt ((double)f_square);
    return f_distance;
}

//////---------->
////// Operator overloads
//////---------->

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::operator= (
        const CMatrix<_CMDATA>& rval) {
    //
    // ARGS
    //  rval            in              right-hand argument to operator
    //
    // DESC
    //  Assignment operator.
    //
    // PRECONDITIONS
    //  o A "pointer" equivalence is performed, i.e. *no* deepcopy. The
    //    pointer to *this matrix's data is simply decremented (and
    //    possibly destroyed) and reset to point to rval's data (which
    //    has its refcount increased.
    // o *this must have the same size as rval
    //
    // HISTORY
    // 20 March 2003
    //        o coreMatrix_[construct/destruct] integration.
    //
    // 26 March 2003
    //         o Templatization.
    //
    // 27 January 2004
    //  o Added refcnt_dec()
    //

    char*        pch_name = "operator=";
    char*        pch_errormsg = "Matrix assignment error. Incompatible dimensions";

    if (!compatible (rval))
        _error (pch_name, pch_errormsg);

    // clean up current value:
    //if (--ps_mat->refcnt == 0)
    if (refcnt_dec()== 0)
        coreMatrix_destruct();

    // assign to new value:
    //rval.ps_mat->refcnt++;              // tell the rval it has another reference
    rval.refcnt_inc();                  // tell the rval it has another reference
    ps_mat = rval.ps_mat;               // point at the rval matrix structure
    return *this;
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::operator - ()
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements unary negation, i.e. each matrix element is
    //  is multiplied by -1.
    //
    // POSTCONDITIONS
    //  o A new matrix containing the negative is created and returned.
    //
    // 26 March 2003
    //         o Templatization.
    //

    CMatrix<_CMDATA>        M_negative (ps_mat->rows, ps_mat->cols, (_CMDATA) 0.0);
    for (int i = 0; i < ps_mat->rows; i++)
        for (int j = 0; j < ps_mat->cols; j++)
            M_negative._mval (i, j) = -_mval (i, j);
    return (M_negative);
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::operator *(
        const CMatrix<_CMDATA>&                 rval) {
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements matrix multiplication, i.e. *this is multiplied
    //  by rval.
    //
    // POSTCONDITIONS
    //  o A new matrix containing the product is created and returned.
    //
    // HISTORY
    // 20 March 2003
    // o Made GSL and CBLAS aware
    //
    // 26 March 2003
    // o Templatization.
    //
    // 25 September 2003
    // o Level 2 BLAS interface is required for matrix-vector operations! This has
    //	 horrible implications, since the API for matrix/vector multiplication is different
    //	 than matrix/matrix. Moreover, it requires *different* data types (a gsl_vector
    //	 type in the input argument as well as the return value). Maintaining such types
    //	 is really a headache.
    //
    //	 The machinery to elegantly accommodate this seems quite extensive at the moment,
    //	 hence the resort to a really ugly hack. If a GSL call is required, but a vector
    //	 argument is detected, the execution thread will "revert" to the original hard
    //	 coded matix multiplication routine.
    //
    // 30 January 2004
    //  o GSL float / int
    //

    const type_info&        Tinfo = typeid(*this);                        // RTTI data

    // check for compatibility
    if (cols_get() != rval.rows_get())
        _error ("Matrix Multiplication Error!",
                 "Columns of first matrix are not equal to rows of second.");
    CMatrix<_CMDATA>        M_result (ps_mat->rows, rval.cols_get(), (_CMDATA) 0.0);
#ifndef GSL_USE
    int                n1rows = rows_get();
    int                n1cols = cols_get();
    int                n2cols = rval.cols_get();
    for (int row = 0; row < n1rows; row++)
        for (int col = 0; col < n2cols; col++) {
            _CMDATA
                sum = 0;
            for (int i = 0; i < n1cols; i++)
                sum += _mval (row, i) * rval._mval (i, col);
            M_result._mval (row, col) = sum;
        }
#else
    if(rval.is_vector() || (Tinfo==Ginfo_intMatrix)) {
	// Ugly re-duplication of code here!
	int                n1rows = rows_get();
	int                n1cols = cols_get();
	int                n2cols = rval.cols_get();
	for (int row = 0; row < n1rows; row++)
	    for (int col = 0; col < n2cols; col++) {
		_CMDATA
		    sum = 0;
		for (int i = 0; i < n1cols; i++)
		    sum += _mval (row, i) * rval._mval (i, col);
		M_result._mval (row, col) = sum;
	    }
    } else {
        if(Tinfo==Ginfo_floatMatrix)
            gsl_blas_sgemm (CblasNoTrans, CblasNoTrans,
                  1.0,
                  this->ps_mat->ps_gsl->matrix_float,
                  rval.ps_mat->ps_gsl->matrix_float,
                  0.0,
                  M_result.ps_mat->ps_gsl->matrix_float);
        if(Tinfo==Ginfo_doubleMatrix)
            gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                  1.0,
                  this->ps_mat->ps_gsl->matrix_double,
                  rval.ps_mat->ps_gsl->matrix_double,
                  0.0,
                  M_result.ps_mat->ps_gsl->matrix_double);
    }
#endif
    return (M_result);
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::operator * (
        const _CMDATA           rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements scalar multiplication, i.e. each matrix element is
    //  is increased by rval.
    //
    // POSTCONDITIONS
    //  o A new matrix containing the scalar product is created and returned.
    //
    // HISTORY
    // 20 March 2003
    // o Made GSL and CBLAS aware
    //
    // 26 March 2003
    // o Templatization.
    //
    // 30 January 2004
    //  o GSL float / int
    //

    const type_info&        Tinfo = typeid(*this);                        // RTTI data

    CMatrix<_CMDATA>        M_result(ps_mat->rows, ps_mat->cols, (_CMDATA) 0.0);

#ifndef GSL_USE
    for (int row = 0; row < ps_mat->rows; row++)
        for (int col = 0; col < ps_mat->cols; col++)
            M_result._mval (row, col) = _mval (row, col) * rval;
#else
    M_result.copy(*this);
    if(Tinfo==Ginfo_intMatrix)
        gsl_matrix_int_scale (  M_result.ps_mat->ps_gsl->matrix_int, rval);
    if(Tinfo==Ginfo_floatMatrix)
        gsl_matrix_float_scale (M_result.ps_mat->ps_gsl->matrix_float, rval);
    if(Tinfo==Ginfo_doubleMatrix)
        gsl_matrix_scale (      M_result.ps_mat->ps_gsl->matrix_double, rval);
#endif
    return (M_result);
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::operator *= (
        const CMatrix<_CMDATA>&                        rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements element by element multiplication i.e. each matrix element
    //  multiplied by its corresponding element in rval.
    //
    // POSTCONDITIONS
    //  o The current matrix data values are overwritten.
    //
    // HISTORY
    // 20 March 2003
    // o Made GSL and CBLAS aware
    //
    // 26 March 2003
    // o Templatization.
    //
    // 30 January 2004
    //  o GSL float / int
    //

    const type_info&        Tinfo = typeid(*this);                        // RTTI data


    // check for compatibility
    char* pch_name = "operator*=";
    if (!compatible (rval))
        _error (pch_name, "Matrices must have identical sizes");

#ifndef GSL_USE
    for (int i = 0; i < ps_mat->rows; i++)
        for (int j = 0; j < ps_mat->cols; j++)
            _mval (i, j) = _mval (i, j) * rval._mval (i, j);
#else
    if(Tinfo==Ginfo_intMatrix)
        gsl_matrix_int_mul_elements(    ps_mat->ps_gsl->matrix_int,
                                (const gsl_matrix_int*)rval.ps_mat->ps_gsl->matrix_int);
    if(Tinfo==Ginfo_floatMatrix)
        gsl_matrix_float_mul_elements(  ps_mat->ps_gsl->matrix_float,
                                (const gsl_matrix_float*)rval.ps_mat->ps_gsl->matrix_float);
    if(Tinfo==Ginfo_doubleMatrix)
        gsl_matrix_mul_elements(        ps_mat->ps_gsl->matrix_double,
                                (const gsl_matrix*)rval.ps_mat->ps_gsl->matrix_double);
#endif

    return *this;
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::operator *= (
        const _CMDATA                                rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements element by constant multiplication i.e. a "scaling"
    //        of matrix values by rval.
    //
    // POSTCONDITIONS
    //  o The current matrix data values are overwritten.
    //
    // HISTORY
    // 20 March 2003
    //  o Made GSL and CBLAS aware
    //
    // 26 March 2003
    //  o Templatization.
    //
    // 30 January 2004
    //  o GSL float / int
    //

    const type_info&        Tinfo = typeid(*this);                        // RTTI data


    // check for compatibility
    char* pch_name = "operator*=";

#ifndef GSL_USE
    for (int i = 0; i < ps_mat->rows; i++)
        for (int j = 0; j < ps_mat->cols; j++)
            _mval (i, j) = _mval (i, j) * rval;
#else
    if(Tinfo==Ginfo_intMatrix)
        gsl_matrix_int_scale(   ps_mat->ps_gsl->matrix_int, rval);
    if(Tinfo==Ginfo_floatMatrix)
        gsl_matrix_float_scale( ps_mat->ps_gsl->matrix_float, rval);
    if(Tinfo==Ginfo_doubleMatrix)
        gsl_matrix_scale(       ps_mat->ps_gsl->matrix_double, rval);
#endif

    return *this;
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::operator /(
        const CMatrix<_CMDATA>&                 rval) {
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements matrix division, i.e. *this is multiplied
    //  by rval.
    //
    // POSTCONDITIONS
    //  o A new matrix containing the product is created and returned.
    //
    // HISTORY
    // 20 March 2003
    // o Made GSL and CBLAS aware
    //
    // 26 March 2003
    // o Templatization.
    //
    // 30 January 2004
    //  o GSL float / int
    //
    // 29 March 2004
    //	o Compatibility checks changed.
    //

    const type_info&        Tinfo = typeid(*this);                        // RTTI data

    // check for compatibility
    if(!compatible(rval))
        _error ("operator/(const CMatrix<_CMDATA>& rval)",
                 "Matrices must have identical dimensions");
    CMatrix<_CMDATA>        M_result (ps_mat->rows, ps_mat->cols, (_CMDATA) 0.0);

#ifndef GSL_USE
    for (int i = 0; i < ps_mat->rows; i++)
        for (int j = 0; j < ps_mat->cols; j++)
            M_result._mval (i, j) = _mval (i, j) / rval._mval (i, j);
#else
    M_result.copy(*this);
    if(Tinfo==Ginfo_intMatrix)
        gsl_matrix_int_div_elements(    M_result.ps_mat->ps_gsl->matrix_int,
                                (const gsl_matrix_int*)rval.ps_mat->ps_gsl->matrix_int);
    if(Tinfo==Ginfo_floatMatrix)
        gsl_matrix_float_div_elements(  M_result.ps_mat->ps_gsl->matrix_float,
                                (const gsl_matrix_float*)rval.ps_mat->ps_gsl->matrix_float);
    if(Tinfo==Ginfo_doubleMatrix)
        gsl_matrix_div_elements(        M_result.ps_mat->ps_gsl->matrix_double,
                                (const gsl_matrix*)rval.ps_mat->ps_gsl->matrix_double);
#endif
    return (M_result);
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::operator / (
        const _CMDATA           rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements scalar division, i.e. each matrix element is
    //  is increased by rval.
    //
    // POSTCONDITIONS
    //  o A new matrix containing the scalar product is created and returned.
    //
    // HISTORY
    // 20 March 2003
    // o Made GSL and CBLAS aware
    //
    // 26 March 2003
    // o Templatization.
    //
    // 30 January 2004
    //  o GSL float / int
    //
    // 29 March 2004
    //	o Compatibility checks changed.
    //

    const type_info&        Tinfo = typeid(*this);                        // RTTI data

    CMatrix<_CMDATA>        M_result(ps_mat->rows, ps_mat->cols, (_CMDATA) 0.0);

#ifndef GSL_USE
    for (int row = 0; row < ps_mat->rows; row++)
        for (int col = 0; col < ps_mat->cols; col++)
            M_result._mval (row, col) = _mval (row, col) / rval;
#else
    M_result.copy(*this);
    if(Tinfo==Ginfo_intMatrix)
       gsl_matrix_int_scale (   M_result.ps_mat->ps_gsl->matrix_int, (_CMDATA) (1/rval));
    if(Tinfo==Ginfo_floatMatrix)
       gsl_matrix_float_scale ( M_result.ps_mat->ps_gsl->matrix_float, (_CMDATA) (1/rval));
    if(Tinfo==Ginfo_doubleMatrix)
       gsl_matrix_scale (       M_result.ps_mat->ps_gsl->matrix_double, (_CMDATA) (1/rval));
#endif
    return (M_result);
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::operator /= (
        const CMatrix<_CMDATA>&                        rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements element by element division i.e. each matrix element
    //  multiplied by its corresponding element in rval.
    //
    // POSTCONDITIONS
    //  o The current matrix data values are overwritten.
    //
    // HISTORY
    // 20 March 2003
    //  o Made GSL and CBLAS aware
    //
    // 26 March 2003
    //  o Templatization.
    //
    // 30 January 2004
    //  o GSL float / int
    //

    const type_info&        Tinfo = typeid(*this);                        // RTTI data


    // check for compatibility
    char* pch_name = "operator/=";
    if (!compatible (rval))
        _error (pch_name, "Matrices must have identical sizes");

#ifndef GSL_USE
    for (int i = 0; i < ps_mat->rows; i++)
        for (int j = 0; j < ps_mat->cols; j++)
            _mval (i, j) = _mval (i, j) / rval._mval (i, j);
#else
    if(Tinfo==Ginfo_intMatrix)
        gsl_matrix_int_div_elements(    ps_mat->ps_gsl->matrix_int,
                                (const gsl_matrix_int*)rval.ps_mat->ps_gsl->matrix_int);
    if(Tinfo==Ginfo_floatMatrix)
        gsl_matrix_float_div_elements(  ps_mat->ps_gsl->matrix_float,
                                (const gsl_matrix_float*)rval.ps_mat->ps_gsl->matrix_float);
    if(Tinfo==Ginfo_doubleMatrix)
        gsl_matrix_div_elements(        ps_mat->ps_gsl->matrix_double,
                                (const gsl_matrix_double*)rval.ps_mat->ps_gsl->matrix_double);
#endif

    return *this;
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::operator /= (
        const _CMDATA                                rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements element by constant division i.e. a "scaling"
    //        of matrix values by 1/rval.
    //
    // POSTCONDITIONS
    //  o The current matrix data values are overwritten.
    //
    // HISTORY
    // 20 March 2003
    //  o Made GSL and CBLAS aware
    //
    // 26 March 2003
    //  o Templatization.
    //
    // 30 January 2004
    //  o GSL float / int
    //

    const type_info&        Tinfo = typeid(*this);                        // RTTI data


    // check for compatibility
    char* pch_name = "operator/=";

#ifndef GSL_USE
    for (int i = 0; i < ps_mat->rows; i++)
        for (int j = 0; j < ps_mat->cols; j++)
            _mval (i, j) = _mval (i, j) / rval;
#else
    if(Tinfo==Ginfo_intMatrix)
        gsl_matrix_int_scale(   ps_mat->ps_gsl->matrix_int, (_CMDATA) (1/rval));
    if(Tinfo==Ginfo_floatMatrix)
        gsl_matrix_float_scale( ps_mat->ps_gsl->matrix_float, (_CMDATA) (1/rval));
    if(Tinfo==Ginfo_doubleMatrix)
        gsl_matrix_scale(       ps_mat->ps_gsl->matrix_double, (_CMDATA) (1/rval));
#endif

    return *this;
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::operator - (
        const CMatrix<_CMDATA>&          rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements matrix subtraction.
    //
    // POSTCONDITIONS
    //  o A new matrix containing the sum is created and returned.
    //
    // WARNING
    //  o Be aware that the returned matrix needs to be assigned
    //    to an already created matrix otherwise memory leaks
    //    may occur.
    //
    // HISTORY
    // 20 March 2003
    //  o GSL integration.
    //
    // 26 March 2003
    //  o Templatization.
    //
    // 30 January 2004
    //  o GSL float / int
    //
    // 10 June 2004
    //	o Caught a float- doubleMatrix bug!
    //
    
    const type_info&        Tinfo = typeid(*this);                        // RTTI data

    // check for compatibility
    char* pch_name = "operator - (const CMatrix<_CMDATA>&rval)";
    if ((ps_mat->rows != rval.rows_get())
        || (ps_mat->cols !=
            rval.cols_get()))_error (pch_name,
                                       "Reason: Incompatible dimensions.");
    CMatrix<_CMDATA>         M_sum (ps_mat->rows, ps_mat->cols, (_CMDATA) 0.0);
#ifndef GSL_USE
    for (int i = 0; i < ps_mat->rows; i++)
        for (int j = 0; j < ps_mat->cols; j++)
            M_sum._mval (i, j) = _mval (i, j) - rval._mval (i, j);
#else
    M_sum.copy(*this);
    if(Tinfo==Ginfo_intMatrix)
        gsl_matrix_int_sub (
            M_sum.ps_mat->ps_gsl->matrix_int,
            rval.ps_mat->ps_gsl->matrix_int
        );
    if(Tinfo==Ginfo_floatMatrix)
        gsl_matrix_float_sub (
            M_sum.ps_mat->ps_gsl->matrix_float,
            rval.ps_mat->ps_gsl->matrix_float
        );
    if(Tinfo==Ginfo_doubleMatrix)
        gsl_matrix_sub (
            M_sum.ps_mat->ps_gsl->matrix_double,
            rval.ps_mat->ps_gsl->matrix_double
        );
#endif
    return (M_sum);
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::operator - (
        const _CMDATA           rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements scalar addition, i.e. each matrix element is
    //  is increased by rval.
    //
    // POSTCONDITIONS
    //  o A new matrix containing the sum is created and returned.
    //
    // HISTORY
    // 20 March 2003
    //  o GSL integration.
    //
    // 26 March 2003
    //  o Templatization.
    //
    // 30 January 2004
    //  o GSL float / int
    //

    const type_info&        Tinfo = typeid(*this);                        // RTTI data

    CMatrix<_CMDATA>        M_sum(ps_mat->rows, ps_mat->cols, (_CMDATA) 0.0);

#ifndef GSL_USE
    for (int i = 0; i < ps_mat->rows; i++)
        for (int j = 0; j < ps_mat->cols; j++)
            M_sum._mval (i, j) = _mval (i, j) - rval;
#else
    M_sum.copy(*this);
    if(Tinfo==Ginfo_intMatrix)
        gsl_matrix_int_add_constant(    M_sum.ps_mat->ps_gsl->matrix_int, -rval);
    if(Tinfo==Ginfo_floatMatrix)
        gsl_matrix_float_add_constant(  M_sum.ps_mat->ps_gsl->matrix_float, -rval);
    if(Tinfo==Ginfo_doubleMatrix)
        gsl_matrix_add_constant(        M_sum.ps_mat->ps_gsl->matrix_double, -rval);
#endif
    return(M_sum);
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::operator -= (
        const CMatrix<_CMDATA>&          rval) {
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements matrix - matrix subtraction.
    //
    // POSTCONDITIONS
    //  o Implements A = A-B matrix type operation, i.e. contents of current
    //          matrix are overwritten.
    //
    // HISTORY
    // 21 March 2003
    //  o GSL integration.
    //
    // 26 March 2003
    //  o Templatization.
    //
    // 30 January 2004
    //  o GSL float / int
    //
    // 10 June 2004
    //	o Caught a float- doubleMatrix bug!
    //

    const type_info&        Tinfo = typeid(*this);                        // RTTI data

#ifndef GSL_USE
    for (int i = 0; i < ps_mat->rows; i++)
        for (int j = 0; j < ps_mat->cols; j++)
            _mval (i, j) = _mval (i, j) - rval._mval(i, j);
#else
    if(Tinfo==Ginfo_intMatrix)
        gsl_matrix_int_sub (
            ps_mat->ps_gsl->matrix_int,
            rval.ps_mat->ps_gsl->matrix_int
        );
    if(Tinfo==Ginfo_floatMatrix)
        gsl_matrix_float_sub (
            ps_mat->ps_gsl->matrix_float,
            rval.ps_mat->ps_gsl->matrix_float
        );
    if(Tinfo==Ginfo_doubleMatrix)
        gsl_matrix_sub (
            ps_mat->ps_gsl->matrix_double,
            rval.ps_mat->ps_gsl->matrix_double
        );
#endif
    return (*this);
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::operator -= (
        const _CMDATA                                rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements "constant" matrix subtraction.
    //
    // POSTCONDITIONS
    //  o Implements A = A-b scalar type operation, i.e. contents of current
    //          matrix are overwritten.
    //
    // HISTORY
    // 21 March 2003
    //  o GSL integration.
    //
    // 26 March 2003
    //  o Templatization.
    //
    // 30 January 2004
    //  o GSL float / int
    //

    const type_info&        Tinfo = typeid(*this);                        // RTTI data

#ifndef GSL_USE
    for (int i = 0; i < ps_mat->rows; i++)
        for (int j = 0; j < ps_mat->cols; j++)
            _mval (i, j) = _mval (i, j) - rval;
#else
    if(Tinfo==Ginfo_intMatrix)
        gsl_matrix_int_add_constant (
            ps_mat->ps_gsl->matrix_int,
            -rval
        );
    if(Tinfo==Ginfo_floatMatrix)
        gsl_matrix_float_add_constant (
            ps_mat->ps_gsl->matrix_float,
            -rval
        );
    if(Tinfo==Ginfo_doubleMatrix)
        gsl_matrix_add_constant (
            ps_mat->ps_gsl->matrix_double,
            -rval
        );

#endif
    return (*this);
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::operator + (
        const CMatrix<_CMDATA>&          rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements matrix addition.
    //
    // POSTCONDITIONS
    //  o A new matrix containing the sum is created and returned.
    //
    // WARNING
    //  o Be aware that the returned matrix needs to be assigned
    //    to an already created matrix otherwise memory leaks
    //    may occur.
    //
    // HISTORY
    // 20 March 2003
    // o GSL integration.
    //
    // 26 March 2003
    // o Templatization.
    //
    // 30 January 2004
    //  o GSL float / int
    //

    const type_info&        Tinfo = typeid(*this);                        // RTTI data

    // check for compatibility
    char* pch_name = "operator + (const CMatrix<_CMDATA>&rval)";
    if ((ps_mat->rows != rval.rows_get())
        || (ps_mat->cols !=
            rval.cols_get()))_error (pch_name,
                                       "Reason: Incompatible dimensions.");
    CMatrix<_CMDATA>         M_sum (ps_mat->rows, ps_mat->cols, (_CMDATA) 0.0);

#ifndef GSL_USE
    for (int i = 0; i < ps_mat->rows; i++)
        for (int j = 0; j < ps_mat->cols; j++)
            M_sum._mval (i, j) = _mval (i, j) + rval._mval (i, j);
#else
    M_sum.copy(*this);
    if(Tinfo==Ginfo_intMatrix)
        gsl_matrix_int_add (
            M_sum.ps_mat->ps_gsl->matrix_int,
            rval.ps_mat->ps_gsl->matrix_int
        );
    if(Tinfo==Ginfo_floatMatrix)
        gsl_matrix_float_add (
            M_sum.ps_mat->ps_gsl->matrix_float,
            rval.ps_mat->ps_gsl->matrix_float
        );
    if(Tinfo==Ginfo_doubleMatrix)
        gsl_matrix_add (
            M_sum.ps_mat->ps_gsl->matrix_double,
            rval.ps_mat->ps_gsl->matrix_double
        );
#endif
    return (M_sum);
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::operator + (
        const _CMDATA           rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements scalar addition, i.e. each matrix element is
    //  is increased by rval.
    //
    // POSTCONDITIONS
    //  o A new matrix containing the sum is created and returned.
    //
    // HISTORY
    // 20 March 2003
    // o GSL integration.
    //
    // 26 March 2003
    // o Templatization.
    //
    // 30 January 2004
    //  o GSL float / int
    //

    const type_info&        Tinfo = typeid(*this);                        // RTTI data

    CMatrix<_CMDATA>        M_sum(ps_mat->rows, ps_mat->cols, (_CMDATA) 0.0);

#ifndef GSL_USE
    for (int i = 0; i < ps_mat->rows; i++)
        for (int j = 0; j < ps_mat->cols; j++)
            M_sum._mval (i, j) = _mval (i, j) + rval;
#else
    M_sum.copy(*this);
    if(Tinfo==Ginfo_intMatrix)
        gsl_matrix_int_add_constant(M_sum.ps_mat->ps_gsl->matrix_int, rval);
    if(Tinfo==Ginfo_floatMatrix)
        gsl_matrix_float_add_constant(M_sum.ps_mat->ps_gsl->matrix_float, rval);
    if(Tinfo==Ginfo_doubleMatrix)
        gsl_matrix_add_constant(M_sum.ps_mat->ps_gsl->matrix_double, rval);
#endif
    return(M_sum);
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::operator += (
        const CMatrix<_CMDATA>&          rval) {
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements matrix - matrix addition.
    //
    // POSTCONDITIONS
    //  o Implements A = A+B matrix type operation, i.e. contents of current
    //          matrix are overwritten.
    //
    // HISTORY
    // 21 March 2003
    // o GSL integration.
    //
    // 26 March 2003
    // o Templatization.
    //
    // 30 January 2004
    //  o GSL float / int
    //

    const type_info&        Tinfo = typeid(*this);                        // RTTI data

#ifndef GSL_USE
    for (int i = 0; i < ps_mat->rows; i++)
        for (int j = 0; j < ps_mat->cols; j++)
            _mval (i, j) = _mval (i, j) + rval._mval(i, j);
#else
    if(Tinfo==Ginfo_intMatrix)
        gsl_matrix_int_add (
            ps_mat->ps_gsl->matrix_int,
            rval.ps_mat->ps_gsl->matrix_int
        );
    if(Tinfo==Ginfo_floatMatrix)
        gsl_matrix_float_add (
            ps_mat->ps_gsl->matrix_float,
            rval.ps_mat->ps_gsl->matrix_float
        );
    if(Tinfo==Ginfo_doubleMatrix)
        gsl_matrix_add (
            ps_mat->ps_gsl->matrix_double,
            rval.ps_mat->ps_gsl->matrix_double
        );
#endif
    return (*this);
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::operator += (
        const _CMDATA                                rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements "constant" matrix addition.
    //
    // POSTCONDITIONS
    //  o Implements A = A-b scalar type operation, i.e. contents of current
    //          matrix are overwritten.
    //
    // HISTORY
    // 21 March 2003
    // o GSL integration.
    //
    // 26 March 2003
    // o Templatization.
    //
    // 30 January 2004
    //  o GSL float / int
    //

    const type_info&        Tinfo = typeid(*this);                        // RTTI data

#ifndef GSL_USE
    for (int i = 0; i < ps_mat->rows; i++)
        for (int j = 0; j < ps_mat->cols; j++)
            _mval (i, j) = _mval (i, j) + rval;
#else
    if(Tinfo==Ginfo_intMatrix)
        gsl_matrix_int_add_constant (
            ps_mat->ps_gsl->matrix_int,
            rval
        );
    if(Tinfo==Ginfo_floatMatrix)
        gsl_matrix_float_add_constant (
            ps_mat->ps_gsl->matrix_float,
            rval
        );
    if(Tinfo==Ginfo_doubleMatrix)
        gsl_matrix_add_constant (
            ps_mat->ps_gsl->matrix_double,
            rval
        );
#endif
    return (*this);
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::operator ! ()
{
    //
    // DESC
    //  Transposes *this matrix.
    //
    // POSTCONDITIONS
    //  o A new matrix containing the result is created and returned.
    //
    // HISTORY
    // 01 April 2003
    //        o Templatization.
    //

    CMatrix        trans (ps_mat->cols, ps_mat->rows, (_CMDATA) 0.0);
    int                row, col;
    for (row = 0; row < ps_mat->rows; row++)
        for (col = 0; col < ps_mat->cols; col++)
            trans._mval (col, row) = _mval (row, col);
    return (trans);
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::operator ~ ()
{
    //
    // DESC
    //  "Vectorizes" matrix, i.e. creates a new column matrix that
    //   contains the data of the original matrix.
    //
    // POSTCONDITIONS
    //  o A new matrix containing the result is created and returned.
    //
    // HISTORY
    // 03 April 2003
    //        o Templatization.
    //

    CMatrix<_CMDATA>        column (rows_get () * cols_get ());
    for (int col = 0; col < cols_get (); col++)
        for (int row = 0; row < rows_get (); row++)
            column._mval (row + col * rows_get (), 0) = _mval (row, col);
    return column;
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::operator >>= (const CMatrix<_CMDATA>&        rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements element by element comparison between two similarly sized
    //  matrices. Each cell of the resultant matrix is either 1 or 0
    //  depending on whether or not the corresponding component cells
    //  are equal or unequal.
    //
    // POSTCONDITIONS
    //  o A new matrix containing the result is created and returned.
    //
    // HISTORY
    // 03 April 2003
    //        o Templatization.
    //    

    CMatrix<_CMDATA>         compared (rows_get(), cols_get());
    for (int row = 0; row < ps_mat->rows; row++)
        for (int col = 0; col < ps_mat->cols; col++)
            compared._mval (row, col) = (_mval (row, col) > rval._mval(row, col) ? 1 : 0);
    return compared;
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::operator >>= (const _CMDATA rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements element by element comparison between two a matrix
    //  and a constant. Each cell of the resultant matrix is either 1 or 0
    //  depending on whether or not the corresponding component cells
    //  are equal or unequal.
    //
    // POSTCONDITIONS
    //  o A new matrix containing the result is created and returned.
    //
    // HISTORY
    // 03 April 2003
    //        o Templatization.
    //

    CMatrix compared (rows_get(), cols_get());
    for (int row = 0; row < ps_mat->rows; row++)
        for (int col = 0; col < ps_mat->cols; col++)
            compared._mval (row, col) = (_mval (row, col) > rval ? 1 : 0);
    return compared;
}

template<typename _CMDATA>
_CMDATA
CMatrix<_CMDATA>::operator++ ()
{
    //
    // DESC
    //  Finds the sum of all the elements in a matrix.
    //
    // POSTCONDITIONS
    //  o The sum (inner sum) is calculated and returned
    //
    // HISTORY
    // 03 April 2003
    //        o Templatization.
    //    

    _CMDATA        sum =         (_CMDATA) 0;
    for (int row = 0; row < ps_mat->rows; row++)
        for (int col = 0; col < ps_mat->cols; col++)
            sum += _mval (row, col);
    return sum;
}

//////---------->
////// Private member functions
//////---------->

template<typename _CMDATA>
void
CMatrix<_CMDATA>::_column_copy (
        CMatrix<_CMDATA>&        mm,
        int                     from_col,
        int                     to_col)
{
    //
    // ARGS
    //  mm              in              matrix from which to copy a column
    //  from_col        in              column in mm to copy
    //  to_col          in              column in *this to place copy
    //
    // DESC
    //  Copies a column from mm to here.
    //
     // 03 April 2003
    //        o Templatization.
    //
       
    if (rows_get () != mm.rows_get ())
        _error ("(Private Function) Copy Column Error!",
                 "Number of rows must be equal.");
    int                r = rows_get();
    for (int row = 0; row < r; row++)
        _mval (row, to_col) = mm._mval (row, from_col);
}

template<typename _CMDATA>
void
CMatrix<_CMDATA>::_column_switch (int col1, int col2)
{
    //
    // ARGS
    //  col1            in              column 1
    //  col2            in              column 2
    //
    // DESC
    //  Switches the contents of column col1 with col2.
    //
    // 03 April 2003
    //        o Templatization.
    //

    CMatrix<_CMDATA>    temp(rows_get());
    int                        r         = rows_get();
    int                        row;
    
    for (row = 0; row < r; row++)
        temp._mval (row, 0) = _mval (row, col1);
    for (row = 0; row < r; row++)
        _mval (row, col1) = _mval (row, col2);
    for (row = 0; row < r; row++)
        _mval (row, col2) = temp._mval (row, 0);
}

template<typename _CMDATA>
void
CMatrix<_CMDATA>::_deepcopy (
        CMatrix<_CMDATA>&        from,
        CMatrix<_CMDATA>&        to)        
{
    //
    // ARGS
    //  from            in              source matrix
    //  to              in              target matrix
    //
    // DESC
    //  Copies the contents of `from' to `to'.
    //
    // 03 April 2003
    //        o Templatization.
    //
        
    if (        from.rows_get() != to.rows_get()
        ||         from.cols_get() != to.cols_get()        )
        _error ("(Private Function) Deep-copy Error!",
                                "Matrices must be of equal dimensions.");
    int                r = from.rows_get();
    int                c = from.cols_get();
    for (int row = 0; row < r; row++)
        for (int col = 0; col < c; col++)
            to._mval (row, col) = from._mval (row, col);
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CMatrix<_CMDATA>::_scale ()
{
    //
    // DESC
    //  Private scale() method.
    //
    // POSTCONDITIONS
    //  o A matrix is created and returned.
    //
    // 03 April 2003
    //        o Templatization.
    //

    _CMDATA        temp;
    if (rows_get () <= 0 || cols_get () <= 0)
        _error ("Matrix Error!",
                 "One of the dimensions has a negative length.");
    if (rows_get () != cols_get ())
        _error ("Scale Error!", "Matrix must be square for scale.");
    CMatrix<_CMDATA>    scale_vector (rows_get ());
    int        r = rows_get ();
    int        c = cols_get ();
    for (int col = 0; col < c; col++) {
        _CMDATA            max = (_CMDATA) 0;
        for (int row = 0; row < r; row++)
            if ((temp = (_CMDATA) fabs ((double)_mval (row, col))) > max)
                max = temp;
        if (max == 0)
            _error ("Scale Error!", "Singular matrix encountered.");
        scale_vector._mval (col, 0) = 1 / max;
    }
    return scale_vector;
}


//////---------->
////// Specializations: GSL_complex
//////---------->

int CMatrix<GSL_complex>::createdCounter        = 0;

void
CMatrix<GSL_complex>::_error (
        char *apch_proc,
        char *apch_msg,
        int   code) {
    //
    // ARGS
    //  apch_proc               in              (optional) function name
    //                                                  in which error has
    //                                                  occured
    //  apch_msg                in              (optional) error message
    //                                                  text
    //  code                    in              (optional) error code
    //
    // DESC
    //  Main error reporting method for the CMatrix class.
    //
    // TODO
    //        o Wrap exception throwing around #ifdefs?
    //
    // HISTORY
    // 12 March 2003
    //        o Added id, rows, and col output
    //

    cerr << "\nFatal error encountered.\n";
    cerr << "\tCMatrix object<GSL_complex> id: " << ps_mat->id << endl;
    cerr << "\tRows: " << ps_mat->rows << ", cols: " << ps_mat->cols << endl;
    cerr << "\tCurrent function: " << "CMatrix::" << apch_proc << "\n";
    cerr << "\t" << apch_msg << "\n";
    cerr << "Throwing an exception to (this) with code " << code << endl;
    throw (this);
}

void
CMatrix <GSL_complex>::coreMatrix_construct(
        int                arows,
        int                acols
) {
    // ARGS
    //        arows                in                        number of rows in matrix
    //        acols                in                        number of cols in matrixtype
    //
    // DESCRIPTION
    //        This method is a consolidation of the main initialization code
    //        that is shared across all the *matrix* constructors.
    //
    // PRECONDITIONS
    //        o Should only be called from a constructor method.
    //
    // POSTCONDITIONS
    //        o A gsl-aware ps_mat structure is created.
    //
    // HISTORY
    // 19 March 2003
    //        o Initial design and coding.
    //

    // create the structure
    try {
        ps_mat = new matstruct;
    } catch (bad_alloc xa) {
        _error("base constructor",
                "could not allocate memory for matstruct (heap exhausted?)");
    }
    ps_mat->rows = arows;
    ps_mat->cols = acols;
    // allocate memory for the structure

    //ps_mat->ps_gsl = (gslwrap*)malloc(sizeof(gslwrap));
    ps_mat->ps_gsl                              = new gslwrap;
    ps_mat->ps_gsl->matrix_int                  = 0x0;
    ps_mat->ps_gsl->matrix_float                = 0x0;
    ps_mat->ps_gsl->matrix_double               = 0x0;
    ps_mat->ps_gsl->matrix_complex              = gsl_matrix_complex_alloc(arows, acols);
    ps_mat->ps_gsl->matrix_complex_float        = 0x0;

    if(!ps_mat->ps_gsl->matrix_complex)
            _error("base constructor", "gsl_matrix_alloc error");
    ps_mat->data = (gsl_complex**) ps_mat->ps_gsl->matrix_complex->data;

    // set first reference/id to this data
    ps_mat->refcnt = 1;
    ps_mat->id           = CMatrix::createdCounter++;
}

CMatrix<GSL_complex>::CMatrix(
        int             mrows,
        int             mcols,
        GSL_complex        initvalue) {
    //
    // ARGS
    //  mrows                   in              number of rows for matrix
    //  mcols                   in              number of cols for matrix
    //  initvalue               in              initial value for matrix
    //                                                  elements
    //
    // DESC
    //  Matrix constructor
    //
    // PRECONDITIONS
    //  o Make sure that mrows and mcols are not negative numbers!
    //
    // HISTORY
    // 16 April 2003
    //        o Template specialization.
    //

    coreMatrix_construct(mrows, mcols);

    for (int i = 0; i < mrows; i++)
        for (int j = 0; j < mcols; j++)
            _mval (i, j) = initvalue;
}

CMatrix<GSL_complex>::CMatrix(
        int                     mrows,
        int                     mcols,
        const GSL_complex*        p_initvalues) {
    //
    // ARGS
    //  mrows                   in              number of rows for matrix
    //  mcols                   in              number of cols for matrix
    //  p_initvalues            in              pointer to block of memory
    //                                                  containing initial
    //                                                  values
    //
    // DESC
    //  Matrix constructor
    //
    // PRECONDITIONS
    //  o Make sure that mrows and mcols are not negative numbers!
    //
    // HISTORY
    // 19 March 2003
    //        o Integration with gsl/coreMatrix_construct
    //

    coreMatrix_construct(mrows, mcols);

    int        c = 0;
    for (int i = 0; i < mrows; i++)
        for (int j = 0; j < mcols; j++)
            _mval (i, j) = p_initvalues[c++];
}

CMatrix<GSL_complex>::CMatrix (
        char            *matfile)
{
    //
    // ARGS
    //  matfile         in              name of file containing the
    //                                          matrix to be parsed
    //
    // DESC
    //  Constructor.
    //
    //  Creates a matrix from the specified filename.
    //
    // PRECONDITIONS
    //  o The information in *matfile must be in the format that
    //    this method expects!
    //
    // POSTCONDITIONS
    //  o Matrix is created.
    //
    // HISTORY
    // 18 April 2003
    //         o Template specialization.
    //

    const int        BUFSIZE = 120;
    FILE*        infile;
    char            buf[BUFSIZE], *cp, *cp2;
    int                rfound = 0, cfound = 0, colonsfound = 0;
    int                rows = 0, cols = 0;

    if ((infile = fopen (matfile, "r")) == 0)
        _error ("File open error! Cannot open file:-", matfile);
    while (fgets (buf, BUFSIZE, infile)) {        // parse file initialization header
        // for each header line remove comments
        if ((cp = strpbrk (buf, "#")) != NULL)        // if comment found,
            *cp = '\0';                // terminate string at comment
        if ((cp = strpbrk (buf, "r")) != NULL)
            if ((cp2 = strpbrk (cp, "=")) != NULL)
                if ((cp = strpbrk (cp2, "0123456789")) != NULL) {
                    rows = atoi (cp);
                    rfound++;
                }
        if ((cp = strpbrk (buf, "c")) != NULL)
            if ((cp2 = strpbrk (cp, "=")) != NULL)
                if ((cp = strpbrk (cp2, "0123456789")) != NULL) {
                    cols = atoi (cp);
                    cfound++;
                }
        if (strstr (buf, ":::::::") != NULL) {
            colonsfound++;
            break;                // out of while loop
        }
    }
    if (!rfound || !cfound || !colonsfound) {
        fprintf (stderr, "%s %s", matfile, nonstandard);
        exit (1);
    }

    coreMatrix_construct(rows, cols);

    for (int row = 0; row < ps_mat->rows; row++)
        for (int col = 0; col < ps_mat->cols; col++) {
            char                pch_real[20], pch_imag[20];
            double                v_real, v_imag;
            fscanf (infile, "%s %s", pch_real, pch_imag);
            v_real        = atof(pch_real);
            v_imag        = atof(pch_imag);
            GSL_complex vc_z(v_real, v_imag);
            _mval (row, col) = vc_z;
            if (ferror (infile))
                _error ("Input file error! Problem accessing file:-",
                         matfile);
        }
    fclose (infile);
}

CMatrix<GSL_complex>::CMatrix(const CMatrix<GSL_complex> & z)
{
    //
    // DESC
    //  Copy constructor - simple pointer based.
    //
    // POSTCONDITIONS
    //  o NB!! NB!! NB!!
    //    The copy constructor merely increases the reference count of the
    //    "source" data, and directs the target to point to the source.
    //    This is *not* a deepcopy!. Although allowing for fast copies between
    //    matrices, it can potentially suffer from problems relating to scope
    //    local variable variable problems. If used in function arguments
    //    or as returns out of functions, then it is not really a problem.
    //
    // HISTORY
    // 27 January 2004
    //  o Added refcnt_inc()
    //

    //z.ps_mat->refcnt++;                 // adding another reference
    z.refcnt_inc();                     // adding another reference
    ps_mat = z.ps_mat;                  // point to the new matstruct
}

void
CMatrix<GSL_complex>::coreMatrix_destruct() {
    //
    // DESC
    //        Explicitly destroys the core data stuctures of a gsl-aware *matrix*
    //        object.
    //
    // HISTORY
    // 20 March 2003
    //        o Initial design and coding.
    //
    // 28 January 2004
    //  o Corrected non-GSL destruction using delete [] properly.
    //

#ifdef GSL_USE
        gsl_matrix_complex_free(ps_mat->ps_gsl->matrix_complex);
        delete ps_mat->ps_gsl;
#else
        for (int y = 0; y < ps_mat->rows; y++)
            delete [] ps_mat->data[y];
        delete [] ps_mat->data;
#endif
        delete ps_mat;
}

CMatrix<GSL_complex>::~CMatrix () {
    //
    // DESC
    //  Destructor.
    //
    // POSTCONDITIONS
    //  Decrements reference count. If counter reaches zero, frees allocated
    //  memory.
    //
    // HISTORY
    // 27 January 2004
    //  o Added refcnt_dec()
    //

    //if (--ps_mat->refcnt == 0)
    if (refcnt_dec()== 0)
             coreMatrix_destruct();
}

GSL_complex &
CMatrix<GSL_complex>::val(
        int         row,
        int         col)
{
    //
    // ARGS
    //        row                in                row to access
    //        col                in                col to access
    //
    // DESC
    //         Element selection: can be used to read or write.
    //
    // PRECONDITIONS
    //        o Range checking can be expensive, particularly in tightly
    //          nested iterative loops. By defining CMATRIX_NORANGE range
    //          checking is disabled. This will result in tremendous performance
    //          boost, but at the risk of a single out of bound access
    //          killing the entire process.
    //
    // HISTORY
    //  3 Decemeber 2001
    //        o Added compiler directives to turn off range checking.
    //

#ifndef CMATRIX_NORANGE
    char *
        pch_name = "val";
    char *
        pch_errormsg = "Addressing error. Index out of range";
    if (!(row >= 0 && row < ps_mat->rows && col >= 0 && col < ps_mat->cols))
        _error (pch_name, pch_errormsg);
    else
#endif
    return (_mval (row, col));
}

GSL_complex &
CMatrix<GSL_complex>::operator()(
        int         row,
        int         col)
{
    //
    // ARGS
    //        row                in                row to access
    //        col                in                col to access
    //
    // DESC
    //         Element selection: can be used to read or write.
    //
    // PRECONDITIONS
    //        o Range checking can be expensive, particularly in tightly
    //          nested iterative loops. By defining CMATRIX_NORANGE range
    //          checking is disabled. This will result in tremendous performance
    //          boost, but at the risk of a single out of bound access
    //          killing the entire process.
    //
    // HISTORY
    //  3 Decemeber 2001
    //        o Added compiler directives to turn off range checking.
    //

#ifndef CMATRIX_NORANGE
    char *
        pch_name = "val";
    char *
        pch_errormsg = "Addressing error. Index out of range";
    if (!(row >= 0 && row < ps_mat->rows && col >= 0 && col < ps_mat->cols))
        _error (pch_name, pch_errormsg);
    else
#endif
    return (_mval (row, col));
}

GSL_complex &
CMatrix<GSL_complex>::val(
        int         element)
{
    //
    // ARGS
    //  element            in                element to access
    //
    // DESC
    //  Element selection: can be used to read or write. This method
    //	is for use with *vectors*.
    //
    // PRECONDITIONS
    //	o Only use with *vectors*, not matrices. 
    //  o Range checking can be expensive, particularly in tightly
    //    nested iterative loops. By defining CMATRIX_NORANGE range
    //    checking is disabled. This will result in tremendous performance
    //    boost, but at the risk of a single out of bound access
    //    killing the entire process.
    //
    // HISTORY
    //  3 Decemeber 2001
    //  o Added compiler directives to turn off range checking.
    //
    // 24 September 2003
    //	o Designed for *vectors*.
    //

#ifndef CMATRIX_NORANGE
    char* pch_name = "val";
    char* pch_errormsg = "Addressing error. Index out of range";
    if ((!is_vectorRow() && !is_vectorColumn()))
        _error (pch_name, pch_errormsg);
#endif
    if(is_vectorRow())
	return (val (0, element));
    if(is_vectorColumn())
	return (val(element, 0));
}

GSL_complex &
CMatrix<GSL_complex>::operator()(
        int         element)
{
    //
    // ARGS
    //  element            in                element to access
    //
    // DESC
    //  Element selection: can be used to read or write. This method
    //	is for use with *vectors*.
    //
    // PRECONDITIONS
    //	o Only use with *vectors*, not matrices. 
    //  o Range checking can be expensive, particularly in tightly
    //    nested iterative loops. By defining CMATRIX_NORANGE range
    //    checking is disabled. This will result in tremendous performance
    //    boost, but at the risk of a single out of bound access
    //    killing the entire process.
    //
    // HISTORY
    //  3 Decemeber 2001
    //  o Added compiler directives to turn off range checking.
    //
    // 24 September 2003
    //	o Designed for *vectors*.
    //

#ifndef CMATRIX_NORANGE
    char* pch_name = "val";
    char* pch_errormsg = "Addressing error. Index out of range";
    if ((!is_vectorRow() && !is_vectorColumn()))
        _error (pch_name, pch_errormsg);
#endif
    if(is_vectorRow())
	return (val (0, element));
    if(is_vectorColumn())
	return (val(element, 0));
}


GSL_complex &        
CMatrix<GSL_complex>::_mval(  
        int        row,
        int        col)  const {
    //
    // ARGS
    //        row                in                row to access
    //        col                in                col to access
    //
    // DESC
    //         Element selection: can be used to read or write.
    //
    // PRECONDITIONS
    //         This routine does *not* use range checking! It is an
    //        interal access routine that is used by class methods
    //        and is not available for "public" use.
    //
    // HISTORY
    // 18 March 2003
    //        o Moved from header file into main source body as part 
    //          of gsl integration.
    //
    
    return( (GSL_complex&)*gsl_matrix_complex_ptr(ps_mat->ps_gsl->matrix_complex, row, col));
}

void
CMatrix<GSL_complex>::print (
        char*                        apch_msg         /*= "ans"       */,
	int                         a_precision         /*= 6           */,
        int                        a_width                /*= 12          */,
        ios::fmtflags           a_userFormat        /*= 0           */,
        int                         a_marginOffset         /*= 0           */,
        char                         ach_left         /*= (char) 0    */,
        int                         a_tabOffset         /*= 0           */,
        bool                         ab_fancy         /*= 1           */) {
    //
    // ARGS
    //         apch_msg                in          message that prepends matrix dump
    //         a_precision             in          the precision of numerical output
    //        a_width                        in        the width of the numerical field
    //        a_userFormat                in        if non-zero, set stream format to arg
    //         a_marginOffset                 in          the amount of tab spaces prepending ch_left
    //         ach_left                in          also somewhat primitive. If non-zero, will
    //                                          print (char) left_char as left most character
    //                                          of each new line. Useful for when the matrix is
    //                                          part of a larger output dump
    //         a_tabOffset                    in          the amount of tab spaces prepending dump
    //         ab_fancy                in          still rather primitive. If false, will print
    //                                                *only* the numbers, nothing else
    //
    // DESC
    //        Print a matrix in a variety of ways. The a_userFormat parameter allows the
    //        user to send extra formatting flags to the underlying stream. Consult a
    //        C++ manual for flags.
    //
    // NOTE
    //         The `marginOffset' and `tabOffset' are used in conjunction
    //         with ch_left as follows:
    //
    //                 [<--- marginOffset --->]<ch_left>[<--- tabOffset --->]
    //
    // HISTORY
    //  18 February 2002
    //        o Changed all c-style printf's to c++-style cout
    //

    ios::fmtflags
            streamFlags        =  cout.flags();        // get current stream flags
    char        pch_tab[64]        = "";
    char        pch_margin[64]        = "";
    char        pch_indent[64]        = "";
    int                i;
    
    for (i = 0; i < a_marginOffset; i++)
        strcat (pch_margin, "\t");

    for (i = 0; i < a_tabOffset; i++)
        strcat (pch_tab, "\t");

    sprintf (pch_indent, "%s%c%s", pch_margin, ach_left, pch_tab);

    if (ab_fancy) {
        if (ps_mat->rows > 1)
            cout << pch_tab << apch_msg << " = [" << endl << pch_indent;
        else
            cout << pch_tab << apch_msg << " = [ ";
    }
    for (int row = 0; row < ps_mat->rows; row++) {
        for (int col = 0; col < ps_mat->cols; col++) {
//            cout.width(a_width);
            cout.precision(a_precision);
            cout.setf(ios::internal);
            cout.setf(ios::skipws);
            if(a_userFormat)
                    cout.flags(a_userFormat);
            cout << _mval(row, col) << "\t";
            //if(ab_fancy) cout << "| ";
        }
        if (ps_mat->rows > 1)
            cout << endl << pch_indent;
    }
    if (ab_fancy)
        if (ps_mat->rows > 1)
            cout << "]" << endl << pch_indent;
        else
            cout << "] ";
    cout.flags(streamFlags);                // Restore stream flags
}

string
CMatrix<GSL_complex>::sprint (        
        char*                        apch_msg         /*= "ans"       */,
        int                         a_precision         /*= 6           */,
        int                        a_width                /*= 12          */,
        ios::fmtflags           a_userFormat        /*= 0           */,
        int                         a_marginOffset         /*= 0           */,
        char                         ach_left         /*= (char) 0    */,
        int                         a_tabOffset         /*= 0           */,
        bool                         ab_fancy         /*= 1           */)
{
    //
    // ARGS
    //         apch_msg                in          message that prepends matrix dump
    //         a_precision             in          the precision of numerical output
    //        a_width                        in        the width of the numerical field
    //        a_userFormat                in        if non-zero, set stream format to arg
    //         a_marginOffset                 in          the amount of tab spaces prepending ch_left
    //         ach_left                in          also somewhat primitive. If non-zero, will
    //                                          print (char) left_char as left most character
    //                                          of each new line. Useful for when the matrix is
    //                                          part of a larger output dump
    //         a_tabOffset                    in          the amount of tab spaces prepending dump
    //         ab_fancy                in          still rather primitive. If false, will print
    //                                                *only* the numbers, nothing else
    // DESC
    //        Print a matrix in a variety of ways. The a_userFormat parameter allows the
    //        user to send extra formatting flags to the underlying stream. Consult a
    //        C++ manual for flags.
    //
    // NOTE
    //  The `marginOffset' and `tabOffset' are used in conjunction
    //  with ch_left as follows:
    //
    //  [<--- marginOffset --->]<ch_left>[<--- tabOffset --->]
    //
    // HISTORY
    // 18 February 2002
    //        o Changed all c-style printf's to c++-style cout
    //
    // 01 April 2003
    //        o Templatization.
    //

    stringstream        sstream("");

    ios::fmtflags            streamFlags        = cout.flags();        // get current stream flags
    char                pch_tab[64]        = "";
    char                pch_margin[64]        = "";
    char                pch_indent[64]        = "";
    int                        i;

    for (i = 0; i < a_marginOffset; i++)
        strcat (pch_margin, "\t");

    for (i = 0; i < a_tabOffset; i++)
        strcat (pch_tab, "\t");

    sprintf (pch_indent, "%s%c%s", pch_margin, ach_left, pch_tab);

    if (ab_fancy) {
        if (ps_mat->rows > 1)
            sstream << pch_tab << apch_msg << " = [" << endl << pch_indent;
        else
            sstream << pch_tab << apch_msg << " = [ ";
    }
    for (int row = 0; row < ps_mat->rows; row++) {
        for (int col = 0; col < ps_mat->cols; col++) {
//            sstream.width(a_width);
            sstream.precision(a_precision);
            sstream.setf(ios::internal);
            sstream.setf(ios::skipws);
            if(a_userFormat)
                    sstream.flags(a_userFormat);
            sstream << GSL_REAL(_mval(row, col)) << " " << GSL_IMAG(_mval(row, col)) << "\t";
        }
        if (ps_mat->rows > 1)
            sstream << endl << pch_indent;
    }
    if (ab_fancy)
        if (ps_mat->rows > 1)
            sstream << "]" << endl << pch_indent;
        else
            sstream << "] ";
    sstream.flags(streamFlags);                // Restore stream flags
    return(sstream.str());
}

void
CMatrix<GSL_complex>::fprint(
        string                        astr_filename,
        char*                        apch_msg        /*= ""          */,
        int                         a_precision         /*= 6           */,
        int                        a_width                /*= 12          */,
        ios::fmtflags           a_userFormat        /*= 0           */
) {
    //
    // ARGS
    //        astr_filename                in                filename to write matrix to
    //         apch_msg                in                  message that prepends matrix dump
    //         a_precision             in                  the precision of numerical output
    //        a_width                        in                the width of the numerical field
    //        a_userFormat                in                if non-zero, set stream format to arg
    //
    // DESC
    //        Print a matrix to file in a variety of ways. The a_userFormat parameter allows
    //        the user to send extra formatting flags to the underlying stream. Consult a
    //        C++ manual for flags.
    //
    // HISTORY
    // 18 February 2002
    //        o Changed all c-style printf's to c++-style cout
    //
    // 01 April 2003
    //        o Templatization.
    //

    char*        pch_proc        = "fprint";

    ofstream        sout(astr_filename.c_str());
    if(!sout)
            _error(pch_proc, "Could not create output file");

    time_t         time_now        = time(NULL);
    string        str_time        = ctime(&time_now);
    // strip trailing \n
    str_time.erase(str_time.length()-1);
    string        str_hostname(getenv("HOSTNAME"));
    string        str_user(getenv("USER"));

    sout << "#"                                        << endl;
    sout << "# Standard CMatrix<GSL_complex> save file."        << endl;
    sout << "#        Created by\t"         << str_user         << endl;
    sout << "#        Date stamp:\t"         << str_time        << endl;
    sout << "#        Machine name:\t"<< str_hostname        << endl;
    if(strlen(apch_msg))
    sout << "#        Matrix name:\t"        << apch_msg        << endl;
    sout << "#"                                        << endl;
    sout << "rows="         << rows_get();
    sout << " columns=" << cols_get()                 << endl;
    sout << ":::::::"                                << endl;

    string        str_data        = sprint("",
                                            a_precision,
                                            a_width,
                                            a_userFormat,
                                            0, 0, 0,
                                            0);

    sout << str_data;
    sout.close();
}

CMatrix<GSL_complex>
CMatrix<GSL_complex>::copy (
        const CMatrix<GSL_complex>&          source)
{
    // ARGS
    //  source                  in              matrix whose contents to copy
    //                                                  into *this
    //
    // DESC
    //  Explicit deepcopy.
    //
    // POSTCONDITIONS.
    //  o The data in source is copied over into *this. Space for this data
    //    is created as in necessary.
    //
    // HISTORY
    // 20 March 2003
    //        o coreMatrix_[construct/destruct] integration.
    //
    // 26 March 2003
    //         o Templatization.
    //

    int                row, col, y;

    if (!compatible (source)) {
            coreMatrix_destruct();
        coreMatrix_construct(source.rows_get(), source.cols_get());
    }
    for (row = 0; row < ps_mat->rows; row++)
        for (col = 0; col < ps_mat->cols; col++)
            _mval (row, col) = source._mval (row, col);
    return *this;
}

CMatrix<GSL_complex>
CMatrix<GSL_complex>::matrix_remove(
        CMatrix<GSL_complex>*&  apM,
        int                     a_trow,
        int                     a_tcol,
        int                     a_rows,
        int                     a_cols) {
    //
    // ARGS
    //  apM             in/out          structure to contain removed
    //                                          sub-matrix
    //  a_trow, a_tcol  in              top left coordinate in base matrix
    //                                          of submatrix
    //  a_rows, a_cols  in              relative to trow and tcol, number
    //                                          of rows and cols to extract
    //                                          into submatrix
    //
    // DESC
    //  Removes a submatrix from a base matrix
    //
    // PRECONDITIONS
    //  o 0 < a_trow < rows_get()
    //  o 0 < a_tcol < cols_get()
    //  o a_trow + a_rows <= rows_get() + 1
    //  o a_tcol + a_cols <= cols_get() + 1
    //
    // POSTCONDITIONS
    //  o Returns submatrix (in name) and pointer to submatrix (in argument list).
    //
    // HISTORY
    // 03 December 2000
    //  o Changed calling parameters to explicitly allow for
    //    holding memory structure.
    //
    // 31 March 2003
    //  o Templatization.
    //
    // 07 September 2003
    //  o GSL_complex
    //

    char*       pch_name = "remove_matrix";
    int         row, col;

    if (a_trow + a_rows > ps_mat->rows + 1 ||
        a_tcol + a_cols > ps_mat->cols + 1)
            _error (pch_name, "Specified target dimensions invalid!");

    if (!apM->compatible (a_rows, a_cols)) {
        delete        apM;
        apM         = new CMatrix<GSL_complex>(a_rows, a_cols, (GSL_complex) 0.0);
    }

    for (row = 0; row < a_rows; row++)
        for (col = 0; col < a_cols; col++) {
            apM->val (row, col) = val (row + a_trow, col + a_tcol);
        }
    return *apM;
}

CMatrix<GSL_complex>
CMatrix<GSL_complex>::matrix_replace(
        int                     trow,
        int                     tcol,
        CMatrix<GSL_complex>&   aM_replacement)
{
    //
    // ARGS
    //  trow, tcol      in      top left coordinate in base matrix
    //                                  where replacement will be added.
    //  aM_replacement  in      replacement matrix
    //
    // DESC
    //  Overwrites part of a matrix with another matrix.
    //
    // 01 April 2003
    //  o Templatization.
    // 07 September 2003
    //  o GSL_complex
    //

    char*       pch_name        = "matrix_replace";
    char*       pch_errorSize   = "Replacement cannot fit into base.";
    char*       pch_errorCoord  = "Specified insert point is out of range.";
    int         row, col;

    if ((trow < 0) || (tcol < 0) ||
        (trow > ps_mat->rows) || (tcol > ps_mat->cols))
            _error (pch_name, pch_errorCoord);
    if ((trow + aM_replacement.rows_get () > ps_mat->rows) ||
        (tcol + aM_replacement.cols_get () > ps_mat->cols))
        _error (pch_name, pch_errorSize);
    for (row = 0; row < aM_replacement.rows_get (); row++) {
        for (col = 0; col < aM_replacement.cols_get (); col++) {
            val (row + trow, col + tcol) = aM_replacement.val (row, col);
        }
    }
    return *this;
}

CMatrix<GSL_complex>
CMatrix<GSL_complex>::replicate(
        int             times,
        e_DOMINANCE     dir) {
    //
    // ARGS
    //  times           in      number of times the matrix is replicated
    //  dir             in      direction of replication
    //
    // DESC
    //  Replicates a matrix
    //
    //  This routine replaces the base matrix with a much larger
    //  copy of itself, made up of (times) instances of the original
    //  base. The replications occurs in either a column- or row-dominant
    //  direction.
    //
    //  In other words, if we have an [N x M] matrix that we wish
    //  to replicate p times in a row-dominant fashion, this routine
    //  replaces the base [N x M] matrix with a new matrix:
    //
    //          [ [ N x M] [N x M] ... [N x M] ]
    //
    //  (The column dominant is simply the transpose of the above)
    //
    // PRECONDITIONS
    //        o *this is destroyed and rebuilt!
    //
    // HISTORY
    // 01 April 2003
    //  o Templatization.
    //


    int                     y, i;
    int                     mrows, mcols;
    CMatrix<GSL_complex>    M(ps_mat->rows, ps_mat->cols, (GSL_complex) 0.0);

    // Make a deepcopy of myself
    M.copy (*this);

    // Track the current matrix size...
    mrows = ps_mat->rows;
    mcols = ps_mat->cols;

    // Then destroy and create a new matrix in its core structures
    coreMatrix_destruct();
    switch (dir) {
    case e_row:
        coreMatrix_construct(mrows, mcols*times);
        break;
    case e_column:
        coreMatrix_construct(mrows*times, mcols);
        break;
    }

    for (i = 0; i < times; i++) {
        switch (dir) {
        case e_row:
            matrix_replace (0, i * mcols, M);
            break;
        case e_column:
            matrix_replace (i * mrows, 0, M);
            break;
        }
    }
    return *this;
}

CMatrix<GSL_complex>
CMatrix<GSL_complex>::zeroPad(
    e_DOMINANCE         dir,
    int                 size) {
    //
    // ARGS
    //  dir             in          the dimension along which to zero pad
    //  size            in          size of padding along a dimension
    //
    // DESC
    //  This method "zeropads" a matrix with 2 submatrices of padding
    //  length 'size' along the direction 'dir'.
    //
    //  Zero padding is necessary before calling fft-type methods on
    //  non-power of two lengthed dimensions.
    //
    // PRECONDITIONS
    //  o None.
    //
    // POSTCONDITIONS
    //  o Non destructive
    //
    // HISTORY
    // 08 September 2003
    //  o Initial design and coding.
    //

    GSL_complex             z_zero(0, 0);
    int                     rows            = rows_get();
    int                     cols            = cols_get();
    int                     i, j;
    int                     trow, tcol;

    CMatrix<GSL_complex>    M(rows, cols, (GSL_complex) 0.0);
    CMatrix<GSL_complex>*   pM_zeroPadded;

    // Make a deepcopy of myself
    M.copy (*this);

    switch (dir) {
    case e_column:
        pM_zeroPadded   = new CMatrix<GSL_complex>(rows, cols+2*size, z_zero);
        trow    = 0;
        tcol    = size;
        break;
    case e_row:
        pM_zeroPadded   = new CMatrix<GSL_complex>(rows+2*size, cols, z_zero);
        trow    = size;
        tcol    = 0;
        break;
    }

    // and insert the original data in the center
    pM_zeroPadded->matrix_replace(trow, tcol, M);
    CMatrix<GSL_complex>    M_padded(*pM_zeroPadded);
    delete  pM_zeroPadded;
    return M_padded;
}

bool
CMatrix<GSL_complex>::is_vector() const
{
    //
    // Checks if object is a vector
    //
    return (is_vectorRow() || is_vectorColumn());
}

bool
CMatrix<GSL_complex>::is_vectorColumn() const
{
    //
    // Checks if object is column vector
    //
    return ((ps_mat->cols == 1) && (ps_mat->rows > 1));
}

bool
CMatrix<GSL_complex>::is_vectorRow() const
{
    //
    // Checks if object is row vector
    //
    return ((ps_mat->rows == 1) && (ps_mat->cols > 1));
}

bool
CMatrix<GSL_complex>::is_square() const
{
    //
    // Checks if matrix is square
    //
    return (ps_mat->cols == ps_mat->rows);
}

int
CMatrix<GSL_complex>::size1D ()		const
{
    //
    // DESC
    //  Returns the 1D size of a matrix, i.e. the product of the rows
    //  and columns
    //
    // HISTORY
    //  14 November 2000
    //  o Initial design and coding.
    //

    return (rows_get () * cols_get ());
}

CMatrix<int>
CMatrix<GSL_complex>::size()		const
{
    //
    // DESC
    //  Returns a pointer to a matrix containing the size of *this.
    //
    // HISTORY
    //  23 September 2000
    //  o Initial design and coding
    //
    //  14 November 2000
    //  o Possible dangling / lost pointer!
    //
    // 31 March 2003
    //         o Templatization.
    //

    CMatrix<int>        iM_size(1, 2);

    // rows
    iM_size.val (0, 0) = rows_get ();
    //cols
    iM_size.val (0, 1) = cols_get ();

    return iM_size;
}

e_MATRIXTYPE
CMatrix<GSL_complex>::matrix_type () {
    //
    // Returns the matrix type
    // of its argument
    //

    if (is_vectorColumn())
        return e_columnVector;
    if (is_vectorRow())
        return e_rowVector;
    if (is_square())
        return e_square;
    return e_matrix;
}

bool
CMatrix<GSL_complex>::compatible (
        int                     dim1,
        int                     dim2) {
    //
    // ARGS
    //  dim1, dim2              in              dimensions under consideration
    //
    // DESC
    //  Checks if base is compatible with passed arguments
    //
    // 26 March 2003
    //         o Templatization.
    //

    if (ps_mat->rows != dim1 || ps_mat->cols != dim2)
        return false;
    return true;
}

bool
CMatrix<GSL_complex>::compatible (
        const CMatrix<GSL_complex>&          M)
{
    //
    // ARGS
    //
    //  M                       in      Matrix under consideration
    //
    // DESC
    //  Checks if base is compatible with passed arguments
    //
    // 26 March 2003
    //         o Templatization.
    //

    int                dim1, dim2;

    dim1         = M.rows_get();
    dim2         = M.cols_get();
    if (ps_mat->rows != dim1 || ps_mat->cols != dim2)
        return false;
    return true;
}

//////---------->
////// Miscellaneous Maths
//////---------->

CMatrix<GSL_complex>
CMatrix<GSL_complex>::inverse ()
{
    //
    // DESC
    //  Finds the inverse of *this matrix.
    //
    // POSTCONDITIONS
    //  o A new matrix containing the result is created and returned.
    //
    // HISTORY
    // 19 April 2003
    //	o Template specialization.
    //
    // 22 July 2004
    //	o gsl_permutation_free() added
    //

    CMatrix<GSL_complex>        zM_inverse(ps_mat->cols, ps_mat->rows);

    if (rows_get () != cols_get ())
        _error ("Inverse Error!", "Matrix must be square for inverse.");

    gsl_permutation*    gp_permute = gsl_permutation_alloc(rows_get());
    int                 s;

    gsl_linalg_complex_LU_decomp(ps_mat->ps_gsl->matrix_complex, gp_permute, &s);
    gsl_linalg_complex_LU_invert(ps_mat->ps_gsl->matrix_complex, gp_permute,
                                    zM_inverse.ps_mat->ps_gsl->matrix_complex);
				    
    gsl_permutation_free(gp_permute);
    return(zM_inverse);
}

CMatrix<GSL_complex>
CMatrix<GSL_complex>::fft1D(
        e_FFTDIR                e_direction     /*= e_forward   */,
        bool                    b_MKLwrap       /*= false       */) {
    //
    // ARGS
    //  e_direction     in              type of FFT operation to perform:
    //                                          - e_forward FFT
    //                                          - e_backward FFT
    //                                          - e_inverse FFT
    //  b_MKLwrap       in              if true, wrap around low level MKL
    //                                          routines directly as opposed
    //                                          to wrapping around the default
    //                                          GSL function calls.
    //
    // DESC
    //  Wrap around GSL (or MKL)-level complex fast Fourier transforms.
    //
    // PRECONDITIONS
    // o Note that if wrapping around MKL routines, the length of the
    //   vector *must* be a power of two.
    //
    // POSTCONDITIONS
    // o Internal contents of vector are *not* preserved!
    // o *this is returned.
    //
    // HISTORY
    // 22 April 2003
    //  o Initial design and coding. Currently, e_direction is ignored
    //    and only a forward FFT is calculated.
    //
    // 03 September 2003
    //  o After long and trying testing, it was decided that the most
    //    effective mechanism for plugging a slow memory leak was to
    //    remove the "inPlace/exPlace" ability. FFT will thus always be
    //    "inPlace".
    //

    char*       pch_proc = "fft1D(e_FFTDIR e_direction)";
    if(!is_vectorRow())
        _error(pch_proc, "Input must be a row vector.", 1);

    int                         i;
    const int                   n       = cols_get();

    if(!b_MKLwrap) {
        gsl_fft_complex_wavetable*  pgsl_wavetable;
        gsl_fft_complex_workspace*  pgsl_workspace;

        // First, create packed array from the complex data set.
        double*     pv_dataPacked = new double[2*n];
        for (i = 0; i < n; i++) {
            pv_dataPacked[2*i]      = GSL_REAL(_mval(0, i));
            pv_dataPacked[2*i+1]    = GSL_IMAG(_mval(0, i));
        }

        // Create the wave tables and workspace
        pgsl_wavetable = gsl_fft_complex_wavetable_alloc (n);
        pgsl_workspace = gsl_fft_complex_workspace_alloc (n);

        // Perform the required operation
        gsl_fft_complex_forward(        pv_dataPacked, 1, n,
                                        pgsl_wavetable,
                                        pgsl_workspace);

        // Pack the results in the return CMatrix<GSL_complex>
        for (i = 0; i < n; i++) {
            GSL_complex     z_data(0, 0);
            GSL_SET_COMPLEX(&z_data, pv_dataPacked[2*i], pv_dataPacked[2*i+1]);
            _mval(0, i)             = z_data;
        }

        gsl_fft_complex_wavetable_free(     pgsl_wavetable);
        gsl_fft_complex_workspace_free(     pgsl_workspace);
        delete                              pv_dataPacked;

    } else {
        int                  direction       = FFTW_FORWARD;

        switch(e_direction) {
            case e_forward:
                direction = FFTW_FORWARD;
                break;
            case e_backward:
                direction = FFTW_BACKWARD;
                break;
            case e_inverse:
                direction = FFTW_BACKWARD;
                break;
         }

         fftw_complex *in, *out;
         fftw_plan p;
         
         in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
         out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
         
         // Fill the input vectors
         for (i = 0; i < n; i++) {
            in[i][0]          = GSL_REAL(_mval(0, i));
            in[i][1]          = GSL_IMAG(_mval(0, i));
         }
         
         p = fftw_plan_dft_1d(n, in, out, direction, FFTW_ESTIMATE);
     
         // Pack the results in the return CMatrix<GSL_complex>
         for (i = 0; i < n; i++) {
             GSL_complex     z_data(0, 0);
             GSL_SET_COMPLEX(&z_data, out[i][0], out[i][1]);
             _mval(0, i)             = z_data;
         }

         fftw_destroy_plan(p);
         fftw_free(in); 
         fftw_free(out);
        
    }
    return                  *this;
}

CMatrix<GSL_complex>
CMatrix<GSL_complex>::fft2D(
        e_FFTDIR                e_direction     /*= e_forward   */,
        bool                    b_MKLwrap       /*= false       */) {
    ///
    /// ARGS
    ///  e_direction     in              type of FFT operation to perform:
    ///                                          - e_forward FFT
    ///                                          - e_backward FFT
    ///                                          - e_inverse FFT
    ///  b_MKLwrap       in              if true, wrap around low level MKL
    ///                                          routines directly as opposed
    ///                                          to wrapping around the default
    ///                                          GSL function calls.
    ///
    /// DESC
    ///  Wrap around GSL (or MKL)-level complex fast Fourier transforms.
    ///
    /// PRECONDITIONS
    /// o Note that if wrapping around MKL routines, the length of the
    ///   rows and columns *must* be a power of two.
    /// o It appears as though the GSL routines do not offer native 2D support.
    ///
    /// POSTCONDITIONS
    /// o Internal matrix contents are overwritten!
    /// o *this is returned.
    ///
    /// HISTORY
    /// 26 August 2003
    ///  o Initial design and coding. Currently, e_direction is ignored
    ///    and only a forward FFT is calculated.
    ///  o Only MKL support is coded.
    ///
    /// 03 September 2003
    ///  o After long and trying testing, it was decided that the most
    ///    effective mechanism for plugging a slow memory leak was to
    ///    remove the "inPlace/exPlace" ability. FFT will thus always be
    ///    "inPlace".
    ///

    char*       pch_proc = "fft2D(e_FFTDIR e_direction)";

    int                         i, j;
    const int                   n       = cols_get();
    const int                   m       = rows_get();

    if(!b_MKLwrap) {
        _error(pch_proc, "2D fft is not yet defined for GSL operations");
    } else {
        
        int                  direction       = FFTW_FORWARD;

        switch(e_direction) {
            case e_forward:
                direction = FFTW_FORWARD;
                break;
            case e_backward:
                direction = FFTW_BACKWARD;
                break;
            case e_inverse:
                direction = FFTW_BACKWARD;
                break;
         }

         fftw_complex *in, *out;
         fftw_plan p;

         in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n * m);
         out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n * m);
         
         for (i = 0; i<m; i++) {
             for(j = 0; j<n; j++) {
                 in[n*i+j][0]          = GSL_REAL(_mval(i, j));
                 in[n*i+j][1]          = GSL_IMAG(_mval(i, j));
             }
         }
   
         p = fftw_plan_dft_2d(n, m, in, out, direction, FFTW_ESTIMATE);

         // Pack the results in the return CMatrix<GSL_complex>
         for (i=0; i<m; i++) {
             for(j=0; j<n; j++) {
                 GSL_complex     z_data(0, 0);
                 GSL_SET_COMPLEX(&z_data, out[n*i+j][0], out[n*i+j][1]);
                 _mval(i, j)             = z_data;
             }
         }
             
         fftw_destroy_plan(p);
         fftw_free(in); 
         fftw_free(out);               
    }
    return                  *this;
}

CMatrix<GSL_complex>
CMatrix<GSL_complex>::ifftshiftNewMem() {
    //
    // DESC
    //  This method implements an ifftshift operation, analogous to the MatLAB
    //  function of the same name.
    //
    //  Note that although "fft" appears in the name, no actual fft is performed;
    //  this method merely "reorganises" the internal data contents, shifting the
    //  center of k-space.
    //
    //  In the case of a matrix, the quadrants are labelled as follows:
    //
    //      [B] [C]
    //      [D] [A]
    //
    //  and are "shifted" to
    //
    //      [A] [D]
    //      [C] [B]
    //
    //  where [A] is swapped with [B], and [C] is swapped with [D].
    //  For vectors, this is [B] [A] and [[B] [A]]' (transpose) with
    //  [A] and [B] swapped.
    //
    // POSTCONDITIONS
    //  o This is a non-destructive operation. An internal matrix is
    //    created locally and returned (using the copy constructor).
    //    Make certain to capture the returned matrix!
    //
    // HISTORY
    // 06 September 2003
    //  o Initial design and coding.
    //

    CMatrix<GSL_complex> M_dispatch(*this);

    M_dispatch = shiftNewMem(-1);

    return(M_dispatch);
}

CMatrix<GSL_complex>
CMatrix<GSL_complex>::fftshiftNewMem() {
    //
    // DESC
    //  This method implements an fftshift operation, analogous to the MatLAB
    //  function of the same name.
    //
    //  Note that although "fft" appears in the name, no actual fft is performed;
    //  this method merely "reorganises" the internal data contents, shifting the
    //  center of k-space.
    //
    //  In the case of a matrix, the quadrants are labelled as follows:
    //
    //      [B] [C]
    //      [D] [A]
    //
    //  and are "shifted" to
    //
    //      [A] [D]
    //      [C] [B]
    //
    //  where [A] is swapped with [B], and [C] is swapped with [D].
    //  For vectors, this is [B] [A] and [[B] [A]]' (transpose) with
    //  [A] and [B] swapped.
    //
    // POSTCONDITIONS
    //  o This is a non-destructive operation. An internal matrix is
    //    created locally and returned (using the copy constructor).
    //    Make certain to capture the returned matrix!
    //
    // HISTORY
    // 06 September 2003
    //  o Initial design and coding.
    //

    CMatrix<GSL_complex> M_dispatch(*this);

    M_dispatch  = shiftNewMem(+1);

    return(M_dispatch);
}

CMatrix<GSL_complex>
CMatrix<GSL_complex>::shiftNewMem(
    int     a_pivotOffset
) {
    //
    // ARGS
    //  a_pivotOffset       in          offset to pivot point
    //                                      about which matrix
    //                                      elements are swapped.
    //
    // DESC
    //  This method performs the actual innards of both the
    //  fftshift and iffshift operation. Conceptually, both
    //  methods are identical, differing only in the pivot
    //  point offset.
    //
    // PRECONDITIONS
    //  o None.
    //
    //
    // POSTCONDITIONS
    //  o This is a non-destructive operation. An internal matrix is
    //    created locally and returned (using the copy constructor).
    //    Make certain to capture the returned matrix!
    //
    // HISTORY
    // 08 September 2003
    //  o Initial design and coding.
    //

    CMatrix<GSL_complex>*   pM_shifted = new CMatrix<GSL_complex>(1, 1);
    pM_shifted->copy(*this);

    if(is_vectorRow()) {
        // shift for row vectors
        int     Atcol, Acols;
        int     Btcol, Bcols;
        int     Pcol;     // Insertion pivot
        int     cols            = cols_get();

        if(isOdd(cols)) {
            Atcol   = (cols+a_pivotOffset)/2;
            Pcol    = Atcol-a_pivotOffset;
        } else {
            Atcol   = cols/2;
            Pcol    = Atcol;
        }
        Acols   = cols-Atcol;
        Bcols   = Atcol;
        Btcol   = 0;

        CMatrix<GSL_complex>* pM_quadA  = new CMatrix<GSL_complex>(1, Acols);
        CMatrix<GSL_complex>* pM_quadB  = new CMatrix<GSL_complex>(1, Bcols);

        matrix_remove(pM_quadA, 0, Atcol, 1, Acols);
        matrix_remove(pM_quadB, 0, Btcol, 1, Bcols);

        pM_shifted->matrix_replace(0,   0,      *pM_quadA);
        pM_shifted->matrix_replace(0,   Pcol,   *pM_quadB);

        delete pM_quadA;
        delete pM_quadB;

    }

    if(is_vectorColumn()) {
        // shift for col vectors
        int     Atrow, Arows;
        int     Btrow, Brows;
        int     Prow;     // Insertion pivot
        int     rows            = rows_get();

        if(isOdd(rows)) {
            Atrow   = (rows+a_pivotOffset)/2;
            Prow    = Atrow-a_pivotOffset;
        } else {
            Atrow   = rows/2;
            Prow    = Atrow;
        }
        Arows   = rows-Atrow;
        Brows   = Atrow;
        Btrow   = 0;

        CMatrix<GSL_complex>* pM_quadA  = new CMatrix<GSL_complex>(Arows, 1);
        CMatrix<GSL_complex>* pM_quadB  = new CMatrix<GSL_complex>(Brows, 1);

        matrix_remove(pM_quadA, Atrow, 0, Arows, 1);
        matrix_remove(pM_quadB, Btrow, 0, Brows, 1);

        pM_shifted->matrix_replace(0,       0,      *pM_quadA);
        pM_shifted->matrix_replace(Prow,    0,      *pM_quadB);

        delete pM_quadA;
        delete pM_quadB;

    }

    if(!is_vector()) {
        // shift for matrices
        //  We need to determine the top left coords and lengths of
        //  each quadrant. This is also dependant on whether or not
        //  the dimension length is even or odd.
        int     Atrow, Atcol, Arows, Acols;
        int     Btrow, Btcol, Brows, Bcols;
        int     Ctrow, Ctcol, Crows, Ccols;
        int     Dtrow, Dtcol, Drows, Dcols;

        int     Prow, Pcol;     // Insertion pivot

        int     rows    = rows_get();
        int     cols    = cols_get();

        Btrow   = 0;
        Btcol   = 0;
        if(isOdd(cols_get())) {
            Atcol   = (cols+a_pivotOffset)/2;
            Pcol    = Atcol-a_pivotOffset;
        } else {
            Atcol   = cols/2;
            Pcol    = Atcol;
        }
        Acols   =  cols-Atcol;
        Bcols   = Atcol;
        Ctcol   = Atcol;
        Ccols   = Acols;
        Dtcol   = 0;
        Dcols   = Bcols;

        if(isOdd(rows_get())) {
            Atrow   = (rows+a_pivotOffset)/2;
            Prow    = Atrow -a_pivotOffset;
        } else {
            Atrow   = rows/2;
            Prow    = Atrow;
        }
        Arows   = rows-Atrow;
        Brows   = Atrow;
        Ctrow   = 0;
        Crows   = Brows;
        Dtrow   = Atrow;
        Drows   = Arows;

        CMatrix<GSL_complex>* pM_quadA  = new CMatrix<GSL_complex>(Arows, Acols);
        CMatrix<GSL_complex>* pM_quadB  = new CMatrix<GSL_complex>(Brows, Bcols);
        CMatrix<GSL_complex>* pM_quadC  = new CMatrix<GSL_complex>(Crows, Ccols);
        CMatrix<GSL_complex>* pM_quadD  = new CMatrix<GSL_complex>(Drows, Dcols);

        matrix_remove(pM_quadA, Atrow, Atcol, Arows, Acols);
        matrix_remove(pM_quadB, Btrow, Btcol, Brows, Bcols);
        matrix_remove(pM_quadC, Ctrow, Ctcol, Crows, Ccols);
        matrix_remove(pM_quadD, Dtrow, Dtcol, Drows, Dcols);

        pM_shifted->matrix_replace(0,       0,      *pM_quadA);
        pM_shifted->matrix_replace(Prow,    Pcol,   *pM_quadB);
        pM_shifted->matrix_replace(0,       Pcol,   *pM_quadD);
        pM_shifted->matrix_replace(Prow,    0,      *pM_quadC);

        delete pM_quadA;
        delete pM_quadB;
        delete pM_quadC;
        delete pM_quadD;

    }
    CMatrix M_shifted(*pM_shifted);
    delete pM_shifted;
    return(M_shifted);
}

CMatrix<GSL_complex>
CMatrix<GSL_complex>::fftshift(
        bool    ab_inPlace
) {
    //
    // ARGS
    //  ab_inPlace              in              if true, branch to the
    //                                                  inPlace() method
    //                                          else, branch to the
    //                                                  newMem() method
    //
    // DESCRIPTION
    //  shift*() dispatching layer. Depending on the value of ab_inPlace,
    //  either an inPlace or new memory operation is performed.
    //
    // HISTORY
    //  22 January 2003
    //  o Initial design and coding.
    //

    if(ab_inPlace) {
        shiftInPlace(+1);
        return *this;
    } else {
        CMatrix<GSL_complex> M_dispatch(*this);
        M_dispatch  = shiftNewMem(+1);
        return(M_dispatch);
    }
}

CMatrix<GSL_complex>
CMatrix<GSL_complex>::ifftshift(
        bool    ab_inPlace
) {
    //
    // ARGS
    //  ab_inPlace              in              if true, branch to the
    //                                                  inPlace() method
    //                                          else, branch to the
    //                                                  newMem() method
    //
    // DESCRIPTION
    //  shift*() dispatching layer. Depending on the value of ab_inPlace,
    //  either an inPlace or new memory operation is performed.
    //
    // HISTORY
    //  22 January 2003
    //  o Initial design and coding.
    //

    if(ab_inPlace) {
        shiftInPlace(-1);
        return *this;
    } else {
        CMatrix<GSL_complex> M_dispatch(*this);
        M_dispatch  = shiftNewMem(-1);
        return(M_dispatch);
    }
}

CMatrix<GSL_complex>
CMatrix<GSL_complex>::fftshiftInPlace()
{
    shiftInPlace(+1);
    return *this;
}

CMatrix<GSL_complex>
CMatrix<GSL_complex>::ifftshiftInPlace()
{
    shiftInPlace(-1);
    return *this;
}

CMatrix<GSL_complex>
CMatrix<GSL_complex>::shiftInPlace(
    int     a_pivotOffset
) {
    //
    // ARGS
    //  a_pivotOffset       in          offset to pivot point
    //                                      about which matrix
    //                                      elements are swapped.
    //
    // DESC
    //  This method performs the actual innards of both the
    //  fftshift and iffshift operation. Conceptually, both
    //  methods are identical, differing only in the pivot
    //  point offset.
    //
    // PRECONDITIONS
    //  o None.
    //
    //
    // POSTCONDITIONS
    //  o This is a destructive (in place) operation. The original matrix
    //    has its internal contents shifted, and is not preserved.
    //
    // HISTORY
    // 20 January 2004
    //  o Initial design and coding.
    //    Usage of the non-destructive shift quickly becomes prohibitive with
    //    when using large volumes. As a matter of fact, the apparent simplicity
    //    of a non-destructive operation has not only a memory penalty, but a
    //    time penalty as well: the time required to create a new block of
    //    memory can be quite large (especially if the operating system starts
    //    to swap).
    //

    //  For odd lengthed dimensions, the "swap" is asymmetrical. An element is
    //  not simply swapped with its "shifted" target, but this target
    //  displaces a different element.
    //
    //  This displaced element in turn "jumps" forward, displacing
    //  another (different) element, which displaces another element
    //  etc. This chain continues until the jumps
    //  terminate again at the original start point.
    //
    CMatrix<int>                M_indexDelta(1,2);      // the "jump" shift

    long int            rows            = rows_get();
    long int            cols            = cols_get();

    // The following keep track of the indices in the matrix index space that
    //  are being shifted. M_startPos begins at (0, 0) and moves through the
    //  indices as governed by M_currentDelta. The shiftCounter counts how
    //  many cells have been shifted. Once this shiftCounter equals the number
    //  of matrix elements, the shift operation terminates.
    CMatrix<int>        M_indexStart(1, 2);             // start position of
                                                        //      current shift;
    CMatrix<int>        M_indexCurrent(1, 2);           // current shift index
    CMatrix<int>        M_indexNext(1, 2);              // next shift index
    CMatrix<int>        M_indexStartDelta(1, 2);        // delta for start index
    long int            shiftCounter    = 0;            // running counter of
                                                        //      shifted cells
    long int            totalElements   = 0;
    long int            juggleSeqLength = 0;
    GSL_complex         cellCurrentValue        = (GSL_complex) 0;
    GSL_complex         cellNextValue           = (GSL_complex) 0;
    GSL_complex         cellJuggleValue         = (GSL_complex) 0;

    // The indexStartDelta is used to update the "head" of a new shift vector in
    //  the matrix index space. The matrix elements are updated *orthogonally*
    //  to the indexDelta direction
    M_indexStartDelta(0)        = isEven(rows);
    M_indexStartDelta(1)        = isEven(cols);
    // If matrix dimensions are odd/odd,
    if(!M_indexStartDelta(0) && !M_indexStartDelta(1)) {
        M_indexStartDelta(0) = -1;      // "up" one row
        M_indexStartDelta(1) =  1;      // "over" one column
    }

    if(isOdd(rows))     M_indexDelta(0)   =  (rows-a_pivotOffset)/2;
        else            M_indexDelta(0)   =  rows/2;
    if(isOdd(cols))     M_indexDelta(1)   =  (cols-a_pivotOffset)/2;
        else            M_indexDelta(1)   =  cols/2;

    // start at (0, 0) in the index space
    M_indexStart(0)         = 0;            M_indexStart(1)   = 0;
    M_indexCurrent(0)       = 0;            M_indexCurrent(1) = 0;

    totalElements   = rows*cols;
    while(++shiftCounter <= totalElements) {
        M_indexNext = M_indexCurrent + M_indexDelta;

        // Check for boundary violation in indices
        if(M_indexNext(0)>=rows)        M_indexNext(0)-=rows;   // rows
        if(M_indexNext(0)<0)            M_indexNext(0)+=rows;
        if(M_indexNext(1)>=cols)        M_indexNext(1)-=cols;   // cols
        if(M_indexNext(1)<0)            M_indexNext(1)+=cols;

        if(!juggleSeqLength++)
            cellJuggleValue    = val(M_indexCurrent(0), M_indexCurrent(1));

        cellNextValue       = val(M_indexNext(0), M_indexNext(1));
        val(M_indexNext(0), M_indexNext(1))     = cellJuggleValue;

        M_indexCurrent      = M_indexNext;
        cellJuggleValue     = cellNextValue;

        if(M_indexNext.equal(M_indexStart)) {
            M_indexStart += M_indexStartDelta;
            // Check for boundary violation in indices
            if(M_indexStart(0)>=rows)   M_indexStart(0)-=rows;  // rows
            if(M_indexStart(0)<0)       M_indexStart(0)+=rows;
            if(M_indexStart(1)>=cols)   M_indexStart(1)-=cols;  // cols
            if(M_indexStart(1)<0)       M_indexStart(1)+=cols;

            if(M_indexStartDelta(0)==1 && M_indexStartDelta(1)==1) {
                // In this case the matrix dimensions are even/even
                //  and the shift is completely symmetrical: i.e.
                //  a straight swap of elements. In this case, we simply
                //  advance the M_indexStart linearly across the matrix
                int row     = 0;
                int col     = 0;
                                                // Since there at least two
                row         = (shiftCounter/2)/cols;    // jumps required to
                col         = (shiftCounter/2)%cols;    // return to start point
                                                // we divide shiftCounter by 2
                M_indexStart(0)     = row;
                M_indexStart(1)     = col;
            }
            M_indexCurrent = M_indexStart;
            juggleSeqLength = 0;
        }

    }
    return *this;
}

CMatrix<GSL_complex>
CMatrix<GSL_complex>::operator = (
        const CMatrix<GSL_complex>& rval) {
    //
    // ARGS
    //  rval            in              right-hand argument to operator
    //
    // DESC
    //  Assignment operator.
    //
    // PRECONDITIONS
    //  o A "pointer" equivalence is performed, i.e. *no* deepcopy. The
    //    pointer to *this matrix's data is simply decremented (and
    //    possibly destroyed) and reset to point to rval's data (which
    //    has its refcount increased.
    // o *this must have the same size as rval
    //
    // HISTORY
    // 17 April 2003
    //         o Template specialization.
    //
    // 27 January 2004
    //  o Added refcnt_inc()
    //


    char*        pch_name = "operator=";
    char*        pch_errormsg = "Matrix assignment error. Incompatible dimensions";

    if (!compatible (rval))
        _error (pch_name, pch_errormsg);

    // clean up current value:
    //if (--ps_mat->refcnt == 0)
    if (refcnt_dec() == 0)
        coreMatrix_destruct();

    // assign to new value:
    //rval.ps_mat->refcnt++;              // tell the rval it has another reference
    rval.refcnt_inc();                  // tell the rval it has another reference
    ps_mat = rval.ps_mat;               // point at the rval matrix structure
    return *this;
}

CMatrix<GSL_complex>
CMatrix<GSL_complex>::operator - ()
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements unary negation, i.e. each matrix element is
    //  is multiplied by -1.
    //
    // POSTCONDITIONS
    //  o A new matrix containing the negative is created and returned.
    //
    // HISTORY
    // 18 April 2003
    //         o Template specialization.
    //

    CMatrix<GSL_complex>        M_negative (ps_mat->rows, ps_mat->cols);
    for (int i = 0; i < ps_mat->rows; i++)
        for (int j = 0; j < ps_mat->cols; j++)
            M_negative._mval (i, j) = -_mval (i, j);
    return (M_negative);
}

CMatrix<GSL_complex>
CMatrix<GSL_complex>::operator *(
        const CMatrix<GSL_complex>&                 rval) {
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements matrix multiplication, i.e. *this is multiplied
    //  by rval.
    //
    // POSTCONDITIONS
    //  o A new matrix containing the product is created and returned.
    //
    // HISTORY
    // 17 April 2003
    //         o Template specialization.
    //

    // check for compatibility
    if (cols_get() != rval.rows_get())
        _error ("Matrix Multiplication Error!",
                 "Columns of first matrix are not equal to rows of second.");
    CMatrix<GSL_complex>        M_result (ps_mat->rows, rval.cols_get());

        gsl_blas_zgemm (CblasNoTrans, CblasNoTrans,
                  (GSL_complex) 1.0,
                  this->ps_mat->ps_gsl->matrix_complex,
                  rval.ps_mat->ps_gsl->matrix_complex,
                  (GSL_complex) 0.0,
                  M_result.ps_mat->ps_gsl->matrix_complex);

    return (M_result);
}

CMatrix<GSL_complex>
CMatrix<GSL_complex>::operator * (
        const GSL_complex&        rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements scalar multiplication, i.e. each matrix element is
    //  is increased by rval.
    //
    // POSTCONDITIONS
    //  o A new matrix containing the scalar product is created and returned.
    //
    // HISTORY
    // 17 April 2003
    //         o Template specialization.
    //
    
    CMatrix<GSL_complex>        M_result(ps_mat->rows, ps_mat->cols);

    M_result.copy(*this);
    gsl_matrix_complex_scale (M_result.ps_mat->ps_gsl->matrix_complex, rval);
    return (M_result);
}

CMatrix<GSL_complex>
CMatrix<GSL_complex>::operator *= (
        const CMatrix<GSL_complex>&                        rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements element by element multiplication i.e. each matrix element
    //  multiplied by its corresponding element in rval.
    //
    // POSTCONDITIONS
    //  o The current matrix data values are overwritten.
    //
    // HISTORY
    // 17 April 2003
    //         o Template specialization.
    //
    
    // check for compatibility
    char* pch_name = "operator*=";
    if (!compatible (rval))
        _error (pch_name, "Matrices must have identical sizes");
    
    gsl_matrix_complex_mul_elements(ps_mat->ps_gsl->matrix_complex,
                                    (const gsl_matrix_complex*)rval.ps_mat->ps_gsl->matrix_complex); 

    return *this;
}

CMatrix<GSL_complex>
CMatrix<GSL_complex>::operator *= (
        const GSL_complex&                                rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements element by constant multiplication i.e. a "scaling"
    //        of matrix values by rval.
    //
    // POSTCONDITIONS
    //  o The current matrix data values are overwritten.
    //
    // HISTORY
    // 17 April 2003
    //         o Template specialization.
    //
    
    // check for compatibility
    char* pch_name = "operator*=";
    
    gsl_matrix_complex_scale(ps_mat->ps_gsl->matrix_complex, rval); 

    return *this;
}

CMatrix<GSL_complex>
CMatrix<GSL_complex>::operator /(
        const CMatrix<GSL_complex>&                rval) {
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements matrix division, i.e. *this is multiplied
    //  by rval.
    //
    // POSTCONDITIONS
    //  o A new matrix containing the product is created and returned.
    //
    // HISTORY
    // 18 April 2003
    //         o Template specialization.
    //
    // 29 March 2004
    //	o Compatibility checks changed.
    //
    
    // check for compatibility
    if (!compatible(rval))
        _error ("operator/(const CMatrix<GSL_complex>& rval)",
                 "Matrices must have identical dimensions");
    CMatrix<GSL_complex>        M_result (ps_mat->rows, ps_mat->cols);

    M_result.copy(*this);
    gsl_matrix_complex_div_elements(M_result.ps_mat->ps_gsl->matrix_complex, 
                            (const gsl_matrix_complex*)rval.ps_mat->ps_gsl->matrix_complex); 
    return (M_result);
}

CMatrix<GSL_complex>
CMatrix<GSL_complex>::operator / (
        const GSL_complex&                        rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements scalar division, i.e. each matrix element is
    //  is increased by rval.
    //
    // POSTCONDITIONS
    //  o A new matrix containing the scalar product is created and returned.
    //
    // HISTORY
    // 18 April 2003
    //         o Template specialization.
    //
    
    CMatrix<GSL_complex>        M_result(ps_mat->rows, ps_mat->cols);
    GSL_complex                        inv;
    
    inv        = rval;
    inv.inverse(); 

    M_result.copy(*this);
    gsl_matrix_complex_scale (M_result.ps_mat->ps_gsl->matrix_complex, inv);                 
    return (M_result);
}

CMatrix<GSL_complex>
CMatrix<GSL_complex>::operator /= (
        const CMatrix<GSL_complex>&                rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements element by element division i.e. each matrix element
    //  multiplied by its corresponding element in rval.
    //
    // POSTCONDITIONS
    //  o The current matrix data values are overwritten.
    //
    // HISTORY
    // 18 April 2003
    //         o Template specialization.
    //
    
    // check for compatibility
    char* pch_name = "operator/=(const CMatrix<GSL_complex>& rval)";
    if (!compatible (rval))
        _error (pch_name, "Matrices must have identical sizes");
    
    gsl_matrix_complex_div_elements(ps_mat->ps_gsl->matrix_complex, 
                                    (const gsl_matrix_complex*)rval.ps_mat->ps_gsl->matrix_complex); 

    return *this;    
}

CMatrix<GSL_complex>
CMatrix<GSL_complex>::operator /= (
        const GSL_complex&                        rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements element by constant division i.e. a "scaling"
    //        of matrix values by 1/rval.
    //
    // POSTCONDITIONS
    //  o The current matrix data values are overwritten.
    //
    // HISTORY
    // 18 April 2003
    //         o Template specialization.
    //
    
    GSL_complex                        inv;
    
    inv        = rval;
    inv.inverse(); 
       
    gsl_matrix_complex_scale(ps_mat->ps_gsl->matrix_complex, inv); 
    
    return *this;    
}

CMatrix<GSL_complex>
CMatrix<GSL_complex>::operator - (
        const CMatrix<GSL_complex>&          rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements matrix subtraction.
    //
    // POSTCONDITIONS
    //  o A new matrix containing the sum is created and returned.
    //
    // WARNING
    //  o Be aware that the returned matrix needs to be assigned
    //    to an already created matrix otherwise memory leaks
    //    may occur.
    //
    // HISTORY
    // 18 April 2003
    //         o Template specialization.
    //
    
    // check for compatibility
    char* pch_name = "operator - (const CMatrix<GSL_complex>&rval)";
    if ((ps_mat->rows != rval.rows_get())
        || (ps_mat->cols !=
            rval.cols_get()))_error (pch_name,
                                       "Reason: Incompatible dimensions.");
    CMatrix<GSL_complex>         M_sum (ps_mat->rows, ps_mat->cols);
    
    M_sum.copy(*this);
    gsl_matrix_complex_sub (
            M_sum.ps_mat->ps_gsl->matrix_complex,
        rval.ps_mat->ps_gsl->matrix_complex
    );
    return (M_sum);
}

CMatrix<GSL_complex>
CMatrix<GSL_complex>::operator - (
        const GSL_complex&                        rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements scalar subtraction, i.e. each matrix element is
    //  is increased by rval.
    //
    // POSTCONDITIONS
    //  o A new matrix containing the sum is created and returned.
    //
    // HISTORY
    // 18 April 2003
    //         o Template specialization.
    //
    
    CMatrix<GSL_complex>        M_sum(ps_mat->rows, ps_mat->cols);
    GSL_complex                        neg;
    neg.dat[0]        = -rval.dat[0];
    neg.dat[1]        = -rval.dat[1];

    M_sum.copy(*this);
    gsl_matrix_complex_add_constant(M_sum.ps_mat->ps_gsl->matrix_complex, neg);
    return(M_sum);
}

CMatrix<GSL_complex>
CMatrix<GSL_complex>::operator -= (
        const CMatrix<GSL_complex>&                rval) {
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements matrix - matrix subtraction.
    //
    // POSTCONDITIONS
    //  o Implements A = A-B matrix type operation, i.e. contents of current
    //          matrix are overwritten.
    //
    // HISTORY
    // 18 April 2003
    //         o Template specialization.
    //
      
    gsl_matrix_complex_sub (
            ps_mat->ps_gsl->matrix_complex,
        rval.ps_mat->ps_gsl->matrix_complex
    );
    return (*this);
}

CMatrix<GSL_complex>
CMatrix<GSL_complex>::operator -= (
        const GSL_complex&                        rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements "constant" matrix subtraction.
    //
    // POSTCONDITIONS
    //  o Implements A = A-b scalar type operation, i.e. contents of current
    //          matrix are overwritten.
    //
    // HISTORY
    // 18 April 2003
    //         o Template specialization.
    //
      
    GSL_complex                        neg;
    neg.dat[0]        = -rval.dat[0];
    neg.dat[1]        = -rval.dat[1];
        
    gsl_matrix_complex_add_constant (
            ps_mat->ps_gsl->matrix_complex,
        neg
    );
    return (*this);
}

CMatrix<GSL_complex>
CMatrix<GSL_complex>::operator + (
        const CMatrix<GSL_complex>&          rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements matrix addition.
    //
    // POSTCONDITIONS
    //  o A new matrix containing the sum is created and returned.
    //
    // WARNING
    //  o Be aware that the returned matrix needs to be assigned
    //    to an already created matrix otherwise memory leaks
    //    may occur.
    //
    // HISTORY
    // 18 April 2003
    //         o Template specialization.
    //

    // check for compatibility
    char* pch_name = "operator + (const CMatrix<GSL_complex>&rval)";
    if ((ps_mat->rows != rval.rows_get())
        || (ps_mat->cols !=
            rval.cols_get()))_error (pch_name,
                                       "Reason: Incompatible dimensions.");
    CMatrix<GSL_complex>         M_sum (ps_mat->rows, ps_mat->cols);
    
    M_sum.copy(*this);
    gsl_matrix_complex_add (
            M_sum.ps_mat->ps_gsl->matrix_complex,
        rval.ps_mat->ps_gsl->matrix_complex
    );
    return (M_sum);
}

CMatrix<GSL_complex>
CMatrix<GSL_complex>::operator + (
        const GSL_complex&                        rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements scalar addition, i.e. each matrix element is
    //  is increased by rval.
    //
    // POSTCONDITIONS
    //  o A new matrix containing the sum is created and returned.
    //
    // HISTORY
    // 18 April 2003
    //         o Template specialization.
    //
    
    CMatrix<GSL_complex>        M_sum(ps_mat->rows, ps_mat->cols);
    
    M_sum.copy(*this);
    gsl_matrix_complex_add_constant(M_sum.ps_mat->ps_gsl->matrix_complex, rval);
    
    return(M_sum);
}

CMatrix<GSL_complex>
CMatrix<GSL_complex>::operator += (
        const CMatrix<GSL_complex>&          rval) {
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements matrix + matrix addition.
    //
    // POSTCONDITIONS
    //  o Implements A = A+B matrix type operation, i.e. contents of current
    //          matrix are overwritten.
    //
    // HISTORY
    // 18 April 2003
    //         o Template specialization.
    //
      
    gsl_matrix_complex_add (
            ps_mat->ps_gsl->matrix_complex,
        rval.ps_mat->ps_gsl->matrix_complex
    );
    return (*this);
}

CMatrix<GSL_complex>
CMatrix<GSL_complex>::operator += (
        const GSL_complex&                        rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements "constant" matrix addition.
    //
    // POSTCONDITIONS
    //  o Implements A = A+b scalar type operation, i.e. contents of current
    //          matrix are overwritten.
    //
    // HISTORY
    // 18 April 2003
    //         o Template specialization.
    //  
    
    gsl_matrix_complex_add_constant (
            ps_mat->ps_gsl->matrix_complex,
        rval
    );
    return (*this);
}

CMatrix<GSL_complex>**
CM_fft3Dl(
        CMatrix<GSL_complex>**  ppzM_input,
        CMatrix<GSL_complex>**& ppzM_output,
        int                     depth,
        e_FFTDIR                e_direction     /*= e_forward   */,
        bool                    b_inPlace       /*= false       */,
        bool                    b_MKLwrap       /*= false       */)
{
    // ARGS
    //  ppzM_input              in              input volume to be FFT'd
    //  ppzM_output             in/out          volume to contain results of
    //                                                  FFT.
    //                                                  If b_inPlace is true,
    //                                                  this will merely point
    //                                                  to pzM_input.
    //  depth                   in              size of the 3rd dimension
    //                                                  of volume. The other
    //                                                  two dimensions are
    //                                                  spec'd by the matrix
    //                                                  slices.
    //  e_direction             in              forward/inverse FFT
    //  b_inPlace               in              if true, change contents of
    //                                                  input vector;
    //                                          if false, conserve input and
    //                                                  create output with
    //                                                  FFT'd volume.
    //  b_MKLwrap               in              if true, wrap around core MKL
    //                                                  routines, else
    //                                                  revert to GSL routines.
    //
    // DESC
    //  Performs a 3D volume FFT by calling a 1D FFT routine across all the
    //  vectors comprising the volume along all three dimensions.
    //
    //  This method calls 1D line FFTs along each orthogonal direction. Compared
    //  with MatLab, this approach is fractionally slower (roughly a factor of 2).
    //
    // PRECONDITIONS
    //  o Volumes are double indirected since they are arrays of pointers to
    //    matrices.
    //  o ppzM_input/output are arrays of matrices (i.e. volumes).
    //  o depth specs the 3rd dimension.
    //  o NO CHECKING IS DONE ON VOLUME BOUNDS!!
    //    Make sure that the volume passed to this function is
    //    correctly specified by the depth parameter!
    //    Also, make sure that each matrix in the array is the
    //    *same* size! This is not checked for!
    //  o ppzM_output:
    //    If b_inPlace is true, and ppzM_output != ppzM_input, an error will
    //    be thrown.
    //    If b_inPlace is false, and ppzM_output is not NULL, an error will
    //    be thrown.
    //
    // POSTCONDITIONS
    //  o A pointer to the transformed volume is returned in function namespace.
    //  o If b_inPlace is true, the passed volume is transformed in place, and
    //    pzM_output = pzM_input.
    //

    char*                       pch_proc = "CMatrix<GSL_complex>** fft3D(...)";
    int                         i, j, k;
    CMatrix<GSL_complex>*       pzM_1D;

    // Some preliminary processing on the input and output volumes.
    if(b_inPlace) {
        if(ppzM_input != ppzM_output)
            CM_error(   pch_proc,
                        "For inPlace 3DFFT the input vol must point to output vol");
    } else {
        if(ppzM_output != NULL)
            CM_error(   pch_proc,
                        "For non-inPlace 3DFFT, output vol must be passed as NULL");
        // Allocate memory for output volume...
        ppzM_output     = new CMatrix<GSL_complex>* [depth];
        for(i=0; i<depth; i++) {
            ppzM_output[i]      = new CMatrix<GSL_complex>(ppzM_input[i]->rows_get(),
                                                           ppzM_input[i]->cols_get());
            // and copy the input volume to the output...
            ppzM_output[i]->copy(*ppzM_input[i]);
        }
    }

    int rows    = ppzM_input[0]->rows_get();
    int cols    = ppzM_input[0]->cols_get();
    // Process all vectors pointing in the `i' (row) dimension direction
    pzM_1D      = new CMatrix<GSL_complex>(1, rows);

    for(k=0; k<depth; k++) {
        for(j=0; j<cols; j++) {
            // Extract the vector
            for(i=0; i<rows; i++)
                pzM_1D->val(0, i) = ppzM_output[k]->val(i, j);
            // 1D FFT - do an inPlace on this dummy vector
            pzM_1D->fft1D(e_direction, b_MKLwrap);
            // Insert back into volume
            for(i=0; i<rows; i++)
                ppzM_output[k]->val(i, j)       = pzM_1D->val(0, i);
        }
    }

    // Process all vectors pointing in the `j' (col) dimension direction
    delete pzM_1D;
    pzM_1D      = new CMatrix<GSL_complex>(1, cols);
    for(k=0; k<depth; k++) {
        for(i=0; i<rows; i++) {
            // Extract the vector
            for(j=0; j<cols; j++)
                pzM_1D->val(0, j) = ppzM_output[k]->val(i, j);
            // 1D FFT - do an inPlace on this dummy vector
            pzM_1D->fft1D(e_direction, b_MKLwrap);
            // Insert back into volume
            for(j=0; j<cols; j++)
                ppzM_output[k]->val(i, j)       = pzM_1D->val(0, j);
        }
    }

    // Process all vectors pointing in the `k' (slice) dimension direction
    delete pzM_1D;
    pzM_1D      = new CMatrix<GSL_complex>(1, depth);
    for(i=0; i<rows; i++) {
        for(j=0; j<cols; j++) {
            // Extract the vector
            for(k=0; k<depth; k++)
                pzM_1D->val(0, k) = ppzM_output[k]->val(i, j);
            // 1D FFT - do an inPlace on this dummy vector
            pzM_1D->fft1D(e_direction, b_MKLwrap);
            // Insert back into volume
            for(k=0; k<depth; k++)
                ppzM_output[k]->val(i, j)       = pzM_1D->val(0, k);
        }
    }

    delete pzM_1D;
    return ppzM_output;
}

CMatrix<GSL_complex>**
CM_fft3D(
        CMatrix<GSL_complex>**  ppzM_input,
        CMatrix<GSL_complex>**& ppzM_output,
        int                     depth,
        e_FFTDIR                e_direction     /*= e_forward   */,
        bool                    b_inPlace       /*= false       */,
        bool                    b_MKLwrap       /*= false       */)
{
    // ARGS
    //  ppzM_input              in              input volume to be FFT'd
    //  ppzM_output             in/out          volume to contain results of
    //                                                  FFT.
    //                                                  If b_inPlace is true,
    //                                                  this will merely point
    //                                                  to pzM_input.
    //  depth                   in              size of the 3rd dimension
    //                                                  of volume. The other
    //                                                  two dimensions are
    //                                                  spec'd by the matrix
    //                                                  slices.
    //  e_direction             in              forward/inverse FFT
    //  b_inPlace               in              if true, change contents of
    //                                                  input vector;
    //                                          if false, conserve input and
    //                                                  create output with
    //                                                  FFT'd volume.
    //  b_MKLwrap               in              if true, wrap around core MKL
    //                                                  routines, else
    //                                                  revert to GSL routines.
    //
    // DESC
    //  Performs a 3D volume FFT by calling a 1D FFT routine across all the
    //  vectors comprising the volume along all three dimensions.
    //
    //  This method calls 2D FFTs in each slice and then 1D FFTs along each line
    //  in the slice direction. It is a factor 2 faster than the FFT3Dl method,
    //  and is thus favourably comparable with MatLab.
    //
    // PRECONDITIONS
    //  o Volumes are double indirected since they are arrays of pointers to
    //    matrices.
    //  o ppzM_input/output are arrays of matrices (i.e. volumes).
    //  o depth specs the 3rd dimension.
    //  o NO CHECKING IS DONE ON VOLUME BOUNDS!!
    //    Make sure that the volume passed to this function is
    //    correctly specified by the depth parameter!
    //    Also, make sure that each matrix in the array is the
    //    *same* size! This is not checked for!
    //  o ppzM_output:
    //    If b_inPlace is true, and ppzM_output != ppzM_input, an error will
    //    be thrown.
    //    If b_inPlace is false, and ppzM_output is not NULL, an error will
    //    be thrown.
    //
    // POSTCONDITIONS
    //  o A pointer to the transformed volume is returned in function namespace.
    //  o If b_inPlace is true, the passed volume is transformed in place, and
    //    pzM_output = pzM_input.
    //

    char*                       pch_proc = "CMatrix<GSL_complex>** fft3D(...)";
    int                         i, j, k;
    CMatrix<GSL_complex>*       pzM_1D;

    // Some preliminary processing on the input and output volumes.
    if(b_inPlace) {
        if(ppzM_input != ppzM_output)
            CM_error(   pch_proc,
                        "For inPlace 3DFFT the input vol must point to output vol");
    } else {
        if(ppzM_output != NULL)
            CM_error(   pch_proc,
                        "For non-inPlace 3DFFT, output vol must be passed as NULL");
        // Allocate memory for output volume...
        ppzM_output     = new CMatrix<GSL_complex>* [depth];
        for(i=0; i<depth; i++) {
            ppzM_output[i]      = new CMatrix<GSL_complex>(ppzM_input[i]->rows_get(),
                                                           ppzM_input[i]->cols_get());
            // and copy the input volume to the output...
            ppzM_output[i]->copy(*ppzM_input[i]);
        }
    }

    int rows    = ppzM_input[0]->rows_get();
    int cols    = ppzM_input[0]->cols_get();

    // Process all matrices orthogonal to the `k' (depth) direction
    CMatrix<GSL_complex>        zM_2D(rows, cols);
    for(k=0; k<depth; k++) {
        zM_2D.copy(*(ppzM_output[k]));
        zM_2D.fft2D(e_direction, b_MKLwrap);
        ppzM_output[k]->copy(zM_2D);
    }

    // Process all vectors pointing in the `k' (slice) dimension direction
    pzM_1D      = new CMatrix<GSL_complex>(1, depth);
    for(i=0; i<rows; i++) {
        for(j=0; j<cols; j++) {
            // Extract the vector
            for(k=0; k<depth; k++)
                pzM_1D->val(0, k) = ppzM_output[k]->val(i, j);
            // 1D FFT - do an inPlace on this dummy vector
            pzM_1D->fft1D(e_direction, b_MKLwrap);
            // Insert back into volume
            for(k=0; k<depth; k++)
                ppzM_output[k]->val(i, j)       = pzM_1D->val(0, k);
        }
    }

    delete pzM_1D;
    return ppzM_output;
}

CMatrix<GSL_complex>**
CM_vol_zeroPad(
        CMatrix<GSL_complex>**  ppzM_input,
        int                     depth,
        CMatrix<GSL_complex>**& ppzM_output,
        e_DOMINANCE             dir,
        int                     size )
{
    // ARGS
    //  ppzM_input              in              input volume to be zeroPadded
    //  depth                   in              size of the 3rd dimension
    //                                                  of volume. The other
    //                                                  two dimensions are
    //                                                  spec'd by the matrix
    //                                                  slices.
    //  ppzM_output             in/out          output volume that is zeroPadded;
    //                                              note that this will
    //                                              always be *larger* than the
    //                                              input volume.
    //  dir                     in              dimension along which to zero pad.
    //  size                    in              depth along the dimension to zero pad.
    //                                              note that this prepended and
    //                                              appended to the dimension.
    //
    // DESC
    //  Performs a 3D volume zero padding operation. Additional memory containing
    //  zeroes is added to the front and back of the dimension specified by
    //  'dir'. The length of this addition is 'size', thus the total volume
    //  increased by a factor of 2*size.
    //
    // PRECONDITIONS
    //  o Volumes are double indirected since they are arrays of pointers to
    //    matrices.
    //  o ppzM_input/output are arrays of matrices (i.e. volumes).
    //
    // POSTCONDITIONS
    //  o A pointer to the transformed volume is returned in function namespace.
    //  o Additional volume is added to the front *and* back of the specified
    //    dimension!
    //

    char*   pch_proc    = "CM_vol_zeroPad(...)";

    int     totalSlices;        //  These define the new dimensions
    int     totalRows;          //  of the created volume.
    int     totalCols;

    int     inputRows   =   ppzM_input[0]->rows_get();
    int     inputCols   =   ppzM_input[0]->cols_get();

    int     i;

    if(ppzM_output != NULL)
        CM_error(   pch_proc,
                    "For vol_zeroPad, output vol must be passed as NULL");

     switch(dir) {
        case e_row:
            totalRows   =   inputRows   + 2*size;
            totalCols   =   inputCols;
            totalSlices =   depth;
            break;
        case e_column:
            totalRows   =   inputRows;
            totalCols   =   inputCols   + 2*size;
            totalSlices =   depth;
            break;
        case e_slice:
            totalRows   =   inputRows;
            totalCols   =   inputCols;
            totalSlices =   depth       + 2*size;
            break;
     }

    // Allocate memory for output volume container
    ppzM_output     = new CMatrix<GSL_complex>* [totalSlices];
    GSL_complex     z_zero(0, 0);
    for(i=0; i<totalSlices; i++) {
        ppzM_output[i]      = new CMatrix<GSL_complex>(    totalRows,
                                                           totalCols,
                                                           z_zero);
    }

    if(dir!=e_slice) {
        for(i=0; i<totalSlices; i++)
            *ppzM_output[i]  = ppzM_input[i]->zeroPad(dir, size);
    } else {
        for(i=size; i<totalSlices-size; i++) {
            ppzM_output[i]->copy(*ppzM_input[i-size]);
        }
    }

    return ppzM_output;

}

CMatrix<GSL_complex>**
CM_coreVolume_construct(
    CMatrix<GSL_complex>**&     ppzM_volume,
    int                         rows,
    int                         cols,
    int                         slices
) {
    //
    // ARGS
    //  ppzM_volume             in          block of memory to construct
    //  rows                    in          number of rows
    //  cols                    in          number of cols
    //  slices                  in          number of slices
    //
    // DESC
    //  Construct a volume block of memory. This is in many ways
    //  a precursor to a class-based definition of CMatrix volumes.
    //
    // PRECONDITIONS
    //  o Make sure that the depth parameter is correct!
    //  o ppzM_volume must be a NULL pointer.
    //
    // POSTCONDITIONS
    //  o The volume pointed to by ppzM_volume is allocated and returned.
    //
    // HISTORY
    // 09 September 2003
    //  o Initial design and coding.
    //

    char*   pch_proc    = "CM_coreVolume_construct";
    int     i;

    if(ppzM_volume != NULL)
        CM_error(   pch_proc,
                    "For construction, input vol must be passed as NULL");

    // Allocate memory for volume container
    ppzM_volume     = new CMatrix<GSL_complex>* [slices];
    GSL_complex     z_zero(0, 0);
    for(i=0; i<slices; i++) {
        ppzM_volume[i]      = new CMatrix<GSL_complex>(    rows,
                                                           cols,
                                                           z_zero);
    }

    return ppzM_volume;

}


void
CM_coreVolume_destruct(
    CMatrix<GSL_complex>**&     ppzM_volume,
    int                         depth
) {
    //
    // ARGS
    //  ppzM_volume             in          block of memory to free
    //  depth                   in          number of slices
    //
    // DESC
    //  Destroy a volume block of memory. This is in many ways
    //  a precursor to a class-based definition of CMatrix volumes.
    //
    // PRECONDITIONS
    //  o Make sure that the depth parameter is correct!
    //
    // POSTCONDITIONS
    //  o The volume pointed to by ppzM_volume is freed.
    //
    // HISTORY
    // 09 September 2003
    //  o Initial design and coding.
    //

    for(int i=0; i<depth; i++) {
        delete ppzM_volume[i];
    }
    delete ppzM_volume;
    ppzM_volume     = NULL;
}

void
CM_vol_print(
    CMatrix<GSL_complex>**      ppzM_volume,
    int                         depth,
    char*                       pch_msg
) {
    //
    // ARGS
    //  ppzM_volume             in          block of memory to free
    //  depth                   in          number of slices
    //  pch_msg                 in          print msg
    //
    // DESC
    //  CMatrix prints a volume.
    //
    // HISTORY
    // 09 September 2003
    //  o Initial design and coding.
    //

    char        pch_txt[16384];


    for(int i=0; i<depth; i++) {
        sprintf(pch_txt, "(slice %d) %s", i, pch_msg);
        ppzM_volume[i]->print(pch_msg);
    }
}

CMatrix<GSL_complex>**
CM_vol_shift(
        CMatrix<GSL_complex>**  ppzM_input,
        CMatrix<GSL_complex>**& ppzM_output,
        int                     depth,
        e_FFTSHIFTDIR           e_dir
) {
    //
    // ARGS
    //  ppzM_input              in              input volume to be ifftshifted
    //  ppzM_output             in/out          output volume that is ifftshifted;
    //  depth                   in              size of the 3rd dimension
    //                                                  of volume. The other
    //                                                  two dimensions are
    //                                                  spec'd by the matrix
    //                                                  slices.
    //  e_dir                   in              the direction in which to fftshift
    //
    // DESC
    //  Performs a volume shift operation - either fftshift or ifftshift.
    //
    // PRECONDITIONS
    //  o Since the *shift() operation is non-destructive, a new memory block is
    //    allocated that will contain the result of the operation. ppzM_output
    //    must be passed as a NULL - on return, it will point to a valid block
    //    of memory.
    //  o depth must be correct!
    //
    // POSTCONDITIONS
    //  o ppzM_output will point to a newly created block of memory housing the
    //    result of the *shift() operation.
    //  o ppzM_output is also returned in function name.
    //
    // HISTORY
    // 11 September 2003
    //  o Initial design and coding.
    //

    char*   pch_proc    = "CM_vol_ifftshift(...)";

    int     i;
    int     rows;
    int     cols;
    int     slices;

    if(ppzM_output != NULL)
        CM_error(   pch_proc,
                    "For vol_zeroPad, output vol must be passed as NULL");

    rows    = ppzM_input[0]->rows_get();
    cols    = ppzM_input[0]->cols_get();
    slices  = depth;

    CM_coreVolume_construct(    ppzM_output,
                                rows,
                                cols,
                                slices);

    for(i=0; i<slices; i++) {
        switch(e_dir) {
            case e_ifftshift:
                *(ppzM_output[i])       = ppzM_input[i]->ifftshift();
                break;
            case e_fftshift:
                *(ppzM_output[i])       = ppzM_input[i]->fftshift();
                break;
        }
    }

    CMatrix<int>                M_list(1, slices);
    CMatrix<int>                M_listifftshifted(1, slices);
    CMatrix<GSL_complex>**      ppzM_origOrder  = new CMatrix<GSL_complex>* [slices];

    // Make a backup of the original output list to conserve the pointers
    for(i=0; i<slices; i++) {
        M_list(0, i)        = i;
        ppzM_origOrder[i]   = ppzM_output[i];
    }

    switch(e_dir) {
        case e_ifftshift:
            M_listifftshifted   = M_list.ifftshift();
            break;
        case e_fftshift:
            M_listifftshifted   = M_list.fftshift();
            break;
    }
    // reorganize the actual slices
    for(i=0; i<slices; i++) {
        ppzM_output[i]  = ppzM_origOrder[M_listifftshifted(0, i)];
    }

    delete [] ppzM_origOrder;
    return ppzM_output;

}

//////---------->
////// Specializations: GSL_complex_float
//////---------->

int CMatrix<GSL_complex_float>::createdCounter        = 0;

void
CMatrix<GSL_complex_float>::_error (
        char *apch_proc,
        char *apch_msg,
        int   code) {
    //
    // ARGS
    //  apch_proc               in              (optional) function name
    //                                                  in which error has
    //                                                  occured
    //  apch_msg                in              (optional) error message
    //                                                  text
    //  code                    in              (optional) error code
    //
    // DESC
    //  Main error reporting method for the CMatrix class.
    //
    // TODO
    //        o Wrap exception throwing around #ifdefs?
    //
    // HISTORY
    // 12 March 2003
    //        o Added id, rows, and col output
    //

    cerr << "\nFatal error encountered.\n";
    cerr << "\tCMatrix object<GSL_complex_float> id: " << ps_mat->id << endl;
    cerr << "\tRows: " << ps_mat->rows << ", cols: " << ps_mat->cols << endl;
    cerr << "\tCurrent function: " << "CMatrix::" << apch_proc << "\n";
    cerr << "\t" << apch_msg << "\n";
    cerr << "Throwing an exception to (this) with code " << code << endl;
    throw (this);
}

void
CMatrix <GSL_complex_float>::coreMatrix_construct(
        int                arows,
        int                acols
) {
    // ARGS
    //        arows                in                        number of rows in matrix
    //        acols                in                        number of cols in matrixtype
    //
    // DESCRIPTION
    //        This method is a consolidation of the main initialization code
    //        that is shared across all the *matrix* constructors.
    //
    // PRECONDITIONS
    //        o Should only be called from a constructor method.
    //
    // POSTCONDITIONS
    //        o A gsl-aware ps_mat structure is created.
    //
    // HISTORY
    // 19 March 2003
    //        o Initial design and coding.
    //
    // 19 November 2003
    //	o complex_float templatisation
    //

    // create the structure
    try {
        ps_mat = new matstruct;
    } catch (bad_alloc xa) {
        _error("base constructor",
                "could not allocate memory for matstruct (heap exhausted?)");
    }
    ps_mat->rows = arows;
    ps_mat->cols = acols;
    // allocate memory for the structure

//    ps_mat->ps_gsl = (gslwrap*)malloc(sizeof(gslwrap));
    ps_mat->ps_gsl                              = new gslwrap;
    ps_mat->ps_gsl->matrix_int                  = 0x0;
    ps_mat->ps_gsl->matrix_float                = 0x0;
    ps_mat->ps_gsl->matrix_double               = 0x0;
    ps_mat->ps_gsl->matrix_complex              = 0x0;
    ps_mat->ps_gsl->matrix_complex_float        = gsl_matrix_complex_float_alloc(arows, acols);

    if(!ps_mat->ps_gsl->matrix_complex_float)
            _error("base constructor", "gsl_matrix_alloc error");
    ps_mat->data = (GSL_complex_float**) ps_mat->ps_gsl->matrix_complex_float->data;

    // set first reference/id to this data
    ps_mat->refcnt = 1;
    ps_mat->id           = CMatrix::createdCounter++;
}

CMatrix<GSL_complex_float>::CMatrix(
        int                     mrows,
	int                     mcols,
        GSL_complex_float       initvalue) {
    //
    // ARGS
    //  mrows                   in              number of rows for matrix
    //  mcols                   in              number of cols for matrix
    //  initvalue               in              initial value for matrix
    //                                                  elements
    //
    // DESC
    //  Matrix constructor
    //
    // PRECONDITIONS
    //  o Make sure that mrows and mcols are not negative numbers!
    //
    // HISTORY
    // 16 April 2003
    //        o Template specialization.
    //
    // 19 November 2003
    //	o complex_float templatisation
    //

    coreMatrix_construct(mrows, mcols);

    for (int i = 0; i < mrows; i++)
        for (int j = 0; j < mcols; j++)
            _mval (i, j) = initvalue;
}

CMatrix<GSL_complex_float>::CMatrix(
        int                             mrows,
        int                             mcols,
        const GSL_complex_float*        p_initvalues) {
    //
    // ARGS
    //  mrows                   in              number of rows for matrix
    //  mcols                   in              number of cols for matrix
    //  p_initvalues            in              pointer to block of memory
    //                                                  containing initial
    //                                                  values
    //
    // DESC
    //  Matrix constructor
    //
    // PRECONDITIONS
    //  o Make sure that mrows and mcols are not negative numbers!
    //
    // HISTORY
    // 19 March 2003
    //        o Integration with gsl/coreMatrix_construct
    //
    // 19 November 2003
    //	o complex_float templatisation
    //

    coreMatrix_construct(mrows, mcols);

    int        c = 0;
    for (int i = 0; i < mrows; i++)
        for (int j = 0; j < mcols; j++)
            _mval (i, j) = p_initvalues[c++];
}

CMatrix<GSL_complex_float>::CMatrix (
        char            *matfile)
{
    //
    // ARGS
    //  matfile         in              name of file containing the
    //                                          matrix to be parsed
    //
    // DESC
    //  Constructor.
    //
    //  Creates a matrix from the specified filename.
    //
    // PRECONDITIONS
    //  o The information in *matfile must be in the format that
    //    this method expects!
    //
    // POSTCONDITIONS
    //  o Matrix is created.
    //
    // HISTORY
    // 18 April 2003
    //         o Template specialization.
    //
    // 19 November 2003
    //	o complex_float templatisation
    //

    const int		BUFSIZE = 120;
    FILE*               infile;
    char                buf[BUFSIZE], *cp, *cp2;
    int                 rfound = 0, cfound = 0, colonsfound = 0;
    int                 rows = 0, cols = 0;

    if ((infile = fopen (matfile, "r")) == 0)
        _error ("File open error! Cannot open file:-", matfile);
    while (fgets (buf, BUFSIZE, infile)) {        // parse file initialization header
        // for each header line remove comments
        if ((cp = strpbrk (buf, "#")) != NULL)        // if comment found,
            *cp = '\0';                // terminate string at comment
        if ((cp = strpbrk (buf, "r")) != NULL)
            if ((cp2 = strpbrk (cp, "=")) != NULL)
                if ((cp = strpbrk (cp2, "0123456789")) != NULL) {
                    rows = atoi (cp);
                    rfound++;
                }
        if ((cp = strpbrk (buf, "c")) != NULL)
            if ((cp2 = strpbrk (cp, "=")) != NULL)
                if ((cp = strpbrk (cp2, "0123456789")) != NULL) {
                    cols = atoi (cp);
                    cfound++;
                }
        if (strstr (buf, ":::::::") != NULL) {
            colonsfound++;
            break;                // out of while loop
        }
    }
    if (!rfound || !cfound || !colonsfound) {
        fprintf (stderr, "%s %s", matfile, nonstandard);
        exit (1);
    }

    coreMatrix_construct(rows, cols);

    for (int row = 0; row < ps_mat->rows; row++)
        for (int col = 0; col < ps_mat->cols; col++) {
            char                pch_real[20], pch_imag[20];
            float               f_real, f_imag;
            fscanf (infile, "%s %s", pch_real, pch_imag);
            f_real        = atof(pch_real);
            f_imag        = atof(pch_imag);
            GSL_complex_float vc_z(f_real, f_imag);
            _mval (row, col) = vc_z;
            if (ferror (infile))
                _error ("Input file error! Problem accessing file:-",
                         matfile);
        }
    fclose (infile);
}

CMatrix<GSL_complex_float>::CMatrix(const CMatrix<GSL_complex_float> & z)
{
    //
    // DESC
    //  Copy constructor - simple pointer based.
    //
    // POSTCONDITIONS
    //  o NB!! NB!! NB!!
    //    The copy constructor merely increases the reference count of the
    //    "source" data, and directs the target to point to the source.
    //    This is *not* a deepcopy!. Although allowing for fast copies between
    //    matrices, it can potentially suffer from problems relating to scope
    //    local variable variable problems. If used in function arguments
    //    or as returns out of functions, then it is not really a problem.
    //
    // HISTORY
    // 27 January 2004
    //  o Added refcnt_inc()
    //

    //z.ps_mat->refcnt++;                 // adding another reference
    z.refcnt_inc();                     // adding another reference
    ps_mat = z.ps_mat;                  // point to the new matstruct
}

void
CMatrix<GSL_complex_float>::coreMatrix_destruct() {
    //
    // DESC
    //        Explicitly destroys the core data stuctures of a gsl-aware *matrix*
    //        object.
    //
    // HISTORY
    // 20 March 2003
    //        o Initial design and coding.
    //
    // 19 November 2003
    //	o complex_float templatisation
    //
    // 28 January 2004
    //  o Corrected non-GSL destruction using delete [] properly.
    //

#ifdef GSL_USE
        gsl_matrix_complex_float_free(ps_mat->ps_gsl->matrix_complex_float);
        delete ps_mat->ps_gsl;
#else
        for (int y = 0; y < ps_mat->rows; y++)
            delete [] ps_mat->data[y];
        delete [] ps_mat->data;
#endif
        delete ps_mat;
}

CMatrix<GSL_complex_float>::~CMatrix () {
    //
    // DESC
    //  Destructor.
    //
    // POSTCONDITIONS
    //  Decrements reference count. If counter reaches zero, frees allocated
    //  memory.
    //
    // HISTORY
    // 27 January 2004
    //  o Added refcnt_dec()
    //

    //if (--ps_mat->refcnt == 0)
    if (refcnt_dec() == 0)
             coreMatrix_destruct();
}

GSL_complex_float &
CMatrix<GSL_complex_float>::val(
        int         row,
        int         col)
{
    //
    // ARGS
    //        row                in                row to access
    //        col                in                col to access
    //
    // DESC
    //         Element selection: can be used to read or write.
    //
    // PRECONDITIONS
    //        o Range checking can be expensive, particularly in tightly
    //          nested iterative loops. By defining CMATRIX_NORANGE range
    //          checking is disabled. This will result in tremendous performance
    //          boost, but at the risk of a single out of bound access
    //          killing the entire process.
    //
    // HISTORY
    //  3 Decemeber 2001
    //        o Added compiler directives to turn off range checking.
    //

#ifndef CMATRIX_NORANGE
    char *
        pch_name = "val";
    char *
        pch_errormsg = "Addressing error. Index out of range";
    if (!(row >= 0 && row < ps_mat->rows && col >= 0 && col < ps_mat->cols))
        _error (pch_name, pch_errormsg);
    else
#endif
    return (_mval (row, col));
}

GSL_complex_float &
CMatrix<GSL_complex_float>::operator()(
        int         row,
        int         col)
{
    //
    // ARGS
    //        row                in                row to access
    //        col                in                col to access
    //
    // DESC
    //         Element selection: can be used to read or write.
    //
    // PRECONDITIONS
    //        o Range checking can be expensive, particularly in tightly
    //          nested iterative loops. By defining CMATRIX_NORANGE range
    //          checking is disabled. This will result in tremendous performance
    //          boost, but at the risk of a single out of bound access
    //          killing the entire process.
    //
    // HISTORY
    //  3 Decemeber 2001
    //        o Added compiler directives to turn off range checking.
    //

#ifndef CMATRIX_NORANGE
    char *
        pch_name = "val";
    char *
        pch_errormsg = "Addressing error. Index out of range";
    if (!(row >= 0 && row < ps_mat->rows && col >= 0 && col < ps_mat->cols))
        _error (pch_name, pch_errormsg);
    else
#endif
    return (_mval (row, col));
}

GSL_complex_float &
CMatrix<GSL_complex_float>::val(
        int         element)
{
    //
    // ARGS
    //  element            in                element to access
    //
    // DESC
    //  Element selection: can be used to read or write. This method
    //	is for use with *vectors*.
    //
    // PRECONDITIONS
    //	o Only use with *vectors*, not matrices. 
    //  o Range checking can be expensive, particularly in tightly
    //    nested iterative loops. By defining CMATRIX_NORANGE range
    //    checking is disabled. This will result in tremendous performance
    //    boost, but at the risk of a single out of bound access
    //    killing the entire process.
    //
    // HISTORY
    //  3 Decemeber 2001
    //  o Added compiler directives to turn off range checking.
    //
    // 24 September 2003
    //	o Designed for *vectors*.
    //

#ifndef CMATRIX_NORANGE
    char* pch_name = "val";
    char* pch_errormsg = "Addressing error. Index out of range";
    if ((!is_vectorRow() && !is_vectorColumn()))
        _error (pch_name, pch_errormsg);
#endif
    if(is_vectorRow())
	return (val (0, element));
    if(is_vectorColumn())
	return (val(element, 0));
}

GSL_complex_float &
CMatrix<GSL_complex_float>::operator()(
        int         element)
{
    //
    // ARGS
    //  element            in                element to access
    //
    // DESC
    //  Element selection: can be used to read or write. This method
    //	is for use with *vectors*.
    //
    // PRECONDITIONS
    //	o Only use with *vectors*, not matrices. 
    //  o Range checking can be expensive, particularly in tightly
    //    nested iterative loops. By defining CMATRIX_NORANGE range
    //    checking is disabled. This will result in tremendous performance
    //    boost, but at the risk of a single out of bound access
    //    killing the entire process.
    //
    // HISTORY
    //  3 Decemeber 2001
    //  o Added compiler directives to turn off range checking.
    //
    // 24 September 2003
    //	o Designed for *vectors*.
    //

#ifndef CMATRIX_NORANGE
    char* pch_name = "val";
    char* pch_errormsg = "Addressing error. Index out of range";
    if ((!is_vectorRow() && !is_vectorColumn()))
        _error (pch_name, pch_errormsg);
#endif
    if(is_vectorRow())
	return (val (0, element));
    if(is_vectorColumn())
	return (val(element, 0));
}


GSL_complex_float &        
CMatrix<GSL_complex_float>::_mval(  
        int        row,
        int        col)  const {
    //
    // ARGS
    //        row                in                row to access
    //        col                in                col to access
    //
    // DESC
    //         Element selection: can be used to read or write.
    //
    // PRECONDITIONS
    //         This routine does *not* use range checking! It is an
    //        interal access routine that is used by class methods
    //        and is not available for "public" use.
    //
    // HISTORY
    // 18 March 2003
    //        o Moved from header file into main source body as part 
    //          of gsl integration.
    //
    // 19 November 2003
    //	o complex_float templatisation
    //
    
    return( (GSL_complex_float&)*gsl_matrix_complex_float_ptr(ps_mat->ps_gsl->matrix_complex_float, row, col));
}

void
CMatrix<GSL_complex_float>::print (
    char*			apch_msg		/*="ans"	*/,
    int                         a_precision             /*= 6           */,
    int                         a_width                 /*= 12          */,
    ios::fmtflags               a_userFormat            /*= 0           */,
    int                         a_marginOffset          /*= 0           */,
    char                        ach_left                /*= (char) 0    */,
    int                         a_tabOffset             /*= 0           */,
    bool                        ab_fancy                /*= 1           */) {
    //
    // ARGS
    //	apch_msg		in			message that prepends matrix dump
    //  a_precision             in                      the precision of numerical output
    //  a_width                 in                      the width of the numerical field
    //  a_userFormat            in                      if non-zero, set stream format to arg
    //  a_marginOffset          in                      the amount of tab spaces prepending ch_left
    //  ach_left                in                      also somewhat primitive. If non-zero, will
    //                                                      print (char) left_char as left most character
    //                                                      of each new line. Useful for when the matrix is
    //                                                      part of a larger output dump
    //  a_tabOffset             in                      the amount of tab spaces prepending dump
    //  ab_fancy                in                      still rather primitive. If false, will print
    //                                                      *only* the numbers, nothing else
    //
    // DESC
    //        Print a matrix in a variety of ways. The a_userFormat parameter allows the
    //        user to send extra formatting flags to the underlying stream. Consult a
    //        C++ manual for flags.
    //
    // NOTE
    //         The `marginOffset' and `tabOffset' are used in conjunction
    //         with ch_left as follows:
    //
    //                 [<--- marginOffset --->]<ch_left>[<--- tabOffset --->]
    //
    // HISTORY
    //  18 February 2002
    //        o Changed all c-style printf's to c++-style cout
    //
    // 19 November 2003
    //	o complex_float templatisation
    //

    ios::fmtflags	streamFlags	= cout.flags();		// get current stream flags
    char                pch_tab[64]     = "";
    char                pch_margin[64]  = "";
    char                pch_indent[64]  = "";
    int                 i;
    
    for (i = 0; i < a_marginOffset; i++)
        strcat (pch_margin, "\t");

    for (i = 0; i < a_tabOffset; i++)
        strcat (pch_tab, "\t");

    sprintf (pch_indent, "%s%c%s", pch_margin, ach_left, pch_tab);

    if (ab_fancy) {
        if (ps_mat->rows > 1)
            cout << pch_tab << apch_msg << " = [" << endl << pch_indent;
        else
            cout << pch_tab << apch_msg << " = [ ";
    }
    for (int row = 0; row < ps_mat->rows; row++) {
        for (int col = 0; col < ps_mat->cols; col++) {
//            cout.width(a_width);
            cout.precision(a_precision);
            cout.setf(ios::internal);
            cout.setf(ios::skipws);
            if(a_userFormat)
                    cout.flags(a_userFormat);
            cout << _mval(row, col) << "\t";
            //if(ab_fancy) cout << "| ";
        }
        if (ps_mat->rows > 1)
            cout << endl << pch_indent;
    }
    if (ab_fancy)
        if (ps_mat->rows > 1)
            cout << "]" << endl << pch_indent;
        else
            cout << "] ";
    cout.flags(streamFlags);                // Restore stream flags
}

string
CMatrix<GSL_complex_float>::sprint (        
    char*			apch_msg		/*="ans"	*/,
    int                         a_precision             /*= 6           */,
    int                         a_width                 /*= 12          */,
    ios::fmtflags               a_userFormat            /*= 0           */,
    int                         a_marginOffset          /*= 0           */,
    char                        ach_left                /*= (char) 0    */,
    int                         a_tabOffset             /*= 0           */,
    bool                        ab_fancy                /*= 1           */) {
    //
    // ARGS
    //	apch_msg		in			message that prepends matrix dump
    //  a_precision             in                      the precision of numerical output
    //  a_width                 in                      the width of the numerical field
    //  a_userFormat            in                      if non-zero, set stream format to arg
    //  a_marginOffset          in                      the amount of tab spaces prepending ch_left
    //  ach_left                in                      also somewhat primitive. If non-zero, will
    //                                                      print (char) left_char as left most character
    //                                                      of each new line. Useful for when the matrix is
    //                                                      part of a larger output dump
    //  a_tabOffset             in                      the amount of tab spaces prepending dump
    //  ab_fancy                in                      still rather primitive. If false, will print
    //                                                      *only* the numbers, nothing else
    //                                                *only* the numbers, nothing else
    // DESC
    //        Print a matrix in a variety of ways. The a_userFormat parameter allows the
    //        user to send extra formatting flags to the underlying stream. Consult a
    //        C++ manual for flags.
    //
    // NOTE
    //  The `marginOffset' and `tabOffset' are used in conjunction
    //  with ch_left as follows:
    //
    //  [<--- marginOffset --->]<ch_left>[<--- tabOffset --->]
    //
    // HISTORY
    // 18 February 2002
    //        o Changed all c-style printf's to c++-style cout
    //
    // 01 April 2003
    //        o Templatization.
    //
    // 19 November 2003
    //	o complex_float templatisation
    //

    stringstream        sstream("");

    ios::fmtflags       streamFlags        = cout.flags();        // get current stream flags
    char                pch_tab[64]        = "";
    char                pch_margin[64]        = "";
    char                pch_indent[64]        = "";
    int                        i;

    for (i = 0; i < a_marginOffset; i++)
        strcat (pch_margin, "\t");

    for (i = 0; i < a_tabOffset; i++)
        strcat (pch_tab, "\t");

    sprintf (pch_indent, "%s%c%s", pch_margin, ach_left, pch_tab);

    if (ab_fancy) {
        if (ps_mat->rows > 1)
            sstream << pch_tab << apch_msg << " = [" << endl << pch_indent;
        else
            sstream << pch_tab << apch_msg << " = [ ";
    }
    for (int row = 0; row < ps_mat->rows; row++) {
        for (int col = 0; col < ps_mat->cols; col++) {
//            sstream.width(a_width);
            sstream.precision(a_precision);
            sstream.setf(ios::internal);
            sstream.setf(ios::skipws);
            if(a_userFormat)
                    sstream.flags(a_userFormat);
            sstream << GSL_REAL(_mval(row, col)) << " " << GSL_IMAG(_mval(row, col)) << "\t";
        }
        if (ps_mat->rows > 1)
            sstream << endl << pch_indent;
    }
    if (ab_fancy)
        if (ps_mat->rows > 1)
            sstream << "]" << endl << pch_indent;
        else
            sstream << "] ";
    sstream.flags(streamFlags);                // Restore stream flags
    return(sstream.str());
}

void
CMatrix<GSL_complex_float>::fprint(
    string                  astr_filename,
    char*                   apch_msg            /*= ""          */,
    int                     a_precision         /*= 6           */,
    int                     a_width             /*= 12          */,
    ios::fmtflags           a_userFormat        /*= 0           */
) {
    //
    // ARGS
    //  astr_filename       in                  filename to write matrix to
    //  apch_msg            in                  message that prepends matrix dump
    //  a_precision         in                  the precision of numerical output
    //  a_width             in                  the width of the numerical field
    //  a_userFormat        in                  if non-zero, set stream format to arg
    //
    // DESC
    //        Print a matrix to file in a variety of ways. The a_userFormat parameter allows
    //        the user to send extra formatting flags to the underlying stream. Consult a
    //        C++ manual for flags.
    //
    // HISTORY
    // 18 February 2002
    //        o Changed all c-style printf's to c++-style cout
    //
    // 01 April 2003
    //        o Templatization.
    //
    // 19 November 2003
    //	o complex_float templatisation
    //

    char*           pch_proc        = "fprint";

    ofstream        sout(astr_filename.c_str());
    if(!sout)
            _error(pch_proc, "Could not create output file");

    time_t         time_now        = time(NULL);
    string        str_time        = ctime(&time_now);
    // strip trailing \n
    str_time.erase(str_time.length()-1);
    string        str_hostname(getenv("HOSTNAME"));
    string        str_user(getenv("USER"));

    sout << "#"                                                     << endl;
    sout << "# Standard CMatrix<GSL_complex_float> save file."      << endl;
    sout << "#        Created by\t"         << str_user             << endl;
    sout << "#        Date stamp:\t"        << str_time             << endl;
    sout << "#        Machine name:\t"      << str_hostname         << endl;
    if(strlen(apch_msg))
    sout << "#        Matrix name:\t"       << apch_msg             << endl;
    sout << "#"                                                     << endl;
    sout << "rows="         << rows_get();
    sout << " columns="     << cols_get()                           << endl;
    sout << ":::::::"                                               << endl;

    string        str_data        = sprint("",
                                            a_precision,
                                            a_width,
                                            a_userFormat,
                                            0, 0, 0,
                                            0);

    sout << str_data;
    sout.close();
}

CMatrix<GSL_complex_float>
CMatrix<GSL_complex_float>::copy (
        const CMatrix<GSL_complex_float>&          source)
{
    // ARGS
    //  source                  in              matrix whose contents to copy
    //                                                  into *this
    //
    // DESC
    //  Explicit deepcopy.
    //
    // POSTCONDITIONS.
    //  o The data in source is copied over into *this. Space for this data
    //    is created as in necessary.
    //
    // HISTORY
    // 20 March 2003
    //        o coreMatrix_[construct/destruct] integration.
    //
    // 26 March 2003
    //         o Templatization.
    //

    int                row, col, y;

    if (!compatible (source)) {
            coreMatrix_destruct();
        coreMatrix_construct(source.rows_get(), source.cols_get());
    }
    for (row = 0; row < ps_mat->rows; row++)
        for (col = 0; col < ps_mat->cols; col++)
            _mval (row, col) = source._mval (row, col);
    return *this;
}

CMatrix<GSL_complex_float>
CMatrix<GSL_complex_float>::matrix_remove(
        CMatrix<GSL_complex_float>*&  apM,
        int                     a_trow,
        int                     a_tcol,
        int                     a_rows,
        int                     a_cols) {
    //
    // ARGS
    //  apM             in/out          structure to contain removed
    //                                          sub-matrix
    //  a_trow, a_tcol  in              top left coordinate in base matrix
    //                                          of submatrix
    //  a_rows, a_cols  in              relative to trow and tcol, number
    //                                          of rows and cols to extract
    //                                          into submatrix
    //
    // DESC
    //  Removes a submatrix from a base matrix
    //
    // PRECONDITIONS
    //  o 0 < a_trow < rows_get()
    //  o 0 < a_tcol < cols_get()
    //  o a_trow + a_rows <= rows_get() + 1
    //  o a_tcol + a_cols <= cols_get() + 1
    //
    // POSTCONDITIONS
    //  o Returns submatrix (in name) and pointer to submatrix (in argument list).
    //
    // HISTORY
    // 03 December 2000
    //  o Changed calling parameters to explicitly allow for
    //    holding memory structure.
    //
    // 31 March 2003
    //  o Templatization.
    //
    // 19 November 2003
    //	o complex_float templatisation
    //

    char*       pch_name = "remove_matrix";
    int         row, col;

    if (a_trow + a_rows > ps_mat->rows + 1 ||
        a_tcol + a_cols > ps_mat->cols + 1)
            _error (pch_name, "Specified target dimensions invalid!");

    if (!apM->compatible (a_rows, a_cols)) {
        delete        apM;
        apM         = new CMatrix<GSL_complex_float>(a_rows, a_cols, (GSL_complex_float) 0.0);
    }

    for (row = 0; row < a_rows; row++)
        for (col = 0; col < a_cols; col++) {
            apM->val (row, col) = val (row + a_trow, col + a_tcol);
        }
    return *apM;
}

CMatrix<GSL_complex_float>
CMatrix<GSL_complex_float>::matrix_replace(
        int                     trow,
        int                     tcol,
        CMatrix<GSL_complex_float>&   aM_replacement)
{
    //
    // ARGS
    //  trow, tcol      in      top left coordinate in base matrix
    //                                  where replacement will be added.
    //  aM_replacement  in      replacement matrix
    //
    // DESC
    //  Overwrites part of a matrix with another matrix.
    //
    // 01 April 2003
    //  o Templatization.
    //
    // 19 November 2003
    //  o GSL_complex_float
    //

    char*       pch_name        = "matrix_replace";
    char*       pch_errorSize   = "Replacement cannot fit into base.";
    char*       pch_errorCoord  = "Specified insert point is out of range.";
    int         row, col;

    if ((trow < 0) || (tcol < 0) ||
        (trow > ps_mat->rows) || (tcol > ps_mat->cols))
            _error (pch_name, pch_errorCoord);
    if ((trow + aM_replacement.rows_get () > ps_mat->rows) ||
        (tcol + aM_replacement.cols_get () > ps_mat->cols))
        _error (pch_name, pch_errorSize);
    for (row = 0; row < aM_replacement.rows_get (); row++) {
        for (col = 0; col < aM_replacement.cols_get (); col++) {
            val (row + trow, col + tcol) = aM_replacement.val (row, col);
        }
    }
    return *this;
}

CMatrix<GSL_complex_float>
CMatrix<GSL_complex_float>::replicate(
        int             times,
        e_DOMINANCE     dir) {
    //
    // ARGS
    //  times           in      number of times the matrix is replicated
    //  dir             in      direction of replication
    //
    // DESC
    //  Replicates a matrix
    //
    //  This routine replaces the base matrix with a much larger
    //  copy of itself, made up of (times) instances of the original
    //  base. The replications occurs in either a column- or row-dominant
    //  direction.
    //
    //  In other words, if we have an [N x M] matrix that we wish
    //  to replicate p times in a row-dominant fashion, this routine
    //  replaces the base [N x M] matrix with a new matrix:
    //
    //          [ [ N x M] [N x M] ... [N x M] ]
    //
    //  (The column dominant is simply the transpose of the above)
    //
    // PRECONDITIONS
    //        o *this is destroyed and rebuilt!
    //
    // HISTORY
    // 01 April 2003
    //  o Templatization.
    //


    int                         y, i;
    int                         mrows, mcols;
    CMatrix<GSL_complex_float>  M(ps_mat->rows, ps_mat->cols, (GSL_complex_float) 0.0);

    // Make a deepcopy of myself
    M.copy (*this);

    // Track the current matrix size...
    mrows = ps_mat->rows;
    mcols = ps_mat->cols;

    // Then destroy and create a new matrix in its core structures
    coreMatrix_destruct();
    switch (dir) {
    case e_row:
        coreMatrix_construct(mrows, mcols*times);
        break;
    case e_column:
        coreMatrix_construct(mrows*times, mcols);
        break;
    }

    for (i = 0; i < times; i++) {
        switch (dir) {
        case e_row:
            matrix_replace (0, i * mcols, M);
            break;
        case e_column:
            matrix_replace (i * mrows, 0, M);
            break;
        }
    }
    return *this;
}

CMatrix<GSL_complex_float>
CMatrix<GSL_complex_float>::zeroPad(
    e_DOMINANCE         dir,
    int                 size) {
    //
    // ARGS
    //  dir             in          the dimension along which to zero pad
    //  size            in          size of padding along a dimension
    //
    // DESC
    //  This method "zeropads" a matrix with 2 submatrices of padding
    //  length 'size' along the direction 'dir'.
    //
    //  Zero padding is necessary before calling fft-type methods on
    //  non-power of two lengthed dimensions.
    //
    // PRECONDITIONS
    //  o None.
    //
    // POSTCONDITIONS
    //  o Non destructive
    //
    // HISTORY
    // 08 September 2003
    //  o Initial design and coding.
    //

    GSL_complex_float           z_zero(0, 0);
    int                         rows            = rows_get();
    int                         cols            = cols_get();
    int                         i, j;
    int                         trow, tcol;

    CMatrix<GSL_complex_float>  M(rows, cols, (GSL_complex_float) 0.0);
    CMatrix<GSL_complex_float>* pM_zeroPadded;

    // Make a deepcopy of myself
    M.copy (*this);

    switch (dir) {
    case e_column:
        pM_zeroPadded   = new CMatrix<GSL_complex_float>(rows, cols+2*size, z_zero);
        trow    = 0;
        tcol    = size;
        break;
    case e_row:
        pM_zeroPadded   = new CMatrix<GSL_complex_float>(rows+2*size, cols, z_zero);
        trow    = size;
        tcol    = 0;
        break;
    }

    // and insert the original data in the center
    pM_zeroPadded->matrix_replace(trow, tcol, M);
    CMatrix<GSL_complex_float>    M_padded(*pM_zeroPadded);
    delete  pM_zeroPadded;
    return M_padded;
}

bool
CMatrix<GSL_complex_float>::is_vector() const
{
    //
    // Checks if object is a vector
    //
    return (is_vectorRow() || is_vectorColumn());
}

bool
CMatrix<GSL_complex_float>::is_vectorColumn() const
{
    //
    // Checks if object is column vector
    //
    return ((ps_mat->cols == 1) && (ps_mat->rows > 1));
}

bool
CMatrix<GSL_complex_float>::is_vectorRow() const
{
    //
    // Checks if object is row vector
    //
    return ((ps_mat->rows == 1) && (ps_mat->cols > 1));
}

bool
CMatrix<GSL_complex_float>::is_square() const
{
    //
    // Checks if matrix is square
    //
    return (ps_mat->cols == ps_mat->rows);
}

int
CMatrix<GSL_complex_float>::size1D ()	const
{
    //
    // DESC
    //  Returns the 1D size of a matrix, i.e. the product of the rows
    //  and columns
    //
    // HISTORY
    //  14 November 2000
    //  o Initial design and coding.
    //

    return (rows_get () * cols_get ());
}

CMatrix<int>
CMatrix<GSL_complex_float>::size()	const
{
    //
    // DESC
    //  Returns a pointer to a matrix containing the size of *this.
    //
    // HISTORY
    //  23 September 2000
    //  o Initial design and coding
    //
    //  14 November 2000
    //  o Possible dangling / lost pointer!
    //
    // 31 March 2003
    //         o Templatization.
    //

    CMatrix<int>        iM_size(1, 2);

    // rows
    iM_size.val (0, 0) = rows_get ();
    //cols
    iM_size.val (0, 1) = cols_get ();

    return iM_size;
}

e_MATRIXTYPE
CMatrix<GSL_complex_float>::matrix_type () {
    //
    // Returns the matrix type
    // of its argument
    //

    if (is_vectorColumn())
        return e_columnVector;
    if (is_vectorRow())
        return e_rowVector;
    if (is_square())
        return e_square;
    return e_matrix;
}

bool
CMatrix<GSL_complex_float>::compatible (
        int                     dim1,
        int                     dim2) {
    //
    // ARGS
    //  dim1, dim2              in              dimensions under consideration
    //
    // DESC
    //  Checks if base is compatible with passed arguments
    //
    // 26 March 2003
    //         o Templatization.
    //

    if (ps_mat->rows != dim1 || ps_mat->cols != dim2)
        return false;
    return true;
}

bool
CMatrix<GSL_complex_float>::compatible (
        const CMatrix<GSL_complex_float>&          M)
{
    //
    // ARGS
    //
    //  M                       in      Matrix under consideration
    //
    // DESC
    //  Checks if base is compatible with passed arguments
    //
    // 26 March 2003
    //         o Templatization.
    //

    int                dim1, dim2;

    dim1         = M.rows_get();
    dim2         = M.cols_get();
    if (ps_mat->rows != dim1 || ps_mat->cols != dim2)
        return false;
    return true;
}

//////---------->
////// Miscellaneous Maths
//////---------->

CMatrix<GSL_complex_float>
CMatrix<GSL_complex_float>::inverse ()
{
    //
    // DESC
    //  Finds the inverse of *this matrix.
    //
    // POSTCONDITIONS
    //  o A new matrix containing the result is created and returned.
    //
    // HISTORY
    // 19 April 2003
    //         o Template specialization.
    //
    // 19 November 2003
    //	o complex_float templatisation
    //
    
    char*	pch_proc	= "CMatrix<GSL_complex_float>::inverse ()";
    _error(pch_proc, "Inverse on complex *float* matrices not defined! Use complex *double* instead.", 1);
}

CMatrix<GSL_complex_float>
CMatrix<GSL_complex_float>::fft1D(
        e_FFTDIR                e_direction     /*= e_forward   */,
        bool                    b_MKLwrap       /*= false       */) {
    //
    // ARGS
    //  e_direction     in              type of FFT operation to perform:
    //                                          - e_forward FFT
    //                                          - e_backward FFT
    //                                          - e_inverse FFT
    //  b_MKLwrap       in              if true, wrap around low level MKL
    //                                          routines directly as opposed
    //                                          to wrapping around the default
    //                                          GSL function calls.
    //
    // DESC
    //  Wrap around GSL (or MKL)-level complex fast Fourier transforms.
    //
    // PRECONDITIONS
    // o Note that if wrapping around MKL routines, the length of the
    //   vector *must* be a power of two.
    //
    // POSTCONDITIONS
    // o Internal contents of vector are *not* preserved!
    // o *this is returned.
    //
    // HISTORY
    // 22 April 2003
    //  o Initial design and coding. Currently, e_direction is ignored
    //    and only a forward FFT is calculated.
    //
    // 03 September 2003
    //  o After long and trying testing, it was decided that the most
    //    effective mechanism for plugging a slow memory leak was to
    //    remove the "inPlace/exPlace" ability. FFT will thus always be
    //    "inPlace".
    //

    char*       pch_proc = "fft1D(e_FFTDIR e_direction)";
    if(!is_vectorRow())
        _error(pch_proc, "Input must be a row vector.", 1);

    int                         i;
    const int                   n       = cols_get();

    if(!b_MKLwrap) {
        gsl_fft_complex_wavetable_float*  pgsl_wavetable;
        gsl_fft_complex_workspace_float*  pgsl_workspace;

        // First, create packed array from the complex data set.
        float*     pv_dataPacked = new float[2*n];
        for (i = 0; i < n; i++) {
            pv_dataPacked[2*i]      = GSL_REAL(_mval(0, i));
            pv_dataPacked[2*i+1]    = GSL_IMAG(_mval(0, i));
        }

        // Create the wave tables and workspace
        pgsl_wavetable = gsl_fft_complex_wavetable_float_alloc (n);
        pgsl_workspace = gsl_fft_complex_workspace_float_alloc (n);

        // Perform the required operation
        gsl_fft_complex_float_forward(  pv_dataPacked, 1, n,
                                        pgsl_wavetable,
                                        pgsl_workspace);

        // Pack the results in the return CMatrix<GSL_complex_float>
        for (i = 0; i < n; i++) {
            GSL_complex_float     z_data(0, 0);
            GSL_SET_COMPLEX(&z_data, pv_dataPacked[2*i], pv_dataPacked[2*i+1]);
            _mval(0, i)             = z_data;
        }

        gsl_fft_complex_wavetable_float_free(     pgsl_wavetable);
        gsl_fft_complex_workspace_float_free(     pgsl_workspace);
        delete                              pv_dataPacked;

    } else {
        // TODO: This is currently using double's, to use floats need to
        // upgrade to newer version of FFTW that allows both float and
        // double to be used in the same program
        int                  direction       = FFTW_FORWARD;

        switch(e_direction) {
            case e_forward:
                direction = FFTW_FORWARD;
                break;
            case e_backward:
                direction = FFTW_BACKWARD;
                break;
            case e_inverse:
                direction = FFTW_BACKWARD;
                break;
         }

         fftw_complex *in, *out;
         fftw_plan p;
         
         in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
         out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
         
         // Fill the input vectors
         for (i = 0; i < n; i++) {
            in[i][0]          = GSL_REAL(_mval(0, i));
            in[i][1]          = GSL_IMAG(_mval(0, i));
         }
         
         p = fftw_plan_dft_1d(n, in, out, direction, FFTW_ESTIMATE);
     
         // Pack the results in the return CMatrix<GSL_complex>
         for (i = 0; i < n; i++) {
             GSL_complex_float     z_data(0, 0);
             GSL_SET_COMPLEX(&z_data, (float)out[i][0], (float)out[i][1]);
             _mval(0, i)             = z_data;
         }

         fftw_destroy_plan(p);
         fftw_free(in); 
         fftw_free(out);

    }
    return                  *this;
}

CMatrix<GSL_complex_float>
CMatrix<GSL_complex_float>::fft2D(
        e_FFTDIR                e_direction     /*= e_forward   */,
        bool                    b_MKLwrap       /*= false       */) {
    //
    // ARGS
    //  e_direction     in              type of FFT operation to perform:
    //                                          - e_forward FFT
    //                                          - e_backward FFT
    //                                          - e_inverse FFT
    //  b_MKLwrap       in              if true, wrap around low level MKL
    //                                          routines directly as opposed
    //                                          to wrapping around the default
    //                                          GSL function calls.
    //
    // DESC
    //  Wrap around GSL (or MKL)-level complex fast Fourier transforms.
    //
    // PRECONDITIONS
    // o Note that if wrapping around MKL routines, the length of the
    //   rows and columns *must* be a power of two.
    // o It appears as though the GSL routines do not offer native 2D support.
    //
    // POSTCONDITIONS
    // o Internal matrix contents are overwritten!
    // o *this is returned.
    //
    // HISTORY
    // 26 August 2003
    //  o Initial design and coding. Currently, e_direction is ignored
    //    and only a forward FFT is calculated.
    //  o Only MKL support is coded.
    //
    // 03 September 2003
    //  o After long and trying testing, it was decided that the most
    //    effective mechanism for plugging a slow memory leak was to
    //    remove the "inPlace/exPlace" ability. FFT will thus always be
    //    "inPlace".
    //
    // 19 November 2003
    //	o complex_float templatisation
    //

    char*       pch_proc = "fft2D(e_FFTDIR e_direction)";

    int                         i, j;
    const int                   n       = cols_get();
    const int                   m       = rows_get();

    if(!b_MKLwrap) {
        _error(pch_proc, "2D fft is not yet defined for GSL operations");
    } else {
        // TODO: This is currently using double's, to use floats need to
        // upgrade to newer version of FFTW that allows both float and
        // double to be used in the same program
        int                  direction       = FFTW_FORWARD;

        switch(e_direction) {
            case e_forward:
                direction = FFTW_FORWARD;
                break;
            case e_backward:
                direction = FFTW_BACKWARD;
                break;
            case e_inverse:
                direction = FFTW_BACKWARD;
                break;
         }

         fftw_complex *in, *out;
         fftw_plan p;

         in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n * m);
         out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n * m);
         
         for (i = 0; i<m; i++) {
             for(j = 0; j<n; j++) {
                 in[n*i+j][0]          = GSL_REAL(_mval(i, j));
                 in[n*i+j][1]          = GSL_IMAG(_mval(i, j));
             }
         }
   
         p = fftw_plan_dft_2d(n, m, in, out, direction, FFTW_ESTIMATE);

         // Pack the results in the return CMatrix<GSL_complex>
         for (i=0; i<m; i++) {
             for(j=0; j<n; j++) {
                 GSL_complex_float    z_data(0, 0);
                 GSL_SET_COMPLEX(&z_data, (float)out[n*i+j][0], (float)out[n*i+j][1]);
                 _mval(i, j)             = z_data;
             }
         }
             
         fftw_destroy_plan(p);
         fftw_free(in); 
         fftw_free(out);             
    }
    return                  *this;
}

CMatrix<GSL_complex_float>
CMatrix<GSL_complex_float>::ifftshiftNewMem() {
    //
    // DESC
    //  This method implements an ifftshift operation, analogous to the MatLAB
    //  function of the same name.
    //
    //  Note that although "fft" appears in the name, no actual fft is performed;
    //  this method merely "reorganises" the internal data contents, shifting the
    //  center of k-space.
    //
    //  In the case of a matrix, the quadrants are labelled as follows:
    //
    //      [B] [C]
    //      [D] [A]
    //
    //  and are "shifted" to
    //
    //      [A] [D]
    //      [C] [B]
    //
    //  where [A] is swapped with [B], and [C] is swapped with [D].
    //  For vectors, this is [B] [A] and [[B] [A]]' (transpose) with
    //  [A] and [B] swapped.
    //
    // POSTCONDITIONS
    //  o This is a non-destructive operation. An internal matrix is
    //    created locally and returned (using the copy constructor).
    //    Make certain to capture the returned matrix!
    //
    // HISTORY
    // 06 September 2003
    //  o Initial design and coding.
    //

    CMatrix<GSL_complex_float> M_dispatch(*this);

    M_dispatch  = shiftNewMem(-1);

    return(M_dispatch);
}

CMatrix<GSL_complex_float>
CMatrix<GSL_complex_float>::fftshiftNewMem() {
    //
    // DESC
    //  This method implements an fftshift operation, analogous to the MatLAB
    //  function of the same name.
    //
    //  Note that although "fft" appears in the name, no actual fft is performed;
    //  this method merely "reorganises" the internal data contents, shifting the
    //  center of k-space.
    //
    //  In the case of a matrix, the quadrants are labelled as follows:
    //
    //      [B] [C]
    //      [D] [A]
    //
    //  and are "shifted" to
    //
    //      [A] [D]
    //      [C] [B]
    //
    //  where [A] is swapped with [B], and [C] is swapped with [D].
    //  For vectors, this is [B] [A] and [[B] [A]]' (transpose) with
    //  [A] and [B] swapped.
    //
    // POSTCONDITIONS
    //  o This is a non-destructive operation. An internal matrix is
    //    created locally and returned (using the copy constructor).
    //    Make certain to capture the returned matrix!
    //
    // HISTORY
    // 06 September 2003
    //  o Initial design and coding.
    //

    CMatrix<GSL_complex_float> M_dispatch(*this);

    M_dispatch  = shiftNewMem(+1);

    return(M_dispatch);
}

CMatrix<GSL_complex_float>
CMatrix<GSL_complex_float>::shiftNewMem(
    int     a_pivotOffset
) {
    //
    // ARGS
    //  a_pivotOffset       in          offset to pivot point
    //                                      about which matrix
    //                                      elements are swapped.
    //
    // DESC
    //  This method performs the actual innards of both the
    //  fftshift and iffshift operation. Conceptually, both
    //  methods are identical, differing only in the pivot
    //  point offset.
    //
    // PRECONDITIONS
    //  o None.
    //
    //
    // POSTCONDITIONS
    //  o This is a non-destructive operation. An internal matrix is
    //    created locally and returned (using the copy constructor).
    //    Make certain to capture the returned matrix!
    //
    // HISTORY
    // 08 September 2003
    //  o Initial design and coding.
    //

    CMatrix<GSL_complex_float>*   pM_shifted = new CMatrix<GSL_complex_float>(1, 1);
    pM_shifted->copy(*this);

    if(is_vectorRow()) {
        // shift for row vectors
        int     Atcol, Acols;
        int     Btcol, Bcols;
        int     Pcol;     // Insertion pivot
        int     cols            = cols_get();

        if(isOdd(cols)) {
            Atcol   = (cols+a_pivotOffset)/2;
            Pcol    = Atcol-a_pivotOffset;
        } else {
            Atcol   = cols/2;
            Pcol    = Atcol;
        }
        Acols   = cols-Atcol;
        Bcols   = Atcol;
        Btcol   = 0;

        CMatrix<GSL_complex_float>* pM_quadA  = new CMatrix<GSL_complex_float>(1, Acols);
        CMatrix<GSL_complex_float>* pM_quadB  = new CMatrix<GSL_complex_float>(1, Bcols);

        matrix_remove(pM_quadA, 0, Atcol, 1, Acols);
        matrix_remove(pM_quadB, 0, Btcol, 1, Bcols);

        pM_shifted->matrix_replace(0,   0,      *pM_quadA);
        pM_shifted->matrix_replace(0,   Pcol,   *pM_quadB);

        delete pM_quadA;
        delete pM_quadB;

    }

    if(is_vectorColumn()) {
        // shift for col vectors
        int     Atrow, Arows;
        int     Btrow, Brows;
        int     Prow;     // Insertion pivot
        int     rows            = rows_get();

        if(isOdd(rows)) {
            Atrow   = (rows+a_pivotOffset)/2;
            Prow    = Atrow-a_pivotOffset;
        } else {
            Atrow   = rows/2;
            Prow    = Atrow;
        }
        Arows   = rows-Atrow;
        Brows   = Atrow;
        Btrow   = 0;

        CMatrix<GSL_complex_float>* pM_quadA  = new CMatrix<GSL_complex_float>(Arows, 1);
        CMatrix<GSL_complex_float>* pM_quadB  = new CMatrix<GSL_complex_float>(Brows, 1);

        matrix_remove(pM_quadA, Atrow, 0, Arows, 1);
        matrix_remove(pM_quadB, Btrow, 0, Brows, 1);

        pM_shifted->matrix_replace(0,       0,      *pM_quadA);
        pM_shifted->matrix_replace(Prow,    0,      *pM_quadB);

        delete pM_quadA;
        delete pM_quadB;

    }

    if(!is_vector()) {
        // shift for matrices
        //  We need to determine the top left coords and lengths of
        //  each quadrant. This is also dependant on whether or not
        //  the dimension length is even or odd.
        int     Atrow, Atcol, Arows, Acols;
        int     Btrow, Btcol, Brows, Bcols;
        int     Ctrow, Ctcol, Crows, Ccols;
        int     Dtrow, Dtcol, Drows, Dcols;

        int     Prow, Pcol;     // Insertion pivot

        int     rows    = rows_get();
        int     cols    = cols_get();

        Btrow   = 0;
        Btcol   = 0;
        if(isOdd(cols_get())) {
            Atcol   = (cols+a_pivotOffset)/2;
            Pcol    = Atcol-a_pivotOffset;
        } else {
            Atcol   = cols/2;
            Pcol    = Atcol;
        }
        Acols   =  cols-Atcol;
        Bcols   = Atcol;
        Ctcol   = Atcol;
        Ccols   = Acols;
        Dtcol   = 0;
        Dcols   = Bcols;

        if(isOdd(rows_get())) {
            Atrow   = (rows+a_pivotOffset)/2;
            Prow    = Atrow -a_pivotOffset;
        } else {
            Atrow   = rows/2;
            Prow    = Atrow;
        }
        Arows   = rows-Atrow;
        Brows   = Atrow;
        Ctrow   = 0;
        Crows   = Brows;
        Dtrow   = Atrow;
        Drows   = Arows;

        CMatrix<GSL_complex_float>* pM_quadA  = new CMatrix<GSL_complex_float>(Arows, Acols);
        CMatrix<GSL_complex_float>* pM_quadB  = new CMatrix<GSL_complex_float>(Brows, Bcols);
        CMatrix<GSL_complex_float>* pM_quadC  = new CMatrix<GSL_complex_float>(Crows, Ccols);
        CMatrix<GSL_complex_float>* pM_quadD  = new CMatrix<GSL_complex_float>(Drows, Dcols);

        matrix_remove(pM_quadA, Atrow, Atcol, Arows, Acols);
        matrix_remove(pM_quadB, Btrow, Btcol, Brows, Bcols);
        matrix_remove(pM_quadC, Ctrow, Ctcol, Crows, Ccols);
        matrix_remove(pM_quadD, Dtrow, Dtcol, Drows, Dcols);

        pM_shifted->matrix_replace(0,       0,      *pM_quadA);
        pM_shifted->matrix_replace(Prow,    Pcol,   *pM_quadB);
        pM_shifted->matrix_replace(0,       Pcol,   *pM_quadD);
        pM_shifted->matrix_replace(Prow,    0,      *pM_quadC);

        delete pM_quadA;
        delete pM_quadB;
        delete pM_quadC;
        delete pM_quadD;

    }
    CMatrix<GSL_complex_float> M_shifted(*pM_shifted);
    delete pM_shifted;
    return(M_shifted);
}

CMatrix<GSL_complex_float>
CMatrix<GSL_complex_float>::fftshift(
        bool    ab_inPlace
) {
    //
    // ARGS
    //  ab_inPlace              in              if true, branch to the
    //                                                  inPlace() method
    //                                          else, branch to the
    //                                                  newMem() method
    //
    // DESCRIPTION
    //  shift*() dispatching layer. Depending on the value of ab_inPlace,
    //  either an inPlace or new memory operation is performed.
    //
    // HISTORY
    //  22 January 2003
    //  o Initial design and coding.
    //

    if(ab_inPlace) {
        shiftInPlace(+1);
        return *this;
    } else {
        CMatrix<GSL_complex_float> M_dispatch(*this);
        M_dispatch  = shiftNewMem(+1);
        return(M_dispatch);
    }
}

CMatrix<GSL_complex_float>
CMatrix<GSL_complex_float>::ifftshift(
        bool    ab_inPlace
) {
    //
    // ARGS
    //  ab_inPlace              in              if true, branch to the
    //                                                  inPlace() method
    //                                          else, branch to the
    //                                                  newMem() method
    //
    // DESCRIPTION
    //  shift*() dispatching layer. Depending on the value of ab_inPlace,
    //  either an inPlace or new memory operation is performed.
    //
    // HISTORY
    //  22 January 2003
    //  o Initial design and coding.
    //

    if(ab_inPlace) {
        shiftInPlace(-1);
        return *this;
    } else {
        CMatrix<GSL_complex_float> M_dispatch(*this);
        M_dispatch  = shiftNewMem(-1);
        return(M_dispatch);
    }
}

CMatrix<GSL_complex_float>
CMatrix<GSL_complex_float>::fftshiftInPlace()
{
    shiftInPlace(+1);
    return *this;
}

CMatrix<GSL_complex_float>
CMatrix<GSL_complex_float>::ifftshiftInPlace()
{
    shiftInPlace(-1);
    return *this;
}

CMatrix<GSL_complex_float>
CMatrix<GSL_complex_float>::shiftInPlace(
    int     a_pivotOffset
) {
    //
    // ARGS
    //  a_pivotOffset       in          offset to pivot point
    //                                      about which matrix
    //                                      elements are swapped.
    //
    // DESC
    //  This method performs the actual innards of both the
    //  fftshift and iffshift operation. Conceptually, both
    //  methods are identical, differing only in the pivot
    //  point offset.
    //
    // PRECONDITIONS
    //  o None.
    //
    //
    // POSTCONDITIONS
    //  o This is a destructive (in place) operation. The original matrix
    //    has its internal contents shifted, and is not preserved.
    //
    // HISTORY
    // 20 January 2004
    //  o Initial design and coding.
    //    Usage of the non-destructive shift quickly becomes prohibitive with
    //    when using large volumes. As a matter of fact, the apparent simplicity
    //    of a non-destructive operation has not only a memory penalty, but a
    //    time penalty as well: the time required to create a new block of
    //    memory can be quite large (especially if the operating system starts
    //    to swap).
    //

    //  For odd lengthed dimensions, the "swap" is asymmetrical. An element is
    //  not simply swapped with its "shifted" target, but this target
    //  displaces a different element.
    //
    //  This displaced element in turn "jumps" forward, displacing
    //  another (different) element, which displaces another element
    //  etc. This chain continues until the jumps
    //  terminate again at the original start point.
    //

    CMatrix<int>                M_indexDelta(1,2);      // the "jump" shift

    long int            rows            = rows_get();
    long int            cols            = cols_get();

    // The following keep track of the indices in the matrix index space that
    //  are being shifted. M_startPos begins at (0, 0) and moves through the
    //  indices as governed by M_currentDelta. The shiftCounter counts how
    //  many cells have been shifted. Once this shiftCounter equals the number
    //  of matrix elements, the shift operation terminates.
    CMatrix<int>        M_indexStart(1, 2);             // start position of
                                                        //      current shift;
    CMatrix<int>        M_indexCurrent(1, 2);           // current shift index
    CMatrix<int>        M_indexNext(1, 2);              // next shift index
    CMatrix<int>        M_indexStartDelta(1, 2);        // delta for start index
    long int            shiftCounter    = 0;            // running counter of
                                                        //      shifted cells
    long int            totalElements   = 0;
    long int            juggleSeqLength = 0;
    GSL_complex_float   cellCurrentValue        = (GSL_complex_float) 0;
    GSL_complex_float   cellNextValue           = (GSL_complex_float) 0;
    GSL_complex_float   cellJuggleValue         = (GSL_complex_float) 0;

    // The indexStartDelta is used to update the "head" of a new shift vector in
    //  the matrix index space. The matrix elements are updated *orthogonally*
    //  to the indexDelta direction
    M_indexStartDelta(0)        = isEven(rows);
    M_indexStartDelta(1)        = isEven(cols);
    // If matrix dimensions are odd/odd,
    if(!M_indexStartDelta(0) && !M_indexStartDelta(1)) {
        M_indexStartDelta(0) = -1;      // "up" one row
        M_indexStartDelta(1) =  1;      // "over" one column
    }

    if(isOdd(rows))     M_indexDelta(0)   =  (rows-a_pivotOffset)/2;
        else            M_indexDelta(0)   =  rows/2;
    if(isOdd(cols))     M_indexDelta(1)   =  (cols-a_pivotOffset)/2;
        else            M_indexDelta(1)   =  cols/2;

    // start at (0, 0) in the index space
    M_indexStart(0)         = 0;            M_indexStart(1)   = 0;
    M_indexCurrent(0)       = 0;            M_indexCurrent(1) = 0;

    totalElements   = rows*cols;
    while(++shiftCounter <= totalElements) {
//        cout << endl << "shiftCounter " << shiftCounter << "\t";
//        cout << "totalElements " << totalElements << endl;
//        cout << "rows " << rows << "\tcols " << cols << endl;
//        M_indexCurrent.print("indexCurrent"); cout << endl;
//        M_indexDelta.print("indexDelta"); cout << endl;

        M_indexNext = M_indexCurrent + M_indexDelta;
//        M_indexNext.print("indexNext"); cout << endl;

        // Check for boundary violation in indices
        if(M_indexNext(0)>=rows)        M_indexNext(0)-=rows;   // rows
        if(M_indexNext(0)<0)            M_indexNext(0)+=rows;
        if(M_indexNext(1)>=cols)        M_indexNext(1)-=cols;   // cols
        if(M_indexNext(1)<0)            M_indexNext(1)+=cols;

        if(!juggleSeqLength++)
            cellJuggleValue    = val(M_indexCurrent(0), M_indexCurrent(1));

//        M_indexNext.print("indexNext-corrected"); cout << endl;
        cellNextValue       = val(M_indexNext(0), M_indexNext(1));
        val(M_indexNext(0), M_indexNext(1))     = cellJuggleValue;

        M_indexCurrent      = M_indexNext;
        cellJuggleValue     = cellNextValue;

        if(M_indexNext.equal(M_indexStart)) {
            M_indexStart += M_indexStartDelta;
            // Check for boundary violation in indices
            if(M_indexStart(0)>=rows)   M_indexStart(0)-=rows;  // rows
            if(M_indexStart(0)<0)       M_indexStart(0)+=rows;
            if(M_indexStart(1)>=cols)   M_indexStart(1)-=cols;  // cols
            if(M_indexStart(1)<0)       M_indexStart(1)+=cols;

            if(M_indexStartDelta(0)==1 && M_indexStartDelta(1)==1) {
                // In this case the matrix dimensions are even/even
                //  and the shift is completely symmetrical: i.e.
                //  a straight swap of elements. In this case, we simply
                //  advance the M_indexStart linearly across the matrix
                int row     = 0;
                int col     = 0;
                                                // Since there at least two
                row         = (shiftCounter/2)/cols;    // jumps required to
                col         = (shiftCounter/2)%cols;    // return to start point
                                                // we divide shiftCounter by 2
                M_indexStart(0)     = row;
                M_indexStart(1)     = col;
            }
            M_indexCurrent = M_indexStart;
            juggleSeqLength = 0;
        }

    }
    return *this;
}

CMatrix<GSL_complex_float>
CMatrix<GSL_complex_float>::operator = (
        const CMatrix<GSL_complex_float>& rval) {
    //
    // ARGS
    //  rval            in              right-hand argument to operator
    //
    // DESC
    //  Assignment operator.
    //
    // PRECONDITIONS
    //  o A "pointer" equivalence is performed, i.e. *no* deepcopy. The
    //    pointer to *this matrix's data is simply decremented (and
    //    possibly destroyed) and reset to point to rval's data (which
    //    has its refcount increased.
    // o *this must have the same size as rval
    //
    // HISTORY
    // 17 April 2003
    //         o Template specialization.
    //
    // HISTORY
    // 27 January 2004
    //  o Added refcnt_inc()
    //

    char*        pch_name = "operator=";
    char*        pch_errormsg = "Matrix assignment error. Incompatible dimensions";

    if (!compatible (rval))
        _error (pch_name, pch_errormsg);

    // clean up current value:
    //if (--ps_mat->refcnt == 0)
    if (refcnt_dec() == 0)
        coreMatrix_destruct();

    // assign to new value:
    //rval.ps_mat->refcnt++;              // tell the rval it has another reference
    rval.refcnt_inc();                  // tell the rval it has another reference
    ps_mat = rval.ps_mat;               // point at the rval matrix structure
    return *this;
}

CMatrix<GSL_complex_float>
CMatrix<GSL_complex_float>::operator - ()
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements unary negation, i.e. each matrix element is
    //  is multiplied by -1.
    //
    // POSTCONDITIONS
    //  o A new matrix containing the negative is created and returned.
    //
    // HISTORY
    // 18 April 2003
    //         o Template specialization.
    //

    CMatrix<GSL_complex_float>        M_negative (ps_mat->rows, ps_mat->cols);
    for (int i = 0; i < ps_mat->rows; i++)
        for (int j = 0; j < ps_mat->cols; j++)
            M_negative._mval (i, j) = -_mval (i, j);
    return (M_negative);
}

CMatrix<GSL_complex_float>
CMatrix<GSL_complex_float>::operator *(
        const CMatrix<GSL_complex_float>&                 rval) {
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements matrix multiplication, i.e. *this is multiplied
    //  by rval.
    //
    // POSTCONDITIONS
    //  o A new matrix containing the product is created and returned.
    //
    // HISTORY
    // 17 April 2003
    //         o Template specialization.
    //
    // 19 November 2003
    //	o complex_float templatisation
    //

    // check for compatibility
    if (cols_get() != rval.rows_get())
        _error ("Matrix Multiplication Error!",
                 "Columns of first matrix are not equal to rows of second.");
    CMatrix<GSL_complex_float>        M_result (ps_mat->rows, rval.cols_get());

        gsl_blas_cgemm (CblasNoTrans, CblasNoTrans,
                  (GSL_complex_float) 1.0,
                  this->ps_mat->ps_gsl->matrix_complex_float,
                  rval.ps_mat->ps_gsl->matrix_complex_float,
                  (GSL_complex_float) 0.0,
                  M_result.ps_mat->ps_gsl->matrix_complex_float);

    return (M_result);
}

CMatrix<GSL_complex_float>
CMatrix<GSL_complex_float>::operator * (
        const GSL_complex_float&        rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements scalar multiplication, i.e. each matrix element is
    //  is increased by rval.
    //
    // POSTCONDITIONS
    //  o A new matrix containing the scalar product is created and returned.
    //
    // HISTORY
    // 17 April 2003
    //         o Template specialization.
    //
    // 19 November 2003
    //	o complex_float templatisation
    //
    
    CMatrix<GSL_complex_float>        M_result(ps_mat->rows, ps_mat->cols);

    M_result.copy(*this);
    gsl_matrix_complex_float_scale (M_result.ps_mat->ps_gsl->matrix_complex_float, rval);
    return (M_result);
}

CMatrix<GSL_complex_float>
CMatrix<GSL_complex_float>::operator *= (
        const CMatrix<GSL_complex_float>&                        rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements element by element multiplication i.e. each matrix element
    //  multiplied by its corresponding element in rval.
    //
    // POSTCONDITIONS
    //  o The current matrix data values are overwritten.
    //
    // HISTORY
    // 17 April 2003
    //         o Template specialization.
    //
    // 19 November 2003
    //	o complex_float templatisation
    //
    
    // check for compatibility
    char* pch_name = "operator*=";
    if (!compatible (rval))
        _error (pch_name, "Matrices must have identical sizes");
    
    gsl_matrix_complex_float_mul_elements(ps_mat->ps_gsl->matrix_complex_float,
                                    (const gsl_matrix_complex_float*)rval.ps_mat->ps_gsl->matrix_complex_float); 

    return *this;
}

CMatrix<GSL_complex_float>
CMatrix<GSL_complex_float>::operator *= (
        const GSL_complex_float&                                rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements element by constant multiplication i.e. a "scaling"
    //        of matrix values by rval.
    //
    // POSTCONDITIONS
    //  o The current matrix data values are overwritten.
    //
    // HISTORY
    // 17 April 2003
    //         o Template specialization.
    //
    // 19 November 2003
    //	o complex_float templatisation
    //
    
    // check for compatibility
    char* pch_name = "operator*=";
    
    gsl_matrix_complex_float_scale(ps_mat->ps_gsl->matrix_complex_float, rval); 

    return *this;
}

CMatrix<GSL_complex_float>
CMatrix<GSL_complex_float>::operator /(
        const CMatrix<GSL_complex_float>&                rval) {
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements matrix division, i.e. *this is multiplied
    //  by rval.
    //
    // POSTCONDITIONS
    //  o A new matrix containing the product is created and returned.
    //
    // HISTORY
    // 18 April 2003
    //         o Template specialization.
    //
    // 19 November 2003
    //	o complex_float templatisation
    //
    // 29 March 2004
    //	o Compatibility checks changed.
    //
    
    // check for compatibility
    if (!compatible(rval))
        _error ("operator/(const CMatrix<GSL_complex_float>& rval)",
                 "Matrices must have identical dimensions");
    CMatrix<GSL_complex_float>        M_result (ps_mat->rows, ps_mat->cols);

    M_result.copy(*this);
    gsl_matrix_complex_float_div_elements(M_result.ps_mat->ps_gsl->matrix_complex_float, 
                            (const gsl_matrix_complex_float*)rval.ps_mat->ps_gsl->matrix_complex_float); 
    return (M_result);
}

CMatrix<GSL_complex_float>
CMatrix<GSL_complex_float>::operator / (
        const GSL_complex_float&                        rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements scalar division, i.e. each matrix element is
    //  is increased by rval.
    //
    // POSTCONDITIONS
    //  o A new matrix containing the scalar product is created and returned.
    //
    // HISTORY
    // 18 April 2003
    //         o Template specialization.
    //
    // 19 November 2003
    //	o complex_float templatisation
    //

    CMatrix<GSL_complex_float>        M_result(ps_mat->rows, ps_mat->cols);
    GSL_complex_float                 inv;
    
    inv        = rval;
    inv.inverse(); 

    M_result.copy(*this);
    gsl_matrix_complex_float_scale (M_result.ps_mat->ps_gsl->matrix_complex_float, inv);                 
    return (M_result);
}

CMatrix<GSL_complex_float>
CMatrix<GSL_complex_float>::operator /= (
        const CMatrix<GSL_complex_float>&                rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements element by element division i.e. each matrix element
    //  multiplied by its corresponding element in rval.
    //
    // POSTCONDITIONS
    //  o The current matrix data values are overwritten.
    //
    // HISTORY
    // 18 April 2003
    //         o Template specialization.
    //
    // 19 November 2003
    //	o complex_float templatisation
    //
    
    // check for compatibility
    char* pch_name = "operator/=(const CMatrix<GSL_complex_float>& rval)";
    if (!compatible (rval))
        _error (pch_name, "Matrices must have identical sizes");
    
    gsl_matrix_complex_float_div_elements(ps_mat->ps_gsl->matrix_complex_float, 
                                    (const gsl_matrix_complex_float*)rval.ps_mat->ps_gsl->matrix_complex_float);

    return *this;    
}

CMatrix<GSL_complex_float>
CMatrix<GSL_complex_float>::operator /= (
        const GSL_complex_float&                        rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements element by constant division i.e. a "scaling"
    //        of matrix values by 1/rval.
    //
    // POSTCONDITIONS
    //  o The current matrix data values are overwritten.
    //
    // HISTORY
    // 18 April 2003
    //         o Template specialization.
    //
    // 19 November 2003
    //	o complex_float templatisation
    //
    
    GSL_complex_float                        inv;
    
    inv        = rval;
    inv.inverse(); 
       
    gsl_matrix_complex_float_scale(ps_mat->ps_gsl->matrix_complex_float, inv); 
    
    return *this;    
}

CMatrix<GSL_complex_float>
CMatrix<GSL_complex_float>::operator - (
        const CMatrix<GSL_complex_float>&          rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements matrix subtraction.
    //
    // POSTCONDITIONS
    //  o A new matrix containing the sum is created and returned.
    //
    // WARNING
    //  o Be aware that the returned matrix needs to be assigned
    //    to an already created matrix otherwise memory leaks
    //    may occur.
    //
    // HISTORY
    // 18 April 2003
    //         o Template specialization.
    //
    // 19 November 2003
    //	o complex_float templatisation
    //
    
    // check for compatibility
    char* pch_name = "operator - (const CMatrix<GSL_complex_float>&rval)";
    if ((ps_mat->rows != rval.rows_get())
        || (ps_mat->cols !=
            rval.cols_get()))_error (pch_name,
                                       "Reason: Incompatible dimensions.");
    CMatrix<GSL_complex_float>         M_sum (ps_mat->rows, ps_mat->cols);
    
    M_sum.copy(*this);
    gsl_matrix_complex_float_sub (
            M_sum.ps_mat->ps_gsl->matrix_complex_float,
        rval.ps_mat->ps_gsl->matrix_complex_float
    );
    return (M_sum);
}

CMatrix<GSL_complex_float>
CMatrix<GSL_complex_float>::operator - (
        const GSL_complex_float&                        rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements scalar subtraction, i.e. each matrix element is
    //  is increased by rval.
    //
    // POSTCONDITIONS
    //  o A new matrix containing the sum is created and returned.
    //
    // HISTORY
    // 18 April 2003
    //         o Template specialization.
    //
    
    CMatrix<GSL_complex_float>        M_sum(ps_mat->rows, ps_mat->cols);
    GSL_complex_float                        neg;
    neg.dat[0]        = -rval.dat[0];
    neg.dat[1]        = -rval.dat[1];

    M_sum.copy(*this);
    gsl_matrix_complex_float_add_constant(M_sum.ps_mat->ps_gsl->matrix_complex_float, neg);
    return(M_sum);
}

CMatrix<GSL_complex_float>
CMatrix<GSL_complex_float>::operator -= (
        const CMatrix<GSL_complex_float>&                rval) {
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements matrix - matrix subtraction.
    //
    // POSTCONDITIONS
    //  o Implements A = A-B matrix type operation, i.e. contents of current
    //          matrix are overwritten.
    //
    // HISTORY
    // 18 April 2003
    //         o Template specialization.
    //
    // 19 November 2003
    //	o complex_float templatisation
    //
      
    gsl_matrix_complex_float_sub (
            ps_mat->ps_gsl->matrix_complex_float,
        rval.ps_mat->ps_gsl->matrix_complex_float
    );
    return (*this);
}

CMatrix<GSL_complex_float>
CMatrix<GSL_complex_float>::operator -= (
        const GSL_complex_float&                        rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements "constant" matrix subtraction.
    //
    // POSTCONDITIONS
    //  o Implements A = A-b scalar type operation, i.e. contents of current
    //          matrix are overwritten.
    //
    // HISTORY
    // 18 April 2003
    //         o Template specialization.
    //
    // 19 November 2003
    //	o complex_float templatisation
    //
      
    GSL_complex_float                        neg;
    neg.dat[0]        = -rval.dat[0];
    neg.dat[1]        = -rval.dat[1];

    gsl_matrix_complex_float_add_constant (
            ps_mat->ps_gsl->matrix_complex_float,
        neg
    );
    return (*this);
}

CMatrix<GSL_complex_float>
CMatrix<GSL_complex_float>::operator + (
        const CMatrix<GSL_complex_float>&          rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements matrix addition.
    //
    // POSTCONDITIONS
    //  o A new matrix containing the sum is created and returned.
    //
    // WARNING
    //  o Be aware that the returned matrix needs to be assigned
    //    to an already created matrix otherwise memory leaks
    //    may occur.
    //
    // HISTORY
    // 18 April 2003
    //         o Template specialization.
    //
    // 19 November 2003
    //	o complex_float templatisation
    //

    // check for compatibility
    char* pch_name = "operator + (const CMatrix<GSL_complex_float>&rval)";
    if ((ps_mat->rows != rval.rows_get())
        || (ps_mat->cols !=
            rval.cols_get()))_error (pch_name,
                                       "Reason: Incompatible dimensions.");
    CMatrix<GSL_complex_float>         M_sum (ps_mat->rows, ps_mat->cols);
    
    M_sum.copy(*this);
    gsl_matrix_complex_float_add (
            M_sum.ps_mat->ps_gsl->matrix_complex_float,
        rval.ps_mat->ps_gsl->matrix_complex_float
    );
    return (M_sum);
}

CMatrix<GSL_complex_float>
CMatrix<GSL_complex_float>::operator + (
        const GSL_complex_float&                        rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements scalar addition, i.e. each matrix element is
    //  is increased by rval.
    //
    // POSTCONDITIONS
    //  o A new matrix containing the sum is created and returned.
    //
    // HISTORY
    // 18 April 2003
    //         o Template specialization.
    //
    // 19 November 2003
    //	o complex_float templatisation
    //
    
    CMatrix<GSL_complex_float>        M_sum(ps_mat->rows, ps_mat->cols);
    
    M_sum.copy(*this);
    gsl_matrix_complex_float_add_constant(M_sum.ps_mat->ps_gsl->matrix_complex_float, rval);
    
    return(M_sum);
}

CMatrix<GSL_complex_float>
CMatrix<GSL_complex_float>::operator += (
        const CMatrix<GSL_complex_float>&          rval) {
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements matrix + matrix addition.
    //
    // POSTCONDITIONS
    //  o Implements A = A+B matrix type operation, i.e. contents of current
    //          matrix are overwritten.
    //
    // HISTORY
    // 18 April 2003
    //         o Template specialization.
    //
    // 19 November 2003
    //	o complex_float templatisation
    //
      
    gsl_matrix_complex_float_add (
            ps_mat->ps_gsl->matrix_complex_float,
        rval.ps_mat->ps_gsl->matrix_complex_float
    );
    return (*this);
}

CMatrix<GSL_complex_float>
CMatrix<GSL_complex_float>::operator += (
        const GSL_complex_float&                        rval)
{
    //
    // ARGS
    //  rval                    in              right hand value
    //
    // DESC
    //  Implements "constant" matrix addition.
    //
    // POSTCONDITIONS
    //  o Implements A = A+b scalar type operation, i.e. contents of current
    //          matrix are overwritten.
    //
    // HISTORY
    // 18 April 2003
    //         o Template specialization.
    //
    // 19 November 2003
    //	o complex_float templatisation
    //
    
    gsl_matrix_complex_float_add_constant (
            ps_mat->ps_gsl->matrix_complex_float,
        rval
    );
    return (*this);
}


//////----------><----------\\\\\\
//////      CVol class      \\\\\\
//////----------><----------\\\\\\

template<typename _CMDATA>
int CVol<_CMDATA>::createdCounter        = 0;

template<typename _CMDATA>
void
CVol<_CMDATA>::_warn (
        char *apch_proc,
        char *apch_msg,
        int   code) const {
    //
    // ARGS
    //  apch_proc               in              (optional) function name
    //                                                  in which error has
    //                                                  occured
    //  apch_msg                in              (optional) error message
    //                                                  text
    //  code                    in              (optional) error code
    //
    // DESC
    //  Main warn reporting method for the CVol class.
    //
    // HISTORY
    // 19 September 2003
    // o Initial design and coding.
    //

    cerr << "\nWARNING encountered.\n";
    cerr << "\tCVol object id: " << ps_vol->id << endl;
    cerr << "\tps_vol:         " << ps_vol     << endl;
    cerr << "\tppzM_volume:    " << ps_vol->ppzM_volume << endl;
    cerr << "\trows:           " << rows_get() << endl;
    cerr << "\tcolumns:        " << cols_get() << endl;
    cerr << "\tslices:         " << slices_get() << endl;
    cerr << "\tCurrent function: " << "CVol::" << apch_proc << "\n";
    cerr << "\t" << apch_msg << "\n";
}

template<typename _CMDATA>
void
CVol<_CMDATA>::_error (
        char *apch_proc,
        char *apch_msg,
        int   code) {
    //
    // ARGS
    //  apch_proc               in              (optional) function name
    //                                                  in which error has
    //                                                  occured
    //  apch_msg                in              (optional) error message
    //                                                  text
    //  code                    in              (optional) error code
    //
    // DESC
    //  Main error reporting method for the CVol class.
    //
    // TODO
    //        o Wrap exception throwing around #ifdefs?
    //
    // HISTORY
    // 19 September 2003
    // o Initial design and coding.
    //

    cerr << "\nFatal error encountered.\n";
    cerr << "\tCVol object id: " << ps_vol->id << endl;
    cerr << "\tps_vol:         " << ps_vol     << endl;
    cerr << "\tppzM_volume:    " << ps_vol->ppzM_volume << endl;
    cerr << "\trows:           " << rows_get() << endl;
    cerr << "\tcolumns:        " << cols_get() << endl;
    cerr << "\tslices:         " << slices_get() << endl;
    cerr << "\tCurrent function: " << "CVol::" << apch_proc << "\n";
    cerr << "\t" << apch_msg << "\n";
    cerr << "Throwing an exception to (this) with code " << code << endl;
    throw (this);
}

//////---------->
////// Constructor, destructor, and miscellaneous initialization
////// routines
//////---------->

template<typename _CMDATA>
void
CVol<_CMDATA>::coreVolume_construct(
        int                 arows,
        int                 acols,
	int                 aslices
) {
    // ARGS
    //        arows                in                        number of rows in matrix
    //        acols                in                        number of cols in matrixtype
    //
    // DESCRIPTION
    //        This method is a consolidation of the main initialization code
    //        that is shared across all the *matrix* constructors.
    //
    // PRECONDITIONS
    //        o Should only be called from a constructor method.
    //
    // POSTCONDITIONS
    //        o A gsl-aware ps_mat structure is created.
    //
    // HISTORY
    // 19 March 2003
    //        o Initial design and coding.
    //

    // create the structure
    try {
        ps_vol = new volstruct;
    } catch (bad_alloc xa) {
        _error("base constructor",
                "could not allocate memory for volstruct (heap exhausted?)");
    }
    // Allocate memory for volume container
    ps_vol->ppzM_volume    = new CMatrix<_CMDATA>* [aslices];
    _CMDATA	zero	   = (_CMDATA) 0;
    for(int i=0; i<aslices; i++) {
	ps_vol->ppzM_volume[i]    = new CMatrix<_CMDATA>(  arows,
                                                           acols,
                                                           zero);
    }

    ps_vol->slices	= aslices;

    // set first reference/id to this data
    ps_vol->refcnt       = 1;
    ps_vol->id           = CVol::createdCounter++;
}

template <typename _CMDATA>
void
CVol<_CMDATA>::coreVolume_destruct() {
    //
    // DESC
    //  Explicitly destroys the core data stuctures of a volume
    //  object.
    //
    // HISTORY
    // 19 September 2003
    //	o Reza's birthday! I miss you, love. Wish I were there in SA with you!
    //  o Initial design and coding.
    //

    for(int i=0; i<ps_vol->slices; i++)
        delete ps_vol->ppzM_volume[i];

    delete [] ps_vol->ppzM_volume;

    delete ps_vol;

    //ps_vol		= NULL;
}

template <typename _CMDATA>
void
CVol<_CMDATA>::coreVolume_replaceWith(
    CVol<_CMDATA>*	pVl) {

    //
    // ARGS
    //	pVl		in		pointer to volume whose core will be used
    //						to replace *this core
    //
    // DESC
    //	This method should only be used in certain rarely occuring high-specific
    //	situations! If used incorrectly it can lead to very hard to trace memory
    //	leaks.
    //
    //  Its purpose is to map the ps_vol pointer of *this to the ps_vol pointer
    //  of pVl.
    //
    //	Why would one use this? As a sort of hack to quickly replace the contents
    //	of one volume with another. Given two volumes A and B, and wanting to
    //  replace A with B, but still keep the main pointer to A intact, one would
    //  coreVolume_destruct A, and then coreVolume_replaceWith(B). B then effectively
    //	ceases to exist, at least from a memory management point of view.
    //
    //	Take care *not* to coreVolume destruct B! This will destruct A as well!
    //
    // POSTCONDITIONS
    //	o this->ps_vol is replaced with pVl->ps_vol. Note that this is a pointer-only
    //	  operation! No memory data is shuffled around, only pointers to this
    //	  data.
    //
    // HISTORY
    // 02 October 2003
    //	o Initial design and coding.
    //

    ps_vol	= 	pVl->ps_vol;

}


template <typename _CMDATA>
CVol<_CMDATA>::~CVol () {
    //
    // DESC
    //  Destructor.
    //
    // POSTCONDITIONS
    //  Decrements reference count. If counter reaches zero, frees allocated
    //  memory.
    //
    // 27 January 2004
    //  o Added refcnt_ methods
    //

    //if (--ps_vol->refcnt == 0) {
    if (refcnt_dec() == 0) {
	coreVolume_destruct();
    }
}

template <typename _CMDATA>
CVol<_CMDATA>::CVol(
        int             mrows,
        int             mcols,
	int             mslices) {
    //
    // ARGS
    //  mrows                   in              number of rows in each slice
    //  mcols                   in              number of cols in each slice
    //  mslices                 in              number of slices in volume
    //
    // DESC
    //  Volume constructor
    //
    // PRECONDITIONS
    //  o Make sure that mrows and mcols are not negative numbers!
    //
    // HISTORY
    // 19 September 2003
    // o Initial design and coding.
    //

    coreVolume_construct(mrows, mcols, mslices);

}

template <typename _CMDATA>
CVol<_CMDATA>::CVol(
	CMatrix<_CMDATA>**	appzM_volume,
	int			aslices
) {
    //
    // ARGS
    //	ppzM_volume		in		pointer to array of CMatrix slices
    //  mslices                 in              number of slices in volume
    //
    // DESC
    //  Volume constructor. This constructor assumes an already allocated block
    //	of memory. It is typically most useful in "wrapping" a class object
    //	around data that was originally created outside of this class definition.
    //
    // PRECONDITIONS
    //  o Make sure that ppzM_volume points to a valid block of memory and that
    //	  slices a correct value.
    //
    // HISTORY
    // 22 September 2003
    //	o Initial design and coding.
    //

    // create the structure
    try {
        ps_vol = new volstruct;
    } catch (bad_alloc xa) {
        _error("base constructor",
                "could not allocate memory for volstruct (heap exhausted?)");
    }

    ps_vol->ppzM_volume		= appzM_volume;
    ps_vol->slices		= aslices;

    // set first reference/id to this data
    ps_vol->refcnt = 1;
    ps_vol->id           = CVol::createdCounter++;


}

template <typename _CMDATA>
CVol<_CMDATA>::CVol(
    const CVol<_CMDATA> & z)
{
    //
    // DESC
    //  Copy constructor - simple pointer based.
    //
    // POSTCONDITIONS
    //  o NB!! NB!! NB!!
    //    The copy constructor merely increases the reference count of the
    //    "source" data, and directs the target to point to the source.
    //    This is *not* a deepcopy!. Although allowing for fast copies between
    //    matrices, it can potentially suffer from problems relating to scope
    //    local variable variable problems. If used in function arguments
    //    or as returns out of functions, then it is not really a problem.
    //
    // 27 January 2004
    //  o Added refcnt_ methods
    //

    //z.ps_vol->refcnt++;                 // adding another reference
    z.refcnt_inc();                     // adding another reference
    ps_vol = z.ps_vol;                  // point to the new matstruct
}

template<typename _CMDATA>
CVol<_CMDATA>
CVol<_CMDATA>::copy (
        const CVol&          source)
{
    // ARGS
    //  source                  in              volume whose contents to copy
    //                                                  into *this
    //
    // DESC
    //  Explicit deepcopy.
    //
    // POSTCONDITIONS.
    //  o The data in source is copied over into *this. Space for this data
    //    is created as in necessary.
    //
    // HISTORY
    // 20 September 2003
    //  o Adaptation from CMatrix code.
    //

    int                row, col, slice;

    if (!compatible (source)) {
        coreVolume_destruct();
        coreVolume_construct(source.rows_get(), source.cols_get(), source.slices_get());
    }
    for (row = 0; row < rows_get(); row++)
        for (col = 0; col < cols_get(); col++)
            for(slice = 0; slice < slices_get(); slice++)
                _mval(row, col, slice) = source._mval (row, col, slice);
    return *this;
}

template <typename _CMDATA>
CVol<_CMDATA>
CVol<_CMDATA>::zeroPad(
        e_DOMINANCE             dir,
        int                     size )
{
    // ARGS
    //  dir                     in              dimension along which to zero pad.
    //  size                    in              depth along the dimension to zero pad.
    //                                              note that this prepended and
    //                                              appended to the dimension.
    //
    // DESC
    //  Performs a 3D volume zero padding operation. Additional memory containing
    //  zeroes is added to the front and back of the dimension specified by
    //  'dir'. The length of this addition is 'size', thus the total volume
    //  increased by a factor of 2*size.
    //
    // PRECONDITIONS
    //  o None.
    //
    // POSTCONDITIONS
    //  o A local volume is constructed, and copy constructed back out of this method.
    //  o Additional volume is added to the front *and* back of the specified
    //    dimension, thus the target dimension increases with 2*size.
    //
    // HISTORY
    // 30 September 2003
    //  o Object encapsulation.
    //

    char*   pch_proc    = "zeroPad(...)";

    CVol<_CMDATA>*      pVl_output;

    int     totalSlices = 0;            // These define the new dimensions
    int     totalRows   = 0;            //      of the created volume.
    int     totalCols   = 0;

    int	    inputRows	= 0;		// These define the dimensions 
    int	    inputCols	= 0;		//	of the current volume.
    int	    inputSlices	= 0;

    inputRows   =       rows_get();
    inputCols   =       cols_get();
    inputSlices	=	slices_get();

    int     i;

     switch(dir) {
        case e_row:
            totalRows   =   inputRows           + 2*size;
            totalCols   =   inputCols;
            totalSlices =   inputSlices;
            break;
        case e_column:
            totalRows   =   inputRows;
            totalCols   =   inputCols           + 2*size;
            totalSlices =   inputSlices;
            break;
        case e_slice:
            totalRows   =   inputRows;
            totalCols   =   inputCols;
            totalSlices =   inputSlices        + 2*size;
            break;
     }

    // Allocate memory for output volume container
    pVl_output  = new CVol<_CMDATA>(totalRows, totalCols, totalSlices);

    if(dir!=e_slice) {
        for(i=0; i<totalSlices; i++)
            pVl_output->slice(i) = slice(i).zeroPad(dir, size);
    } else {
        for(i=size; i<totalSlices-size; i++) {
            pVl_output->slice(i) = slice(i-size);
        }
    }

    CVol<_CMDATA>       V_return(*pVl_output);
    delete              pVl_output;
    return              V_return;
}

template <typename _CMDATA>
CVol<_CMDATA>
CVol<_CMDATA>::fftshift
(
        bool    ab_inPlace
) {
    //
    // ARGS
    //  ab_inPlace              in              if true, branch to the
    //                                                  inPlace() method
    //                                          else, branch to the
    //                                                  newMem() method
    //
    // DESCRIPTION
    //  shift*() dispatching layer. Depending on the value of ab_inPlace,
    //  either an inPlace or new memory operation is performed.
    //
    // HISTORY
    //  22 January 2003
    //  o Initial design and coding.
    //

    if(ab_inPlace) {
        shiftInPlace(+1);
        return *this;
    } else {
        CVol<_CMDATA> V_dispatch(*this);
//        CVol<_CMDATA> V_dispatch(1, 1, 1);
        V_dispatch  = shiftNewMem(+1);
        return(V_dispatch);
    }
}

template <typename _CMDATA>
CVol<_CMDATA>
CVol<_CMDATA>::ifftshift(
        bool    ab_inPlace
) {
    //
    // ARGS
    //  ab_inPlace              in              if true, branch to the
    //                                                  inPlace() method
    //                                          else, branch to the
    //                                                  newMem() method
    //
    // DESCRIPTION
    //  shift*() dispatching layer. Depending on the value of ab_inPlace,
    //  either an inPlace or new memory operation is performed.
    //
    // HISTORY
    //  22 January 2003
    //  o Initial design and coding.
    //

    if(ab_inPlace) {
        shiftInPlace(-1);
        return *this;
    } else {
        CVol<_CMDATA> V_dispatch(*this);
        V_dispatch  = shiftNewMem(-1);
        return(V_dispatch);
    }
}

template <typename _CMDATA>
CVol<_CMDATA>
CVol<_CMDATA>::shiftInPlace(
        int                     a_pivotOffset
) {
    //
    // ARGS
    //  a_pivotOffset           in              pivot offset, i.e. the direction
    //                                                  in which to shift
    //
    // DESC
    //  Performs a volume shift operation - either fftshift or ifftshift - depending
    //  on the a_pivotOffset.
    //
    // PRECONDITIONS
    //
    // POSTCONDITIONS
    //  o Destructive shift operation.
    //
    // HISTORY
    //  21 January 2003
    //  o Initial design and coding.
    //

    char*   pch_proc    = "shiftInPlace(...)";

    int     i;
    int     rows;
    int     cols;
    int     slices;

    int         j;
    int         k=0;
    rows    = rows_get();
    cols    = cols_get();
    slices  = slices_get();

    // Shift the contents of each volume "slice"
    for(i=0; i<slices; i++) {
        slice(i).shiftInPlace(a_pivotOffset);
    }

    // Now, also shift the slice indices themselves:
    CMatrix<int>                M_list(1, slices);
    CMatrix<_CMDATA>**          ppzM_origOrder  = new CMatrix<_CMDATA>* [slices];
    CMatrix<_CMDATA>**          ppzM_order      = ppzM_volume_get();

    // ... make a backup of the original output list to conserve the pointers
    for(i=0; i<slices; i++) {
        M_list(0, i)        = i;
        ppzM_origOrder[i]   = ppzM_order[i];
    }

    // ... shift the slice indices
    M_list.shiftInPlace(a_pivotOffset);

    // ... reorganize the actual slices
    for(i=0; i<slices; i++) {
        ppzM_order[i]  = ppzM_origOrder[M_list(0, i)];
    }

    delete [] ppzM_origOrder;
    return *this;
}

template <typename _CMDATA>
CVol<_CMDATA>
CVol<_CMDATA>::shift(
        e_FFTSHIFTDIR           e_dir
) {
    //
    // ARGS
    //  e_dir                   in              the direction in which to fftshift
    //
    // DESC
    //  Performs a volume shift operation - either fftshift or ifftshift.
    //
    // PRECONDITIONS
    //  o Since the *shift() operation is non-destructive, a new memory block is
    //    allocated that will contain the result of the operation.
    //  o All internal shift calls are assumed to be non-inPlace.
    //
    // POSTCONDITIONS
    //  o A locally constructed volume will be copy constructed out of this method.
    //
    // HISTORY
    // 11 September 2003
    //  o Initial design and coding.
    //
    // 30 September 2003
    //  o Object encapsulation.
    //
    // 21 January 2003
    //  o Updated for changes in underlying *shift() API
    //

    char*       pch_proc    = "shift(...)";

    int         i;
    int         rows;
    int         cols;
    int         slices;

    bool        b_inPlace       = false;

    rows    = rows_get();
    cols    = cols_get();
    slices  = slices_get();

    CVol<_CMDATA>       V_output(rows, cols, slices);

    for(i=0; i<slices; i++) {
        switch(e_dir) {
            case e_ifftshift:
                V_output.slice(i)       = slice(i).ifftshift(b_inPlace);
                break;
            case e_fftshift:
                V_output.slice(i)       = slice(i).fftshift(b_inPlace);
                break;
        }
    }

    CMatrix<int>                M_list(1, slices);
    CMatrix<int>                M_listShifted(1, slices);
    CMatrix<_CMDATA>**          ppzM_origOrder  = new CMatrix<_CMDATA>* [slices];
    CMatrix<_CMDATA>**          ppzM_order      = V_output.ppzM_volume_get();

    // Make a backup of the original output list to conserve the pointers
    for(i=0; i<slices; i++) {
        M_list(0, i)        = i;
        ppzM_origOrder[i]   = ppzM_order[i];
    }

    switch(e_dir) {
        case e_ifftshift:
            M_listShifted   = M_list.ifftshift(b_inPlace);
            break;
        case e_fftshift:
            M_listShifted   = M_list.fftshift(b_inPlace);
            break;
    }
    // reorganize the actual slices
    for(i=0; i<slices; i++) {
        ppzM_order[i]  = ppzM_origOrder[M_listShifted(0, i)];
    }

    delete [] ppzM_origOrder;
    return V_output;
}

template <typename _CMDATA>
CVol<_CMDATA>
CVol<_CMDATA>::shiftNewMem(
        int                     a_pivotOffset
) {
    //
    // ARGS
    //  a_pivotOffset           in              pivot offset, i.e. the direction
    //                                                  in which to shift
    //
    // DESC
    //  Performs a volume shift operation - either fftshift or ifftshift - depending
    //  on the a_pivotOffset.
    //
    // PRECONDITIONS
    //
    // POSTCONDITIONS
    //  o Destructive shift operation.
    //
    // HISTORY
    //  21 January 2003
    //  o Initial design and coding.
    //

    char*   pch_proc    = "shiftNewMem(...)";

    int     i;
    int     rows;
    int     cols;
    int     slices;

    rows    = rows_get();
    cols    = cols_get();
    slices  = slices_get();

    CVol<_CMDATA>       V_output(rows, cols, slices);

    for(i=0; i<slices; i++)
        V_output.slice(i)       = slice(i).shiftNewMem(a_pivotOffset);

    CMatrix<int>                M_list(1, slices);
    CMatrix<int>                M_listShifted(1, slices);
    CMatrix<_CMDATA>**          ppzM_origOrder  = new CMatrix<_CMDATA>* [slices];
    CMatrix<_CMDATA>**          ppzM_order      = V_output.ppzM_volume_get();

    // Make a backup of the original output list to conserve the pointers
    for(i=0; i<slices; i++) {
        M_list(0, i)        = i;
        ppzM_origOrder[i]   = ppzM_order[i];
    }

    M_listShifted       = M_list.shiftNewMem(a_pivotOffset);

    // reorganize the actual slices
    for(i=0; i<slices; i++) {
        ppzM_order[i]  = ppzM_origOrder[M_listShifted(0, i)];
    }

    delete [] ppzM_origOrder;
    return V_output;
}

template <typename _CMDATA>
void
CVol<_CMDATA>::saveBinary(
    string		astr_fileName)
{
    //
    // ARGS
    //	astr_fileName		in 		file name to save
    //
    // DESC
    //	Performs a simple binary mode save of volume.
    //
    // HISTORY
    // 20 October 2003
    //	o Initial design and coding.
    //

    ofstream	fout(astr_fileName.c_str(), ios_base::out | ios_base::binary | ios_base::trunc);
    if(!fout.is_open())
	_error("saveBinary", "Could not create file for saving", 1);

    int	rows = rows_get();
    int	cols = cols_get();
    int	slices = slices_get();

    fout.write((const char*)&rows, sizeof(int));
    fout.write((const char*)&cols, sizeof(int));
    fout.write((const char*)&slices, sizeof(int));
    
    for(int k=0; k<slices; k++)
	for(int i=0; i<rows; i++)
	    for(int j=0; j<cols; j++)
		fout.write((const char*)&val(i, j, k), sizeof(_CMDATA));
        
    fout.close();
}

template <typename _CMDATA>
CVol<_CMDATA>
CVol<_CMDATA>::loadBinary(	
    string		astr_fileName) 
{
    //
    // ARGS
    //	astr_fileName		in 		file name to load
    //
    // DESC
    //	Performs a simple binary mode load of volume.
    //
    // HISTORY
    // 20 October 2003
    //	o Initial design and coding.
    //
    
    ifstream	fin(astr_fileName.c_str(), ios_base::in | ios_base::binary);
    if(!fin.is_open())
	_error("loadBinary", (char*)
	       (string("Could not open target file ") + astr_fileName).c_str(), 1);

    int	rows;
    int	cols;
    int	slices;

    fin.read((char*)&rows, sizeof(int));
    fin.read((char*)&cols, sizeof(int));
    fin.read((char*)&slices, sizeof(int));

    if (!compatible (rows, cols, slices)) {
        coreVolume_destruct();
        coreVolume_construct(rows, cols, slices);
    }
    
    for(int k=0; k<slices; k++)
	for(int i=0; i<rows; i++)
	    for(int j=0; j<cols; j++)
		fin.read((char*)&val(i, j, k), sizeof(_CMDATA));
        
    fin.close();
}

template <typename _CMDATA>
void
CVol<_CMDATA>::print (
    char*                apch_msg        /*= "ans"       */,
    int                  a_precision     /*= 6           */,
    int                  a_width         /*= 12          */,
    ios::fmtflags        a_userFormat    /*= 0           */,
    int                  a_marginOffset  /*= 0           */,
    char                 ach_left        /*= (char) 0    */,
    int                  a_tabOffset     /*= 0           */,
    bool                 ab_fancy        /*= 1           */) {
    //
    // ARGS
    //  apch_msg                in        message that prepends matrix dump
    //  a_precision             in        the precision of numerical output
    //  a_width                 in        the width of the numerical field
    //  a_userFormat            in        if non-zero, set stream format to arg
    //  a_marginOffset          in        the amount of tab spaces prepending ch_left
    //  ach_left                in        also somewhat primitive. If non-zero, will
    //                                          print (char) left_char as left most character
    //                                          of each new line. Useful for when the matrix is
    //                                          part of a larger output dump
    //  a_tabOffset             in          the amount of tab spaces prepending dump
    //  ab_fancy                in          still rather primitive. If false, will print
    //                                                *only* the numbers, nothing else
    //
    // DESC
    //
    //	Simple wrapper around CMatrix::print method. For each slice of volume,
    //	calls the print() method for that slice.
    //
    // NOTE
    //         The `marginOffset' and `tabOffset' are used in conjunction
    //         with ch_left as follows:
    //
    //                 [<--- marginOffset --->]<ch_left>[<--- tabOffset --->]
    //
    // HISTORY
    // 22 September 2003
    //	o Initial design and coding.
    //

    for(int i=0; i<ps_vol->slices; i++) {
	cout << apch_msg << " slice: " << i << endl;
	ps_vol->ppzM_volume[i]->print(
	    apch_msg        /*= "ans"       */,
	    a_precision     /*= 6           */,
	    a_width         /*= 12          */,
	    a_userFormat    /*= 0           */,
	    a_marginOffset  /*= 0           */,
	    ach_left        /*= (char) 0    */,
	    a_tabOffset     /*= 0           */,
	    ab_fancy        /*= 1           */);
    }

}

template<typename _CMDATA>
_CMDATA &
CVol<_CMDATA>::val(
        int         row,
        int         col,
	int         slice)
{
    //
    // ARGS
    //        row                in                row to access
    //        col                in                col to access
    //        slice              in                slice to access
    //
    // DESC
    //         Element selection: can be used to read or write.
    //
    // PRECONDITIONS
    //        o Range checking can be expensive, particularly in tightly
    //          nested iterative loops. By defining CMATRIX_NORANGE range
    //          checking is disabled. This will result in tremendous performance
    //          boost, but at the risk of a single out of bound access
    //          killing the entire process.
    //
    // HISTORY
    // 22 September 2003
    //	o Adaptation from CMatrix.
    //

#ifndef CMATRIX_NORANGE
    char *
        pch_name = "val";
    char *
        pch_errormsg = "Addressing error. Index out of range";

    if(!(slice>=0 && slice<ps_vol->slices))
        _error (pch_name, pch_errormsg);
    else
#endif
    return (ps_vol->ppzM_volume[slice]->val(row, col));
}

template<typename _CMDATA>
_CMDATA &
CVol<_CMDATA>::operator()(
        int         row,
        int         col,
	int         slice)
{
    //
    // ARGS
    //        row                in                row to access
    //        col                in                col to access
    //        slice              in                slice to access
    //
    // DESC
    //         Element selection: can be used to read or write.
    //
    // PRECONDITIONS
    //        o Range checking can be expensive, particularly in tightly
    //          nested iterative loops. By defining CMATRIX_NORANGE range
    //          checking is disabled. This will result in tremendous performance
    //          boost, but at the risk of a single out of bound access
    //          killing the entire process.
    //
    // HISTORY
    // 22 September 2003
    //	o Adaptation from CMatrix.
    //

#ifndef CMATRIX_NORANGE
    char *
        pch_name = "val";
    char *
        pch_errormsg = "Addressing error. Index out of range";
    if(!(slice>=0 && slice<ps_vol->slices))
        _error (pch_name, pch_errormsg);
    else
#endif
    return (ps_vol->ppzM_volume[slice]->val(row, col));
}


//////---------->
////// Access and set internal entire slices
//////---------->
template<typename _CMDATA>
CMatrix<_CMDATA> &
CVol<_CMDATA>::slice(
        int         slice
) {
    //
    // ARGS
    //        slice              in                slice (in slice direction) to access
    //
    // DESC
    //         slice selection: can be used to read or write.
    //
    // PRECONDITIONS
    //  o Use this method to set/get individual "slices" (where a slice is a plane normal
    //     to the slice direction).
    //
    // HISTORY
    // 30 September 2003
    //  o Initial design and coding.
    //

#ifndef CMATRIX_NORANGE
    char *
        pch_name = "val";
    char *
        pch_errormsg = "Addressing error. Slice out of range";
    if (!(slice>= 0 && slice < slices_get()))
        _error (pch_name, pch_errormsg);
    else
#endif
    return (*ps_vol->ppzM_volume[slice]);
}

template<typename _CMDATA>
CMatrix<_CMDATA> &
CVol<_CMDATA>::operator()(
        int         slice)
{
    //
    // ARGS
    //        slice              in                slice to access
    //
    // DESC
    //         slice selection: can be used to read or write.
    //
    // PRECONDITIONS
    //  o Use this method to set/get individual "slices" (where a slice is a plane normal
    //     to the slice direction).
    //
    // HISTORY
    // 30 September 2003
    //  o Initial design and coding.
    //

#ifndef CMATRIX_NORANGE
    char *
        pch_name = "val";
    char *
        pch_errormsg = "Addressing error. Index out of range";
    if(!(slice>=0 && slice<ps_vol->slices))
        _error (pch_name, pch_errormsg);
    else
#endif
    return (*ps_vol->ppzM_volume[slice]);
}

template<typename _CMDATA>
CMatrix<_CMDATA>
CVol<_CMDATA>::plane_get(
        int	                slice,
        e_DOMINANCE             e_dir	/* = e_slice*/) {
    //
    // ARGS
    //	slice			in		slice index to extract
    //	e_dir			in		"direction" along which to extract slice
    //
    // DESC
    //	This method removes a "plane", i.e. a single matrix, from a volume.
    //	By specifying an optional "direction", the plane can be orientated
    //	normal to this direction, allowing extraction normal to any orthogonal
    //  direction.
    //
    // PRECONDITIONS
    //	o None.
    //  o This method is called "plane" specifically to indicate that a matrix
    //    normal to any of the orthogonal directions can be returned.
    //
    // POSTCONDITIONS
    //	o A plane is extracted (created) and returned.
    //	o The row/col sense of returned slices maintain context in that the
    //	  "row" corresponds to the "lower" enumeration of the non-slice directions,
    //	  and the "col" corresponds to the "upper" enumeration of the non-slice
    //	  directions.
    //	o Be sure to "capture" the slice so that it can eventually be destroyed!
    //
    // HISTORY
    // 30 September 2003
    //	o Initial design and coding.
    //

    CMatrix<_CMDATA>*	pM;
    int			i, j, k;

    switch(e_dir) {
	case e_row:
	    pM		= new CMatrix<_CMDATA>(cols_get(), slices_get());
	    for(i=0; i<cols_get(); i++)
		for(j=0; j<slices_get(); j++)
		    pM->val(i, j)	= val(slice, i, j);
	break;
	case e_column:
	    pM		= new CMatrix<_CMDATA>(rows_get(), slices_get());
	    for(i=0; i<rows_get(); i++)
		for(j=0; j<slices_get(); j++)
		    pM->val(i, j)	= val(i, slice, j);
	break;
	case e_slice:
	    pM		= new CMatrix<_CMDATA>(rows_get(), cols_get());
	    for(i=0; i<rows_get(); i++)
		for(j=0; j<cols_get(); j++)
		    pM->val(i, j)	= val(i, j, slice);
	break;
    }
    CMatrix<_CMDATA>	M_ret(*pM);
    delete pM;
    return M_ret;
}

template<typename _CMDATA>
void
CVol<_CMDATA>::plane_set(
    CMatrix<_CMDATA>&	M_replacement,
    int	                slice,
    e_DOMINANCE         e_dir	/* = e_slice*/) {
    //
    // ARGS
    //	M_replacement		in		the replacement matrix
    //	slice			in		slice index to insert
    //	e_dir			in		"direction" along which to insert slice
    //
    // DESC
    //	This method replaces a "slice", in a volume.
    //	By specifying an optional "direction", the slice can be orientated
    //	normal to this direction, allowing replacement of any arbitrary
    //	matrix.
    //
    // PRECONDITIONS
    //	o None.
    //  o This method is called "plane" specifically to indicate that a matrix
    //    normal to any of the orthogonal directions can be replaced.
    //
    // POSTCONDITIONS
    //	o A slice is replaced.
    //	o The row/col sense of returned slices maintain context in that the
    //	  "row" corresponds to the "lower" enumeration of the non-slice directions,
    //	  and the "col" corresponds to the "upper" enumeration of the non-slice
    //	  directions.
    //	o Logical inverse of slice_get(...)
    //
    // HISTORY
    // 30 September 2003
    //	o Initial design and coding.
    //

    char*               pch_proc = "slice_set(...)";
    int			i, j, k;

    switch(e_dir) {
	case e_row:
	    if(!(M_replacement.compatible(cols_get(), slices_get())))
		_error(pch_proc, "Replacement matrix will not fit in row direction.", 1);
	    for(i=0; i<cols_get(); i++)
		for(j=0; j<slices_get(); j++)
		    val(slice, i, j) = M_replacement(i, j);
	break;
	case e_column:
	    if(!(M_replacement.compatible(rows_get(), slices_get())))
		_error(pch_proc, "Replacement matrix will not fit in column direction.", 1);
	    for(i=0; i<rows_get(); i++)
		for(j=0; j<slices_get(); j++)
		    val(i, slice, j) = M_replacement(i, j);
	break;
	case e_slice:
	    if(!(M_replacement.compatible(rows_get(), cols_get())))
		_error(pch_proc, "Replacement matrix will not fit in slice direction.", 1);
	    for(i=0; i<rows_get(); i++)
		for(j=0; j<cols_get(); j++)
		    val(i, j, slice) = M_replacement(i, j);
	break;
    }
}


//////---------->
////// Fourier transforms
//////---------->

// GSL_complex (double) FFT
template<typename _CMDATA>
CVol<_CMDATA>
CVol<_CMDATA>::fft3D(
        e_FFTDIR                e_direction     /*= e_forward   */,
        bool                    b_MKLwrap       /*= false       */)
{
    // ARGS
    //  e_direction             in              forward/inverse FFT
    //  b_MKLwrap               in              if true, wrap around core MKL
    //                                                  routines, else
    //                                                  revert to GSL routines.
    //
    // DESC
    //  Performs a 3D volume FFT by calling a 1D FFT routine across all the
    //  vectors comprising the volume along all three dimensions.
    //
    //  This method calls 2D FFTs in each slice and then 1D FFTs along each line
    //  in the slice direction. It is a factor 2 faster than the FFT3Dl method,
    //  and is thus favourably comparable with MatLab.
    //
    // PRECONDITIONS
    //
    // POSTCONDITIONS
    //  o *this is returned.
    //	o Note that this is a *destructive* operation! Current volume is replaced
    //	  with its FFT!
    //

    char*                       pch_proc = "CVol<_CMDATA>::fft3D(...)";
    int                         i, j, k;
    CMatrix<_CMDATA>*       pzM_1D;


    int rows    = rows_get();
    int cols    = cols_get();

    // Process all matrices orthogonal to the `k' (depth) direction
    CMatrix<_CMDATA>        zM_2D(rows, cols);
    for(k=0; k<slices_get(); k++) {
        zM_2D = slice(k);
        zM_2D.fft2D(e_direction, b_MKLwrap);
	slice(k).copy(zM_2D);
    }

    // Process all vectors pointing in the `k' (slice) dimension direction
    pzM_1D      = new CMatrix<_CMDATA>(1, slices_get());
    for(i=0; i<rows; i++) {
        for(j=0; j<cols; j++) {
            // Extract the vector
            for(k=0; k<slices_get(); k++)
                pzM_1D->val(0, k) = val(i, j, k);
            // 1D FFT - do an inPlace on this dummy vector
            pzM_1D->fft1D(e_direction, b_MKLwrap);
            // Insert back into volume
            for(k=0; k<slices_get(); k++)
                val(i, j, k)       = pzM_1D->val(0, k);
        }
    }
    delete pzM_1D;
    return *this;
}


#if(0)
CVol<GSL_complex_float>
CVol<GSL_complex_float>::fft3D(
        e_FFTDIR                e_direction     /*= e_forward   */,
        bool                    b_MKLwrap       /*= false       */)
{
    // ARGS
    //  e_direction             in              forward/inverse FFT
    //  b_MKLwrap               in              if true, wrap around core MKL
    //                                                  routines, else
    //                                                  revert to GSL routines.
    //
    // DESC
    //  Performs a 3D volume FFT by calling a 1D FFT routine across all the
    //  vectors comprising the volume along all three dimensions.
    //
    //  This method calls 2D FFTs in each slice and then 1D FFTs along each line
    //  in the slice direction. It is a factor 2 faster than the FFT3Dl method,
    //  and is thus favourably comparable with MatLab.
    //
    // PRECONDITIONS
    //
    // POSTCONDITIONS
    //  o *this is returned.
    //	o Note that this is a *destructive* operation! Current volume is replaced
    //	  with its FFT!
    //

    char*                       pch_proc = "CVol<_CMDATA>::fft3D(...)";
    int                         i, j, k;
    CMatrix<GSL_complex_float>* pzM_1D;


    int rows    = rows_get();
    int cols    = cols_get();

    // Process all matrices orthogonal to the `k' (depth) direction
    CMatrix<GSL_complex_float>        zM_2D(rows, cols);
    for(k=0; k<slices_get(); k++) {
        zM_2D = slice(k);
        zM_2D.fft2D(e_direction, b_MKLwrap);
	slice(k).copy(zM_2D);
    }

    // Process all vectors pointing in the `k' (slice) dimension direction
    pzM_1D      = new CMatrix<GSL_complex_float>(1, slices_get());
    for(i=0; i<rows; i++) {
        for(j=0; j<cols; j++) {
            // Extract the vector
            for(k=0; k<slices_get(); k++)
                pzM_1D->val(0, k) = val(i, j, k);
            // 1D FFT - do an inPlace on this dummy vector
            pzM_1D->fft1D(e_direction, b_MKLwrap);
            // Insert back into volume
            for(k=0; k<slices_get(); k++)
                val(i, j, k)       = pzM_1D->val(0, k);
        }
    }
    delete pzM_1D;
    return *this;
}
#endif

///////---------->
////// Operator overloads
//////---------->

template<typename _CMDATA>
CVol<_CMDATA>
CVol<_CMDATA>::operator = (
        const CVol<_CMDATA>& rval) {
    //
    // ARGS
    //  rval            in              right-hand argument to operator
    //
    // DESC
    //  Assignment operator.
    //
    // PRECONDITIONS
    //  o A "pointer" equivalence is performed, i.e. *no* deepcopy. The
    //    pointer to *this matrix's data is simply decremented (and
    //    possibly destroyed) and reset to point to rval's data (which
    //    has its refcount increased.
    // o *this must have the same size as rval
    //
    // HISTORY
    // 20 March 2003
    //  o coreMatrix_[construct/destruct] integration.
    //
    // 26 March 2003
    //  o Templatization.
    //
    // 30 September 2003
    //	o CVol conversion.
    //
    // 27 January 2004
    //  o Added refcnt_ methods
    //

    char*        pch_name = "operator=";
    char*        pch_errormsg = "Volume assignment error. Incompatible dimensions";

    if (!compatible (rval)) {
	rval._warn(pch_name, "rval");
        _error(pch_name, pch_errormsg);
    }

    // clean up current value:
    //if (--ps_vol->refcnt == 0)
    if (refcnt_dec() == 0) {
        coreVolume_destruct();
        //for(int i=0; i<ps_vol->slices; i++)
        //    delete ps_vol->ppzM_volume[i];

        //delete [] ps_vol->ppzM_volume;

        //delete ps_vol;
    }

    // assign base structure to new value:
    //rval.ps_vol->refcnt++;      // tell the rval it has another reference
    rval.refcnt_inc();          // tell the rval it has another reference
    ps_vol = rval.ps_vol;       // point at the rval matrix structure


    return *this;
}

template<typename _CMDATA>
bool
CVol<_CMDATA>::compatible (
    const CVol&                 V)
{
    //
    // ARGS
    //
    //  V                       in      Volume under consideration
    //
    // DESC
    //  Checks if base is compatible with passed arguments
    //
    // 26 March 2003
    //  o Templatization.
    //
    // 30 September 2003
    //  o Expanded to CVol class.
    //

    int         dim1, dim2, dim3;

    dim1        = V.rows_get();
    dim2        = V.cols_get();
    dim3	= V.slices_get();
    if (rows_get() != dim1 || cols_get() != dim2 || slices_get() != dim3)
        return false;
    return true;
}

template<typename _CMDATA>
bool
CVol<_CMDATA>::compatible (
    const int   rows,
    const int   cols,
    const int   slices)
{
    //
    // ARGS
    //
    //  rows, cols, slices     in      Volume under consideration
    //
    // DESC
    //  Checks if base is compatible with passed arguments
    //
    // 26 March 2003
    //  o Templatization.
    //
    // 30 September 2003
    //  o Expanded to CVol class.
    //

    if (rows_get() != rows || cols_get() != cols || slices_get() != slices)
        return false;
    return true;
}

template<typename _CMDATA>
_CMDATA &
CVol<_CMDATA>::_mval(
        int             row,
        int             col,
        int             slice)  const {
    //
    // ARGS
    //        row                in                row to access
    //        col                in                col to access
    //        slice              in                slice to access
    //
    // DESC
    //         Element selection: can be used to read or write.
    //
    // PRECONDITIONS
    //         This routine does *not* use range checking! It is an
    //        interal access routine that is used by class methods
    //        and is not available for "public" use.
    //
    // HISTORY
    // 30 September 2003
    //  o Adapted from CMatrix code.
    //

    return(ps_vol->ppzM_volume[slice]->val(row, col));
}

//////----------><----------\\\\\\
//////      CVol4D class    \\\\\\
//////----------><----------\\\\\\

template<typename _CMDATA>
int CVol4D<_CMDATA>::createdCounter     = 0;

template<typename _CMDATA>
void
CVol4D<_CMDATA>::_warn (
        char *apch_proc,
        char *apch_msg,
        int   code) const {
    //
    // ARGS
    //  apch_proc               in              (optional) function name
    //                                                  in which error has
    //                                                  occured
    //  apch_msg                in              (optional) error message
    //                                                  text
    //  code                    in              (optional) error code
    //
    // DESC
    //  Main warn reporting method for the CVol class.
    //
    // HISTORY
    // 19 September 2003
    // o Initial design and coding.
    //

    cerr << "\nWARNING encountered.\n";
    cerr << "\tCVol4D object id: "      << mps_vol4D->m_id      << endl;
    cerr << "\tps_vol4D:         "      << mps_vol4D            << endl;
    cerr << "\tppVol3D:          "      << mps_vol4D->mppVol3D  << endl;
    cerr << "\trows:             "      << rows_get()           << endl;
    cerr << "\tcolumns:          "      << cols_get()           << endl;
    cerr << "\tslices:           "      << slices_get()         << endl;
    cerr << "\tdim4:             "      << dim4_get()           << endl;
    cerr << "\tCurrent function: " << "CVol4D::" << apch_proc << "\n";
    cerr << "\t" << apch_msg << "\n";
}

template<typename _CMDATA>
void
CVol4D<_CMDATA>::_error (
        char *apch_proc,
        char *apch_msg,
        int   code) {
    //
    // ARGS
    //  apch_proc               in              (optional) function name
    //                                                  in which error has
    //                                                  occured
    //  apch_msg                in              (optional) error message
    //                                                  text
    //  code                    in              (optional) error code
    //
    // DESC
    //  Main error reporting method for the CVol class.
    //
    // TODO
    //        o Wrap exception throwing around #ifdefs?
    //
    // HISTORY
    // 19 September 2003
    // o Initial design and coding.
    //

    cerr << "\nFatal error encountered.\n";
    cerr << "\tCVol4D object id: "      << mps_vol4D->m_id      << endl;
    cerr << "\tps_vol4D:         "      << mps_vol4D            << endl;
    cerr << "\tppVol3D:          "      << mps_vol4D->mppVol3D  << endl;
    cerr << "\trows:             "      << rows_get()           << endl;
    cerr << "\tcolumns:          "      << cols_get()           << endl;
    cerr << "\tslices:           "      << slices_get()         << endl;
    cerr << "\tdim4:             "      << dim4_get()           << endl;
    cerr << "\tCurrent function: " << "CVol4D::" << apch_proc << "\n";
    cerr << "Throwing an exception to (this) with code " << code << endl;
    throw (this);
}

//////---------->
////// Constructor, destructor, and miscellaneous initialization
////// routines
//////---------->

template<typename _CMDATA>
void
CVol4D<_CMDATA>::coreVolume_construct(
        int                     arows,
        int                     acols,
	int                     aslices,
        int                     adim4
) {
    // ARGS
    //  arows                   in              number of rows in core matrix
    //  acols                   in              number of cols in core matrix
    //  aslices                 in              number of slices in volume
    //  adim4                   in              number of volumes in structure
    //
    // DESCRIPTION
    //        This method is a consolidation of the main initialization code
    //        that is shared across all the *matrix* constructors.
    //
    // PRECONDITIONS
    //        o Should only be called from a constructor method.
    //
    // POSTCONDITIONS
    //        o A gsl-aware ps_mat structure is created.
    //
    // HISTORY
    // 19 March 2003
    //        o Initial design and coding.
    //

    // create the structure
    try {
        mps_vol4D = new vol4Dstruct;
    } catch (bad_alloc xa) {
        _error("base constructor",
                "could not allocate memory for vol4Dstruct (heap exhausted?)");
    }
    // Allocate memory for volume container
    mps_vol4D->mppVol3D = new CVol<_CMDATA>* [adim4];
    for(int i=0; i<adim4; i++) {
        mps_vol4D->mppVol3D[i]  = new CVol<_CMDATA>(    arows,
                                                        acols,
                                                        aslices);
    }

    mps_vol4D->m_dim4   = adim4;

    // set first reference/id to this data
    mps_vol4D->m_refcnt = 1;
    mps_vol4D->m_id     = CVol4D::createdCounter++;
}

template <typename _CMDATA>
void
CVol4D<_CMDATA>::coreVolume_destruct() {
    //
    // DESC
    //  Explicitly destroys the core data stuctures of a volume
    //  object.
    //
    // HISTORY
    // 19 September 2003
    //	o Reza's birthday! I miss you, love. Wish I were there in SA with you!
    //  o Initial design and coding.
    //

    for(int i=0; i<mps_vol4D->m_dim4; i++)
        delete mps_vol4D->mppVol3D[i];

    delete [] mps_vol4D->mppVol3D;

    delete mps_vol4D;

}

template <typename _CMDATA>
CVol4D<_CMDATA>::~CVol4D () {
    //
    // DESC
    //  Destructor.
    //
    // POSTCONDITIONS
    //  Decrements reference count. If counter reaches zero, frees allocated
    //  memory.
    //
    // 27 January 2004
    //  o Added refcnt_ methods
    //

    //if (--ps_vol->refcnt == 0) {
    if (refcnt_dec() == 0) {
	coreVolume_destruct();
    }
}

template <typename _CMDATA>
CVol4D<_CMDATA>::CVol4D(
        int             mrows,
        int             mcols,
	int             mslices,
        int             mdim4) {
    //
    // ARGS
    //  mrows                   in              number of rows in each slice
    //  mcols                   in              number of cols in each slice
    //  mslices                 in              number of slices in volume
    //  mdim4                   in              dim4 size
    //
    // DESC
    //  Volume 4D constructor
    //
    // PRECONDITIONS
    //  o Make sure that mrows and mcols are not negative numbers!
    //
    // HISTORY
    // 19 September 2003
    // o Initial design and coding.
    //

    coreVolume_construct(mrows, mcols, mslices, mdim4);

}

template <typename _CMDATA>
CVol4D<_CMDATA>::CVol4D(
    const CVol4D<_CMDATA> & z)
{
    //
    // DESC
    //  Copy constructor - simple pointer based.
    //
    // POSTCONDITIONS
    //  o NB!! NB!! NB!!
    //    The copy constructor merely increases the reference count of the
    //    "source" data, and directs the target to point to the source.
    //    This is *not* a deepcopy!. Although allowing for fast copies between
    //    matrices, it can potentially suffer from problems relating to scope
    //    local variable variable problems. If used in function arguments
    //    or as returns out of functions, then it is not really a problem.
    //
    // 27 January 2004
    //  o Added refcnt_ methods
    //

    //z.ps_vol->refcnt++;               // adding another reference
    z.refcnt_inc();                     // adding another reference
    mps_vol4D = z.mps_vol4D;            // point to the new matstruct
}

template<typename _CMDATA>
CVol4D<_CMDATA>
CVol4D<_CMDATA>::copy (
        const CVol4D&          source)
{
    // ARGS
    //  source                  in              volume whose contents to copy
    //                                                  into *this
    //
    // DESC
    //  Explicit deepcopy.
    //
    // POSTCONDITIONS.
    //  o The data in source is copied over into *this. Space for this data
    //    is created as in necessary.
    //
    // HISTORY
    // 20 September 2003
    //  o Adaptation from CMatrix code.
    //

    int                row, col, slice, dim4;

    if (!compatible (source)) {
        coreVolume_destruct();
        coreVolume_construct(   source.rows_get(),      source.cols_get(),
                                source.slices_get(),    source.dim4_get());
    }
    for (row = 0; row < rows_get(); row++)
        for (col = 0; col < cols_get(); col++)
            for(slice = 0; slice < slices_get(); slice++)
                for(dim4 = 0; dim4 < dim4_get(); dim4++)
                    _mval(row, col, slice, dim4) = source._mval (row, col, slice, dim4);
    return *this;
}

template <typename _CMDATA>
void
CVol4D<_CMDATA>::saveBinary(
    string		astr_fileName)
{
    //
    // ARGS
    //	astr_fileName		in 		file name to save
    //
    // DESC
    //	Performs a simple binary mode save of volume.
    //
    // HISTORY
    // 20 October 2003
    //	o Initial design and coding.
    //

    ofstream	fout(astr_fileName.c_str(), ios_base::out | ios_base::binary | ios_base::trunc);
    if(!fout.is_open())
	_error("saveBinary", "Could not create file for saving", 1);

    int	rows    = rows_get();
    int	cols    = cols_get();
    int	slices  = slices_get();
    int dim4    = dim4_get();

    fout.write((const char*)&rows,      sizeof(int));
    fout.write((const char*)&cols,      sizeof(int));
    fout.write((const char*)&slices,    sizeof(int));
    fout.write((const char*)&dim4,      sizeof(int));

    for(int l=0; l<dim4; l++)
        for(int k=0; k<slices; k++)
	    for(int i=0; i<rows; i++)
	        for(int j=0; j<cols; j++)
		    fout.write((const char*)&val(i, j, k, l), sizeof(_CMDATA));

    fout.close();
}

template <typename _CMDATA>
CVol4D<_CMDATA>
CVol4D<_CMDATA>::loadBinary(
    string		astr_fileName)
{
    //
    // ARGS
    //	astr_fileName		in 		file name to load
    //
    // DESC
    //	Performs a simple binary mode load of volume.
    //
    // HISTORY
    // 20 October 2003
    //	o Initial design and coding.
    //
    
    ifstream	fin(astr_fileName.c_str(), ios_base::in | ios_base::binary);
    if(!fin.is_open())
	_error("loadBinary", (char*)
	       (string("Could not open target file ") + astr_fileName).c_str(), 1);

    int	rows;
    int	cols;
    int	slices;
    int dim4;

    fin.read((char*)&rows,      sizeof(int));
    fin.read((char*)&cols,      sizeof(int));
    fin.read((char*)&slices,    sizeof(int));
    fin.read((char*)&dim4,      sizeof(int));

    if (!compatible (rows, cols, slices, dim4)) {
        coreVolume_destruct();
        coreVolume_construct(rows, cols, slices, dim4);
    }

    for(int l=0; l<dim4; l++)
        for(int k=0; k<slices; k++)
	    for(int i=0; i<rows; i++)
	        for(int j=0; j<cols; j++)
		    fin.read((char*)&val(i, j, k, l), sizeof(_CMDATA));

    fin.close();
}

template <typename _CMDATA>
void
CVol4D<_CMDATA>::print (
    char*                apch_msg        /*= "ans"       */,
    int                  a_precision     /*= 6           */,
    int                  a_width         /*= 12          */,
    ios::fmtflags        a_userFormat    /*= 0           */,
    int                  a_marginOffset  /*= 0           */,
    char                 ach_left        /*= (char) 0    */,
    int                  a_tabOffset     /*= 0           */,
    bool                 ab_fancy        /*= 1           */) {
    //
    // ARGS
    //  apch_msg                in        message that prepends matrix dump
    //  a_precision             in        the precision of numerical output
    //  a_width                 in        the width of the numerical field
    //  a_userFormat            in        if non-zero, set stream format to arg
    //  a_marginOffset          in        the amount of tab spaces prepending ch_left
    //  ach_left                in        also somewhat primitive. If non-zero, will
    //                                          print (char) left_char as left most character
    //                                          of each new line. Useful for when the matrix is
    //                                          part of a larger output dump
    //  a_tabOffset             in          the amount of tab spaces prepending dump
    //  ab_fancy                in          still rather primitive. If false, will print
    //                                                *only* the numbers, nothing else
    //
    // DESC
    //
    //	Simple wrapper around CMatrix::print method. For each slice of volume,
    //	calls the print() method for that slice.
    //
    // NOTE
    //         The `marginOffset' and `tabOffset' are used in conjunction
    //         with ch_left as follows:
    //
    //                 [<--- marginOffset --->]<ch_left>[<--- tabOffset --->]
    //
    // HISTORY
    // 22 September 2003
    //	o Initial design and coding.
    //

    for(int i=0; i<dim4_get(); i++) {
	cout << apch_msg << " dim4: " << i << endl;
	mps_vol4D->mppVol3D[i]->print(
	    apch_msg        /*= "ans"       */,
	    a_precision     /*= 6           */,
	    a_width         /*= 12          */,
	    a_userFormat    /*= 0           */,
	    a_marginOffset  /*= 0           */,
	    ach_left        /*= (char) 0    */,
	    a_tabOffset     /*= 0           */,
	    ab_fancy        /*= 1           */);
    }

}

template<typename _CMDATA>
_CMDATA &
CVol4D<_CMDATA>::val(
        int             row,
        int             col,
	int             slice,
        int             dim4)
{
    //
    // ARGS
    //  row             in              row to access
    //  col             in              col to access
    //  slice           in              slice to access
    //  dim4            in              dim4 to access
    //
    // DESC
    //         Element selection: can be used to read or write.
    //
    // PRECONDITIONS
    //        o Range checking can be expensive, particularly in tightly
    //          nested iterative loops. By defining CMATRIX_NORANGE range
    //          checking is disabled. This will result in tremendous performance
    //          boost, but at the risk of a single out of bound access
    //          killing the entire process.
    //
    // HISTORY
    // 22 September 2003
    //	o Adaptation from CMatrix.
    //

#ifndef CMATRIX_NORANGE
    char *
        pch_name = "val";
    char *
        pch_errormsg = "Addressing error. Index out of range";

    if(!(dim4>=0 && dim4<dim4_get()))
        _error (pch_name, pch_errormsg);
    else
#endif
    return (mps_vol4D->mppVol3D[dim4]->val(row, col, slice));
}

template<typename _CMDATA>
_CMDATA &
CVol4D<_CMDATA>::operator()(
        int             row,
        int             col,
	int             slice,
        int             dim4)
{
    //
    // ARGS
    //  row             in              row to access
    //  col             in              col to access
    //  slice           in              slice to access
    //  dim4            in              dim4 to access
    //
    // DESC
    //         Element selection: can be used to read or write.
    //
    // PRECONDITIONS
    //        o Range checking can be expensive, particularly in tightly
    //          nested iterative loops. By defining CMATRIX_NORANGE range
    //          checking is disabled. This will result in tremendous performance
    //          boost, but at the risk of a single out of bound access
    //          killing the entire process.
    //
    // HISTORY
    // 22 September 2003
    //	o Adaptation from CMatrix.
    //

#ifndef CMATRIX_NORANGE
    char *
        pch_name = "val";
    char *
        pch_errormsg = "Addressing error. Index out of range";
    if(!(dim4>=0 && dim4<dim4_get()))
        _error (pch_name, pch_errormsg);
    else
#endif
    return (mps_vol4D->mppVol3D[dim4]->val(row, col, slice));
}


//////---------->
////// Access and set internal entire slices
//////---------->
template<typename _CMDATA>
CVol<_CMDATA> &
CVol4D<_CMDATA>::val4(
        int             dim4
) {
    //
    // ARGS
    //        dim4              in                volume (in dim4 direction) to access
    //
    // DESC
    //         volume selection: can be used to read or write.
    //
    // PRECONDITIONS
    //  o Use this method to set/get individual "volumes"
    //
    // HISTORY
    // 04 February 2004
    //  o Adaptation for new system
    //

#ifndef CMATRIX_NORANGE
    char *
        pch_name = "val4";
    char *
        pch_errormsg = "Addressing error. dim4 out of range";
    if (!(dim4>= 0 && dim4 < dim4_get()))
        _error (pch_name, pch_errormsg);
    else
#endif
    return (*mps_vol4D->mppVol3D[dim4]);
}

template<typename _CMDATA>
CVol<_CMDATA> &
CVol4D<_CMDATA>::operator()(
        int             dim4
) {
    //
    // ARGS
    //        dim4              in                volume (in dim4 direction) to access
    //
    // DESC
    //         volume selection: can be used to read or write.
    //
    // PRECONDITIONS
    //  o Use this method to set/get individual "volumes"
    //
    // HISTORY
    // 04 February 2004
    //  o Adaptation for new system
    //

#ifndef CMATRIX_NORANGE
    char *
        pch_name = "operator()(int)";
    char *
        pch_errormsg = "Addressing error. dim4 out of range";
    if (!(dim4>= 0 && dim4 < dim4_get()))
        _error (pch_name, pch_errormsg);
    else
#endif
    return (*mps_vol4D->mppVol3D[dim4]);
}

template<typename _CMDATA>
CVol<_CMDATA> &
CVol4D<_CMDATA>::vol3D(
        int             dim4
) {
    //
    // ARGS
    //        dim4              in                volume (in dim4 direction) to access
    //
    // DESC
    //         volume selection: can be used to read or write.
    //
    // PRECONDITIONS
    //  o Use this method to set/get individual "volumes"
    //
    // HISTORY
    // 04 February 2004
    //  o Adaptation for new system
    //

#ifndef CMATRIX_NORANGE
    char *
        pch_name = "operator()(int)";
    char *
        pch_errormsg = "Addressing error. dim4 out of range";
    if (!(dim4>= 0 && dim4 < dim4_get()))
        _error (pch_name, pch_errormsg);
    else
#endif
    return (*mps_vol4D->mppVol3D[dim4]);
}

///////---------->
////// Operator overloads
//////---------->

template<typename _CMDATA>
CVol4D<_CMDATA>
CVol4D<_CMDATA>::operator = (
        const CVol4D<_CMDATA>& rval) {
    //
    // ARGS
    //  rval            in              right-hand argument to operator
    //
    // DESC
    //  Assignment operator.
    //
    // PRECONDITIONS
    //  o A "pointer" equivalence is performed, i.e. *no* deepcopy. The
    //    pointer to *this matrix's data is simply decremented (and
    //    possibly destroyed) and reset to point to rval's data (which
    //    has its refcount increased.
    // o *this must have the same size as rval
    //
    // HISTORY
    // 04 February 2004
    //  o Vol4D extension
    //

    char*        pch_name = "operator=";
    char*        pch_errormsg = "Volume assignment error. Incompatible dimensions";

    if (!compatible (rval)) {
	rval._warn(pch_name, "rval");
        _error(pch_name, pch_errormsg);
    }

    // clean up current value:
    //if (--ps_vol->refcnt == 0)
    if (refcnt_dec() == 0) {
        coreVolume_destruct();
    }

    // assign base structure to new value:
    //rval.ps_vol->refcnt++;      // tell the rval it has another reference
    rval.refcnt_inc();          // tell the rval it has another reference
    mps_vol4D = rval.mps_vol4D; // point at the rval matrix structure


    return *this;
}

template<typename _CMDATA>
bool
CVol4D<_CMDATA>::compatible (
    const CVol4D<_CMDATA>&      V)
{
    //
    // ARGS
    //
    //  V                       in      Volume under consideration
    //
    // DESC
    //  Checks if base is compatible with passed arguments
    //
    // HISTORY
    // 04 February 2004
    //  o Vol4D extension
    //

    int         dim1, dim2, dim3, dim4;

    dim1        = V.rows_get();
    dim2        = V.cols_get();
    dim3	= V.slices_get();
    dim4        = V.dim4_get();
    if (rows_get() != dim1 || cols_get() != dim2 || slices_get() != dim3 || dim4_get() != dim4)
        return false;
    return true;
}

template<typename _CMDATA>
bool
CVol4D<_CMDATA>::compatible (
    const int   rows,
    const int   cols,
    const int   slices,
    const int   dim4)
{
    //
    // ARGS
    //
    //  rows, cols, slices, dim4        in      Volume under consideration
    //
    // DESC
    //  Checks if base is compatible with passed arguments
    //
    // HISTORY
    // 04 February 2004
    //  o Vol4D extension
    //

    if (rows_get() != rows || cols_get() != cols || slices_get() != slices || dim4_get() != dim4)
        return false;
    return true;
}

template<typename _CMDATA>
_CMDATA &
CVol4D<_CMDATA>::_mval(
        int             row,
        int             col,
        int             slice,
        int             dim4)  const {
    //
    // ARGS
    //  row                     in              row to access
    //  col                     in              col to access
    //  slice                   in              slice to access
    //  dim4                    in              dmi4 to access
    //
    // DESC
    //         Element selection: can be used to read or write.
    //
    // PRECONDITIONS
    //         This routine does *not* use range checking! It is an
    //        interal access routine that is used by class methods
    //        and is not available for "public" use.
    //
    // HISTORY
    // 30 September 2003
    //  o Adapted from CMatrix code.
    //

    return(mps_vol4D->mppVol3D[dim4]->val(row, col, slice));
}

//////----------><----------\\\\\\
//////      CVol5D class    \\\\\\
//////----------><----------\\\\\\

template<typename _CMDATA>
int CVol5D<_CMDATA>::createdCounter     = 0;

template<typename _CMDATA>
void
CVol5D<_CMDATA>::_warn (
        char *apch_proc,
        char *apch_msg,
        int   code) const {
    //
    // ARGS
    //  apch_proc               in              (optional) function name
    //                                                  in which error has
    //                                                  occured
    //  apch_msg                in              (optional) error message
    //                                                  text
    //  code                    in              (optional) error code
    //
    // DESC
    //  Main warn reporting method for the CVol class.
    //
    // HISTORY
    // 05 February 2004
    //  o Initial adaptiation from Vol4D
    //

    cerr << "\nWARNING encountered.\n";
    cerr << "\tCVol5D object id: "      << mps_vol5D->m_id      << endl;
    cerr << "\tps_vol5D:         "      << mps_vol5D            << endl;
    cerr << "\tppVol4D:          "      << mps_vol5D->mppVol4D  << endl;
    cerr << "\trows:             "      << rows_get()           << endl;
    cerr << "\tcolumns:          "      << cols_get()           << endl;
    cerr << "\tslices:           "      << slices_get()         << endl;
    cerr << "\tdim4:             "      << dim4_get()           << endl;
    cerr << "\tdim5:             "      << dim5_get()           << endl;
    cerr << "\tCurrent function: " << "CVol5D::" << apch_proc << "\n";
    cerr << "\t" << apch_msg << "\n";
}

template<typename _CMDATA>
void
CVol5D<_CMDATA>::_error (
        char *apch_proc,
        char *apch_msg,
        int   code) {
    //
    // ARGS
    //  apch_proc               in              (optional) function name
    //                                                  in which error has
    //                                                  occured
    //  apch_msg                in              (optional) error message
    //                                                  text
    //  code                    in              (optional) error code
    //
    // DESC
    //  Main error reporting method for the CVol class.
    //
    // TODO
    //        o Wrap exception throwing around #ifdefs?
    //
    // HISTORY
    // 05 February 2004
    //  o Initial adaptiation from Vol4D
    //

    cerr << "\nFatal error encountered.\n";
    cerr << "\tCVol5D object id: "      << mps_vol5D->m_id      << endl;
    cerr << "\tps_vol5D:         "      << mps_vol5D            << endl;
    cerr << "\tppVol4D:          "      << mps_vol5D->mppVol4D  << endl;
    cerr << "\trows:             "      << rows_get()           << endl;
    cerr << "\tcolumns:          "      << cols_get()           << endl;
    cerr << "\tslices:           "      << slices_get()         << endl;
    cerr << "\tdim4:             "      << dim4_get()           << endl;
    cerr << "\tdim5:             "      << dim5_get()           << endl;
    cerr << "\tCurrent function: " << "CVol5D::" << apch_proc << "\n";
    cerr << "Throwing an exception to (this) with code " << code << endl;
    throw (this);
}

//////---------->
////// Constructor, destructor, and miscellaneous initialization
////// routines
//////---------->

template<typename _CMDATA>
void
CVol5D<_CMDATA>::coreVolume_construct(
        int                     arows,
        int                     acols,
	int                     aslices,
        int                     adim4,
        int                     adim5
) {
    // ARGS
    //  arows                   in              number of rows in core matrix
    //  acols                   in              number of cols in core matrix
    //  aslices                 in              number of slices in volume
    //  adim4                   in              number of volumes in structure
    //  adim5                   in              number of volume vectors in
    //                                                  structure
    //
    // DESCRIPTION
    //        This method is a consolidation of the main initialization code
    //        that is shared across all the *matrix* constructors.
    //
    // PRECONDITIONS
    //        o Should only be called from a constructor method.
    //
    // POSTCONDITIONS
    //        o A gsl-aware ps_mat structure is created.
    //
    // HISTORY
    // 05 February 2004
    //  o Initial adaptiation from Vol4D
    //

    // create the structure
    try {
        mps_vol5D = new vol5Dstruct;
    } catch (bad_alloc xa) {
        _error("base constructor",
                "could not allocate memory for vol4Dstruct (heap exhausted?)");
    }
    // Allocate memory for volume container
    mps_vol5D->mppVol4D = new CVol4D<_CMDATA>* [adim5];
    for(int i=0; i<adim5; i++) {
        mps_vol5D->mppVol4D[i]  = new CVol4D<_CMDATA>(  arows,
                                                        acols,
                                                        aslices,
                                                        adim4);
    }

    mps_vol5D->m_dim5   = adim5;

    // set first reference/id to this data
    mps_vol5D->m_refcnt = 1;
    mps_vol5D->m_id     = CVol5D::createdCounter++;
}

template <typename _CMDATA>
void
CVol5D<_CMDATA>::coreVolume_destruct() {
    //
    // DESC
    //  Explicitly destroys the core data stuctures of a volume
    //  object.
    //
    // HISTORY
    // 05 February 2004
    //  o Initial adaptiation from Vol4D
    //

    for(int i=0; i<mps_vol5D->m_dim5; i++)
        delete mps_vol5D->mppVol4D[i];

    delete [] mps_vol5D->mppVol4D;

    delete mps_vol5D;

}

template <typename _CMDATA>
CVol5D<_CMDATA>::~CVol5D () {
    //
    // DESC
    //  Destructor.
    //
    // POSTCONDITIONS
    //  Decrements reference count. If counter reaches zero, frees allocated
    //  memory.
    //
    // HISTORY
    // 05 February 2004
    //  o Initial adaptiation from Vol4D
    //

    //if (--ps_vol->refcnt == 0) {
    if (refcnt_dec() == 0) {
	coreVolume_destruct();
    }
}

template <typename _CMDATA>
CVol5D<_CMDATA>::CVol5D(
        int             mrows,
        int             mcols,
	int             mslices,
        int             mdim4,
        int             mdim5) {
    //
    // ARGS
    //  mrows                   in              number of rows in each slice
    //  mcols                   in              number of cols in each slice
    //  mslices                 in              number of slices in volume
    //  mdim4                   in              dim4 size
    //
    // DESC
    //  Volume 4D constructor
    //
    // PRECONDITIONS
    //  o Make sure that mrows and mcols are not negative numbers!
    //
    // HISTORY
    // 05 February 2004
    //  o Initial adaptiation from Vol4D
    //

    coreVolume_construct(mrows, mcols, mslices, mdim4, mdim5);

}

template <typename _CMDATA>
CVol5D<_CMDATA>::CVol5D(
    const CVol5D<_CMDATA> & z)
{
    //
    // DESC
    //  Copy constructor - simple pointer based.
    //
    // POSTCONDITIONS
    //  o NB!! NB!! NB!!
    //    The copy constructor merely increases the reference count of the
    //    "source" data, and directs the target to point to the source.
    //    This is *not* a deepcopy!. Although allowing for fast copies between
    //    matrices, it can potentially suffer from problems relating to scope
    //    local variable variable problems. If used in function arguments
    //    or as returns out of functions, then it is not really a problem.
    //
    // HISTORY
    // 05 February 2004
    //  o Initial adaptiation from Vol4D
    //

    //z.ps_vol->refcnt++;               // adding another reference
    z.refcnt_inc();                     // adding another reference
    mps_vol5D = z.mps_vol5D;            // point to the new matstruct
}

template<typename _CMDATA>
CVol5D<_CMDATA>
CVol5D<_CMDATA>::copy (
        const CVol5D&          source)
{
    // ARGS
    //  source                  in              volume whose contents to copy
    //                                                  into *this
    //
    // DESC
    //  Explicit deepcopy.
    //
    // POSTCONDITIONS.
    //  o The data in source is copied over into *this. Space for this data
    //    is created as in necessary.
    //
    // HISTORY
    // 05 February 2004
    //  o Initial adaptiation from Vol4D
    //

    int                 row, col, slice, dim4, dim5;

    if (!compatible (source)) {
        coreVolume_destruct();
        coreVolume_construct(   source.rows_get(),      source.cols_get(),
                                source.slices_get(),    source.dim4_get(),
                                source.dim5_get());
    }
        
    for (row = 0; row < rows_get(); row++)
        for (col = 0; col < cols_get(); col++)
            for(slice = 0; slice < slices_get(); slice++)
                for(dim4 = 0; dim4 < dim4_get(); dim4++)
                    for(dim5 = 0; dim5 < dim5_get(); dim5++) {
                        _mval(row, col, slice, dim4, dim5) =
                                source._mval (row, col, slice, dim4, dim5);
                     }
    return *this;
}

template <typename _CMDATA>
void
CVol5D<_CMDATA>::saveBinary(
    string		astr_fileName)
{
    //
    // ARGS
    //	astr_fileName		in 		file name to save
    //
    // DESC
    //	Performs a simple binary mode save of volume.
    //
    // HISTORY
    // 05 February 2004
    //  o Initial adaptiation from Vol4D
    //

    ofstream	fout(astr_fileName.c_str(), ios_base::out | ios_base::binary | ios_base::trunc);
    if(!fout.is_open())
	_error("saveBinary", "Could not create file for saving", 1);

    int	rows    = rows_get();
    int	cols    = cols_get();
    int	slices  = slices_get();
    int dim4    = dim4_get();
    int dim5    = dim5_get();

    fout.write((const char*)&rows,      sizeof(int));
    fout.write((const char*)&cols,      sizeof(int));
    fout.write((const char*)&slices,    sizeof(int));
    fout.write((const char*)&dim4,      sizeof(int));
    fout.write((const char*)&dim5,      sizeof(int));

    for(int m=0; m<dim5; m++)
        for(int l=0; l<dim4; l++)
            for(int k=0; k<slices; k++)
	        for(int i=0; i<rows; i++)
	            for(int j=0; j<cols; j++)
		        fout.write((const char*)&val(i, j, k, l, m), sizeof(_CMDATA));

    fout.close();
}

template <typename _CMDATA>
CVol5D<_CMDATA>
CVol5D<_CMDATA>::loadBinary(
    string		astr_fileName)
{
    //
    // ARGS
    //	astr_fileName		in 		file name to load
    //
    // DESC
    //	Performs a simple binary mode load of volume.
    //
    // HISTORY
    // 05 February 2004
    //  o Initial adaptiation from Vol4D
    //

    ifstream	fin(astr_fileName.c_str(), ios_base::in | ios_base::binary);
    if(!fin.is_open())
	_error("loadBinary", (char*)
	       (string("Could not open target file ") + astr_fileName).c_str(), 1);

    int	rows;
    int	cols;
    int	slices;
    int dim4;
    int dim5;

    fin.read((char*)&rows,      sizeof(int));
    fin.read((char*)&cols,      sizeof(int));
    fin.read((char*)&slices,    sizeof(int));
    fin.read((char*)&dim4,      sizeof(int));
    fin.read((char*)&dim5,      sizeof(int));

    if (!compatible (rows, cols, slices, dim4, dim5)) {
        coreVolume_destruct();
        coreVolume_construct(rows, cols, slices, dim4, dim5);
    }

    for(int m=0; m<dim5; m++)
        for(int l=0; l<dim4; l++)
            for(int k=0; k<slices; k++)
	        for(int i=0; i<rows; i++)
	            for(int j=0; j<cols; j++)
		        fin.read((char*)&val(i, j, k, l, m), sizeof(_CMDATA));

    fin.close();
}

template <typename _CMDATA>
void
CVol5D<_CMDATA>::print (
    char*                apch_msg        /*= "ans"       */,
    int                  a_precision     /*= 6           */,
    int                  a_width         /*= 12          */,
    ios::fmtflags        a_userFormat    /*= 0           */,
    int                  a_marginOffset  /*= 0           */,
    char                 ach_left        /*= (char) 0    */,
    int                  a_tabOffset     /*= 0           */,
    bool                 ab_fancy        /*= 1           */) {
    //
    // ARGS
    //  apch_msg                in        message that prepends matrix dump
    //  a_precision             in        the precision of numerical output
    //  a_width                 in        the width of the numerical field
    //  a_userFormat            in        if non-zero, set stream format to arg
    //  a_marginOffset          in        the amount of tab spaces prepending ch_left
    //  ach_left                in        also somewhat primitive. If non-zero, will
    //                                          print (char) left_char as left most character
    //                                          of each new line. Useful for when the matrix is
    //                                          part of a larger output dump
    //  a_tabOffset             in          the amount of tab spaces prepending dump
    //  ab_fancy                in          still rather primitive. If false, will print
    //                                                *only* the numbers, nothing else
    //
    // DESC
    //
    //	Simple wrapper around CMatrix::print method. For each slice of volume,
    //	calls the print() method for that slice.
    //
    // NOTE
    //         The `marginOffset' and `tabOffset' are used in conjunction
    //         with ch_left as follows:
    //
    //                 [<--- marginOffset --->]<ch_left>[<--- tabOffset --->]
    //
    // HISTORY
    // 05 February 2004
    //  o Initial adaptiation from Vol4D
    //

    for(int i=0; i<dim5_get(); i++) {
	cout << apch_msg << " dim5: " << i << endl;
	mps_vol5D->mppVol4D[i]->print(
	    apch_msg        /*= "ans"       */,
	    a_precision     /*= 6           */,
	    a_width         /*= 12          */,
	    a_userFormat    /*= 0           */,
	    a_marginOffset  /*= 0           */,
	    ach_left        /*= (char) 0    */,
	    a_tabOffset     /*= 0           */,
	    ab_fancy        /*= 1           */);
    }

}

template<typename _CMDATA>
_CMDATA &
CVol5D<_CMDATA>::val(
        int             row,
        int             col,
	int             slice,
        int             dim4,
        int             dim5)
{
    //
    // ARGS
    //  row             in              row to access
    //  col             in              col to access
    //  slice           in              slice to access
    //  dim4            in              dim4 to access
    //  dim5            in              dim5 to access
    //
    // DESC
    //         Element selection: can be used to read or write.
    //
    // PRECONDITIONS
    //        o Range checking can be expensive, particularly in tightly
    //          nested iterative loops. By defining CMATRIX_NORANGE range
    //          checking is disabled. This will result in tremendous performance
    //          boost, but at the risk of a single out of bound access
    //          killing the entire process.
    //
    // HISTORY
    // 05 February 2004
    //  o Initial adaptiation from Vol4D
    //

#ifndef CMATRIX_NORANGE
    char *
        pch_name = "val";
    char *
        pch_errormsg = "Addressing error. Index out of range";

    if(!(dim5>=0 && dim5<dim5_get()))
        _error (pch_name, pch_errormsg);
    else
#endif
    return (mps_vol5D->mppVol4D[dim5]->val(row, col, slice, dim4));
}

template<typename _CMDATA>
_CMDATA &
CVol5D<_CMDATA>::operator()(
        int             row,
        int             col,
	int             slice,
        int             dim4,
        int             dim5)
{
    //
    // ARGS
    //  row             in              row to access
    //  col             in              col to access
    //  slice           in              slice to access
    //  dim4            in              dim4 to access
    //  dim5            in              dim5 to access
    //
    // DESC
    //         Element selection: can be used to read or write.
    //
    // PRECONDITIONS
    //        o Range checking can be expensive, particularly in tightly
    //          nested iterative loops. By defining CMATRIX_NORANGE range
    //          checking is disabled. This will result in tremendous performance
    //          boost, but at the risk of a single out of bound access
    //          killing the entire process.
    //
    // HISTORY
    // 05 February 2004
    //  o Initial adaptiation from Vol4D
    //

#ifndef CMATRIX_NORANGE
    char *
        pch_name = "val";
    char *
        pch_errormsg = "Addressing error. Index out of range";
    if(!(dim5>=0 && dim5<dim5_get()))
        _error (pch_name, pch_errormsg);
    else
#endif
    return (mps_vol5D->mppVol4D[dim5]->val(row, col, slice, dim4));
}


//////---------->
////// Access and set internal entire slices
//////---------->
template<typename _CMDATA>
CVol<_CMDATA> &
CVol5D<_CMDATA>::operator()(
        int             dim4,
        int             dim5
) {
    //
    // ARGS
    //        dim4              in                volume (in dim4 direction) to access
    //        dim5              in                volume (in dim5 direction) to access
    //
    // DESC
    //         volume selection: can be used to read or write.
    //
    // PRECONDITIONS
    //  o Use this method to set/get individual "volumes"
    //
    // HISTORY
    // 05 February 2004
    //  o Initial adaptiation from Vol4D
    //

#ifndef CMATRIX_NORANGE
    char *
        pch_name = "operator()(int)";
    char *
        pch_errormsg = "Addressing error. dim4 out of range";
    if (!(dim5>= 0 && dim5 < dim5_get()))
        _error (pch_name, pch_errormsg);
    else
#endif
    return (mps_vol5D->mppVol4D[dim5]->vol3D(dim4));
}

template<typename _CMDATA>
CVol<_CMDATA> &
CVol5D<_CMDATA>::vol3D(
        int             dim4,
        int             dim5
) {
    //
    // ARGS
    //        dim4              in                volume (in dim4 direction) to access
    //        dim5              in                volume (in dim5 direction) to access
    //
    // DESC
    //         volume selection: can be used to read or write.
    //
    // PRECONDITIONS
    //  o Use this method to set/get individual "volumes"
    //
    // HISTORY
    // 05 February 2004
    //  o Initial adaptiation from Vol4D
    //

#ifndef CMATRIX_NORANGE
    char *
        pch_name = "operator()(int)";
    char *
        pch_errormsg = "Addressing error. dim4 out of range";
    if (!(dim5>= 0 && dim5 < dim5_get()))
        _error (pch_name, pch_errormsg);
    else
#endif
    return (mps_vol5D->mppVol4D[dim5]->vol3D(dim4));
}

///////---------->
////// Operator overloads
//////---------->

template<typename _CMDATA>
CVol5D<_CMDATA>
CVol5D<_CMDATA>::operator = (
        const CVol5D<_CMDATA>& rval) {
    //
    // ARGS
    //  rval            in              right-hand argument to operator
    //
    // DESC
    //  Assignment operator.
    //
    // PRECONDITIONS
    //  o A "pointer" equivalence is performed, i.e. *no* deepcopy. The
    //    pointer to *this matrix's data is simply decremented (and
    //    possibly destroyed) and reset to point to rval's data (which
    //    has its refcount increased.
    // o *this must have the same size as rval
    //
    // HISTORY
    // 05 February 2004
    //  o Vol4D extension
    //

    char*        pch_name = "operator=";
    char*        pch_errormsg = "Volume assignment error. Incompatible dimensions";

    if (!compatible (rval)) {
	rval._warn(pch_name, "rval");
        _error(pch_name, pch_errormsg);
    }

    // clean up current value:
    //if (--ps_vol->refcnt == 0)
    if (refcnt_dec() == 0) {
        coreVolume_destruct();
    }

    // assign base structure to new value:
    //rval.ps_vol->refcnt++;      // tell the rval it has another reference
    rval.refcnt_inc();          // tell the rval it has another reference
    mps_vol5D = rval.mps_vol5D; // point at the rval matrix structure

    return *this;
}

template<typename _CMDATA>
bool
CVol5D<_CMDATA>::compatible (
    const CVol5D<_CMDATA>&      V)
{
    //
    // ARGS
    //
    //  V                       in      Volume under consideration
    //
    // DESC
    //  Checks if base is compatible with passed arguments
    //
    // HISTORY
    // 05 February 2004
    //  o Vol5D extension
    //

    int         dim1, dim2, dim3, dim4, dim5;

    dim1        = V.rows_get();
    dim2        = V.cols_get();
    dim3	= V.slices_get();
    dim4        = V.dim4_get();
    dim5        = V.dim5_get();
    if (rows_get()      != dim1 || cols_get()   != dim2 ||
        slices_get()    != dim3 || dim4_get()   != dim4 ||
        dim5_get()      != dim5)
        return false;
    return true;
}

template<typename _CMDATA>
bool
CVol5D<_CMDATA>::compatible (
    const int   rows,
    const int   cols,
    const int   slices,
    const int   dim4,
    const int   dim5)
{
    //
    // ARGS
    //
    //  rows, cols, slices, dim4, dim5  in      Volume under consideration
    //
    // DESC
    //  Checks if base is compatible with passed arguments
    //
    // HISTORY
    // 05 February 2004
    //  o Vol4D extension
    //

    if (rows_get()      != rows         ||      cols_get()      != cols ||
        slices_get()    != slices       ||      dim4_get()      != dim4 ||
        dim5_get()      != dim5)
        return false;
    return true;
}

template<typename _CMDATA>
_CMDATA &
CVol5D<_CMDATA>::_mval(
        int             row,
        int             col,
        int             slice,
        int             dim4,
        int             dim5)  const {
    //
    // ARGS
    //  row                     in              row to access
    //  col                     in              col to access
    //  slice                   in              slice to access
    //  dim4                    in              dim4 to access
    //  dim5                    in              dim5 to access
    //
    // DESC
    //         Element selection: can be used to read or write.
    //
    // PRECONDITIONS
    //         This routine does *not* use range checking! It is an
    //        interal access routine that is used by class methods
    //        and is not available for "public" use.
    //
    // HISTORY
    // 05 February 2004
    //  o Initial adaptiation from Vol4D
    //

    return(mps_vol5D->mppVol4D[dim5]->val(row, col, slice, dim4));
}


// Explicit defines for particular data types
//	Note that GSL interfaces are only defined for
//	a subset of possible data types, most notably:
//	<complex>{_float} double
//
template class CMatrix<GSL_complex>;
template class CMatrix<GSL_complex_float>;
template class CMatrix<double>;
template class CMatrix<float>;
template class CMatrix<int>;

template class CVol<GSL_complex>;
template class CVol<GSL_complex_float>;
template class CVol<double>;
template class CVol<float>;
template class CVol<int>;

template class CVol4D<GSL_complex_float>;
template class CVol4D<int>;

template class CVol5D<GSL_complex_float>;
template class CVol5D<int>;

//template CMatrix< complex<double> >;

//typedef complex<double> complex_double;
//template CMatrix<complex_double>;

