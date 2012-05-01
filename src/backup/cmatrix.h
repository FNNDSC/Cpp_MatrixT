/***************************************************************************
                          cmatrix.h  -  description
                             -------------------
    begin                : Thu Feb 10 2000
    copyright            : (C) 2000 by Rudolph Pienaar
    email                : pienaar@bme.ri.ccf.org
                           rudolph@nmr.mgh.harvard.edu
 ***************************************************************************/
//
//                                  NOTE
// Portions of this library are based on source code orginally written by
// Bruce Eckel. These portions are copyright either Osborne/McGraw-Hill or
// copyright Bruce Eckel. See http://www.mindview.net
//
// Additional contributors include:
// Antonio      (alrl1@alu.um.es) - 2002
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
//  cmatrix.h       $Id: cmatrix.h,v 1.2 2004/06/22 18:30:55 rudolph Exp $
//
// DESCRIPTION
//
//  `cmatrix.h' contains the class declarations for the associated
//  source cmatrix.cpp.
//
// HISTORY
// o see `README' for detailed history information.
//

#ifndef __CMATRIX_H__
#define __CMATRIX_H__

//#define CMATRIX_NORANGE 1

#include <config.h>

#define GSL_USE
#ifdef GSL_USE
#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_matrix_float>
//#include <gsl/gsl_matrix_complex_double>
//#include <gsl/gsl_matrix_complex_float>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics.h>

#include <gsl/gsl_complex.h>
//#include <gsl/gsl_complex_float.h>
#include <gsl/gsl_complex_math.h>

#include <gsl/gsl_linalg.h>

#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_fft_complex_float.h>

typedef gsl_matrix      gsl_matrix_double;

typedef struct gslstruct {
    gsl_matrix_int*             matrix_int;             // int matrix
    gsl_matrix_float*           matrix_float;           // float matrix
    gsl_matrix_double*          matrix_double;          // double (default) matrix
    gsl_matrix_complex*         matrix_complex;         // double complex matrix
    gsl_matrix_complex_float*   matrix_complex_float;   // float complex matrix
} gslwrap;

struct GSL_complex : gsl_complex {

    GSL_complex(double v) {
        dat[0] = v;
        dat[1] = 0.0;
    };

    GSL_complex(double av_real, double av_imag) {
        dat[0] = av_real;
        dat[1] = av_imag;
    };

    GSL_complex() {
        dat[0] = 0.0;
        dat[1] = 0.0;
    };

    GSL_complex	operator=(const double&		rval) {
        dat[0]	= rval;
        dat[1]	= 0.0;
    };

    GSL_complex	operator=(const GSL_complex&	rval) {
        dat[0]	= GSL_REAL(rval);
        dat[1]	= GSL_IMAG(rval);
    };

    GSL_complex operator-() {
        GSL_complex gz_neg;
        gz_neg.dat[0]   = -dat[0];
        gz_neg.dat[1]   = -dat[1];
        return gz_neg;
    };

    GSL_complex	inverse() {
        double  v_realsq = dat[0] * dat[0];
        double  v_imagsq = dat[1] * dat[1];

        dat[0]  =  dat[0] / (v_realsq + v_imagsq);
        dat[1]  = -dat[1] / (v_realsq + v_imagsq);
        return *this;
    };

    friend std::ostream& operator<< (std::ostream& o, const GSL_complex& rval);

};

struct GSL_complex_float : gsl_complex_float {

    GSL_complex_float(float f) {
        dat[0] = f;
        dat[1] = 0.0;
    };

    GSL_complex_float(double af_real, double af_imag) {
        dat[0] = af_real;
        dat[1] = af_imag;
    };

    GSL_complex_float() {
        dat[0] = 0.0;
        dat[1] = 0.0;
    };

    GSL_complex_float	operator=(const float&		rval) {
        dat[0]	= rval;
        dat[1]	= 0.0;
    };

    GSL_complex_float	operator=(const GSL_complex_float&	rval) {
        dat[0]	= GSL_REAL(rval);
        dat[1]	= GSL_IMAG(rval);
    };

    GSL_complex_float operator-() {
        GSL_complex_float gz_neg;
        gz_neg.dat[0]   = -dat[0];
        gz_neg.dat[1]   = -dat[1];
        return gz_neg;
    };

    GSL_complex_float	inverse() {
        float  f_realsq = dat[0] * dat[0];
        float  f_imagsq = dat[1] * dat[1];

        dat[0]  =  dat[0] / (f_realsq + f_imagsq);
        dat[1]  = -dat[1] / (f_realsq + f_imagsq);
        return *this;
    };

    friend std::ostream& operator<< (std::ostream& o, const GSL_complex_float& rval);

};


#endif // GSL_USE

#include <string>
#include <typeinfo>

typedef enum {
    e_matrix, e_rowVector, e_columnVector, e_square
} e_MATRIXTYPE;

typedef enum {
    e_row, e_column, e_slice
} e_DOMINANCE;
// The e_slice enum is meaningful only in the context of volumes of
//  matrices

typedef enum {
    e_MAGNITUDE, e_MAXIMUM, e_SUM, e_MEANSTD
} e_NORMALISATIONTYPE;

typedef enum {
    e_ascending, e_descending
} e_SORTTYPE;

typedef enum {
        e_forward, e_backward, e_inverse, e_forward2, e_inverse2
} e_FFTDIR;

// FFT direction enums:
// For the MKL:
//      e_forward:      isign = -1
//      e_backward:     isign = -1 (e_backward is a non-scaled GSL parameter)
//      e_inverse       isign = +1
//      e_forward2      isign = -2 (input normal order, output bit-reversed)
//      e_inverse2      isign = +2 (input bit-reversed, output normal order)

typedef enum {
        e_ifftshift,
        e_fftshift
} e_FFTSHIFTDIR;

//////----------------><----------------\\\\\\
//////      CMatrix<_CMDATA> class      \\\\\\
//////----------------><----------------\\\\\\

template<typename _CMDATA>
class CMatrix {

    // Data members - declared as protected so that derived classes can
    // also access information

 protected:
    static int createdCounter;      // Class static member keeps track
                                    //  of absolute number of created
                                    //  matrix objects

    struct matstruct {
        _CMDATA     **data;         // pointer to actual matrix data
        int         rows, cols;     // rows and columns in data structure
        int         refcnt;         // number of times a given instance
                                    //  has been referenced. Used in
                                    //  copy and operator= calls
        int         id;             // internal id number of matrix

#ifdef GSL_USE
        gslwrap*    ps_gsl;         // pointer to gsl wrapper
#endif

    } *ps_mat;

    // Public member functions

 public:

//////---------->
////// Constructor, destructor, and miscellaneous initialization
////// routines
//////---------->

    int         refcnt_inc()    const
                        {return ++(ps_mat->refcnt);};
    int         refcnt_dec()    const
                        {return --(ps_mat->refcnt);};

    void    coreMatrix_construct(
                int             arows,
                int             acols
            );
    CMatrix(    int             mrows           = 1,
                int             mcols           = 1,
                _CMDATA         initvalue       = (_CMDATA) 0.0
            );
    CMatrix(    int             mrows,
                int             mcols,
                const _CMDATA*  data
            );
    CMatrix(    int                 mrows,
                int                 mcols,
                CMatrix<_CMDATA>&   V,
                e_DOMINANCE         dominance
            );
    CMatrix(    CMatrix<_CMDATA>&       M,
                e_DOMINANCE             dominance,
                e_MATRIXTYPE            type
            );
    CMatrix(    char*           flag,
                int             dimension
            );
    //  create an Identity matrix (flag == "I")

    CMatrix(    char*           matfile);
    //  create a matrix from the contents of file "matfile"

    CMatrix(    CMatrix<_CMDATA>& x);
    //  copy initializer
    ~CMatrix();
    void 	coreMatrix_destruct();
    //  destructor

    CMatrix<_CMDATA>    reinitialise(	const _CMDATA*  apf_data,
                                        int             a_elements = 0
    );
    //  overwrite internal contents of matrix with the data
    //          in memory block apf_data. The a_elements is
    //          an (optional) but recommended additional check
    //          on size of memory block.

    int                 positionVectorise(CMatrix<_CMDATA>& aM_dim);
    //  given a vector of dimensions in aM_dim, implements a
    //          N->1 dimensional positional mapping.

    void        randomize(      _CMDATA	 lower,
                                _CMDATA	 upper,
                                bool     b_asInt = false
                );
    //  randomly fill a matrix with values between the
    //          lower and upper ranges. If b_asInt is true,
    //          then fill matrix with integers.

    void        randomly_fill(  _CMDATA	af_fillValue,
                                _CMDATA	af_density
                );
    //  fill a matrix with af_fillValue at density af_density.

    int         hardlimit(      _CMDATA   low,
                                _CMDATA   high,
                                _CMDATA   tolow,
                                _CMDATA   tohigh
                );
    //  a hardlimiting filter. Filter any matrix values less than
    //          "low" to "tolow" and higher than "high" to "tohigh"


    void        quantize(       _CMDATA   f_grain);
    //  quantize a vector with the given f_grain.

    bool        quantizeIndex(  CMatrix & M_lookup,
                                CMatrix *&pM_index
                );
    //  quantize a vector into a lookup table and return the
    //          lookup index

    void        quantize_linearlyShifted(
                                        _CMDATA af_start,
                                        _CMDATA af_stop

                );
    //  quantize a vector over the given range in a linear manner

    void        quantize_linearly(      _CMDATA af_start,
                                        _CMDATA af_stop

                );
    //  quantize a vector over the given range in a linear manner

    void	delta_add(      _CMDATA   af_delta,
                                _CMDATA   af_lower,
                                _CMDATA   af_upper);
    //  adds af_delta to matrix cells that fall between the ranges defined
    //          by af_lower to af_upper


//////---------->
////// Access and set internal values
//////---------->
#ifdef GSL_USE
    gslwrap*    ps_gsl_get()
                    const {return ps_mat->ps_gsl;};
#endif
    int         rows_get()
                                const {return ps_mat->rows;};
    int         cols_get()
                                const {return ps_mat->cols;};
    _CMDATA**   data_get()
                                const {return ps_mat->data;};

    _CMDATA*    dataAsArray_get(    _CMDATA*    apf_CMDATA,
                                    int         a_size      = 0);

    //  element selection: these methods are used to either read values from
    //          or write values to the matrix
    _CMDATA &	val(            int         row,
                                int         col);
    _CMDATA &	operator() (    int         row,
                                int         col);
    _CMDATA &	val(            int         element);
    _CMDATA &	operator() (    int         element);

//////---------->
////// Output routines
//////---------->

    void    print(  char*           apch_msg        = "ans",
                    int             a_precision     = 6,
                    int             a_width         = 12,
                    ios::fmtflags   a_userFormat    = ios::internal,
                    int             a_marginOffset  = 0,
                    char            ach_left        = (char) 0,
                    int             a_tabOffset     = 0,
                    bool            ab_fancy        = 1);
    //  print a matrix to stdout

    string  sprint (
                    char*           apch_msg        = "ans",
                    int             a_precision     = 6,
                    int             a_width         = 12,
                    ios::fmtflags   a_userFormat    = ios::internal,
                    int             a_marginOffset  = 0,
                    char            ach_left        = (char) 0,
                    int             a_tabOffset     = 0,
                    bool            ab_fancy        = 1);
    //  print a matrix to a string    V_RO	= V_SS.cros


    void    fprint(
                    string          astr_filename,
                    char*           apch_msg        = "",
                    int             a_precision     = 6,
                    int             a_width         = 12,
                    ios::fmtflags   a_userFormat    = ios::internal
    );


//////---------->
////// Matrix information routines
//////---------->
    bool            is_outOfBounds( CMatrix<_CMDATA>&	aM_boundaryUpper,
                                    CMatrix<_CMDATA>&	aM_boundaryLower,
                                    bool&               ab_violate
                    );
    //  checks if matrix exceeds upper or lower boundaries

    //  checks if matrix exceeds upper or lower boundaries

    bool            is_vector()         const;

    bool            is_vectorColumn()   const;

    bool            is_vectorRow()      const;

    bool            is_square()         const;

    CMatrix<int>    size()		const;
    //  returns the matrix size, ie [rows cols] as a vector

    int             size1D()		const;
    //  returns the 1D size of the matrix, i.e. the
    //          rows multiplied by the cols (rows X cols)

    e_MATRIXTYPE    matrix_type();
    //  returns the matrix type (square, rowVector, colVector)

    bool            compatible( int     dim1,
                                int     dim2);
    //  checks if *this is compatible with the passed dimensions

    bool            compatible( const CMatrix& M);
    //  checks if *this is compatible with passed matrix

//////---------->
////// Matrix structural manipulation routines
//////---------->
    CMatrix<_CMDATA>	copy(   const   CMatrix<_CMDATA>&	source);
    //  explicit deepcopy

    CMatrix<_CMDATA>    matrix_remove(  CMatrix<_CMDATA>*&  apM,
                                        int                 a_trow,
                                        int                 a_tcol,
                                        int                 a_rows,
                                        int                 a_cols);
    CMatrix<_CMDATA>    matrix_replace( int                 a_trow,
                                        int                 a_tcol,
                                        CMatrix<_CMDATA>&   aM_replacement);
    CMatrix<_CMDATA>    row_remove(     CMatrix<_CMDATA>*&  apV,
                                        int                 a_rownum);
    CMatrix<_CMDATA>    row_insert(     CMatrix<_CMDATA>&	aV,
                                        int                 a_rownum);
    CMatrix<_CMDATA>    row_replace(    CMatrix<_CMDATA>&	aV,
                                        int                 a_rownum);
    CMatrix<_CMDATA>    col_remove(     CMatrix<_CMDATA>*&  apV,
                                        int                 colnum);
    CMatrix<_CMDATA>    col_insert(     CMatrix<_CMDATA>&   aV,
                                        int                 a_colnum);
    CMatrix<_CMDATA>    col_replace(    CMatrix<_CMDATA>&   aV,
                                        int                 tcol);
    CMatrix<_CMDATA>    diagonal_replace(CMatrix<_CMDATA>&	v,
                                        bool                ab_dominant = true);
    //  replace the diagonal of *this with `v'. Defaults to dominant
    //          diagonal (i.e. top left to bottom right)

    CMatrix<_CMDATA>	diagonalise(	_CMDATA         af_offdiag      = 0.0);
    //  creates an NxN matrix with diagonal (N=*this dimension) with
    //          *this as diagonal and off_diagonal value `f_offdiag'

    CMatrix<_CMDATA>    replicate(      int             times,
                                        e_DOMINANCE     dir);
    CMatrix<_CMDATA>    sort(           e_SORTTYPE      e_sorttype  = e_ascending,
                                        e_DOMINANCE     e_dominance = e_row);

    CMatrix<_CMDATA>    zeroPad(        e_DOMINANCE         dir,
                                        int                 size);

    CMatrix<_CMDATA>	shiftIndex(	CMatrix<_CMDATA>&   aM_index,
					int		    a_pivotOffset);
    CMatrix<_CMDATA>	ifftshiftIndex( CMatrix<_CMDATA>&   aM_index);
    CMatrix<_CMDATA>	fftshiftIndex(  CMatrix<_CMDATA>&   aM_index);
    
    CMatrix<_CMDATA>    ifftshift(      bool            ab_inPlace = true);
    CMatrix<_CMDATA>    fftshift(       bool            ab_inPlace = true);
    CMatrix<_CMDATA>    ifftshiftNewMem();
    CMatrix<_CMDATA>    fftshiftNewMem();
    CMatrix<_CMDATA>    ifftshiftInPlace();         // Simple front ends to
    CMatrix<_CMDATA>    fftshiftInPlace();          //  shiftInPlace() method
    CMatrix<_CMDATA>    shiftNewMem(    int a_pivotOffset);
    CMatrix<_CMDATA>    shiftInPlace(   int a_pivotOffset);

//////---------->
////// Search/replace/fill routines
//////---------->
    CMatrix findAll(        _CMDATA     what,
                            int&        occurences);
    //  searches for all instances of `what' and returns their locations
    //          in a vector

    int         findAndReplace( _CMDATA     target,
                            _CMDATA     source);
    //  find all instances of `target' and replace with `source'

    bool        find_quantized( _CMDATA             af_what,
                            CMatrix<_CMDATA>*   pM_index);
    //  find the position corresponding to the quantum of `af_what'

    void	perimeter_set(  _CMDATA     af_value,
                            int         offset          = 0);
    //  sets the perimeter of the matrix to af_value


//////---------->
////// Mathematical considerations
//////---------->

    CMatrix<_CMDATA>        functionApply(  double (*func)(double));
    //  apply a function to each element

//////---------->
////// Simple statistics
//////---------->

    _CMDATA             mean();
    //  mean of *all* elements

    _CMDATA             mean(   e_DOMINANCE     ae_type,
                                int             a_index);
    //  mean of a particular row or column

    CMatrix             mean(   e_DOMINANCE     ae_type);
    //  mean value on a per row/col basis

    _CMDATA             std(    bool            ab_sample = false);
    //  std dev of *all* elements (with/without sampling)

    _CMDATA             std(    e_DOMINANCE     ae_type,
                                int             a_index,
                                bool            ab_sample = false);
    //  std dev of particular row/column (with/without sampling)

    CMatrix<_CMDATA>    std(    e_DOMINANCE     ae_type,
                                bool            ab_sample = false);
    //  std dev value on a per row/col basis

    _CMDATA		sum(	bool            ab_abs = false);
    //  finds the sum of all elements

    _CMDATA		sum(	e_DOMINANCE	ae_type,
    				int		a_index,
    				bool		ab_abs = false);
    //  sum of a particular row or column

    CMatrix<_CMDATA>	sum(	e_DOMINANCE	ae_type,
    				bool		ab_abs = false);
    //  sum value on a per row/col basis

    _CMDATA		min();
    //  find the minimum

    _CMDATA		min(	e_DOMINANCE	ae_type,
    				int		a_index);
    //  min of a particular row or column

    CMatrix<_CMDATA>	min(	e_DOMINANCE	ae_type);
    //  min value on a per row/col basis

    _CMDATA		min_filter(     CMatrix<_CMDATA>& 	M_mask);
    //  find the min (filtering out M_mask)

    _CMDATA		min_epsilon(    _CMDATA   af_center,
                        	        _CMDATA   af_range);
    //  find the epsilon min

    _CMDATA		min_find(       CMatrix<_CMDATA>& 	where);
    //  find the location of the minimum

    _CMDATA		minabs();
    //  find the absolute minimum

    _CMDATA		max();
    //  find the maximum

    _CMDATA		max(	e_DOMINANCE	ae_type,
    				int		a_index);
    //  max of a particular row or column

    CMatrix<_CMDATA>	max(	e_DOMINANCE	ae_type);
    //  max value on a per row/col basis

    _CMDATA		max_filter(     CMatrix<_CMDATA>& 	M_mask);
    //  find the max (filtering out M_mask)

    _CMDATA		max_epsilon(	_CMDATA   af_center,
                                _CMDATA   af_range);
    // find the epsilon max

    _CMDATA		max_find(       CMatrix<_CMDATA>& 		where);
    //  find the location of maximum

    _CMDATA		maxabs();
    //  find the absolute maximum

//////---------->
////// Miscellaneous Maths
//////---------->

    _CMDATA		dot(            const CMatrix<_CMDATA>& 	rval);
    //  dot product

    CMatrix<_CMDATA>	cross(		const CMatrix<_CMDATA>& 	rval);
    //  cross product
        
    _CMDATA		sqrdDistance(	CMatrix<_CMDATA>& 		M_to);
    //  determines the squared distance between *this and M_to

    _CMDATA		distance(	CMatrix<_CMDATA>& 		M_to);
    //  determines the distance between *this and M_to

    _CMDATA		determinant();
    //  find the determinant
    CMatrix		inverse();
    //  find the invese

    CMatrix<_CMDATA>   	normalise(      e_NORMALISATIONTYPE     ae_type = e_MAGNITUDE,
					bool                    ab_inPlace = false);
    //  normalise contents of matrix

    CMatrix<_CMDATA>	normalise(      e_NORMALISATIONTYPE     ae_type,
                                	e_DOMINANCE             ae_dominance);
    //  normalise contents of matrix in eithera row or column dominant manner

    void        	flip(           _CMDATA   af_prob,
                        	        _CMDATA   af_false      = 0.0,
                                	_CMDATA   af_true       = 1.0);
    //  flips (logical inverse) matrix contents

    _CMDATA		mag();
    //  magnitude of matrix

    _CMDATA		mag(		e_DOMINANCE	ae_type,
    					int		a_index);
    //  mag of a particular row or column

    CMatrix<_CMDATA>	mag(		e_DOMINANCE	ae_type);
    //  mag value on a per row/col basis

    CMatrix		abs(		bool		ab_inPlace = false);
    //  absolute value

    _CMDATA		variance();
    //  statistical variance

    CMatrix		lu_decompose(   CMatrix<_CMDATA>&	indx,
                                	int&       		d);
    //  lower/upper decomposition

    void		lu_back_subst(  CMatrix<_CMDATA>&	indx,
                                	CMatrix<_CMDATA>&        b);
    //  uses L-U decomp for matrix inverse

    CMatrix<_CMDATA>     scale(  	const _CMDATA		rval);
    //  scales a matrix

    CMatrix<_CMDATA>     scale(  	const CMatrix<_CMDATA>&	M);
    //  does element by element multiplication


    _CMDATA		innerProd();
    //  finds the inner product (product of all cells)

    _CMDATA		innerSum();
    //  finds the inner sum (sum of all cells) - identical to certain calls
    //          of mag()

    _CMDATA		rms(		CMatrix<_CMDATA>&	rval);
    //  rms between *this and rval

    int			b10_convertTo(	int			radix);
    //  convert a row vector (of given radix) to a base 10 number

    CMatrix<_CMDATA>	b10_convertFrom(int             num,
                                	int             radix,
                                	int             forcelength     = 0);
    //  convert from number in base 10 to number of base radix, with optional length

    bool        	equal(		CMatrix<_CMDATA>&	aM_B);
    //  am I (numerically) equal to B?

    _CMDATA		distanceTo(	CMatrix<_CMDATA>&	aM_B);
    //  distance from *this to B

//////---------->
////// Fourier transforms
//////---------->

    CMatrix<_CMDATA>	    nextPowerOf2();

    CMatrix<_CMDATA>        fft1D(
                                        e_FFTDIR        e_direction     = e_forward,
                                        bool            b_MKLwrap       = false
    );
    CMatrix<_CMDATA>        fft2D(
                                        e_FFTDIR        e_direction     = e_forward,
                                        bool            b_MKLwrap       = false
    );

//////---------->
////// Operator overloads
//////---------->
    CMatrix<_CMDATA>    operator=(      const CMatrix<_CMDATA>&     rval);

    CMatrix<_CMDATA>    operator-();
    //  unary negation

    CMatrix<_CMDATA>    operator*(      const CMatrix<_CMDATA>&     rval);
    CMatrix<_CMDATA>    operator*(      const _CMDATA               rval);
    CMatrix<_CMDATA>    operator*=(     const CMatrix<_CMDATA>&     rval);
    CMatrix<_CMDATA>    operator*=(     const _CMDATA               rval);
    //  element by element multiplication

    CMatrix<_CMDATA>    operator/(      const CMatrix<_CMDATA>&     rval);
    CMatrix<_CMDATA>    operator/(      const _CMDATA               rval);
    CMatrix<_CMDATA>    operator/=(     const CMatrix<_CMDATA>&     rval);
    CMatrix<_CMDATA>    operator/=(     const _CMDATA               rval);
    //  element by element division

    CMatrix<_CMDATA>    operator-(      const CMatrix<_CMDATA>&     rval);
    CMatrix<_CMDATA>    operator-(      const _CMDATA               rval);
    CMatrix<_CMDATA>    operator-=(     const CMatrix<_CMDATA>&     rval);
    CMatrix<_CMDATA>    operator-=(     const _CMDATA               rval);

    CMatrix<_CMDATA>    operator+(      const CMatrix<_CMDATA>&     rval);
    CMatrix<_CMDATA>    operator+(      const _CMDATA               rval);
    CMatrix<_CMDATA>    operator+=(     const CMatrix<_CMDATA>&     rval);
    CMatrix<_CMDATA>    operator+=(     const _CMDATA               rval);

    CMatrix<_CMDATA>    operator!();
    //  transpose a matrix

    CMatrix<_CMDATA>    operator>>=(    const CMatrix<_CMDATA>&     rval);
    //  element by element comparison with rval
    CMatrix<_CMDATA>    operator>>=(    const _CMDATA               rval);
    //  element by element comparison with rval

    CMatrix<_CMDATA>    operator~();
    //  writes the elements of rval into a vector

    _CMDATA             operator++();
    //  inner sum

//////---------->
////// Private member functions
//////---------->
 private:

    void	_error(         char*       apch_proc,
                                char*       apch_msg,
                                int         code            = 1);
    void	_column_switch( int         col1,
                                int         col2);
    void	_column_copy(   CMatrix&        m,
                                int         from_col,
                                int         to_col);
    CMatrix	_scale();
    void	_deepcopy(      CMatrix&        from,
                                CMatrix&        to);
    //  deep copy matrix

    _CMDATA &	_mval(  int        row,
                        int        col) const;
    //  used by matrix functions that know they aren't exceeding the
    //          boundaries
};

//////---------->
////// Non Class functions
//////---------->

template<typename _CMDATA>
CMatrix<_CMDATA>        leastSquaresError_find(
                                CMatrix<_CMDATA>&        aM_A,
                                CMatrix<_CMDATA>&        aM_b);
//      return the least square error between A and b

template<typename _CMDATA>
CMatrix<_CMDATA>	linearSystem_solve(
                                CMatrix<_CMDATA>&        aM_A,
                                CMatrix<_CMDATA>&        aM_b);
//      return x, the solution to the linear system Ax = b,


//////---------->
////// Specializations
//////---------->

//////------------------><------------------\\\\\\
//////      CMatrix<GSL_complex> class      \\\\\\
//////------------------><------------------\\\\\\

template<>
class CMatrix<GSL_complex> {
 protected:
    static int createdCounter;  // Class static member keeps track
                                //  of absolute number of created
                                //  matrix objects

    struct matstruct {
    gsl_complex     **data;         // pointer to actual matrix data
    int             rows, cols;     // rows and columns in data structure
    int             refcnt;         // number of times a given instance
                                    //  has been referenced. Used in
                                    //  copy and operator= calls
    int             id;             // internal id number of matrix

    gslwrap*        ps_gsl;         // pointer to gsl wrapper

    } *ps_mat;
 public:

//////---------->
////// Constructor, destructor, and miscellaneous initialization
////// routines
//////---------->

    int         refcnt_inc()    const
                        {return ++(ps_mat->refcnt);};
    int         refcnt_dec()    const
                        {return --(ps_mat->refcnt);};

    void	coreMatrix_construct(
                int                 arows,
                int                 acols
            );
    CMatrix(    int                 mrows       = 1,
                int                 mcols       = 1,
                GSL_complex         initvalue   = GSL_complex(0.0)
            );
    CMatrix(    int                 mrows,
                int                 mcols,
                const GSL_complex*  data
            );
    CMatrix(    char*               matfile);

    CMatrix(    CMatrix<GSL_complex>& x);
    //  copy initializer
    ~CMatrix();
    void    coreMatrix_destruct();
    //  destructor

//////---------->
////// Output routines
//////---------->

    void    print(      char*           apch_msg        = "ans",
                        int             a_precision     = 6,
                        int             a_width         = 12,
                        ios::fmtflags   a_userFormat	= ios::internal,
                        int             a_marginOffset  = 0,
                        char            ach_left        = (char) 0,
                        int             a_tabOffset     = 0,
                        bool            ab_fancy        = 1);

    string  sprint (
                        char*           apch_msg        = "ans",
                        int             a_precision     = 6,
                        int             a_width         = 12,
                        ios::fmtflags   a_userFormat    = ios::internal,
                        int             a_marginOffset  = 0,
                        char            ach_left        = (char) 0,
                        int             a_tabOffset     = 0,
                        bool            ab_fancy        = 1);
    //  print a matrix to a string

    void	fprint(
                        string          astr_filename,
                        char*           apch_msg        = "",
                        int             a_precision     = 6,
                        int             a_width         = 12,
                        ios::fmtflags   a_userFormat    = ios::internal
    );

//////---------->
////// Access and set internal values
//////---------->
    gsl_complex**       data_get()
                                const {return ps_mat->data;};

    int                 rows_get()
                                const {return ps_mat->rows;};
    int                 cols_get()
                                const {return ps_mat->cols;};
    GSL_complex &	val(            int         row,
                                        int         col);
    GSL_complex &	operator() (    int         row,
                                        int         col);
    GSL_complex &	val(            int         element);
    GSL_complex &	operator() (    int         element);
    
//////---------->
////// Matrix structural manipulation routines
//////---------->
    CMatrix<GSL_complex>    copy(   const   CMatrix<GSL_complex>&	source);
    //  explicit deepcopy
    CMatrix<GSL_complex>    matrix_remove(
                                        CMatrix<GSL_complex>*&  apM,
                                        int                 a_trow,
                                        int                 a_tcol,
                                        int                 a_rows,
                                        int                 a_cols);
    CMatrix<GSL_complex>    matrix_replace(
                                        int                 a_trow,
                                        int                 a_tcol,
                                        CMatrix<GSL_complex>&   aM_replacement);
    CMatrix<GSL_complex>    replicate(  int                 times,
                                        e_DOMINANCE         dir);

    CMatrix<GSL_complex>    zeroPad(    e_DOMINANCE         dir,
                                        int                 size);


//////---------->
////// Matrix information routines
//////---------->
    bool            is_vector()         const;

    bool            is_vectorColumn()   const;

    bool            is_vectorRow()      const;

    bool            is_square()         const;

    CMatrix<int>    size()		const;
    //  returns the matrix size, ie [rows cols] as a vector

    int             size1D()		const;
    //  returns the 1D size of the matrix, i.e. the
    //          rows multiplied by the cols (rows X cols)

    e_MATRIXTYPE    matrix_type();
    //  returns the matrix type (square, rowVector, colVector)

    bool            compatible(         int     dim1,
                                        int     dim2);
    //  checks if *this is compatible with the passed dimensions

    bool            compatible(const CMatrix<GSL_complex>& M);
    //  checks if *this is compatible with passed matrix


//////---------->
////// Miscellaneous Maths
//////---------->

    CMatrix<GSL_complex>        inverse();

//////---------->
////// Fourier transforms
//////---------->

    CMatrix<GSL_complex>        fft1D(
                                        e_FFTDIR        e_direction     = e_forward,
                                        bool            b_MKLwrap       = false
    );
    CMatrix<GSL_complex>        fft2D(
                                        e_FFTDIR        e_direction     = e_forward,
                                        bool            b_MKLwrap       = false
    );

    CMatrix<GSL_complex>        ifftshift(      bool            ab_inPlace = true);
    CMatrix<GSL_complex>        fftshift(       bool            ab_inPlace = true);
    CMatrix<GSL_complex>        ifftshiftNewMem();
    CMatrix<GSL_complex>        fftshiftNewMem();
    CMatrix<GSL_complex>        ifftshiftInPlace();         // Simple front ends to
    CMatrix<GSL_complex>        fftshiftInPlace();          //  shiftInPlace() method
    CMatrix<GSL_complex>        shiftNewMem(    int a_pivotOffset);
    CMatrix<GSL_complex>        shiftInPlace(   int a_pivotOffset);


//////---------->
////// Operator overloads
//////---------->
    CMatrix<GSL_complex>    operator=(      const CMatrix<GSL_complex>&     rval);
    CMatrix<GSL_complex>    operator-();

    CMatrix<GSL_complex>    operator*(      const CMatrix<GSL_complex>&     rval);
    CMatrix<GSL_complex>    operator*(      const GSL_complex&              rval);
    CMatrix<GSL_complex>    operator*=(     const CMatrix<GSL_complex>&     rval);
    CMatrix<GSL_complex>    operator*=(     const GSL_complex&              rval);

    CMatrix<GSL_complex>    operator/(      const CMatrix<GSL_complex>&     rval);
    CMatrix<GSL_complex>    operator/(      const GSL_complex&              rval);
    CMatrix<GSL_complex>    operator/=(     const CMatrix<GSL_complex>&     rval);
    CMatrix<GSL_complex>    operator/=(     const GSL_complex&              rval);
    //  element by element division

    CMatrix<GSL_complex>    operator-(      const CMatrix<GSL_complex>&     rval);
    CMatrix<GSL_complex>    operator-(      const GSL_complex&              rval);
    CMatrix<GSL_complex>    operator-=(     const CMatrix<GSL_complex>&     rval);
    CMatrix<GSL_complex>    operator-=(     const GSL_complex&              rval);

    CMatrix<GSL_complex>    operator+(      const CMatrix<GSL_complex>&     rval);
    CMatrix<GSL_complex>    operator+(      const GSL_complex&              rval);
    CMatrix<GSL_complex>    operator+=(     const CMatrix<GSL_complex>&     rval);
    CMatrix<GSL_complex>    operator+=(     const GSL_complex&              rval);

//////---------->
////// Private member functions
//////---------->
 private:
    GSL_complex &	_mval(  int        row,
                                int        col) const;

    void    _error(         char*       apch_proc,
                            char*       apch_msg,
                            int         code            = 1);

};


//////---------->
////// Non Class CMatrix<complex> functions
//////---------->

//      Perform a 3D FFT
CMatrix<GSL_complex>**
CM_fft3Dl(
        CMatrix<GSL_complex>**  ppzM_input,
        CMatrix<GSL_complex>**& ppzM_output,
        int                     depth,
        e_FFTDIR                e_direction     = e_forward,
        bool                    b_inPlace       = false,
        bool                    b_MKLwrap       = false
);

CMatrix<GSL_complex>**
CM_fft3D(
        CMatrix<GSL_complex>**  ppzM_input,
        CMatrix<GSL_complex>**& ppzM_output,
        int                     depth,
        e_FFTDIR                e_direction     = e_forward,
        bool                    b_inPlace       = false,
        bool                    b_MKLwrap       = false
);


//      Perform a 3D zeroPad operation along a single dimension
CMatrix<GSL_complex>**
CM_vol_zeroPad(
        CMatrix<GSL_complex>**  ppzM_input,
        int                     depth,
        CMatrix<GSL_complex>**& ppzM_output,
        e_DOMINANCE             dir,
        int                     size
);

void    CM_coreVolume_destruct(
                        CMatrix<GSL_complex>**&     ppzM_volume,
                        int                         depth
);

CMatrix<GSL_complex>**
CM_coreVolume_construct(
                        CMatrix<GSL_complex>**&     ppzM_volume,
                        int                         rows,
                        int                         cols,
                        int                         slices
);

void    CM_vol_print(
                        CMatrix<GSL_complex>**      ppzM_volume,
                        int                         depth,
                        char*                       pch_msg
);

CMatrix<GSL_complex>**
CM_vol_shift(
        CMatrix<GSL_complex>**  ppzM_input,
        CMatrix<GSL_complex>**& ppzM_output,
        int                     depth,
        e_FFTSHIFTDIR           e_dir
);

//////------------------><------------------\\\\\\
//////   CMatrix<GSL_complex_float> class   \\\\\\
//////------------------><------------------\\\\\\

template<>
class CMatrix<GSL_complex_float> {
 protected:
    static int createdCounter;  // Class static member keeps track
                                //  of absolute number of created
                                //  matrix objects

    struct matstruct {
    GSL_complex_float       **data;         // pointer to actual matrix data
    int                     rows, cols;     // rows and columns in data structure
    int                     refcnt;         // number of times a given instance
                                            //  has been referenced. Used in
                                            //  copy and operator= calls
    int                     id;             // internal id number of matrix

    gslwrap*                ps_gsl;         // pointer to gsl wrapper

    } *ps_mat;
 public:

//////---------->
////// Constructor, destructor, and miscellaneous initialization
////// routines
//////---------->

    int         refcnt_inc()    const
                        {return ++(ps_mat->refcnt);};
    int         refcnt_dec()    const
                        {return --(ps_mat->refcnt);};

    void	coreMatrix_construct(
                int                 arows,
                int                 acols
            );
    CMatrix(    int                 mrows       = 1,
                int                 mcols       = 1,
                GSL_complex_float         initvalue   = GSL_complex_float(0.0)
            );
    CMatrix(    int                 mrows,
                int                 mcols,
                const GSL_complex_float*  data
            );
    CMatrix(    char*               matfile);

    CMatrix(    CMatrix<GSL_complex_float>& x);
    //  copy initializer
    ~CMatrix();
    void    coreMatrix_destruct();
    //  destructor

//////---------->
////// Output routines
//////---------->

    void    print(      char*           apch_msg        = "ans",
                        int             a_precision     = 6,
                        int             a_width         = 12,
                        ios::fmtflags   a_userFormat	= ios::internal,
                        int             a_marginOffset  = 0,
                        char            ach_left        = (char) 0,
                        int             a_tabOffset     = 0,
                        bool            ab_fancy        = 1);

    string  sprint (
                        char*           apch_msg        = "ans",
                        int             a_precision     = 6,
                        int             a_width         = 12,
                        ios::fmtflags   a_userFormat    = ios::internal,
                        int             a_marginOffset  = 0,
                        char            ach_left        = (char) 0,
                        int             a_tabOffset     = 0,
                        bool            ab_fancy        = 1);
    //  print a matrix to a string

    void	fprint(
                        string          astr_filename,
                        char*           apch_msg        = "",
                        int             a_precision     = 6,
                        int             a_width         = 12,
                        ios::fmtflags   a_userFormat    = ios::internal
    );

//////---------->
////// Access and set internal values
//////---------->
    GSL_complex_float**       data_get()
                                const {return ps_mat->data;};

    int                 rows_get()
                                const {return ps_mat->rows;};
    int                 cols_get()
                                const {return ps_mat->cols;};
    GSL_complex_float &	val(            int         row,
                                        int         col);
    GSL_complex_float &	operator() (    int         row,
                                        int         col);
    GSL_complex_float &	val(            int         element);
    GSL_complex_float &	operator() (    int         element);

//////---------->
////// Matrix structural manipulation routines
//////---------->
    CMatrix<GSL_complex_float>    copy(   const   CMatrix<GSL_complex_float>&	source);
    //  explicit deepcopy
    CMatrix<GSL_complex_float>    matrix_remove(
                                        CMatrix<GSL_complex_float>*&  apM,
                                        int                 a_trow,
                                        int                 a_tcol,
                                        int                 a_rows,
                                        int                 a_cols);
    CMatrix<GSL_complex_float>    matrix_replace(
                                        int                 a_trow,
                                        int                 a_tcol,
                                        CMatrix<GSL_complex_float>&   aM_replacement);
    CMatrix<GSL_complex_float>    replicate(  int                 times,
                                        e_DOMINANCE         dir);

    CMatrix<GSL_complex_float>    zeroPad(    e_DOMINANCE         dir,
                                        int                 size);


//////---------->
////// Matrix information routines
//////---------->
    bool            is_vector()         const;

    bool            is_vectorColumn()   const;

    bool            is_vectorRow()      const;

    bool            is_square()         const;

    CMatrix<int>    size()		const;	
    //  returns the matrix size, ie [rows cols] as a vector

    int             size1D()		const;
    //  returns the 1D size of the matrix, i.e. the
    //          rows multiplied by the cols (rows X cols)

    e_MATRIXTYPE    matrix_type();
    //  returns the matrix type (square, rowVector, colVector)

    bool            compatible(         int     dim1,
                                        int     dim2);
    //  checks if *this is compatible with the passed dimensions

    bool            compatible(const CMatrix<GSL_complex_float>& M);
    //  checks if *this is compatible with passed matrix


//////---------->
////// Miscellaneous Maths
//////---------->

    CMatrix<GSL_complex_float>        inverse();

//////---------->
////// Fourier transforms
//////---------->

    CMatrix<GSL_complex_float>        fft1D(
                                        e_FFTDIR        e_direction     = e_forward,
                                        bool            b_MKLwrap       = false
    );
    CMatrix<GSL_complex_float>        fft2D(
                                        e_FFTDIR        e_direction     = e_forward,
                                        bool            b_MKLwrap       = false
    );

    CMatrix<GSL_complex_float>        ifftshift(      bool            ab_inPlace = true);
    CMatrix<GSL_complex_float>        fftshift(       bool            ab_inPlace = true);
    CMatrix<GSL_complex_float>        ifftshiftNewMem();
    CMatrix<GSL_complex_float>        fftshiftNewMem();
    CMatrix<GSL_complex_float>        ifftshiftInPlace();         // Simple front ends to
    CMatrix<GSL_complex_float>        fftshiftInPlace();          //  shiftInPlace() method
    CMatrix<GSL_complex_float>        shiftNewMem(    int a_pivotOffset);
    CMatrix<GSL_complex_float>        shiftInPlace(   int a_pivotOffset);


//////---------->
////// Operator overloads
//////---------->
    CMatrix<GSL_complex_float>    operator=(      const CMatrix<GSL_complex_float>&     rval);
    CMatrix<GSL_complex_float>    operator-();

    CMatrix<GSL_complex_float>    operator*(      const CMatrix<GSL_complex_float>&     rval);
    CMatrix<GSL_complex_float>    operator*(      const GSL_complex_float&              rval);
    CMatrix<GSL_complex_float>    operator*=(     const CMatrix<GSL_complex_float>&     rval);
    CMatrix<GSL_complex_float>    operator*=(     const GSL_complex_float&              rval);

    CMatrix<GSL_complex_float>    operator/(      const CMatrix<GSL_complex_float>&     rval);
    CMatrix<GSL_complex_float>    operator/(      const GSL_complex_float&              rval);
    CMatrix<GSL_complex_float>    operator/=(     const CMatrix<GSL_complex_float>&     rval);
    CMatrix<GSL_complex_float>    operator/=(     const GSL_complex_float&              rval);
    //  element by element division

    CMatrix<GSL_complex_float>    operator-(      const CMatrix<GSL_complex_float>&     rval);
    CMatrix<GSL_complex_float>    operator-(      const GSL_complex_float&              rval);
    CMatrix<GSL_complex_float>    operator-=(     const CMatrix<GSL_complex_float>&     rval);
    CMatrix<GSL_complex_float>    operator-=(     const GSL_complex_float&              rval);

    CMatrix<GSL_complex_float>    operator+(      const CMatrix<GSL_complex_float>&     rval);
    CMatrix<GSL_complex_float>    operator+(      const GSL_complex_float&              rval);
    CMatrix<GSL_complex_float>    operator+=(     const CMatrix<GSL_complex_float>&     rval);
    CMatrix<GSL_complex_float>    operator+=(     const GSL_complex_float&              rval);

//////---------->
////// Private member functions
//////---------->
 private:
    GSL_complex_float &	_mval(  int        row,
                                int        col) const;

    void    _error(         char*       apch_proc,
                            char*       apch_msg,
                            int         code            = 1);

};


//////----------><----------\\\\\\
//////      CVol class      \\\\\\
//////----------><----------\\\\\\

template<typename _CMDATA>
class CVol {
    // Data members - declared as protected so that derived classes can
    // also access information

 protected:
    static int createdCounter;      // Class static member keeps track
                                    //  of absolute number of created
                                    //  matrix objects
    struct volstruct {
	CMatrix<_CMDATA>**	ppzM_volume;
                                    // pointer to volume structure - this
	                            //	is an array of CMatrix slices.
	int         slices;	    // "depth" dimension of volume. The other
	                            //	two dimensions are implicitly defined
	                            //	in the CMatrix objects comprising
	                            //	individual slices.

	int         refcnt;         // number of times a given instance
                                    //  has been referenced. Used in
                                    //  copy and operator= calls
        int         id;             // internal id number of matrix


    } *ps_vol;

    // Public member functions

 public:

//////---------->
////// Constructor, destructor, and miscellaneous initialization
////// routines
//////---------->

    int         refcnt_inc()    const
        {
            int k = 0;

//            for(k=0; k<ps_vol->slices; k++)
//                ps_vol->ppzM_volume[k]->refcnt_inc();

            return ++(ps_vol->refcnt);
        };

    int         refcnt_dec()    const
        {
            int k = 0;

//            for(k=0; k<ps_vol->slices; k++)
//                ps_vol->ppzM_volume[k]->refcnt_dec();

            return --(ps_vol->refcnt);
        };

    void    coreVolume_construct(
                int             arows,
                int             acols,
		int		aslices
            );
    void    coreVolume_destruct();

    void    coreVolume_replaceWith(	CVol<_CMDATA>*	pVl);

    CVol(       int             mrows           = 1,
                int             mcols           = 1,
		int		mslices		= 1
    );

    CVol(
	CMatrix<_CMDATA>**	ppzM_volume,
	int			slices
    );

    CVol(    CVol<_CMDATA>& x);
    //  copy initializer

    ~CVol();
    void 	coreMatrix_destruct();
    //  destructor


    CVol<_CMDATA>       zeroPad(
                e_DOMINANCE             dir,
                int                     size
    );
    // perform a 3D zeroPad operation along a single dimension

    CVol<_CMDATA>        ifftshift(      bool            ab_inPlace = true);
    CVol<_CMDATA>        fftshift(       bool            ab_inPlace = true);

    CVol<_CMDATA>       shiftNewMem(    int a_pivotOffset);
    CVol<_CMDATA>       shiftInPlace(   int a_pivotOffset);

    CVol<_CMDATA>       shift(
                        e_FFTSHIFTDIR           e_dir
                        );


//////---------->
////// Output routines
//////---------->

    void    print(      char*           apch_msg        = "ans",
                        int             a_precision     = 6,
                        int             a_width         = 12,
                        ios::fmtflags   a_userFormat	= ios::internal,
                        int             a_marginOffset  = 0,
                        char            ach_left        = (char) 0,
                        int             a_tabOffset     = 0,
                        bool            ab_fancy        = 1);

    void    saveBinary(	string		astr_fileName);
    CVol<_CMDATA>
	    loadBinary(	string		astr_fileName);


//////---------->
////// Volume information routines
//////---------->
    bool            compatible(const    CVol<_CMDATA>&  M);
    //  checks if *this is compatible with passed volume

    bool            compatible(const int    rows,
			       const int    cols,
			       const int    slices);
    //  checks if *this is compatible with passed dimensions

//////---------->
////// Access and set internal values
//////---------->

    CVol<_CMDATA>	copy(   const   CVol<_CMDATA>&	source);
    //  explicit deepcopy

    CMatrix<_CMDATA>**  ppzM_volume_get()
                                const {return ps_vol->ppzM_volume;};


    int		slices_get()
	const {return ps_vol->slices;};

    // The rows/cols_get() methods only "query" the first slice
    //	of the volume. Note that the implicit assumption is that
    //	of a "rectangular" space - i.e. each slice has the same
    //	dimensions.
    int         rows_get(int slice=0)
	const {return ps_vol->ppzM_volume[slice]->rows_get();};
    int         cols_get(int slice=0)
	const {return ps_vol->ppzM_volume[slice]->cols_get();};

    //  element selection: these methods are used to either read values from
    //          or write values to the volume
    _CMDATA &	val(                    int         row,
                                        int         col,
				        int         slice);
    _CMDATA &	operator() (            int         row,
                                        int         col,
				        int         slice);
//////---------->
////// Access and set internal entire slices
//////---------->
    CMatrix<_CMDATA>	plane_get(
                                int	                slice,
				e_DOMINANCE	        e_dir	        = e_slice);
    void		plane_set(
                                CMatrix<_CMDATA>&	M_replacement,
				int 		        slice,
				e_DOMINANCE	        e_dir	        = e_slice);
    CMatrix<_CMDATA> &  slice(          int             slice);
    CMatrix<_CMDATA> &  operator() (    int             slice);

//////---------->
////// Fourier transforms
//////---------->
    CVol<_CMDATA>       fft3D(
        e_FFTDIR                e_direction     = e_forward,
        bool                    b_MKLwrap       = false
    );


//////---------->
////// Operator overloads
//////---------->
    CVol<_CMDATA>       operator=(      const CVol<_CMDATA>&     rval);



//////---------->
////// warn /debug functions
//////---------->
    void	_warn(          char*   apch_proc,
                                char*   apch_msg,
                                int     code            = 1) const;

//////---------->
////// Private member functions
//////---------->
 private:

    _CMDATA &	_mval(          int     row,
                                int     col,
                                int     slice) const;
    void	_error(         char*   apch_proc,
                                char*   apch_msg,
                                int     code            = 1);

};

//////----------><----------\\\\\\
//////      CVol4D class    \\\\\\
//////----------><----------\\\\\\

template<typename _CMDATA>
class CVol4D {

 protected:
    static int createdCounter;          // Class static member keeps track
                                        //  of absolute number of created
                                        //  matrix objects
    struct vol4Dstruct {
	CVol<_CMDATA>**
                mppVol3D;               // pointer to array of CVols
	int     m_dim4;	                // the "size" of the 4th volume
                                        //  dimension

	int     m_refcnt;               // number of times a given instance
                                        //  has been referenced. Used in
                                        //  copy and operator= calls
        int     m_id;                   // internal id number of matrix


    } *mps_vol4D;

    // Public member functions

 public:

//////---------->
////// Constructor, destructor, and miscellaneous initialization
////// routines
//////---------->

    int         refcnt_inc()    const
        {
            return ++(mps_vol4D->m_refcnt);
        };

    int         refcnt_dec()    const
        {
            return --(mps_vol4D->m_refcnt);
        };

    void    coreVolume_construct(
                int             arows,
                int             acols,
		int		aslices,
                int             adim4
            );
    void    coreVolume_destruct();

    CVol4D(     int             mrows           = 1,
                int             mcols           = 1,
		int		mslices		= 1,
                int             mdim4           = 1
    );

    CVol4D(     CVol4D<_CMDATA>& x);
    //  copy initializer

    ~CVol4D();
    //  destructor

//////---------->
////// Ouput/Save/load routines
//////---------->

    void    print(      char*           apch_msg        = "ans",
                        int             a_precision     = 6,
                        int             a_width         = 12,
                        ios::fmtflags   a_userFormat	= ios::internal,
                        int             a_marginOffset  = 0,
                        char            ach_left        = (char) 0,
                        int             a_tabOffset     = 0,
                        bool            ab_fancy        = 1);

    void    saveBinary(	string		astr_fileName);
    CVol4D<_CMDATA>
	    loadBinary(	string		astr_fileName);


//////---------->
////// Volume information routines
//////---------->
    bool            compatible(const    CVol4D<_CMDATA>&  V4D);
    //  checks if *this is compatible with passed volume

    bool            compatible(const int        rows,
			       const int        cols,
			       const int        slices,
                               const int        dim4);
    //  checks if *this is compatible with passed dimensions

//////---------->
////// Access and set internal values
//////---------->

    CVol4D<_CMDATA>	copy(   const   CVol4D<_CMDATA>&	source);
    //  explicit deepcopy

    CVol<_CMDATA>**     ppVol3D_get()
                                const {return mps_vol4D->mppVol3D;};


    int		        dim4_get()
	                        const {return mps_vol4D->m_dim4;};

    // The rows/cols_get() methods only "query" the first slice
    //	of the volume. Note that the implicit assumption is that
    //	of a "rectangular" space - i.e. each slice has the same
    //	dimensions.
    int         rows_get(       int slice       = 0,
                                int dim4        = 0)
	const {return mps_vol4D->mppVol3D[dim4]->rows_get(slice);};
    int         cols_get(       int     slice   = 0,
                                int     dim4    = 0)
	const {return mps_vol4D->mppVol3D[dim4]->cols_get(slice);};
    int         slices_get(     int     dim4    = 0)
        const {return mps_vol4D->mppVol3D[dim4]->slices_get();};

    //  element selection: these methods are used to either read values from
    //          or write values to the volume
    _CMDATA &	val(                    int             row,
                                        int             col,
				        int             slice,
                                        int             dim4);
    _CMDATA &	operator() (            int             row,
                                        int             col,
				        int             slice,
                                        int             dim4);
//////---------->
////// Access and set internal entire vol 3Ds
//////---------->

    CVol<_CMDATA> &     val4(           int             dim4);
    CVol<_CMDATA> &  operator() (       int             dim4);
    CVol<_CMDATA> &     vol3D(          int             dim4);

//////---------->
////// Operator overloads
//////---------->
    CVol4D<_CMDATA>     operator=(      const CVol4D<_CMDATA>&  rval);



//////---------->
////// warn /debug functions
//////---------->
    void	_warn(          char*   apch_proc,
                                char*   apch_msg,
                                int     code            = 1) const;

//////---------->
////// Private member functions
//////---------->
 private:

    _CMDATA &	_mval(          int     row,
                                int     col,
                                int     slice,
                                int     dim4) const;
    void	_error(         char*   apch_proc,
                                char*   apch_msg,
                                int     code            = 1);

};

//////----------><----------\\\\\\
//////      CVol5D class    \\\\\\
//////----------><----------\\\\\\

template<typename _CMDATA>
class CVol5D {

 protected:
    static int createdCounter;          // Class static member keeps track
                                        //  of absolute number of created
                                        //  matrix objects
    struct vol5Dstruct {
	CVol4D<_CMDATA>**
                mppVol4D;               // pointer to array of CVols
	int     m_dim5;	                // the "size" of the 4th volume
                                        //  dimension

	int     m_refcnt;               // number of times a given instance
                                        //  has been referenced. Used in
                                        //  copy and operator= calls
        int     m_id;                   // internal id number of matrix


    } *mps_vol5D;

    // Public member functions

 public:

//////---------->
////// Constructor, destructor, and miscellaneous initialization
////// routines
//////---------->

    int         refcnt_inc()    const
        {
            return ++(mps_vol5D->m_refcnt);
        };

    int         refcnt_dec()    const
        {
            return --(mps_vol5D->m_refcnt);
        };

    void    coreVolume_construct(
                int             arows,
                int             acols,
		int		aslices,
                int             adim4,
                int             adim5
            );
    void    coreVolume_destruct();

    CVol5D(     int             mrows           = 1,
                int             mcols           = 1,
		int		mslices		= 1,
                int             mdim4           = 1,
                int             mdim5           = 1
    );

    CVol5D(     CVol5D<_CMDATA>& x);
    //  copy initializer

    ~CVol5D();
    //  destructor

//////---------->
////// Ouput/Save/load routines
//////---------->

    void    print(      char*           apch_msg        = "ans",
                        int             a_precision     = 6,
                        int             a_width         = 12,
                        ios::fmtflags   a_userFormat	= ios::internal,
                        int             a_marginOffset  = 0,
                        char            ach_left        = (char) 0,
                        int             a_tabOffset     = 0,
                        bool            ab_fancy        = 1);

    void    saveBinary(	string		astr_fileName);
    CVol5D<_CMDATA>
	    loadBinary(	string		astr_fileName);


//////---------->
////// Volume information routines
//////---------->
    bool            compatible(const    CVol5D<_CMDATA>&  V5D);
    //  checks if *this is compatible with passed volume

    bool            compatible(const int        rows,
			       const int        cols,
			       const int        slices,
                               const int        dim4,
                               const int        dim5);
    //  checks if *this is compatible with passed dimensions

//////---------->
////// Access and set internal values
//////---------->

    CVol5D<_CMDATA>	copy(   const   CVol5D<_CMDATA>&	source);
    //  explicit deepcopy

    CVol4D<_CMDATA>**   ppVol4D_get()
                                const {return mps_vol5D->mppVol4D;};


    int		        dim5_get()
	                        const {return mps_vol5D->m_dim5;};

    // The rows/cols_get() methods only "query" the first slice
    //	of the volume. Note that the implicit assumption is that
    //	of a "rectangular" space - i.e. each slice has the same
    //	dimensions.
    int         rows_get(       int slice       = 0,
                                int dim4        = 0,
                                int dim5        = 0)
	const {return mps_vol5D->mppVol4D[dim5]->rows_get(slice, dim4);};
    int         cols_get(       int     slice   = 0,
                                int     dim4    = 0,
                                int     dim5    = 0)
	const {return mps_vol5D->mppVol4D[dim5]->cols_get(slice, dim4);};
    int         slices_get(     int     dim4    = 0,
                                int     dim5    = 0)
        const {return mps_vol5D->mppVol4D[dim5]->slices_get(dim4);};
    int         dim4_get(       int     dim5    = 0)
        const {return mps_vol5D->mppVol4D[dim5]->dim4_get();};

    //  element selection: these methods are used to either read values from
    //          or write values to the volume
    _CMDATA &	val(                    int             row,
                                        int             col,
				        int             slice,
                                        int             dim4,
                                        int             dim5);
    _CMDATA &	operator() (            int             row,
                                        int             col,
				        int             slice,
                                        int             dim4,
                                        int             dim5);
//////---------->
////// Access and set internal entire vol 3Ds
//////---------->

    CVol<_CMDATA> &  operator() (       int             dim4,
                                        int             dim5);
    CVol<_CMDATA> &     vol3D(          int             dim4,
                                        int             dim5);

//////---------->
////// Operator overloads
//////---------->
    CVol5D<_CMDATA>     operator=(      const CVol5D<_CMDATA>&  rval);



//////---------->
////// warn /debug functions
//////---------->
    void	_warn(          char*   apch_proc,
                                char*   apch_msg,
                                int     code            = 1) const;

//////---------->
////// Private member functions
//////---------->
 private:

    _CMDATA &	_mval(          int     row,
                                int     col,
                                int     slice,
                                int     dim4,
                                int     dim5) const;
    void	_error(         char*   apch_proc,
                                char*   apch_msg,
                                int     code            = 1);

};

#endif //__CMATRIX_H__






