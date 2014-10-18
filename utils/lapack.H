/*******************************************************************************
 *    FILE: lapack.H
 *    DESC: Function declarations for some LAPACK routines.
 *    DATE: Thu 25 Apr 2013 10:51 AM EDT
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef __lapack_H__INCLUDED__
#define __lapack_H__INCLUDED__

// #include "mkl_scalapack.h"

extern "C" {
    // DLAMCH determines double precision machine parameters.
    double dlamch_ (char* CMACH);

    // DGBEQU computes row and column scalings intended to equilibrate an
    // M-by-N band matrix A and reduce its condition number.  R returns the
    // row scale factors and C the column scale factors, chosen to try to
    // make the largest element in each row and column of the matrix B with
    // elements B(i,j)=R(i)*A(i,j)*C(j) have absolute value 1.
    void dgbequ_(int* M, int* N, int* KL, int* KU, double* AB, int* LDAB,
                 double* R, double* C, double* ROWCND, double* COLCND, double* AMAX, int* INFO);

    // DSYGV computes all the eigenvalues, and optionally, the eigenvectors
    // of a real generalized symmetric-definite eigenproblem, of the form
    // A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.
    // Here A and B are assumed to be symmetric and B is also
    // positive definite.
    void dsygv_(int* ITYPE, char* JOBZ, char* UPLO, int* N,
                double* A, int* LDA, double* B, int* LDB,
                double* W, double* WORK, int* LWORK, int* INFO);

    // DSYGVD computes all the eigenvalues, and optionally, the eigenvectors
    // of a real generalized symmetric-definite eigenproblem, of the form
    // A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and
    // B are assumed to be symmetric and B is also positive definite.
    // If eigenvectors are desired, it uses a divide and conquer algorithm.

    // The divide and conquer algorithm makes very mild assumptions about
    // floating point arithmetic. It will work on machines with a guard
    // digit in add/subtract, or on those binary machines without guard
    // digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
    // Cray-2. It could conceivably fail on hexadecimal or decimal machines
    // without guard digits, but we know of none.
    void dsygvd_(int* ITYPE, char* JOBZ, char* UPLO, int* N,
                 double* A, int* LDA, double* B, int* LDB,
                 double* W, double* WORK, int* LWORK, int* IWORK, int* LIWORK, int* INFO);

    // DSYGVX computes selected eigenvalues, and optionally, eigenvectors
    // of a real generalized symmetric-definite eigenproblem, of the form
    // A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A
    // and B are assumed to be symmetric and B is also positive definite.
    // Eigenvalues and eigenvectors can be selected by specifying either a
    // range of values or a range of indices for the desired eigenvalues.
    void dsygvx_(int* ITYPE, char* JOBZ, char* UPLO, int* N, double* A, int* LDA, double* B, int* LDB,
                 double* VL, double* VU, int* IL, int* IU, double* ABSTOL, int* M, double* W, double* Z,
                 int* LDZ, double* WORK, int* LWORK, int* IWORK, int* IFAIL, int* INFO);

    // DSPGV computes all the eigenvalues and, optionally, the eigenvectors
    // of a real generalized symmetric-definite eigenproblem, of the form
    // A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.
    // Here A and B are assumed to be symmetric, stored in packed format,
    // and B is also positive definite.
    void dspgv_(int* ITYPE, char* JOBZ, char* UPLO, int* N, double* AP, double* BP, double* W, double* Z,
                int* LDZ, double* WORK, int* INFO);

    // DSPGVD computes all the eigenvalues, and optionally, the eigenvectors
    // of a real generalized symmetric-definite eigenproblem, of the form
    // A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and
    // B are assumed to be symmetric, stored in packed format, and B is also
    // positive definite.
    // If eigenvectors are desired, it uses a divide and conquer algorithm.

    // The divide and conquer algorithm makes very mild assumptions about
    // floating point arithmetic. It will work on machines with a guard
    // digit in add/subtract, or on those binary machines without guard
    // digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
    // Cray-2. It could conceivably fail on hexadecimal or decimal machines
    // without guard digits, but we know of none.
    void dspgvd_(int* ITYPE, char* JOBZ, char* UPLO, int* N, double* AP, double* BP, double* W, double* Z,
                 int* LDZ, double* WORK, int* LWORK, int* IWORK, int* LIWORK, int* INFO);

    // DSPGVX computes selected eigenvalues, and optionally, eigenvectors
    // of a real generalized symmetric-definite eigenproblem, of the form
    // A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A
    // and B are assumed to be symmetric, stored in packed storage, and B
    // is also positive definite.  Eigenvalues and eigenvectors can be
    // selected by specifying either a range of values or a range of indices
    // for the desired eigenvalues.
    void dspgvx_(int* ITYPE, char* JOBZ, char* RANGE, char* UPLO, int* N, double* AP, double* BP, double* VL,
                 double* VU, int* IL, int* IU, double* ABSTOL, int* M, double* W, double* Z, int* LDZ,
                 double* WORK, int* IWORK, int* IFAIL, int* INFO);

    // DSBGV computes all the eigenvalues, and optionally, the eigenvectors
    // of a real generalized symmetric-definite banded eigenproblem, of
    // the form A*x=(lambda)*B*x. Here A and B are assumed to be symmetric
    // and banded, and B is also positive definite.
    void dsbgv_(char* JOBZ, char* UPLO, int* N, int* KA, int* KB, double* AB, int* LDAB, double* BB,
                int* LDBB, double* W, double* Z, int* LDZ, double* WORK, int* INFO);

    // DSBGVD computes all the eigenvalues, and optionally, the eigenvectors
    // of a real generalized symmetric-definite banded eigenproblem, of the
    // form A*x=(lambda)*B*x.  Here A and B are assumed to be symmetric and
    // banded, and B is also positive definite.  If eigenvectors are
    // desired, it uses a divide and conquer algorithm.

    // The divide and conquer algorithm makes very mild assumptions about
    // floating point arithmetic. It will work on machines with a guard
    // digit in add/subtract, or on those binary machines without guard
    // digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
    // Cray-2. It could conceivably fail on hexadecimal or decimal machines
    // without guard digits, but we know of none.
    void dsbgvd_(char* JOBZ, char* UPLO, int* N, int* KA, int* KB, double* AB, int* LDAB,
                 double* BB, int* LDBB, double* W, double* Z, int* LDZ,
                 double* WORK, int* LWORK, int* IWORK, int* LIWORK, int* INFO);

    // DSBGVX computes selected eigenvalues, and optionally, eigenvectors
    // of a real generalized symmetric-definite banded eigenproblem, of
    // the form A*x=(lambda)*B*x.  Here A and B are assumed to be symmetric
    // and banded, and B is also positive definite.  Eigenvalues and
    // eigenvectors can be selected by specifying either all eigenvalues,
    // a range of values or a range of indices for the desired eigenvalues.
    void dsbgvx_ (char* JOBZ, char* RANGE, char* UPLO, int* N, int* KA, int* KB, double* AB, int* LDAB, double* BB,
                  int* LDBB, double* Q, int* LDQ, double* VL, double* VU, int* IL, int* IU, double* ABSTOL, int* M,
                  double* W, double* Z, int* LDZ, double* WORK, int* IWORK, int* IFAIL, int* INFO);
}


#endif //!__lapack_H__INCLUDED__
