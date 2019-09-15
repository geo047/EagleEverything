/*
    -- MAGMA (version 2.5.1-alpha1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2019

       @generated from include/magma_zlapack.h, normal z -> d, Fri May 10 14:21:38 2019
*/

#ifndef MAGMA_DLAPACK_H
#define MAGMA_DLAPACK_H

#include "magma_types.h"
#include "magma_mangling.h"

#define MAGMA_REAL

#ifdef __cplusplus
extern "C" {
#endif

/* ////////////////////////////////////////////////////////////////////////////
   -- BLAS and LAPACK functions (alphabetical order)
*/


#define lapackf77_dlacpy   FORTRAN_NAME( dlacpy, DLACPY )


void   lapackf77_dlacpy( const char *uplo,
                         const magma_int_t *m, const magma_int_t *n,
                         const double *A, const magma_int_t *lda,
                         double *B, const magma_int_t *ldb );


/*
 * BLAS functions (alphabetical order)
 */
magma_int_t blasf77_idamax(
                     const magma_int_t *n,
                     const double *x, const magma_int_t *incx );

void blasf77_daxpy(  const magma_int_t *n,
                     const double *alpha,
                     const double *x, const magma_int_t *incx,
                           double *y, const magma_int_t *incy );

void blasf77_dcopy(  const magma_int_t *n,
                     const double *x, const magma_int_t *incx,
                           double *y, const magma_int_t *incy );

void blasf77_dgemm(  const char *transa, const char *transb,
                     const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                     const double *alpha,
                     const double *A, const magma_int_t *lda,
                     const double *B, const magma_int_t *ldb,
                     const double *beta,
                           double *C, const magma_int_t *ldc );

void blasf77_dgemv(  const char *transa,
                     const magma_int_t *m, const magma_int_t *n,
                     const double *alpha,
                     const double *A, const magma_int_t *lda,
                     const double *x, const magma_int_t *incx,
                     const double *beta,
                           double *y, const magma_int_t *incy );

void blasf77_dger(  const magma_int_t *m, const magma_int_t *n,
                     const double *alpha,
                     const double *x, const magma_int_t *incx,
                     const double *y, const magma_int_t *incy,
                           double *A, const magma_int_t *lda );

void blasf77_dger(  const magma_int_t *m, const magma_int_t *n,
                     const double *alpha,
                     const double *x, const magma_int_t *incx,
                     const double *y, const magma_int_t *incy,
                           double *A, const magma_int_t *lda );

void blasf77_dsymm(  const char *side, const char *uplo,
                     const magma_int_t *m, const magma_int_t *n,
                     const double *alpha,
                     const double *A, const magma_int_t *lda,
                     const double *B, const magma_int_t *ldb,
                     const double *beta,
                           double *C, const magma_int_t *ldc );

void blasf77_dsymv(  const char *uplo,
                     const magma_int_t *n,
                     const double *alpha,
                     const double *A, const magma_int_t *lda,
                     const double *x, const magma_int_t *incx,
                     const double *beta,
                           double *y, const magma_int_t *incy );

void blasf77_dsyr(   const char *uplo,
                     const magma_int_t *n,
                     const double *alpha,
                     const double *x, const magma_int_t *incx,
                           double *A, const magma_int_t *lda );

void blasf77_dsyr2(  const char *uplo,
                     const magma_int_t *n,
                     const double *alpha,
                     const double *x, const magma_int_t *incx,
                     const double *y, const magma_int_t *incy,
                           double *A, const magma_int_t *lda );

void blasf77_dsyr2k( const char *uplo, const char *trans,
                     const magma_int_t *n, const magma_int_t *k,
                     const double *alpha,
                     const double *A, const magma_int_t *lda,
                     const double *B, const magma_int_t *ldb,
                     const double *beta,
                           double *C, const magma_int_t *ldc );

void blasf77_dsyrk(  const char *uplo, const char *trans,
                     const magma_int_t *n, const magma_int_t *k,
                     const double *alpha,
                     const double *A, const magma_int_t *lda,
                     const double *beta,
                           double *C, const magma_int_t *ldc );

void blasf77_dscal(  const magma_int_t *n,
                     const double *alpha,
                           double *x, const magma_int_t *incx );

void blasf77_dscal( const magma_int_t *n,
                     const double *alpha,
                           double *x, const magma_int_t *incx );

void blasf77_dswap(  const magma_int_t *n,
                     double *x, const magma_int_t *incx,
                     double *y, const magma_int_t *incy );

/* real-symmetric (non-symmetric) routines */
void blasf77_dsymm(  const char *side, const char *uplo,
                     const magma_int_t *m, const magma_int_t *n,
                     const double *alpha,
                     const double *A, const magma_int_t *lda,
                     const double *B, const magma_int_t *ldb,
                     const double *beta,
                           double *C, const magma_int_t *ldc );

void blasf77_dsyr2k( const char *uplo, const char *trans,
                     const magma_int_t *n, const magma_int_t *k,
                     const double *alpha,
                     const double *A, const magma_int_t *lda,
                     const double *B, const magma_int_t *ldb,
                     const double *beta,
                           double *C, const magma_int_t *ldc );

void blasf77_dsyrk(  const char *uplo, const char *trans,
                     const magma_int_t *n, const magma_int_t *k,
                     const double *alpha,
                     const double *A, const magma_int_t *lda,
                     const double *beta,
                           double *C, const magma_int_t *ldc );

void blasf77_drotg(  double *ca, const double *cb,
                     double *c, double *s );

void blasf77_drot(   const magma_int_t *n,
                     double *x, const magma_int_t *incx,
                     double *y, const magma_int_t *incy,
                     const double *c, const double *s );

void blasf77_drot(  const magma_int_t *n,
                     double *x, const magma_int_t *incx,
                     double *y, const magma_int_t *incy,
                     const double *c, const double *s );

void blasf77_dtrmm(  const char *side, const char *uplo, const char *transa, const char *diag,
                     const magma_int_t *m, const magma_int_t *n,
                     const double *alpha,
                     const double *A, const magma_int_t *lda,
                           double *B, const magma_int_t *ldb );

void blasf77_dtrmv(  const char *uplo, const char *transa, const char *diag,
                     const magma_int_t *n,
                     const double *A, const magma_int_t *lda,
                           double *x, const magma_int_t *incx );

void blasf77_dtrsm(  const char *side, const char *uplo, const char *transa, const char *diag,
                     const magma_int_t *m, const magma_int_t *n,
                     const double *alpha,
                     const double *A, const magma_int_t *lda,
                           double *B, const magma_int_t *ldb );

void blasf77_dtrsv(  const char *uplo, const char *transa, const char *diag,
                     const magma_int_t *n,
                     const double *A, const magma_int_t *lda,
                           double *x, const magma_int_t *incx );

/* ////////////////////////////////////////////////////////////////////////////
 -- MAGMA wrappers around BLAS functions (alphabetical order)
    The Fortran interface for these is not portable, so we
    provide a C interface identical to the Fortran interface.
*/

double magma_cblas_dasum(
    magma_int_t n,
    const double *x, magma_int_t incx );

double magma_cblas_dnrm2(
    magma_int_t n,
    const double *x, magma_int_t incx );

double magma_cblas_ddot(
    magma_int_t n,
    const double *x, magma_int_t incx,
    const double *y, magma_int_t incy );

double magma_cblas_ddot(
    magma_int_t n,
    const double *x, magma_int_t incx,
    const double *y, magma_int_t incy );


/*
 * LAPACK functions (alphabetical order)
 */
#ifdef MAGMA_REAL
void   lapackf77_dbdsdc( const char *uplo, const char *compq,
                         const magma_int_t *n,
                         double *d, double *e,
                         double *U,  const magma_int_t *ldu,
                         double *VT, const magma_int_t *ldvt,
                         double *Q, magma_int_t *IQ,
                         double *work, magma_int_t *iwork,
                         magma_int_t *info );
#endif  // MAGMA_REAL

void   lapackf77_dbdsqr( const char *uplo,
                         const magma_int_t *n, const magma_int_t *ncvt, const magma_int_t *nru,  const magma_int_t *ncc,
                         double *d, double *e,
                         double *Vt, const magma_int_t *ldvt,
                         double *U, const magma_int_t *ldu,
                         double *C, const magma_int_t *ldc,
                         double *work,
                         magma_int_t *info );

void   lapackf77_dgebak( const char *job, const char *side,
                         const magma_int_t *n,
                         const magma_int_t *ilo, const magma_int_t *ihi,
                         const double *scale, const magma_int_t *m,
                         double *V, const magma_int_t *ldv,
                         magma_int_t *info );

void   lapackf77_dgebal( const char *job,
                         const magma_int_t *n,
                         double *A, const magma_int_t *lda,
                         magma_int_t *ilo, magma_int_t *ihi,
                         double *scale,
                         magma_int_t *info );

void   lapackf77_dgebd2( const magma_int_t *m, const magma_int_t *n,
                         double *A, const magma_int_t *lda,
                         double *d, double *e,
                         double *tauq,
                         double *taup,
                         double *work,
                         magma_int_t *info );

void   lapackf77_dgebrd( const magma_int_t *m, const magma_int_t *n,
                         double *A, const magma_int_t *lda,
                         double *d, double *e,
                         double *tauq,
                         double *taup,
                         double *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_dgbbrd( const char *vect, const magma_int_t *m,
                         const magma_int_t *n, const magma_int_t *ncc,
                         const magma_int_t *kl, const magma_int_t *ku,
                         double *Ab, const magma_int_t *ldab,
                         double *d, double *e,
                         double *Q, const magma_int_t *ldq,
                         double *PT, const magma_int_t *ldpt,
                         double *C, const magma_int_t *ldc,
                         double *work,
                         #ifdef MAGMA_COMPLEX
                         double *rwork,
                         #endif
                         magma_int_t *info );

void   lapackf77_dgbsv( const magma_int_t *n, 
                        const magma_int_t *kl, const magma_int_t *ku, 
                        const magma_int_t *nrhs,
                        double *ab, const magma_int_t *ldab, 
                        magma_int_t *ipiv, 
                        double *B, const magma_int_t *ldb, 
                        magma_int_t *info );

void   lapackf77_dgeev(  const char *jobvl, const char *jobvr,
                         const magma_int_t *n,
                         double *A,    const magma_int_t *lda,
                         #ifdef MAGMA_COMPLEX
                         double *w,
                         #else
                         double *wr, double *wi,
                         #endif
                         double *Vl,   const magma_int_t *ldvl,
                         double *Vr,   const magma_int_t *ldvr,
                         double *work, const magma_int_t *lwork,
                         #ifdef MAGMA_COMPLEX
                         double *rwork,
                         #endif
                         magma_int_t *info );

void   lapackf77_dgehd2( const magma_int_t *n,
                         const magma_int_t *ilo, const magma_int_t *ihi,
                         double *A, const magma_int_t *lda,
                         double *tau,
                         double *work,
                         magma_int_t *info );

void   lapackf77_dgehrd( const magma_int_t *n,
                         const magma_int_t *ilo, const magma_int_t *ihi,
                         double *A, const magma_int_t *lda,
                         double *tau,
                         double *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_dgelqf( const magma_int_t *m, const magma_int_t *n,
                         double *A, const magma_int_t *lda,
                         double *tau,
                         double *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_dgels(  const char *trans,
                         const magma_int_t *m, const magma_int_t *n, const magma_int_t *nrhs,
                         double *A, const magma_int_t *lda,
                         double *B, const magma_int_t *ldb,
                         double *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_dgeqlf( const magma_int_t *m, const magma_int_t *n,
                         double *A, const magma_int_t *lda,
                         double *tau,
                         double *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_dgeqp3( const magma_int_t *m, const magma_int_t *n,
                         double *A, const magma_int_t *lda,
                         magma_int_t *jpvt,
                         double *tau,
                         double *work, const magma_int_t *lwork,
                         #ifdef MAGMA_COMPLEX
                         double *rwork,
                         #endif
                         magma_int_t *info );

void   lapackf77_dgeqrf( const magma_int_t *m, const magma_int_t *n,
                         double *A, const magma_int_t *lda,
                         double *tau,
                         double *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_dgerqf( const magma_int_t *m, const magma_int_t *n, 
                         double *A, const magma_int_t *lda,
                         double *tau, 
                         double *work, const magma_int_t *lwork, 
                         magma_int_t *info);

void   lapackf77_dgesdd( const char *jobz,
                         const magma_int_t *m, const magma_int_t *n,
                         double *A, const magma_int_t *lda,
                         double *s,
                         double *U,  const magma_int_t *ldu,
                         double *Vt, const magma_int_t *ldvt,
                         double *work, const magma_int_t *lwork,
                         #ifdef MAGMA_COMPLEX
                         double *rwork,
                         #endif
                         magma_int_t *iwork,
                         magma_int_t *info );

void   lapackf77_dgesv(  const magma_int_t *n, const magma_int_t *nrhs,
                         double *A, const magma_int_t *lda,
                         magma_int_t *ipiv,
                         double *B,  const magma_int_t *ldb,
                         magma_int_t *info );

void   lapackf77_dgesvd( const char *jobu, const char *jobvt,
                         const magma_int_t *m, const magma_int_t *n,
                         double *A, const magma_int_t *lda,
                         double *s,
                         double *U,  const magma_int_t *ldu,
                         double *Vt, const magma_int_t *ldvt,
                         double *work, const magma_int_t *lwork,
                         #ifdef MAGMA_COMPLEX
                         double *rwork,
                         #endif
                         magma_int_t *info );

void   lapackf77_dgetrf( const magma_int_t *m, const magma_int_t *n,
                         double *A, const magma_int_t *lda,
                         magma_int_t *ipiv,
                         magma_int_t *info );

void   lapackf77_dgetri( const magma_int_t *n,
                         double *A, const magma_int_t *lda,
                         const magma_int_t *ipiv,
                         double *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_dgetrs( const char *trans,
                         const magma_int_t *n, const magma_int_t *nrhs,
                         const double *A, const magma_int_t *lda,
                         const magma_int_t *ipiv,
                         double *B, const magma_int_t *ldb,
                         magma_int_t *info );

void   lapackf77_dgglse( magma_int_t *m, magma_int_t *n, magma_int_t *p,
                         double *A, magma_int_t *lda,
                         double *B, magma_int_t *ldb,
                         double *c, double *d, 
                         double *x,
                         double *work, magma_int_t *lwork,
                         magma_int_t *info);

void   lapackf77_dggrqf( magma_int_t *m, magma_int_t *p, magma_int_t *n,
                         double *A, magma_int_t *lda,
                         double *tauA, double *B,
                         magma_int_t *ldb, double *tauB,
                         double *work, magma_int_t *lwork, 
                         magma_int_t *info);

void   lapackf77_dsytf2( const char *uplo, const magma_int_t *n,
                         double *A, const magma_int_t *lda,
                         magma_int_t *ipiv,
                         magma_int_t *info );

void   lapackf77_dsytrs( const char *uplo,
                         const magma_int_t *n, const magma_int_t *nrhs,
                         const double *A, const magma_int_t *lda,
                         const magma_int_t *ipiv,
                         double *B, const magma_int_t *ldb,
                         magma_int_t *info );

void   lapackf77_dsbtrd( const char *vect, const char *uplo,
                         const magma_int_t *n, const magma_int_t *kd,
                         double *Ab, const magma_int_t *ldab,
                         double *d, double *e,
                         double *Q, const magma_int_t *ldq,
                         double *work,
                         magma_int_t *info );

void   lapackf77_dsyev(  const char *jobz, const char *uplo,
                         const magma_int_t *n,
                         double *A, const magma_int_t *lda,
                         double *w,
                         double *work, const magma_int_t *lwork,
                         #ifdef MAGMA_COMPLEX
                         double *rwork,
                         #endif
                         magma_int_t *info );

void   lapackf77_dsyevd( const char *jobz, const char *uplo,
                         const magma_int_t *n,
                         double *A, const magma_int_t *lda,
                         double *w,
                         double *work, const magma_int_t *lwork,
                         #ifdef MAGMA_COMPLEX
                         double *rwork, const magma_int_t *lrwork,
                         #endif
                         magma_int_t *iwork, const magma_int_t *liwork,
                         magma_int_t *info );

void   lapackf77_dsyevr( const char *jobz, const char *range, const char *uplo,
                         const magma_int_t *n,
                         double *A, const magma_int_t *lda,
                         const double *vl, const double *vu,
                         const magma_int_t *il, const magma_int_t *iu,
                         const double *abstol,
                         magma_int_t *m, double *w,
                         double *Z, const magma_int_t *ldz,
                         magma_int_t *isuppz,
                         double *work, const magma_int_t *lwork,
                         #ifdef MAGMA_COMPLEX
                         double *rwork, const magma_int_t *lrwork,
                         #endif
                         magma_int_t *iwork, const magma_int_t *liwork,
                         magma_int_t *info);

void   lapackf77_dsyevx( const char *jobz, const char *range, const char *uplo,
                         const magma_int_t *n,
                         double *A, const magma_int_t *lda,
                         const double *vl, const double *vu,
                         const magma_int_t *il, const magma_int_t *iu,
                         const double *abstol,
                         magma_int_t *m, double *w,
                         double *Z, const magma_int_t *ldz,
                         double *work, const magma_int_t *lwork,
                         #ifdef MAGMA_COMPLEX
                         double *rwork,
                         #endif
                         magma_int_t *iwork, magma_int_t *ifail,
                         magma_int_t *info);

void   lapackf77_dsygs2( const magma_int_t *itype, const char *uplo,
                         const magma_int_t *n,
                         double *A, const magma_int_t *lda,
                         const double *B, const magma_int_t *ldb,
                         magma_int_t *info );

void   lapackf77_dsygst( const magma_int_t *itype, const char *uplo,
                         const magma_int_t *n,
                         double *A, const magma_int_t *lda,
                         const double *B, const magma_int_t *ldb,
                         magma_int_t *info );

void   lapackf77_dsygvd( const magma_int_t *itype, const char *jobz, const char *uplo,
                         const magma_int_t *n,
                         double *A, const magma_int_t *lda,
                         double *B, const magma_int_t *ldb,
                         double *w,
                         double *work, const magma_int_t *lwork,
                         #ifdef MAGMA_COMPLEX
                         double *rwork, const magma_int_t *lrwork,
                         #endif
                         magma_int_t *iwork, const magma_int_t *liwork,
                         magma_int_t *info );

void   lapackf77_dsysv( const char *uplo,
                        const magma_int_t *n, const magma_int_t *nrhs,
                        double *A, const magma_int_t *lda, magma_int_t *ipiv,
                        double *B, const magma_int_t *ldb,
                        double *work, const magma_int_t *lwork,
                        magma_int_t *info );

void   lapackf77_dsytd2( const char *uplo,
                         const magma_int_t *n,
                         double *A, const magma_int_t *lda,
                         double *d, double *e,
                         double *tau,
                         magma_int_t *info );

void   lapackf77_dsytrd( const char *uplo,
                         const magma_int_t *n,
                         double *A, const magma_int_t *lda,
                         double *d, double *e,
                         double *tau,
                         double *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_dsytrf( const char *uplo,
                         const magma_int_t *n,
                         double *A, const magma_int_t *lda,
                         magma_int_t *ipiv,
                         double *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_dhseqr( const char *job, const char *compz,
                         const magma_int_t *n,
                         const magma_int_t *ilo, const magma_int_t *ihi,
                         double *H, const magma_int_t *ldh,
                         #ifdef MAGMA_COMPLEX
                         double *w,
                         #else
                         double *wr, double *wi,
                         #endif
                         double *Z, const magma_int_t *ldz,
                         double *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_dlabrd( const magma_int_t *m, const magma_int_t *n, const magma_int_t *nb,
                         double *A, const magma_int_t *lda,
                         double *d, double *e,
                         double *tauq,
                         double *taup,
                         double *X, const magma_int_t *ldx,
                         double *Y, const magma_int_t *ldy );

#ifdef MAGMA_COMPLEX
void   lapackf77_dlacgv( const magma_int_t *n,
                         double *x, const magma_int_t *incx );
#endif

#ifdef MAGMA_COMPLEX
void   lapackf77_dlacp2( const char *uplo,
                         const magma_int_t *m, const magma_int_t *n,
                         const double *A, const magma_int_t *lda,
                         double *B, const magma_int_t *ldb );
#endif






#ifdef MAGMA_COMPLEX
void   lapackf77_dlacrm( const magma_int_t *m, const magma_int_t *n,
                         const double *A, const magma_int_t *lda,
                         const double             *B, const magma_int_t *ldb,
                         double       *C, const magma_int_t *ldc,
                         double *rwork );
#endif

#ifdef MAGMA_COMPLEX
void   lapackf77_dladiv( double *ret_val,
                         const double *x,
                         const double *y );
#else // MAGMA_REAL
void   lapackf77_dladiv( const double *a, const double *b,
                         const double *c, const double *d,
                         double *p, double *q );
#endif

void   lapackf77_dlasyf( const char *uplo,
                         const magma_int_t *n, const magma_int_t *nb,
                         magma_int_t *kb,
                         double *A, const magma_int_t *lda,
                         magma_int_t *ipiv,
                         double *work, const magma_int_t *ldwork,
                         magma_int_t *info );

double lapackf77_dlange( const char *norm,
                         const magma_int_t *m, const magma_int_t *n,
                         const double *A, const magma_int_t *lda,
                         double *work );

double lapackf77_dlansy( const char *norm, const char *uplo,
                         const magma_int_t *n,
                         const double *A, const magma_int_t *lda,
                         double *work );

double lapackf77_dlanst( const char *norm, const magma_int_t *n,
                         const double *d, const double *e );

double lapackf77_dlansy( const char *norm, const char *uplo,
                         const magma_int_t *n,
                         const double *A, const magma_int_t *lda,
                         double *work );

double lapackf77_dlantr( const char *norm, const char *uplo, const char *diag,
                         const magma_int_t *m, const magma_int_t *n,
                         const double *A, const magma_int_t *lda,
                         double *work );

void   lapackf77_dlaqp2( const magma_int_t *m, const magma_int_t *n, const magma_int_t *offset,
                         double *A, const magma_int_t *lda,
                         magma_int_t *jpvt,
                         double *tau,
                         double *vn1, double *vn2,
                         double *work );

#ifdef MAGMA_COMPLEX
void   lapackf77_dlarcm( const magma_int_t *m, const magma_int_t *n,
                         const double             *A, const magma_int_t *lda,
                         const double *B, const magma_int_t *ldb,
                         double       *C, const magma_int_t *ldc,
                         double *rwork );
#endif

void   lapackf77_dlarf(  const char *side, const magma_int_t *m, const magma_int_t *n,
                         const double *v, const magma_int_t *incv,
                         const double *tau,
                         double *C, const magma_int_t *ldc,
                         double *work );

void   lapackf77_dlarfb( const char *side, const char *trans, const char *direct, const char *storev,
                         const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         const double *V, const magma_int_t *ldv,
                         const double *T, const magma_int_t *ldt,
                         double *C, const magma_int_t *ldc,
                         double *work, const magma_int_t *ldwork );

void   lapackf77_dlarfg( const magma_int_t *n,
                         double *alpha,
                         double *x, const magma_int_t *incx,
                         double *tau );

void   lapackf77_dlarft( const char *direct, const char *storev,
                         const magma_int_t *n, const magma_int_t *k,
                         const double *V, const magma_int_t *ldv,
                         const double *tau,
                         double *T, const magma_int_t *ldt );

void   lapackf77_dlarfx( const char *side, const magma_int_t *m, const magma_int_t *n,
                         const double *V,
                         const double *tau,
                         double *C, const magma_int_t *ldc,
                         double *work );

void   lapackf77_dlarnv( const magma_int_t *idist, magma_int_t *iseed, const magma_int_t *n,
                         double *x );

void   lapackf77_dlartg( const double *f,
                         const double *g,
                         double *cs,
                         double *sn,
                         double *r );

void   lapackf77_dlascl( const char *type,
                         const magma_int_t *kl, const magma_int_t *ku,
                         const double *cfrom,
                         const double *cto,
                         const magma_int_t *m, const magma_int_t *n,
                         double *A, const magma_int_t *lda,
                         magma_int_t *info );

void   lapackf77_dlaset( const char *uplo,
                         const magma_int_t *m, const magma_int_t *n,
                         const double *alpha,
                         const double *beta,
                         double *A, const magma_int_t *lda );

void   lapackf77_dlaswp( const magma_int_t *n,
                         double *A, const magma_int_t *lda,
                         const magma_int_t *k1, const magma_int_t *k2,
                         const magma_int_t *ipiv,
                         const magma_int_t *incx );

void   lapackf77_dlatrd( const char *uplo,
                         const magma_int_t *n, const magma_int_t *nb,
                         double *A, const magma_int_t *lda,
                         double *e,
                         double *tau,
                         double *work, const magma_int_t *ldwork );

void   lapackf77_dlatrs( const char *uplo, const char *trans, const char *diag,
                         const char *normin,
                         const magma_int_t *n,
                         const double *A, const magma_int_t *lda,
                         double *x, double *scale,
                         double *cnorm,
                         magma_int_t *info );

void   lapackf77_dlauum( const char *uplo,
                         const magma_int_t *n,
                         double *A, const magma_int_t *lda,
                         magma_int_t *info );

void   lapackf77_dlavsy( const char *uplo, const char *trans, const char *diag,
                         magma_int_t *n, magma_int_t *nrhs,
                         double *A, magma_int_t *lda,
                         magma_int_t *ipiv,
                         double *B, magma_int_t *ldb,
                         magma_int_t *info );

void   lapackf77_dposv(  const char *uplo,
                         const magma_int_t *n, const magma_int_t *nrhs,
                         double *A, const magma_int_t *lda,
                         double *B,  const magma_int_t *ldb,
                         magma_int_t *info );

void   lapackf77_dpotrf( const char *uplo,
                         const magma_int_t *n,
                         double *A, const magma_int_t *lda,
                         magma_int_t *info );

void   lapackf77_dpotri( const char *uplo,
                         const magma_int_t *n,
                         double *A, const magma_int_t *lda,
                         magma_int_t *info );

void   lapackf77_dpotrs( const char *uplo,
                         const magma_int_t *n, const magma_int_t *nrhs,
                         const double *A, const magma_int_t *lda,
                         double *B, const magma_int_t *ldb,
                         magma_int_t *info );

void   lapackf77_dstedc( const char *compz,
                         const magma_int_t *n,
                         double *d, double *e,
                         double *Z, const magma_int_t *ldz,
                         double *work, const magma_int_t *lwork,
                         #ifdef MAGMA_COMPLEX
                         double *rwork, const magma_int_t *lrwork,
                         #endif
                         magma_int_t *iwork, const magma_int_t *liwork,
                         magma_int_t *info );

void   lapackf77_dstein( const magma_int_t *n,
                         const double *d, const double *e,
                         const magma_int_t *m,
                         const double *w,
                         const magma_int_t *iblock,
                         const magma_int_t *isplit,
                         double *Z, const magma_int_t *ldz,
                         double *work, magma_int_t *iwork, magma_int_t *ifailv,
                         magma_int_t *info );

void   lapackf77_dstemr( const char *jobz, const char *range,
                         const magma_int_t *n,
                         double *d, double *e,
                         const double *vl, const double *vu,
                         const magma_int_t *il, const magma_int_t *iu,
                         magma_int_t *m,
                         double *w,
                         double *Z, const magma_int_t *ldz,
                         const magma_int_t *nzc, magma_int_t *isuppz, magma_int_t *tryrac,
                         double *work, const magma_int_t *lwork,
                         magma_int_t *iwork, const magma_int_t *liwork,
                         magma_int_t *info );

void   lapackf77_dsteqr( const char *compz,
                         const magma_int_t *n,
                         double *d, double *e,
                         double *Z, const magma_int_t *ldz,
                         double *work,
                         magma_int_t *info );

#ifdef MAGMA_COMPLEX
void   lapackf77_dsymv(  const char *uplo,
                         const magma_int_t *n,
                         const double *alpha,
                         const double *A, const magma_int_t *lda,
                         const double *x, const magma_int_t *incx,
                         const double *beta,
                               double *y, const magma_int_t *incy );

void   lapackf77_dsyr(   const char *uplo,
                         const magma_int_t *n,
                         const double *alpha,
                         const double *x, const magma_int_t *incx,
                               double *A, const magma_int_t *lda );

void   lapackf77_dsysv(  const char *uplo,
                         const magma_int_t *n, const magma_int_t *nrhs,
                         double *A, const magma_int_t *lda, magma_int_t *ipiv,
                         double *B, const magma_int_t *ldb,
                         double *work, const magma_int_t *lwork,
                         magma_int_t *info );

#endif  // MAGMA_COMPLEX

void   lapackf77_dtrevc( const char *side, const char *howmny,
                         // select is [in] for real; [in,out] for real
                         #ifdef MAGMA_COMPLEX
                         const
                         #endif
                         magma_int_t *select, const magma_int_t *n,
                         // T is modified but restored in real; const for real
                         #ifdef MAGMA_REAL
                         const
                         #endif
                         double *T,  const magma_int_t *ldt,
                         double *Vl, const magma_int_t *ldvl,
                         double *Vr, const magma_int_t *ldvr,
                         const magma_int_t *mm, magma_int_t *m,
                         double *work,
                         #ifdef MAGMA_COMPLEX
                         double *rwork,
                         #endif
                         magma_int_t *info );

void   lapackf77_dtrevc3( const char *side, const char *howmny,
                          magma_int_t *select, const magma_int_t *n,
                          double *T,  const magma_int_t *ldt,
                          double *VL, const magma_int_t *ldvl,
                          double *VR, const magma_int_t *ldvr,
                          const magma_int_t *mm,
                          const magma_int_t *mout,
                          double *work, const magma_int_t *lwork,
                          #ifdef MAGMA_COMPLEX
                          double *rwork,
                          #endif
                          magma_int_t *info );

void   lapackf77_dtrtri( const char *uplo, const char *diag,
                         const magma_int_t *n,
                         double *A, const magma_int_t *lda,
                         magma_int_t *info );

void   lapackf77_dorg2r( const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         double *A, const magma_int_t *lda,
                         const double *tau,
                         double *work,
                         magma_int_t *info );

void   lapackf77_dorgbr( const char *vect,
                         const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         double *A, const magma_int_t *lda,
                         const double *tau,
                         double *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_dorghr( const magma_int_t *n,
                         const magma_int_t *ilo, const magma_int_t *ihi,
                         double *A, const magma_int_t *lda,
                         const double *tau,
                         double *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_dorglq( const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         double *A, const magma_int_t *lda,
                         const double *tau,
                         double *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_dorgql( const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         double *A, const magma_int_t *lda,
                         const double *tau,
                         double *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_dorgqr( const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         double *A, const magma_int_t *lda,
                         const double *tau,
                         double *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_dorgtr( const char *uplo,
                         const magma_int_t *n,
                         double *A, const magma_int_t *lda,
                         const double *tau,
                         double *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_dorm2r( const char *side, const char *trans,
                         const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         const double *A, const magma_int_t *lda,
                         const double *tau,
                         double *C, const magma_int_t *ldc,
                         double *work,
                         magma_int_t *info );

void   lapackf77_dormbr( const char *vect, const char *side, const char *trans,
                         const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         const double *A, const magma_int_t *lda,
                         const double *tau,
                         double *C, const magma_int_t *ldc,
                         double *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_dormlq( const char *side, const char *trans,
                         const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         const double *A, const magma_int_t *lda,
                         const double *tau,
                         double *C, const magma_int_t *ldc,
                         double *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_dormql( const char *side, const char *trans,
                         const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         const double *A, const magma_int_t *lda,
                         const double *tau,
                         double *C, const magma_int_t *ldc,
                         double *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_dormqr( const char *side, const char *trans,
                         const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         const double *A, const magma_int_t *lda,
                         const double *tau,
                         double *C, const magma_int_t *ldc,
                         double *work, const magma_int_t *lwork,
                         magma_int_t *info );

void   lapackf77_dormrq( const char *side, const char *trans, 
                         magma_int_t *m, magma_int_t *n, magma_int_t *k, 
                         double *A, magma_int_t *lda,
                         double *tau, 
                         double *C, magma_int_t *ldc,
                         double *work, magma_int_t *lwork, 
                         magma_int_t *info );

void   lapackf77_dormtr( const char *side, const char *uplo, const char *trans,
                         const magma_int_t *m, const magma_int_t *n,
                         const double *A, const magma_int_t *lda,
                         const double *tau,
                         double *C, const magma_int_t *ldc,
                         double *work, const magma_int_t *lwork,
                         magma_int_t *info );

/*
 * Real precision extras
 */
void   lapackf77_dstebz( const char *range, const char *order,
                         const magma_int_t *n,
                         const double *vl, const double *vu,
                         const magma_int_t *il, const magma_int_t *iu,
                         const double *abstol,
                         const double *d, const double *e,
                         magma_int_t *m, magma_int_t *nsplit,
                         double *w,
                         magma_int_t *iblock, magma_int_t *isplit,
                         double *work,
                         magma_int_t *iwork,
                         magma_int_t *info );

void   lapackf77_dlaln2( const magma_int_t *ltrans,
                         const magma_int_t *na, const magma_int_t *nw,
                         const double *smin, const double *ca,
                         const double *a,  const magma_int_t *lda,
                         const double *d1, const double *d2,
                         const double *b,  const magma_int_t *ldb,
                         const double *wr, const double *wi,
                         double *x, const magma_int_t *ldx,
                         double *scale, double *xnorm,
                         magma_int_t *info );

double lapackf77_dlamc3( const double *a, const double *b );

void   lapackf77_dlamrg( const magma_int_t *n1, const magma_int_t *n2,
                         const double *a,
                         const magma_int_t *dtrd1, const magma_int_t *dtrd2,
                         magma_int_t *index );

double lapackf77_dlapy3( const double *x, const double *y, const double *z );

void   lapackf77_dlaed2( magma_int_t *k, const magma_int_t *n, const magma_int_t *n1,
                         double *d,
                         double *q, const magma_int_t *ldq,
                         magma_int_t *indxq,
                         double *rho, const double *z,
                         double *dlamda, double *w, double *q2,
                         magma_int_t *indx, magma_int_t *indxc, magma_int_t *indxp,
                         magma_int_t *coltyp,
                         magma_int_t *info);

void   lapackf77_dlaed4( const magma_int_t *n, const magma_int_t *i,
                         const double *d,
                         const double *z,
                         double *delta,
                         const double *rho,
                         double *dlam,
                         magma_int_t *info );

void   lapackf77_dlasrt( const char *id, const magma_int_t *n, double *d,
                         magma_int_t *info );

/*
 * Testing functions
 */
void   lapackf77_dbdt01( const magma_int_t *m, const magma_int_t *n, const magma_int_t *kd,
                         double *A, const magma_int_t *lda,
                         double *Q, const magma_int_t *ldq,
                         double *d, double *e,
                         double *Pt, const magma_int_t *ldpt,
                         double *work,
                         #ifdef MAGMA_COMPLEX
                         double *rwork,
                         #endif
                         double *resid );

void   lapackf77_dget22( const char *transa, const char *transe, const char *transw, const magma_int_t *n,
                         double *A, const magma_int_t *lda,
                         double *E, const magma_int_t *lde,
                         #ifdef MAGMA_COMPLEX
                         double *w,
                         #else
                         double *wr,
                         double *wi,
                         #endif
                         double *work,
                         #ifdef MAGMA_COMPLEX
                         double *rwork,
                         #endif
                         double *result );

void   lapackf77_dsyt21( const magma_int_t *itype, const char *uplo,
                         const magma_int_t *n, const magma_int_t *kband,
                         double *A, const magma_int_t *lda,
                         double *d, double *e,
                         double *U, const magma_int_t *ldu,
                         double *V, const magma_int_t *ldv,
                         double *tau,
                         double *work,
                         #ifdef MAGMA_COMPLEX
                         double *rwork,
                         #endif
                         double *result );

void   lapackf77_dhst01( const magma_int_t *n, const magma_int_t *ilo, const magma_int_t *ihi,
                         double *A, const magma_int_t *lda,
                         double *H, const magma_int_t *ldh,
                         double *Q, const magma_int_t *ldq,
                         double *work, const magma_int_t *lwork,
                         #ifdef MAGMA_COMPLEX
                         double *rwork,
                         #endif
                         double *result );

void   lapackf77_dstt21( const magma_int_t *n, const magma_int_t *kband,
                         double *AD,
                         double *AE,
                         double *SD,
                         double *SE,
                         double *U, const magma_int_t *ldu,
                         double *work,
                         #ifdef MAGMA_COMPLEX
                         double *rwork,
                         #endif
                         double *result );

void   lapackf77_dort01( const char *rowcol, const magma_int_t *m, const magma_int_t *n,
                         double *U, const magma_int_t *ldu,
                         double *work, const magma_int_t *lwork,
                         #ifdef MAGMA_COMPLEX
                         double *rwork,
                         #endif
                         double *resid );

void   lapackf77_dlarfy( const char *uplo, const magma_int_t *n,
                         double *V, const magma_int_t *incv,
                         double *tau,
                         double *C, const magma_int_t *ldc,
                         double *work );

double lapackf77_dqpt01( const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         double *A,
                         double *Af, const magma_int_t *lda,
                         double *tau, magma_int_t *jpvt,
                         double *work, const magma_int_t *lwork );

void   lapackf77_dqrt02( const magma_int_t *m, const magma_int_t *n, const magma_int_t *k,
                         double *A,
                         double *AF,
                         double *Q,
                         double *R, const magma_int_t *lda,
                         double *tau,
                         double *work, const magma_int_t *lwork,
                         double *rwork,
                         double *result );

void   lapackf77_dlatms( const magma_int_t *m, const magma_int_t *n,
                         const char *dist, magma_int_t *iseed, const char *sym,
                         double *d,
                         const magma_int_t *mode, const double *cond,
                         const double *dmax,
                         const magma_int_t *kl, const magma_int_t *ku, const char *pack,
                         double *A, const magma_int_t *lda,
                         double *work,
                         magma_int_t *info );

#ifdef __cplusplus
}
#endif

#undef MAGMA_REAL

#endif /* MAGMA_DLAPACK_H */
