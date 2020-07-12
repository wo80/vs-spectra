#pragma once

#include "common.h"

#define EXPORT __declspec(dllexport)

#ifdef __cplusplus
extern "C"
{
#endif

    // Method parameters:
    //
    // [in]       which = requested spectrum
    // [in]           k = number of eigenvalues
    // [in]         ncv = number of basis vectors
    // [in]       maxit = maximum number of iterations
    // [in]         tol = tolerance
    // [in]       sigma = sigma for shifted mode
    // [in]           A = csc matrix
    // [in,out]  result = eigenvalues/vectors storage

    // Real symmetric generalized problem - regular mode
    EXPORT int spectra_di_sg(int which, int k, int ncv, int maxit, double tol,
        spectra_spmat *A, spectra_spmat *B, spectra_result *result);

    // Real symmetric generalized problem - shift-and-invert mode (standard, buckling or Caley)
    EXPORT int spectra_di_sg_shift(int which, char mode, int k, int ncv, int maxit, double tol, double sigma,
        spectra_spmat *A, spectra_spmat *B, spectra_result *result);

#ifdef __cplusplus
}
#endif