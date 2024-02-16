#pragma warning (disable : 4250)

#include "spectra_di_ns.h"

#include <stdexcept>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/GenEigsRealShiftSolver.h>
#include <Spectra/GenEigsComplexShiftSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/MatOp/SparseGenRealShiftSolve.h>

using namespace Spectra;

using MapVec = Eigen::Map<Eigen::VectorXcd>;
using MapMat = Eigen::Map<Eigen::MatrixXcd>;
using SpMatrix = Eigen::SparseMatrix<double>;

template <typename TypeA, typename Solver>
int solve_di_ns(const TypeA& A, Solver& eigs, int maxit, double tol, SortRule selection, spectra_result* result)
{
    eigs.init();

    int nconv = eigs.compute(selection, maxit, tol, SortRule::SmallestMagn);

    CompInfo info = eigs.info();

    if (info == CompInfo::Successful)
    {
        auto evals = static_cast<Eigen::dcomplex*>(result->eigval);
        MapVec v(evals, nconv, 1);
        v.noalias() = eigs.eigenvalues();

        if (result->eigvec != NULL)
        {
            auto evecs = static_cast<Eigen::dcomplex*>(result->eigvec);
            MapMat m(evecs, A.rows(), nconv);
            m.noalias() = eigs.eigenvectors();
        }
    }

    result->iterations = eigs.num_iterations();
    result->info = static_cast<int>(info);

    return nconv;
}
int spectra_di_ns(int which, int k, int ncv, int maxit, double tol,
    spectra_spmat* A, spectra_result* result)
{
    try
    {
        SortRule selection;

        SELECTION_RULE_FROM_INT(which)

        // We are going to calculate the eigenvalues of M
        Eigen::Map<const SpMatrix> M(A->n, A->n, A->nnz, A->p, A->i, (double*)A->x);

        SparseGenMatProd<double> op(M);
        GenEigsSolver<SparseGenMatProd<double>> eigs(op, k, ncv);

        return solve_di_ns(M, eigs, maxit, tol, selection, result);
    }
    catch (std::invalid_argument& e)
    {
        result->info = -1002;
    }
    catch (std::exception& e)
    {
        result->info = -1001;
    }
    catch (...)
    {
        result->info = -1000;
    }

    return 0;
}

int spectra_di_ns_shift(int which, int k, int ncv, int maxit, double tol, double sigma,
    spectra_spmat* A, spectra_result* result)
{
    try
    {
        SortRule selection;

        SELECTION_RULE_FROM_INT(which)

        // We are going to calculate the eigenvalues of M
        Eigen::Map<const SpMatrix> M(A->n, A->n, A->nnz, A->p, A->i, (double*)A->x);

        SparseGenRealShiftSolve<double> op(M);
        GenEigsRealShiftSolver<SparseGenRealShiftSolve<double>> eigs(op, k, ncv, sigma);

        return solve_di_ns(M, eigs, maxit, tol, selection, result);
    }
    catch (std::invalid_argument& e)
    {
        result->info = -1002;
    }
    catch (std::exception& e)
    {
        result->info = -1001;
    }
    catch (...)
    {
        result->info = -1000;
    }

    return 0;
}
