#pragma warning (disable : 4250)

#include "spectra_di_ss.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>

using namespace Spectra;

using MapVec = Eigen::Map<Eigen::VectorXd>;
using MapMat = Eigen::Map<Eigen::MatrixXd>;
using SpMatrix = Eigen::SparseMatrix<double>;

template <typename TypeA, typename Solver>
int solve_di_ss(const TypeA& A, Solver& eigs, int maxit, double tol, SortRule selection, spectra_result* result)
{
    eigs.init();

    int nconv = eigs.compute(selection, maxit, tol, SortRule::SmallestAlge);

    CompInfo info = eigs.info();

    if (info == CompInfo::Successful)
    {
        auto evals = static_cast<double*>(result->eigval);
        MapVec v(evals, nconv, 1);
        v.noalias() = eigs.eigenvalues();

        if (result->eigvec != NULL)
        {
            auto evecs = static_cast<double*>(result->eigvec);
            MapMat m(evecs, A.rows(), nconv);
            m.noalias() = eigs.eigenvectors();
        }
    }

    result->iterations = eigs.num_iterations();
    result->info = static_cast<int>(info);

    return nconv;
}

int spectra_di_ss(int which, int k, int ncv, int maxit, double tol,
    spectra_spmat *A, spectra_result *result)
{
    try
    {
        SortRule selection;

        SELECTION_RULE_FROM_INT(which)

        // We are going to calculate the eigenvalues of M
        Eigen::Map<const SpMatrix> M(A->n, A->n, A->nnz, A->p, A->i, (double *)A->x);

        SparseSymMatProd<double> op(M);
        SymEigsSolver<double, SparseSymMatProd<double>> eigs(op, k, ncv);

        return solve_di_ss(M, eigs, maxit, tol, selection, result);
    }
    catch (std::exception& e)
    {
        result->info = -1001;
    }
    catch(...)
    {
        result->info = -1000;
    }

    return 0;
}

int spectra_di_ss_shift(int which,  int k, int ncv, int maxit, double tol, double sigma,
                 spectra_spmat *A, spectra_result *result)
{
    try
    {
        SortRule selection;

        SELECTION_RULE_FROM_INT(which)

        // We are going to calculate the eigenvalues of M
        Eigen::Map<const SpMatrix> M(A->n, A->n, A->nnz, A->p, A->i, (double *)A->x);

        SparseSymShiftSolve<double> op(M);
        SymEigsShiftSolver<double, SparseSymShiftSolve<double>> eigs(op, k, ncv, sigma);

        return solve_di_ss(M, eigs, maxit, tol, selection, result);
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

#undef AR_TYPE
