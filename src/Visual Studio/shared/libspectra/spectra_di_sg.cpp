#pragma warning (disable : 4250)

#include "spectra_di_sg.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Spectra/SymGEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/MatOp/SparseCholesky.h>

using namespace Spectra;

using MapVec = Eigen::Map<Eigen::VectorXd>;
using MapMat = Eigen::Map<Eigen::MatrixXd>;
using SpMatrix = Eigen::SparseMatrix<double>;

int spectra_di_sg(int which, int k, int ncv, int maxit, double tol,
	spectra_spmat *A, spectra_spmat *B, spectra_result *result)
{
	CompInfo info;
	SortRule selection;

	int nconv = 0;
	bool retvec = result->eigvec != NULL;

	try
	{
		SELECTION_RULE_FROM_INT(which)

		// We are going to calculate the eigenvalues of M
		Eigen::Map<const SpMatrix> M(A->n, A->n, A->nnz, A->p, A->i, (double *)A->x);
		Eigen::Map<const SpMatrix> N(A->n, A->n, B->nnz, B->p, B->i, (double *)B->x);

		SparseSymMatProd<double> op(M);
		SparseCholesky<double> Bop(N);

		auto evecs = reinterpret_cast<double *>(result->eigvec);
		auto evals = reinterpret_cast<double *>(result->eigval);

		SymGEigsSolver<double, SparseSymMatProd<double>, SparseCholesky<double>, GEigsMode::Cholesky> eigs(op, Bop, k, ncv);

		eigs.init();
		nconv = eigs.compute(selection, maxit, tol);
		info = eigs.info();

		if (info == CompInfo::Successful)
		{
			MapVec v(evals, nconv, 1);
			v.noalias() = eigs.eigenvalues();

			if (retvec)
			{
				MapMat m(evecs, A->n, nconv);
				m.noalias() = eigs.eigenvectors();
			}
		}

		result->iterations = eigs.num_iterations();
		result->info = static_cast<int>(info);
	}
	catch (std::exception& e)
	{
		result->info = -1001;
	}
	catch (...)
	{
		result->info = -1000;
	}

	return nconv;
}

/*
int spectra_di_sg_shift(int which, int k, int ncv, int maxit, double tol, double sigma,
	spectra_spmat *A, spectra_result *result)
{
	int nconv = 0;
	CompInfo info;
	SortRule selection;

	bool retvec = result->eigvec != NULL;

	int n = A->n;

	try
	{
		SELECTION_RULE_FROM_INT(which)

		// We are going to calculate the eigenvalues of M
		Eigen::Map<const SpMatrix> M(n, n, A->nnz, A->p, A->i, (double *)A->x);

		SparseSymShiftSolve<double> op(M);

		double *evals = (double *)result->eigval;
		double *evecs = (double *)result->eigvec;

		EIG_CODE_GENERATOR(REAL_SHIFT, SparseSymShiftSolve<double>)

			result->iterations = niter;
		result->info = static_cast<int>(info);
	}
	catch (...)
	{
		result->info = -1000;
	}

	return nconv;
}
//*/

#undef AR_TYPE
