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

int spectra_di_ss(int which, int k, int ncv, int maxit, double tol,
	spectra_spmat *A, spectra_result *result)
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

		auto evecs = reinterpret_cast<double *>(result->eigvec);
		auto evals = reinterpret_cast<double *>(result->eigval);

		SparseSymMatProd<double> op(M);
		SymEigsSolver<double, SparseSymMatProd<double>> eigs(op, k, ncv);

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
	catch(...)
	{
		result->info = -1000;
	}

	return nconv;
}

int spectra_di_ss_shift(int which,  int k, int ncv, int maxit, double tol, double sigma,
				 spectra_spmat *A, spectra_result *result)
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

		auto evecs = reinterpret_cast<double *>(result->eigvec);
		auto evals = reinterpret_cast<double *>(result->eigval);

		SparseSymShiftSolve<double> op(M);
		SymEigsShiftSolver<double, SparseSymShiftSolve<double>> eigs(op, k, ncv, sigma);

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
	catch (...)
	{
		result->info = -1000;
	}

	return nconv;
}

#undef AR_TYPE
