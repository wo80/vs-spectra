#pragma warning (disable : 4250)

#include "spectra_di_ns.h"

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

/************************ Macros to generate code ************************

#define EIG_CODE_COMPLEX_SHIFT(RULE, OPTYPE)                                   \
GenEigsComplexShiftSolver<double, RULE, OPTYPE> eigs(&op, k, ncv, sigmar, sigmai); \
EIG_COMMON_CODE

/************************ Macros to generate code ************************/

int spectra_di_ns(int which, int k, int ncv, int maxit, double tol,
	spectra_spmat* A, spectra_result* result)
{
	CompInfo info;
	SortRule selection;

	int nconv = 0;
	bool retvec = result->eigvec != NULL;

	try
	{
		SELECTION_RULE_FROM_INT(which)

		// We are going to calculate the eigenvalues of M
		Eigen::Map<const SpMatrix> M(A->n, A->n, A->nnz, A->p, A->i, (double*)A->x);

		// Construct matrix operation object using the wrapper class SparseSymMatProd
		//SparseGenMatProd<double> op(M);

		auto evecs = reinterpret_cast<std::complex<double>*>(result->eigvec);
		auto evals = reinterpret_cast<std::complex<double>*>(result->eigval);

		SparseGenMatProd<double> op(M);
		GenEigsSolver<double, SparseGenMatProd<double>> eigs(op, k, ncv);

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

int spectra_di_ns_shift(int which, int k, int ncv, int maxit, double tol, double sigma,
	spectra_spmat* A, spectra_result* result)
{
	CompInfo info;
	SortRule selection;

	int nconv = 0;
	bool retvec = result->eigvec != NULL;

	try
	{
		SELECTION_RULE_FROM_INT(which)

		// We are going to calculate the eigenvalues of M
		Eigen::Map<const SpMatrix> M(A->n, A->n, A->nnz, A->p, A->i, (double*)A->x);

		auto evecs = reinterpret_cast<std::complex<double>*>(result->eigvec);
		auto evals = reinterpret_cast<std::complex<double>*>(result->eigval);

		SparseGenRealShiftSolve<double> op(M);
		GenEigsRealShiftSolver<double, SparseGenRealShiftSolve<double>> eigs(op, k, ncv, sigma);

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
