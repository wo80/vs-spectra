#pragma warning (disable : 4250)

#include "spectra_di_sg.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Spectra/SymGEigsSolver.h>
#include <Spectra/SymGEigsShiftSolver.h>
#include <Spectra/MatOp/SymShiftInvert.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/MatOp/SparseCholesky.h>

using namespace Spectra;

using MapVec = Eigen::Map<Eigen::VectorXd>;
using MapMat = Eigen::Map<Eigen::MatrixXd>;
using SpMatrix = Eigen::SparseMatrix<double>;

template <typename TypeA, typename TypeB, typename Solver>
int solve_di_sg(const TypeA& A, const TypeB& B, Solver& eigs, int maxit, double tol, SortRule selection, spectra_result* result)
{
	eigs.init();

	int nconv = eigs.compute(selection, maxit, tol);

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

int spectra_di_sg(int which, int k, int ncv, int maxit, double tol,
	spectra_spmat* A, spectra_spmat* B, spectra_result* result)
{
	try
	{
		SortRule selection;

		SELECTION_RULE_FROM_INT(which)

		// We are going to calculate the eigenvalues of M
		Eigen::Map<const SpMatrix> M(A->n, A->n, A->nnz, A->p, A->i, (double*)A->x);
		Eigen::Map<const SpMatrix> N(A->n, A->n, B->nnz, B->p, B->i, (double*)B->x);

		SparseSymMatProd<double> op(M);
		SparseCholesky<double> Bop(N);

		SymGEigsSolver<double, SparseSymMatProd<double>, SparseCholesky<double>, GEigsMode::Cholesky> eigs(op, Bop, k, ncv);

		return solve_di_sg(M, N, eigs, maxit, tol, selection, result);
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

int spectra_di_sg_shift(int which, char mode, int k, int ncv, int maxit, double tol, double sigma,
	spectra_spmat* A, spectra_spmat* B, spectra_result* result)
{
	try
	{
		SortRule selection;

		SELECTION_RULE_FROM_INT(which)

		// We are going to calculate the eigenvalues of M
		Eigen::Map<const SpMatrix> M(A->n, A->n, A->nnz, A->p, A->i, (double*)A->x);
		Eigen::Map<const SpMatrix> N(A->n, A->n, B->nnz, B->p, B->i, (double*)B->x);

		using OpType = SymShiftInvert<double, Eigen::Sparse, Eigen::Sparse>;
		using BOpType = SparseSymMatProd<double>;

		OpType op(M, N);

		if (mode == 'S')
		{
			BOpType Bop(N);
			SymGEigsShiftSolver<double, OpType, BOpType, GEigsMode::ShiftInvert> eigs(op, Bop, k, ncv, sigma);

			return solve_di_sg(M, N, eigs, maxit, tol, selection, result);
		}
		else if (mode == 'B')
		{
			BOpType Bop(M);
			SymGEigsShiftSolver<double, OpType, BOpType, GEigsMode::Buckling> eigs(op, Bop, k, ncv, sigma);

			return solve_di_sg(M, N, eigs, maxit, tol, selection, result);
		}
		else if (mode == 'C')
		{
			BOpType Bop(N);
			SymGEigsShiftSolver<double, OpType, BOpType, GEigsMode::Cayley> eigs(op, Bop, k, ncv, sigma);

			return solve_di_sg(M, N, eigs, maxit, tol, selection, result);
		}

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
