#include "method/HypreInterface.hpp"
#include <algorithm>
#include <iostream>

using std::fill_n;

#pragma comment(lib, "msmpi.lib")

HypreSolver::HypreSolver()
{
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	first = true;
}
HypreSolver::~HypreSolver()
{
	delete[] x_values, x_sol;

	HYPRE_IJMatrixDestroy(ij_matrix);
	HYPRE_IJVectorDestroy(b);
	HYPRE_IJVectorDestroy(x);
	HYPRE_EuclidDestroy(precond);
	HYPRE_ParCSRBiCGSTABDestroy(solver);
}
void HypreSolver::Init(const int vecSize, const double _absTol, const double _relTol, const double _dropTol)
{
	ilower = 0;		iupper = vecSize - 1;
	nrows = vecSize;
	absTol = _absTol;
	relTol = _relTol;		dropTol = _dropTol;

	res = HYPRE_ParCSRBiCGSTABCreate(MPI_COMM_WORLD, &solver);
	res = HYPRE_EuclidCreate(MPI_COMM_WORLD, &precond);
	res = HYPRE_ParCSRBiCGSTABSetPrintLevel(solver, 10);
	res = HYPRE_ParCSRBiCGSTABSetLogging(solver, 10);

	res = -1;
	res = HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &ij_matrix);
	res = HYPRE_IJMatrixSetObjectType(ij_matrix, HYPRE_PARCSR);
	res = HYPRE_IJMatrixSetPrintLevel(ij_matrix, 10);
	//res = HYPRE_IJMatrixInitialize(ij_matrix);

	res = HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &b);
	res = HYPRE_IJVectorSetPrintLevel(b, 10);
	res = HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
	//res = HYPRE_IJVectorInitialize(b);

	res = HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &x);
	res = HYPRE_IJVectorSetPrintLevel(x, 10);
	res = HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
	//res = HYPRE_IJVectorInitialize(x);

	x_values = new double[nrows];
	x_sol = new double[nrows];
	fill_n(x_values, nrows, 0.0);
}
void HypreSolver::Assemble(int* cols, const int* ind_i, const int* ind_j, const double* a, const int counter, const int* ind_rhs, const double* rhs)
{
	res = -1;

	res = HYPRE_IJMatrixInitialize(ij_matrix);
	res = HYPRE_IJMatrixSetValues(ij_matrix, nrows, cols, ind_rhs, ind_j, a);
	res = HYPRE_IJMatrixAssemble(ij_matrix);
	res = HYPRE_IJMatrixGetObject(ij_matrix, (void **)&parcsr_matrix);

	res = HYPRE_IJVectorInitialize(b);
	res = HYPRE_IJVectorSetValues(b, nrows, ind_rhs, rhs);
	res = HYPRE_IJVectorAssemble(b);
	res = HYPRE_IJVectorGetObject(b, (void **)&par_b);

	HYPRE_IJMatrixPrint(ij_matrix, "snaps/mat_out.mtx");
	HYPRE_IJVectorPrint(b, "snaps/rhs_out.data");

	res = HYPRE_IJVectorInitialize(x);
	res = HYPRE_IJVectorSetValues(x, nrows, ind_rhs, x_values);	
	res = HYPRE_IJVectorAssemble(x);
	res = HYPRE_IJVectorGetObject(x, (void **)&par_x);

	rows = const_cast<int*>(ind_rhs);
}
void HypreSolver::Solve()
{
	res = -1;

	//res = HYPRE_EuclidSetILUT(precond, 1.E-10);
	//res = HYPRE_ParCSRBiCGSTABSetPrecond(solver, (HYPRE_PtrToParSolverFcn)HYPRE_EuclidSolve,
	//	(HYPRE_PtrToParSolverFcn)HYPRE_EuclidSetup, precond);

	if (first)
	{
		res = HYPRE_EuclidSetLevel(precond, 0);
		res = HYPRE_ParCSRBiCGSTABSetPrecond(solver, (HYPRE_PtrToParSolverFcn)HYPRE_EuclidSolve,
			(HYPRE_PtrToParSolverFcn)HYPRE_EuclidSetup, precond);
		HYPRE_Solver precon_gotten;
		HYPRE_ParCSRBiCGSTABGetPrecond(solver, &precon_gotten);
		if (precon_gotten != precond)
			std::cout << "HYPRE_GMRESGetPrecond got bad precon" << std::endl;

		res = HYPRE_ParCSRBiCGSTABSetMaxIter(solver, dropTol);
		//res = HYPRE_ParCSRBiCGSTABSetAbsoluteTol(solver, absTol);
		res = HYPRE_ParCSRBiCGSTABSetTol(solver, relTol);

		res = HYPRE_ParCSRBiCGSTABSetup(solver, parcsr_matrix, par_b, par_x);
		first = false;
	}
	res = HYPRE_ParCSRBiCGSTABSolve(solver, parcsr_matrix, par_b, par_x);

	res = HYPRE_ParCSRBiCGSTABGetNumIterations(solver, &itersNum);
	res = HYPRE_ParCSRBiCGSTABGetFinalRelativeResidualNorm(solver, &final_res);
	HYPRE_IJVectorPrint(x, "snaps/x_out.data");
}
const HypreSolver::Vector& HypreSolver::getSolution()
{
	HYPRE_IJVectorGetValues(x, nrows, rows, x_sol);
	return x_sol;
}