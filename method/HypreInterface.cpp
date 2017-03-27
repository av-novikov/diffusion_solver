#include "method/HypreInterface.hpp"
#include <algorithm>

using std::fill_n;

HypreSolver::HypreSolver()
{
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
}
HypreSolver::~HypreSolver()
{
	delete[] x_values, x_sol;

	HYPRE_ParCSRPCGDestroy(solver);
	HYPRE_EuclidDestroy(precond);
}
void HypreSolver::Init(const int vecSize)
{
	ilower = 0;		iupper = vecSize - 1;
	nrows = vecSize;

	HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &ij_matrix);
	HYPRE_IJMatrixSetObjectType(ij_matrix, HYPRE_PARCSR);
	HYPRE_IJMatrixInitialize(ij_matrix);

	HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &b);
	HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
	HYPRE_IJVectorInitialize(b);

	HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &x);
	HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
	HYPRE_IJVectorInitialize(x);

	HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);
	HYPRE_PCGSetMaxIter(solver, 100);
	HYPRE_PCGSetTol(solver, 1e-8);
	HYPRE_PCGSetTwoNorm(solver, 1);
	HYPRE_PCGSetPrintLevel(solver, 2);
	HYPRE_PCGSetLogging(solver, 1);

	HYPRE_EuclidCreate(MPI_COMM_WORLD, &precond);
	HYPRE_EuclidSetILUT(precond, 0.01);
	HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn)HYPRE_EuclidSolve,
		(HYPRE_PtrToSolverFcn)HYPRE_EuclidSetup, precond);
	
	x_values = new double[nrows];
	x_sol = new double[nrows];
	fill_n(x_values, nrows, 0.0);
}
void HypreSolver::Assemble(const int* cols, const int* ind_j, const double* a, const int counter, const int* ind_rhs, const double* rhs)
{
	auto cols_hypre = const_cast<HYPRE_Int*>(cols);
	HYPRE_IJMatrixSetValues(ij_matrix, nrows, cols_hypre, ind_rhs, ind_j, a);
	HYPRE_IJMatrixAssemble(ij_matrix);
	HYPRE_IJMatrixGetObject(ij_matrix, (void **)&parcsr_matrix);

	HYPRE_IJVectorSetValues(b, nrows, ind_rhs, rhs);
	HYPRE_IJVectorAssemble(b);
	HYPRE_IJVectorGetObject(b, (void **)&par_b);

	HYPRE_IJVectorSetValues(x, nrows, ind_rhs, x_values);	
	HYPRE_IJVectorAssemble(x);
	HYPRE_IJVectorGetObject(x, (void **)&par_x);

	rows = const_cast<int*>(ind_rhs);
}
void HypreSolver::Solve()
{
	HYPRE_ParCSRPCGSetup(solver, parcsr_matrix, par_b, par_x);
	HYPRE_ParCSRPCGSolve(solver, parcsr_matrix, par_b, par_x);
}
const HypreSolver::Vector& HypreSolver::getSolution()
{
	HYPRE_IJVectorGetValues(x, nrows, rows, x_sol);
	return x_sol;
}