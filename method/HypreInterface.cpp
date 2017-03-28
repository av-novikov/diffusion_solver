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

	HYPRE_IJMatrixDestroy(ij_matrix);
	HYPRE_IJVectorDestroy(b);
	HYPRE_IJVectorDestroy(x);

	HYPRE_ParCSRPCGDestroy(solver);
	HYPRE_EuclidDestroy(precond);
}
void HypreSolver::Init(const int vecSize, const double _relTol, const double _dropTol)
{
	ilower = 0;		iupper = vecSize - 1;
	nrows = vecSize;
	relTol = _relTol;		dropTol = _dropTol;

	HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &ij_matrix);
	HYPRE_IJMatrixSetObjectType(ij_matrix, HYPRE_PARCSR);
	HYPRE_IJMatrixInitialize(ij_matrix);

	HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &b);
	HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
	HYPRE_IJVectorInitialize(b);

	HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &x);
	HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
	HYPRE_IJVectorInitialize(x);

	HYPRE_ParCSRBiCGSTABCreate(MPI_COMM_WORLD, &solver);
	HYPRE_BiCGSTABSetMaxIter(solver, 200);
	HYPRE_BiCGSTABSetTol(solver, relTol);
	HYPRE_BiCGSTABSetPrintLevel(solver, 2);
	HYPRE_BiCGSTABSetLogging(solver, 1);

	HYPRE_EuclidCreate(MPI_COMM_WORLD, &precond);
	HYPRE_EuclidSetILUT(precond, dropTol);
	HYPRE_BiCGSTABSetPrecond(solver, (HYPRE_PtrToSolverFcn)HYPRE_EuclidSolve,
		(HYPRE_PtrToSolverFcn)HYPRE_EuclidSetup, precond);
	
	x_values = new double[nrows];
	x_sol = new double[nrows];
	fill_n(x_values, nrows, 0.0);
}
void HypreSolver::Assemble(int* cols, const int* ind_j, const double* a, const int counter, const int* ind_rhs, const double* rhs)
{
	HYPRE_IJMatrixSetValues(ij_matrix, nrows, cols, ind_rhs, ind_j, a);
	HYPRE_IJMatrixAssemble(ij_matrix);
	HYPRE_IJMatrixGetObject(ij_matrix, (void **)&parcsr_matrix);

	HYPRE_IJVectorSetValues(b, nrows, ind_rhs, rhs);
	HYPRE_IJVectorAssemble(b);
	HYPRE_IJVectorGetObject(b, (void **)&par_b);

	//HYPRE_IJMatrixPrint(ij_matrix, "snaps/mat_out.mtx");
	//HYPRE_IJVectorPrint(b, "snaps/rhs_out.data");

	HYPRE_IJVectorSetValues(x, nrows, ind_rhs, x_values);	
	HYPRE_IJVectorAssemble(x);
	HYPRE_IJVectorGetObject(x, (void **)&par_x);

	rows = const_cast<int*>(ind_rhs);
}
void HypreSolver::Solve()
{

	HYPRE_ParCSRBiCGSTABSetup(solver, parcsr_matrix, par_b, par_x);
	HYPRE_ParCSRBiCGSTABSolve(solver, parcsr_matrix, par_b, par_x);

	//HYPRE_IJVectorPrint(x, "snaps/x_out.data");
}
const HypreSolver::Vector& HypreSolver::getSolution()
{
	HYPRE_IJVectorGetValues(x, nrows, rows, x_sol);
	return x_sol;
}