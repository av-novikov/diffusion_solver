#include "method/ParalutionInterface.h"

#include <fstream>
#include <iostream>

using namespace paralution;
using std::ifstream;
using std::cout;
using std::endl;

ParSolver::ParSolver() : resHistoryFile("snaps/resHistory.dat"), max_iter(500)
{
	isAssembled = false;
	isPrecondBuilt = false;
    abs_tol = 1.E-15;
    rel_tol = 1.E-10;
    div_crit = 1.E+8;
	//gmres.Init(1.E-15, 1.E-10, 1E+12, max_iter);
	//bicgstab.Init(1.E-15, 1.E-10, 1E+12, max_iter);
}
ParSolver::~ParSolver()
{
}
void ParSolver::Init(const int vecSize, const double absTol, const double relTol, const double dropTol)
{
    abs_tol = absTol;
    rel_tol = relTol;
    div_crit = dropTol;
	matSize = vecSize;
	x.Allocate("x", vecSize);
}
void ParSolver::Assemble(const int* ind_i, const int* ind_j, const double* a, const int counter, const int* ind_rhs, const double* rhs)
{
	Mat.Zeros();
	Rhs.Zeros();
	x.Zeros();

	if (isAssembled)
	{
		Mat.AssembleUpdate(a);
		Rhs.Assemble(ind_rhs, rhs, matSize, "rhs");
	} else {
		Mat.Assemble(ind_i, ind_j, a, counter, "A", matSize, matSize);
		Rhs.Assemble(ind_rhs, rhs, matSize, "rhs");

		Mat.MoveToAccelerator();
		Rhs.MoveToAccelerator();
		x.MoveToAccelerator();
	}	
}
void ParSolver::Solve()
{
	//SolveGMRES();
	SolveBiCGStab();
		
	x.MoveToHost();
}
void ParSolver::Solve(const PRECOND key)
{
	//SolveGMRES();
	if (key == PRECOND::ILU_SERIOUS)
		SolveBiCGStab();
	else if (key == PRECOND::ILU_SIMPLE)
		SolveBiCGStab_Simple(0);
	else if (key == PRECOND::ILU_GMRES)
		SolveGMRES();

	x.MoveToHost();
}
void ParSolver::Solve(const PRECOND key, bool isHarder, const int degree)
{
	//SolveGMRES();
	if (key == PRECOND::ILU_SERIOUS)
		SolveBiCGStab();
	else if (key == PRECOND::ILU_SIMPLE)
		SolveBiCGStab_Simple(degree, isHarder);
	else if (key == PRECOND::ILU_GMRES)
		SolveGMRES();

	x.MoveToHost();
}
void ParSolver::SolveBiCGStab()
{
	bicgstab.SetOperator(Mat);
	p_ilut.Set(1.E-20, 200);
	//p.Set(1);
	bicgstab.SetPreconditioner(p_ilut);	
	bicgstab.Build();
	//isAssembled = true;

	bicgstab.Init(abs_tol, rel_tol, div_crit, max_iter);
	Mat.info();

    bicgstab.RecordResidualHistory();
    bicgstab.Solve(Rhs, &x);
    bicgstab.RecordHistory(resHistoryFile);
    getResiduals();
    iter_num = bicgstab.GetIterationCount();
    final_res = bicgstab.GetCurrentResidual();
    status = static_cast<RETURN_TYPE>(bicgstab.GetSolverStatus());
	//writeSystem();
	bicgstab.Clear();
}
void ParSolver::SolveBiCGStab_Simple(const int degree, bool isHarder)
{
	bicgstab.SetOperator(Mat);
	if (!isHarder)
	{
		p.Set(degree);
		bicgstab.SetPreconditioner(p);
	}
	else
	{
        p_ilut.Set(1.E-12, 20);
		bicgstab.SetPreconditioner(p_ilut);
	}
	
	bicgstab.Build();
	//isAssembled = true;

	bicgstab.Init(abs_tol, rel_tol, div_crit, max_iter);
	Mat.info();

	bicgstab.RecordResidualHistory();
    bicgstab.Solve(Rhs, &x);
    bicgstab.RecordHistory(resHistoryFile);
    getResiduals();
    iter_num = bicgstab.GetIterationCount();
    final_res = bicgstab.GetCurrentResidual();
	status = static_cast<RETURN_TYPE>(bicgstab.GetSolverStatus());
	//writeSystem();
	bicgstab.Clear();
}
void ParSolver::SolveGMRES()
{
	gmres.SetOperator(Mat);
	p_ilut.Set(1.E-20, 100);
	//p.Set(3);
	gmres.SetPreconditioner(p_ilut);
	gmres.Build();
	//isAssembled = true;

	gmres.Init(abs_tol, rel_tol, div_crit, max_iter);
	Mat.info();

    gmres.RecordResidualHistory();
    gmres.Solve(Rhs, &x);
    gmres.RecordHistory(resHistoryFile);
    getResiduals();
    iter_num = gmres.GetIterationCount();
    final_res = gmres.GetCurrentResidual();
    status = static_cast<RETURN_TYPE>(gmres.GetSolverStatus());
	//writeSystem();
	gmres.Clear();
}
void ParSolver::getResiduals()
{
	double tmp;
	int i = 0;

	ifstream file;
	file.open(resHistoryFile, ifstream::in);

	file >> init_res;
	while ( !file.eof() )
	{
		file >> tmp;
		i++;
	}
	final_res = tmp;
	iter_num = i;

	file.close();
}