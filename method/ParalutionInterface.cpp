#include "method/ParalutionInterface.h"

using namespace paralution;

ParSolver::ParSolver()
{
	isAssembled = false;
	ls.Init(1.E-7, 1.E-6, 1E+4, 20000);
}

ParSolver::~ParSolver()
{
}

void ParSolver::Init(int vecSize)
{
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

const paralution::LocalVector<double>& ParSolver::getSolution()
{
	return x;
}

void ParSolver::Solve()
{
	double tick, tack;
	
//	Mat.WriteFileMTX("snaps/mat.mtx");
//	Rhs.WriteFileASCII("snaps/rhs.dat");

	ls.SetOperator(Mat);
	p.Set(2);
	ls.SetPreconditioner(p);
	ls.Build();
	isAssembled = true;
	
	ls.Init(1.E-7, 1.E-6, 1E+4, 5000);
	Mat.info();

	//ls.RecordResidualHistory();
	//std::cout << ls.GetCurrentResidual() << std::endl;
	tick = paralution_time();
	ls.Solve(Rhs, &x);
	tack = paralution_time();
	//ls.RecordHistory("snaps/history.dat");
	//x.WriteFileASCII("snaps/x.dat");
	//std::cout << ls.GetCurrentResidual() << std::endl;
	std::cout << "Solver execution:" << (tack - tick) / 1000000 << " sec" << std::endl << std::endl;

	x.MoveToHost();

	ls.Clear();
}