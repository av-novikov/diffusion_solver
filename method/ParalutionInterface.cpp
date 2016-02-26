#include "method/ParalutionInterface.h"

using namespace paralution;

ParSolver::ParSolver()
{
	isAssembled = false;
	ls.Init(1.E-6, 1.E-5, 1E+4, 100000);
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
		AssembleUpdate(a, ind_rhs, rhs);
	else {
		Mat.Assemble(ind_i, ind_j, a, counter, "A", matSize, matSize);
		Rhs.Assemble(ind_rhs, rhs, matSize, "rhs");

		Mat.MoveToAccelerator();
		Rhs.MoveToAccelerator();
		x.MoveToAccelerator();
	}	
}

void ParSolver::AssembleUpdate(const double* a, const int* ind_rhs, const double* rhs)
{
	Mat.AssembleUpdate(a);
	Rhs.Assemble(ind_rhs, rhs, matSize, "rhs");
}

const paralution::LocalVector<double>& ParSolver::getSolution()
{
	return x;
}

void ParSolver::Solve()
{
	double tick, tack;
	
	//Mat.WriteFileMTX("snaps/mat.mtx");
	//Rhs.WriteFileASCII("snaps/rhs.dat");

	/*if (isAssembled)
	{
		ls.ResetOperator(Mat);
		//p.Set(1);
		//ls.SetPreconditioner(p);
	}
	else
	{*/
		ls.SetOperator(Mat);
		p.Set(1);
		ls.SetPreconditioner(p);
		ls.Build();
		isAssembled = true;
	//}

	Mat.info();

	tick = paralution_time();
	ls.Solve(Rhs, &x);
	tack = paralution_time();

	//x.WriteFileASCII("snaps/x.dat");

	std::cout << "Solver execution:" << (tack - tick) / 1000000 << " sec" << std::endl;

	x.MoveToHost();

	ls.Clear();
}