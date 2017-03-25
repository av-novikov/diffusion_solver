#include "method/ParalutionInterface.h"

#include <fstream>
#include <iostream>

using namespace paralution;
using std::ifstream;
using std::cout;
using std::endl;

ParSolver::ParSolver() : resHistoryFile("snaps/resHistory.dat")
{
	isAssembled = false;
	isPrecondBuilt = false;
	gmres.Init(1.E-22, 1.E-15, 1E+4, 10000);
	bicgstab.Init(1.E-22, 1.E-15, 1E+4, 10000);
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
	//SolveGMRES();
	SolveBiCGStab();
		
	x.MoveToHost();
}
void ParSolver::SolveBiCGStab()
{
	bicgstab.SetOperator(Mat);
	p.Set(1.E-12, 200);
	bicgstab.SetPreconditioner(p);
	bicgstab.Build();
	isAssembled = true;

	bicgstab.Init(1.E-22, 1.E-15, 1E+4, 10000);
	Mat.info();

	//bicgstab.RecordResidualHistory();
	bicgstab.Solve(Rhs, &x);
	//bicgstab.RecordHistory(resHistoryFile);
	//writeSystem();

	//getResiduals();
	//cout << "Initial residual: " << initRes << endl;
	//cout << "Final residual: " << finalRes << endl;
	//cout << "Number of iterations: " << iterNum << endl << endl;

	bicgstab.Clear();
}
void ParSolver::SolveGMRES()
{
	gmres.SetOperator(Mat);
	p.Set(1.E-20, 100);
	gmres.SetPreconditioner(p);
	gmres.Build();
	isAssembled = true;

	gmres.Init(1.E-22, 1.E-15, 1E+4, 10000);
	Mat.info();

	//gmres.RecordResidualHistory();
	gmres.Solve(Rhs, &x);
	//gmres.RecordHistory(resHistoryFile);
	//writeSystem();


	//getResiduals();
	//cout << "Initial residual: " << initRes << endl;
	//cout << "Final residual: " << finalRes << endl;
	//cout << "Number of iterations: " << iterNum << endl << endl;

	gmres.Clear();
}
void ParSolver::getResiduals()
{
	double tmp;
	int i = 0;

	ifstream file;
	file.open(resHistoryFile, ifstream::in);

	file >> initRes;
	while ( !file.eof() )
	{
		file >> tmp;
		i++;
	}
	finalRes = tmp;
	iterNum = i;

	file.close();
}