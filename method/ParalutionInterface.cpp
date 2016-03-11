#include "method/ParalutionInterface.h"

#include <fstream>

using namespace paralution;
using std::ifstream;
using std::string;

ParSolver::ParSolver() : resHistoryFile("snaps/resHistory.dat"), gmresMaxIter(300), bicgstabMaxIter(1000)
{
	isAssembled = false;
	isGmres = true;
	gmres.Init(1.E-15, 1.E-8, 1E+4, 20000);
	bicgstab.Init(1.E-15, 1.E-8, 1E+4, 20000);
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

void ParSolver::Solve(bool setGmres)
{
	double initRes, finalRes;
	int iterNum;

	if(setGmres == true)
		isGmres = setGmres;

	if (isGmres)
	{
		SolveGMRES();

		getResiduals(&initRes, &finalRes, &iterNum);
		std::cout << "Initial residual: " << initRes << std::endl;
		std::cout << "Final residual: " << finalRes << std::endl;
		std::cout << "Number of iterations: " << iterNum << std::endl << std::endl;

		if ( (iterNum >= 3 * gmresMaxIter / 4) && (iterNum < gmresMaxIter) )
		{
			isGmres = false;
		}
		else if (iterNum >= gmresMaxIter)
		{
			isGmres = false;

			x.Zeros();
			x.MoveToAccelerator();

			SolveBiCGStab();

			getResiduals(&initRes, &finalRes, &iterNum);
			std::cout << "Initial residual: " << initRes << std::endl;
			std::cout << "Final residual: " << finalRes << std::endl;
			std::cout << "Number of iterations: " << iterNum << std::endl << std::endl;
		}
	}
	else
	{
		SolveBiCGStab();

		getResiduals(&initRes, &finalRes, &iterNum);
		std::cout << "Initial residual: " << initRes << std::endl;
		std::cout << "Final residual: " << finalRes << std::endl;
		std::cout << "Number of iterations: " << iterNum << std::endl << std::endl;
	}

	x.MoveToHost();
}

void ParSolver::SolveGMRES()
{
	gmres.SetOperator(Mat);
	p.Set(1);
	gmres.SetPreconditioner(p);
	gmres.Build();
	isAssembled = true;

	gmres.Init(1.E-15, 1.E-8, 1E+4, gmresMaxIter);
	//Mat.info();

	gmres.RecordResidualHistory();
	gmres.Solve(Rhs, &x);
	gmres.RecordHistory(resHistoryFile);

	gmres.Clear();
}

void ParSolver::SolveBiCGStab()
{
	bicgstab.SetOperator(Mat);
	p.Set(1);
	bicgstab.SetPreconditioner(p);
	bicgstab.Build();
	isAssembled = true;

	bicgstab.Init(1.E-15, 1.E-8, 1E+4, bicgstabMaxIter);
	//Mat.info();

	bicgstab.RecordResidualHistory();
	bicgstab.Solve(Rhs, &x);
	bicgstab.RecordHistory(resHistoryFile);

	bicgstab.Clear();
}


void ParSolver::getResiduals(double* const initRes, double* const finalRes, int* const iterNum)
{
	ifstream file;
	file.open(resHistoryFile.c_str(), ifstream::in);

	file >> *initRes;

	double tmp;
	int i = 0;
	while (!file.eof())
	{
		file >> tmp;
		i++;
	}
	*finalRes = tmp;
	*iterNum = i;

	file.close();
}