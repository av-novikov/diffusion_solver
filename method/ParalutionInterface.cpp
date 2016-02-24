#include "method/ParalutionInterface.h"

using namespace paralution;

ParSolver::ParSolver()
{
}

ParSolver::~ParSolver()
{
}

void ParSolver::Solve()
{
	double tick, tack;
	
//	Mat.WriteFileMTX("snaps/mat.mtx");
//	Rhs.WriteFileASCII("snaps/rhs.dat");
	ls.Build();

	//Mat.info();

	tick = paralution_time();
	ls.Solve(Rhs, &x);
	tack = paralution_time();

	//x.WriteFileASCII("snaps/x.dat");

	std::cout << "Solver execution:" << (tack - tick) / 1000000 << " sec" << std::endl;
}