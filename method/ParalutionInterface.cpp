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
	BiCGStab<LocalMatrix<double>, LocalVector<double>, double > ls;

	double tick, tack;

	ls.SetOperator(Mat);
	Mat.WriteFileMTX("snaps/mat.mtx");
	//std::cout << Mat.get_ncol() << " " << Mat.get_nrow() << std::endl;
	Rhs.WriteFileASCII("snaps/rhs.dat");
	//std::cout << Rhs.get_size() << std::endl;
	ls.Build();

	//Mat.info();

	tick = paralution_time();
	ls.Solve(Rhs, &x);
	tack = paralution_time();

	std::cout << "Solver execution:" << (tack - tick) / 1000000 << " sec" << std::endl;
}