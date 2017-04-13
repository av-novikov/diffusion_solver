#include "model/Acid/2d/Acid2dSolver.hpp"

using namespace acid2d;

Acid2dSolver::Acid2dSolver(Acid2d* _model) : basic2d::Basic2dSolver<Acid2d>(_model)
{
	//P.open("snaps/P.dat", ofstream::out);
	//S.open("snaps/S.dat", ofstream::out);
	//qcells.open("snaps/q_cells.dat", ofstream::out);

	chop_mult = 0.1;
	max_sat_change = 0.1;
}
Acid2dSolver::~Acid2dSolver()
{
	//P.close();
	//S.close();
	//qcells.close();
}