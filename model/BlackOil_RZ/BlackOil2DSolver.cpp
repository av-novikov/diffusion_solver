#include "model/BlackOil_RZ/BlackOil2DSolver.hpp"

using namespace std;
using namespace blackoil_rz;

BlackOil2dSolver::BlackOil2dSolver(BlackOil_RZ* _model) : basic2d::Basic2dSolver<BlackOil_RZ>(_model)
{
	P.open("snaps/P.dat", ofstream::out);
	S.open("snaps/S.dat", ofstream::out);
	qcells.open("snaps/q_cells.dat", ofstream::out);
}
BlackOil2dSolver::~BlackOil2dSolver()
{
	P.close();
	S.close();
	qcells.close();
}
void BlackOil2dSolver::writeData()
{
	double p = 0.0, s_w = 0.0, s_o = 0.0;

	qcells << cur_t * t_dim / 3600.0;

	map<int, double>::iterator it;
	for (it = model->Qcell.begin(); it != model->Qcell.end(); ++it)
	{
		p += model->cells[it->first].u_next.p * model->P_dim;
		s_w += model->cells[it->first].u_next.s_w;
		s_o += model->cells[it->first].u_next.s_o;
		if (model->leftBoundIsRate)
			qcells << "\t" << it->second * model->Q_dim * 86400.0;
		else
			qcells << "\t" << model->getRate(it->first) * model->Q_dim * 86400.0;
	}

	P << cur_t * t_dim / 3600.0 << 
		"\t" << p / (double)(model->Qcell.size()) << endl;
	S << cur_t * t_dim / 3600.0 << 
		"\t" << s_w / (double)(model->Qcell.size()) << 
		"\t" << s_o / (double)(model->Qcell.size()) <<
		"\t" << 1.0 - (s_w + s_o) / (double)(model->Qcell.size()) << endl;

	qcells << endl;
}
void BlackOil2dSolver::solveStep()
{
	int cellIdx, varIdx;
	double err_newton = 1.0;
	double averPresPrev = averValue(0);
	double averSatPrev = averValue(1);
	double averPres, averSat;
	double dAverPres = 1.0, dAverSat = 1.0;

	iterations = 0;
	while (err_newton > 1.e-4 && (dAverSat > 1.e-10 || dAverPres > 1.e-10) && iterations < 9)
	{
		copyIterLayer();

		Solve(model->cellsNum_r + 1, 2 * (model->cellsNum_z + 2), PRES);
		construction_from_fz(model->cellsNum_r + 2, 2 * (model->cellsNum_z + 2), PRES);
		model->solveP_bub();

		err_newton = convergance(cellIdx, varIdx);

		averPres = averValue(0);					averSat = averValue(1);
		dAverPres = fabs(averPres - averPresPrev);	dAverSat = fabs(averSat - averSatPrev);
		averPresPrev = averPres;					averSatPrev = averSat;

		iterations++;
	}

	cout << "Newton Iterations = " << iterations << endl;
}
void BlackOil2dSolver::construction_from_fz(int N, int n, int key)
{}
void BlackOil2dSolver::MiddleAppr(int current, int MZ, int key)
{}
void BlackOil2dSolver::LeftBoundAppr(int MZ, int key)
{}
void BlackOil2dSolver::RightBoundAppr(int MZ, int key)
{}