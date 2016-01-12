#include "model/Oil1D_NIT/Oil1DNITSolver.h"

using namespace std;
using namespace oil1D_NIT;

Oil1DNITSolver::Oil1DNITSolver(Oil1D_NIT* _model) : AbstractSolver<Oil1D_NIT>(_model)
{
	Initialize(model->cellsNum_r+2, 1);

	plot_Pdyn.open("snaps/P_dyn.dat", ofstream::out);

	t_dim = model->t_dim;

	Tt = model->period[model->period.size()-1];
}

Oil1DNITSolver::~Oil1DNITSolver()
{
	plot_Pdyn.close();
}

void Oil1DNITSolver::writeData()
{
	plot_Pdyn << cur_t * t_dim / 3600.0 << "\t" << model->cells[idx1].u_next.p << endl;
}

void Oil1DNITSolver::control()
{	
	writeData();

	if(cur_t >= model->period[curTimePeriod])
	{
		curTimePeriod++;
		model->ht = model->ht_min;
		model->setPeriod(curTimePeriod);
	}

	if(model->ht <= model->ht_max && iterations < 6)
		model->ht = model->ht * 1.5;
	else if(iterations > 6 && model->ht > model->ht_min)
		model->ht = model->ht / 1.5;

	if(cur_t + model->ht > model->period[curTimePeriod])
		model->ht = model->period[curTimePeriod] - cur_t;

	cur_t += model->ht;
}

void Oil1DNITSolver::doNextStep()
{
	int cellIdx, varIdx;
	double err_newton = 1.0;
	double averPresPrev = averValue(0);
	double averPres;
	double dAverPres = 1.0;
	
	iterations = 0;
	while( err_newton > 1.e-4 && dAverPres > 1.e-4 && iterations < 8 )
	{	
		copyIterLayer();

		Solve(model->cellsNum_r+1, 1, PRES);
		construction_from_fz(model->cellsNum_r+2, 1, PRES);

		err_newton = convergance(cellIdx, varIdx);

		averPres = averValue(0);
		dAverPres = fabs(averPres - averPresPrev);
		averPresPrev = averPres;

		//cout << "BadPresValue[" << cellIdx  << "]: " << model->cells[cellIdx].u_next.p << endl;

		iterations++;
	}
}

void Oil1DNITSolver::construction_from_fz(int N, int n, int key)
{
	vector<Cell>::iterator it;
	if (key == PRES)
		for(int i = 0; i < N; i++)
			 model->cells[i].u_next.p += fz[i][1];
}

void Oil1DNITSolver::LeftBoundAppr(int MZ, int key)
{
	C[0][0] = model->solve_left_dp();
	B[0][0] = model->solve_left_dp_beta();
	A[0][0] = 0.0;
	RightSide[0][0] = -model->solve_left();

	construction_bz(MZ, 2);
}

void Oil1DNITSolver::RightBoundAppr(int MZ, int key)
{
	C[0][0] = 0.0;
	B[0][0] = model->solve_right_dp_beta();
	A[0][0] = model->solve_right_dp();
	RightSide[0][0] = -model->solve_right();

	construction_bz(MZ, 1);
}

void Oil1DNITSolver::MiddleAppr(int current, int MZ, int key)
{
	C[0][0] = model->solve_eq_dp_beta(current, current-1);
	B[0][0] = model->solve_eq_dp(current);
	A[0][0] = model->solve_eq_dp_beta(current, current+1);
	RightSide[0][0] = -model->solve_eq(current);

	construction_bz(MZ, 2);
}