#include "model/Gas1D/Gas1DSolver.h"

using namespace gas1D;
using std::ofstream;
using std::cout;
using std::endl;

Gas1DSolver::Gas1DSolver(Gas1D* _model) : AbstractSolver<Gas1D>(_model)
{
	Initialize(model->cellsNum_r+2, 1);

	plot_Pdyn.open("snaps/P_dyn.dat", ofstream::out);

	t_dim = model->t_dim;

	Tt = model->period[model->period.size()-1];
}

Gas1DSolver::~Gas1DSolver()
{
	plot_Pdyn.close();
}

void Gas1DSolver::writeData()
{
	plot_Pdyn << cur_t * t_dim / 3600.0 << "\t" << model->cells[idx1].u_next.p << endl;
}

void Gas1DSolver::fill()
{
}

void Gas1DSolver::control()
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

void Gas1DSolver::start()
{
	int counter = 0;
	iterations = 8;

	model->setPeriod(curTimePeriod);
	while(cur_t < Tt)
	{
		control();
		model->snapshot(counter++);
		doNextStep();
		copyTimeLayer();
	}
}

void Gas1DSolver::doNextStep()
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

	cout << "Newton Iterations = " << iterations << endl;
}

void Gas1DSolver::construct_solution()
{
}

void Gas1DSolver::construction_from_fz(int N, int n, int key)
{
	vector<Cell>::iterator it;
	if (key == PRES)
		for(int i = 0; i < N; i++)
			 model->cells[i].u_next.p += fz[i][1];
}

void Gas1DSolver::LeftBoundAppr(int MZ, int key)
{
	C[0][0] = model->solve_eqLeft_dp();
	B[0][0] = model->solve_eqLeft_dp_beta();
	A[0][0] = 0.0;
	RightSide[0][0] = -model->solve_eqLeft();

	construction_bz(MZ, 2);
}

void Gas1DSolver::RightBoundAppr(int MZ, int key)
{
	C[0][0] = 0.0;
	B[0][0] = model->solve_eqRight_dp_beta();
	A[0][0] = model->solve_eqRight_dp();
	RightSide[0][0] = -model->solve_eqRight();

	construction_bz(MZ, 1);
}

void Gas1DSolver::MiddleAppr(int current, int MZ, int key)
{
	C[0][0] = model->solve_eq_dp_beta(current, current-1);
	B[0][0] = model->solve_eq_dp(current);
	A[0][0] = model->solve_eq_dp_beta(current, current+1);
	RightSide[0][0] = -model->solve_eq(current);

	construction_bz(MZ, 2);
}