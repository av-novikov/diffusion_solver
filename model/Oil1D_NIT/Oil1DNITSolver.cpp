#include "model/Oil1D_NIT/Oil1DNITSolver.h"

using namespace std;
using namespace oil1D_NIT;

Oil1DNITSolver::Oil1DNITSolver(Oil1D_NIT* _model) : AbstractSolver<Oil1D_NIT>(_model)
{
	Initialize(model->cellsNum_r+2, 1);

	plot_Tdyn.open("snaps/T_dyn.dat", ofstream::out);
	plot_Pdyn.open("snaps/P_dyn.dat", ofstream::out);
	plot_qcells.open("snaps/q_cells.dat", ofstream::out);

	t_dim = model->t_dim;
	T_dim = model->T_dim;

	Tt = model->period[model->period.size()-1];
}

Oil1DNITSolver::~Oil1DNITSolver()
{
	plot_Tdyn.close();
	plot_Pdyn.close();
	plot_qcells.close();
}

void Oil1DNITSolver::writeData()
{
	plot_Pdyn << cur_t * t_dim / 3600.0 << "\t" << model->cells[idx1].u_next.p << endl;
	plot_Tdyn << cur_t * t_dim / 3600.0 << "\t" << model->cells[idx1].u_next.t * T_dim << endl;
	plot_qcells << cur_t * t_dim;
	if( model->leftBoundIsRate )
		plot_qcells << "\t" << model->Qcell[0] * model->Q_dim * 86400.0 << endl;
	else
		plot_qcells << "\t" << model->getRate() * model->Q_dim * 86400.0 << endl;
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

	Solve(model->cellsNum_r+1, 1, TEMP);
	construction_from_fz(model->cellsNum_r+2, 1, TEMP);
}

void Oil1DNITSolver::construction_from_fz(int N, int n, int key)
{
	vector<Cell>::iterator it;
	if (key == PRES)
		for(int i = 0; i < N; i++)
			model->cells[i].u_next.p += fz[i][1];
	else if(key == TEMP)
		for(int i = 0; i < N; i++)
			model->cells[i].u_next.t += fz[i][1];
}

void Oil1DNITSolver::LeftBoundAppr(int MZ, int key)
{
	if(key == PRES)
	{
		C[0][0] = model->solve_eqLeft_dp();
		B[0][0] = model->solve_eqLeft_dp_beta();
		A[0][0] = 0.0;
		RightSide[0][0] = -model->solve_eqLeft();
	} 
	else if(key == TEMP)
	{
		Cell& curr = model->cells[0];
		Cell& nebr = model->cells[1];
			
		C[0][0] = 1.0;
		B[0][0] = (curr.r - nebr.r) / (model->cells[2].r - nebr.r) - 1.0;
		A[0][0] = -(curr.r - nebr.r) / (model->cells[2].r - nebr.r);
		RightSide[0][0] = 0.0;
	}

	construction_bz(MZ, 2);
}

void Oil1DNITSolver::RightBoundAppr(int MZ, int key)
{
	if( key == PRES )
	{
		C[0][0] = 0.0;
		B[0][0] = model->solve_eqRight_dp_beta();
		A[0][0] = model->solve_eqRight_dp();
		RightSide[0][0] = -model->solve_eqRight();
	}
	else if(key == TEMP)
	{
		C[0][0] = 0.0;
		B[0][0] = 0.0;
		A[0][0] = 1.0;
		RightSide[0][0] = model->props_sk[0].t_init;
	}

	construction_bz(MZ, 1);
}

void Oil1DNITSolver::MiddleAppr(int current, int MZ, int key)
{
	if(key == PRES)
	{
		C[0][0] = model->solve_eq_dp_beta(current, current-1);
		B[0][0] = model->solve_eq_dp(current);
		A[0][0] = model->solve_eq_dp_beta(current, current+1);
		RightSide[0][0] = -model->solve_eq(current);
	}
	else if(key == TEMP)
	{
		Cell& cell_prev = model->cells[current-1];
		Cell& cell = model->cells[current];
		Cell& cell_next = model->cells[current+1];

		C[0][0] = -2.0 * ( max(model->getA(cell, NEXT), 0.0) +
							model->getLambda(cell, cell_prev) * (cell.r - cell.hr / 2.0) / cell.r / cell.hr ) / (cell.hr + cell_prev.hr);
		A[0][0] = 2.0 * ( min(model->getA(cell, NEXT), 0.0) -
							model->getLambda(cell, cell_next) * (cell.r + cell.hr / 2.0) / cell.r / cell.hr ) / (cell.hr + cell_next.hr);
		B[0][0] = model->getCn(cell) / model->ht - C[0][0] - A[0][0];
		
		RightSide[0][0] = model->getAd(cell) * (cell.u_next.p - cell.u_prev.p) / model->ht -
							model->getJT(cell, NEXT) * model->getNablaP(cell, NEXT);
	}
	construction_bz(MZ, 2);
}