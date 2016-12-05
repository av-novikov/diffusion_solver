#include "model/GasOil_RZ/GasOil2DSolver.h"

#include <cassert>

#include "adolc/drivers/drivers.h"
#include "adolc/adolc.h"  

using namespace std;
using namespace gasOil_rz;

GasOil2DSolver::GasOil2DSolver(GasOil_RZ* _model) : AbstractSolver<GasOil_RZ>(_model)
{
	Initialize(model->cellsNum_r+2, 2*(model->cellsNum_z+2));

	plot_Pdyn.open("snaps/P_dyn.dat", ofstream::out);
	plot_Sdyn.open("snaps/S_dyn.dat", ofstream::out);
	plot_qcells.open("snaps/q_cells.dat", ofstream::out);

	t_dim = model->t_dim;

	n = model->Qcell.size();
	dq.Initialize(n);
	q.Initialize(n);

	dpdq.Initialize(n, n-1);
	mat.Initialize(n-1, n-1);
	b.Initialize(n-1);

	Tt = model->period[model->period.size()-1];

	jac = new double*[Variable::size-1];
	for (int i = 0; i < Variable::size-1; i++)
		jac[i] = new double[stencil * Variable::size];
}
GasOil2DSolver::~GasOil2DSolver()
{
	plot_Pdyn.close();
	plot_Sdyn.close();
	plot_qcells.close();

	for (int i = 0; i < Variable::size-1; i++)
		delete[] jac[i];
	delete[] jac;
}
void GasOil2DSolver::writeData()
{
	double p = 0.0, s = 0.0;

	plot_qcells << cur_t * t_dim / 3600.0;

	map<int,double>::iterator it;
	for(it = model->Qcell.begin(); it != model->Qcell.end(); ++it)
	{
		p += model->cells[it->first].u_next.p * model->P_dim;
		s += model->cells[it->first].u_next.s;
		if( model->leftBoundIsRate )
			plot_qcells << "\t" << it->second * model->Q_dim * 86400.0;
		else
			plot_qcells << "\t" << model->getRate(it->first) * model->Q_dim * 86400.0;
	}

	plot_Pdyn << cur_t * t_dim / 3600.0 << "\t" << p / (double)(model->Qcell.size()) << endl;
	plot_Sdyn << cur_t * t_dim / 3600.0 << "\t" << s / (double)(model->Qcell.size()) << endl;

	plot_qcells << endl;
}
void GasOil2DSolver::control()
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
void GasOil2DSolver::doNextStep()
{
	solveStep();

	if(n > 1 && model->Q_sum > EQUALITY_TOLERANCE)
	{
		double H0 = fabs( model->solveH() );
		if( H0 > 0.1 )
		{
			printWellRates();
			fillq();

			double mult = 0.9;
			double H = H0;			

			while(H > H0 / 50.0 || H > 0.05)
			{
				solveDq(mult);

				int i = 0;
				map<int,double>::iterator it = model->Qcell.begin();

				while(it != model->Qcell.end())
				{
					q[i] += mult * dq[i];
					it->second = q[i];
					i++;	++it;
				}
				
				solveStep();
				printWellRates();

				H = fabs( model->solveH() );
			}
		}
	}
}
void GasOil2DSolver::fillq()
{
	int i = 0;
	map<int,double>::iterator it = model->Qcell.begin();
	while(it != model->Qcell.end())
	{
		q[i++] = it->second;
		++it;
	}
}
void GasOil2DSolver::fillDq()
{
	for(int i = 0; i < n; i++)
		dq[i] = 0.0;
}
void GasOil2DSolver::solveDq(double mult)
{
	fillDq();
	filldPdQ(mult);
	solveStep();
	solveSystem();

	int i = 0;
	map<int,double>::iterator it;
	for(it = model->Qcell.begin(); it != model->Qcell.end(); ++it)
	{
		cout << "dq[" << it->first << "] = " << dq[i++] * model->Q_dim * 86400.0 << endl;
	}
	cout << endl;
}
void GasOil2DSolver::solveSystem()
{
	double s = 0.0, p1, p2;
	map<int,double>::iterator it;

	for(int i = 0; i < n-1; i++)
	{
		for(int j = 0; j < n-1; j++)
		{
			s = 0.0;
			for(int k = 0; k < n-1; k++)
				s += ( dpdq[k+1][j] - dpdq[k][j] ) * ( dpdq[k+1][i] - dpdq[k][i] );
			mat[i][j] = s;
		}
		
		s = 0.0;
		it = model->Qcell.begin();
		for(int k = 0; k < n-1; k++)
		{
			p1 = model->cells[ it->first ].u_next.p;
			p2 = model->cells[ (++it)->first ].u_next.p;
			s += ( p2 - p1 ) * ( dpdq[k+1][i] - dpdq[k][i] );
		}
		b[i] = -s;
	}

	MC_LU rateSystem(mat, b);
	rateSystem.LU_Solve();

	s = 0.0;
	for(int i = 0; i < n-1; i++)
	{
		dq[i+1] = rateSystem.ptResult[i];
		s += rateSystem.ptResult[i];
	}
	dq[0] = - s;
}
void GasOil2DSolver::filldPdQ(double mult)
{
	double p1, p2, ratio;
	ratio = mult * 0.001 / (double)(n);
	
	int i = 0, j = 0;
	map<int,double>::iterator it0 = model->Qcell.begin();
	map<int,double>::iterator it1 = model->Qcell.begin();
	map<int,double>::iterator it2 = it1;	++it2;
	while(it1 != model->Qcell.end())
	{
		j = 0;
		it2 = model->Qcell.begin();		++it2;
		while(it2 != model->Qcell.end())
		{
			model->setRateDeviation(it2->first, -ratio);
			model->setRateDeviation(it0->first, ratio);
			solveStep();
			p1 = model->cells[ it1->first ].u_next.p;

			model->setRateDeviation(it2->first, 2.0 * ratio);
			model->setRateDeviation(it0->first, -2.0 * ratio);
			solveStep();
			p2 = model->cells[ it1->first ].u_next.p;

			model->setRateDeviation(it2->first, -ratio);
			model->setRateDeviation(it0->first, ratio);

			dpdq[i][j++] = (p2 - p1) / ( 2.0 * ratio * model->Q_sum);

			++it2;
		}
		i++;
		++it1;
	}
}

void GasOil2DSolver::solveStep()
{
	int cellIdx, varIdx;
	double err_newton = 1.0;
	double averPresPrev = averValue(0);
	double averSatPrev = averValue(1);
	double averPres, averSat;
	double dAverPres = 1.0, dAverSat = 1.0;
	
	iterations = 0;
	while( err_newton > 1.e-4 && ( dAverSat > 1.e-10 || dAverPres > 1.e-10) && iterations < 9 )
	{	
		copyIterLayer();

		Solve(model->cellsNum_r+1, 2*(model->cellsNum_z+2), PRES);
		construction_from_fz(model->cellsNum_r+2, 2*(model->cellsNum_z+2), PRES);
		model->solveP_bub();
 
		err_newton = convergance(cellIdx, varIdx);

		averPres = averValue(0);					averSat = averValue(1);
		dAverPres = fabs(averPres - averPresPrev);	dAverSat = fabs(averSat - averSatPrev);
		averPresPrev = averPres;					averSatPrev = averSat;

		iterations++;
	}

	cout << "Newton Iterations = " << iterations << endl;
}
void GasOil2DSolver::construction_from_fz(int N, int n, int key)
{
	vector<Cell>::iterator it;
	if (key == PRES)
	{
		for(int i = 0; i < N; i++)
		{
			for(int j = 0; j < model->cellsNum_z+2; j++)
			{
				Var2phase& var = model->cells[i*(model->cellsNum_z+2) + j].u_next;
				var.p += fz[i][2*j+1];
				if(var.SATUR)
				{
					var.s += fz[i][2 * j + 2];
				}
				else
				{
					var.p_bub += fz[i][2 * j + 2];
				}
			}
		}
	}
}

void GasOil2DSolver::LeftBoundAppr(int MZ, int key)
{
	for(int i = 0; i < MZ; i++)
	{
		for (int j = 0 ; j < MZ; j++)
		{
			C[i][j] = 0.0;
			B[i][j] = 0.0;
			A[i][j] = 0.0;
		}
		C[i][i] = 1.0;

		const Variable next = model->cells[int(i / (Variable::size - 1))].u_next;
		const Variable nebr = model->cells[int(i / (Variable::size - 1)) + model->cellsNum_z + 2].u_next;

		if (i % (Variable::size - 1) == 1)
		{
			if (nebr.SATUR == next.SATUR)
				B[i][i] = -1.0;
			
			RightSide[i][0] = -next.values[1 + !next.SATUR] + nebr.values[1 + !nebr.SATUR];
		} else
		{
			B[i][i] = -1.0;
			RightSide[i][0] = -next.p + nebr.p;
		}
	}

	map<int,double>::iterator it;
	int idx;
	if(key == PRES)
	{
		for(it = model->Qcell.begin(); it != model->Qcell.end(); ++it)
		{
			idx = (Variable::size-1) * it->first;

			model->setVariables(it->first);

			model->solve_eqLeft(it->first);
			for (int i = 0; i < Variable::size-1; i++)
				RightSide[idx + i][0] = -model->y[i];

			jacobian(left, Variable::size-1, Variable::size * Lstencil, model->x, jac);
			for (int i = 0; i < Variable::size-1; i++)
			{
				// i - equation index
				//if(next.SATUR)
				C[idx + i][idx] = jac[i][0];
				B[idx + i][idx] = jac[i][Variable::size];
				A[idx + i][idx] = jac[i][2 * Variable::size];

				if(model->cells[it->first].u_next.SATUR)
					C[idx + i][idx + 1] = jac[i][1];
				else
					C[idx + i][idx + 1] = jac[i][2];
				if(model->cells[it->first + model->cellsNum_z + 2].u_next.SATUR)
					B[idx + i][idx + 1] = jac[i][Variable::size + 1];
				else
					B[idx + i][idx + 1] = jac[i][Variable::size + 2];
				if(model->cells[it->first + 2 * model->cellsNum_z + 4].u_next.SATUR)
					A[idx + i][idx + 1] = jac[i][2 * Variable::size + 1];
				else
					A[idx + i][idx + 1] = jac[i][2 * Variable::size + 2];
			}
		}
	}

	construction_bz(MZ, 2);
}
void GasOil2DSolver::RightBoundAppr(int MZ, int key)
{
	for(int i = 0; i < MZ; i++)
	{
		for (int j = 0 ; j < MZ; j++)
		{
			C[i][j] = 0.0;
			B[i][j] = 0.0;
			A[i][j] = 0.0;
		}
		RightSide[i][0] = 0.0;
	}

	if(key == PRES)
	{
		int idx = 0;
		
		for(int cell_idx = (model->cellsNum_r+1)*(model->cellsNum_z+2); cell_idx < (model->cellsNum_r+2)*(model->cellsNum_z+2); cell_idx++)
		{
			model->setVariables(cell_idx);

			model->solve_eqRight(cell_idx);
			for (int i = 0; i < Variable::size-1; i++)
				RightSide[idx + i][0] = -model->y[i];


			jacobian(right, Variable::size-1, Variable::size * Rstencil, model->x, jac);
			for (int i = 0; i < Variable::size-1; i++)
			{
				// i - equation index
				A[idx + i][idx] = jac[i][0];
				B[idx + i][idx] = jac[i][Variable::size];

				if (model->cells[cell_idx].u_next.SATUR)
					A[idx + i][idx + 1] = jac[i][1];
				else
					A[idx + i][idx + 1] = jac[i][2];
				if (model->cells[cell_idx - model->cellsNum_z - 2].u_next.SATUR)
					B[idx + i][idx + 1] = jac[i][Variable::size + 1];
				else
					B[idx + i][idx + 1] = jac[i][Variable::size + 2];
			}

			idx += Variable::size-1;
		}
	}

	construction_bz(MZ,1);
}
void GasOil2DSolver::MiddleAppr(int current, int MZ, int key)
{
	for(int i = 0; i < MZ; i++)
	{
		for(int j = 0; j < MZ; j++)
		{
			A[i][j] = 0.0;
			B[i][j] = 0.0;
			C[i][j] = 0.0;
		}
		RightSide[i][0] = 0.0;
	}

	if(key == PRES)
	{
		int idx = 0;
		int cell_idx = current * (model->cellsNum_z+2);
		
		// Top cell
		model->setVariables(cell_idx);

		model->solve_eqVertical(cell_idx);
		for (int i = 0; i < Variable::size-1; i++)
			RightSide[idx + i][0] = -model->y[i];

		jacobian(vertical, Variable::size-1, Variable::size * Vstencil, model->x, jac);
		for (int i = 0; i < Variable::size-1; i++)
		{
			// i - equation index
			B[idx + i][idx] = jac[i][0];
			B[idx + i][idx + Variable::size - 1] = jac[i][Variable::size];

			if(model->cells[cell_idx].u_next.SATUR)
				B[idx + i][idx + 1] = jac[i][1];
			else
				B[idx + i][idx + 1] = jac[i][2];
			if (model->cells[cell_idx + 1].u_next.SATUR)
				B[idx + i][idx + 1 + Variable::size - 1] = jac[i][Variable::size + 1];
			else
				B[idx + i][idx + 1 + Variable::size - 1] = jac[i][Variable::size + 2];
		}

		idx += Variable::size-1;

		// Middle cells
		for(cell_idx = current * (model->cellsNum_z+2) + 1; cell_idx < (current+1) * (model->cellsNum_z+2) - 1; cell_idx++)
		{
			model->setVariables(cell_idx);

			model->solve_eqMiddle(cell_idx);
			for (int i = 0; i < Variable::size-1; i++)
				RightSide[idx + i][0] = -model->y[i];

			jacobian(mid, Variable::size-1, Variable::size * stencil, model->x, jac);
			for (int i = 0; i < Variable::size-1; i++)
			{
				// i - equation index
				B[idx + i][idx] = jac[i][0];
				C[idx + i][idx] = jac[i][Variable::size];
				A[idx + i][idx] = jac[i][Variable::size * 2];
				B[idx + i][idx - Variable::size + 1] = jac[i][Variable::size * 3];
				B[idx + i][idx + Variable::size - 1] = jac[i][Variable::size * 4];

				if (model->cells[cell_idx].u_next.SATUR)
					B[idx + i][idx + 1] = jac[i][1];
				else
					B[idx + i][idx + 1] = jac[i][2];
				if (model->cells[cell_idx - model->cellsNum_z - 2].u_next.SATUR)
					C[idx + i][idx + 1] = jac[i][Variable::size + 1];
				else 
					C[idx + i][idx + 1] = jac[i][Variable::size + 2];
				if (model->cells[cell_idx + model->cellsNum_z + 2].u_next.SATUR)
					A[idx + i][idx + 1] = jac[i][Variable::size * 2 + 1];
				else
					A[idx + i][idx + 1] = jac[i][Variable::size * 2 + 2];
				if (model->cells[cell_idx - 1].u_next.SATUR)
					B[idx + i][idx + 1 - Variable::size + 1] = jac[i][Variable::size * 3 + 1];
				else 
					B[idx + i][idx + 1 - Variable::size + 1] = jac[i][Variable::size * 3 + 2];
				if (model->cells[cell_idx + 1].u_next.SATUR)
					B[idx + i][idx + 1 + Variable::size - 1] = jac[i][Variable::size * 4 + 1];
				else
					B[idx + i][idx + 1 + Variable::size - 1] = jac[i][Variable::size * 4 + 2];
			}

			idx += Variable::size - 1;
		}

		// Bottom cell
		model->setVariables(cell_idx);

		model->solve_eqVertical(cell_idx);
		for (int i = 0; i < Variable::size - 1; i++)
			RightSide[idx + i][0] = -model->y[i];

		jacobian(vertical, Variable::size-1, Variable::size * Vstencil, model->x, jac);
		for (int i = 0; i < Variable::size-1; i++)
		{
			// i - equation index
			B[idx + i][idx] = jac[i][0];
			B[idx + i][idx - Variable::size + 1] = jac[i][Variable::size];

			if (model->cells[cell_idx].u_next.SATUR)
				B[idx + i][idx + 1] = jac[i][1];
			else
				B[idx + i][idx + 1] = jac[i][2];
			if (model->cells[cell_idx - 1].u_next.SATUR)
				B[idx + i][idx + 1 - Variable::size + 1] = jac[i][Variable::size + 1];
			else
				B[idx + i][idx + 1 - Variable::size + 1] = jac[i][Variable::size + 2];
		}
	}

	construction_bz(MZ,2);
}