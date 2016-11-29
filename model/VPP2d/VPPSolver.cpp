#include "model/VPP2d/VPPSolver.hpp"

#include "adolc/drivers/drivers.h"

using namespace vpp2d;
using namespace std;

VPPSolver::VPPSolver(VPP2d* _model) : AbstractSolver<VPP2d>(_model)
{
	Initialize(model->cellsNum_r + 2, Variable::size * (model->cellsNum_z + 2));

	plot_P.open("snaps/P.dat", ofstream::out);
	plot_S.open("snaps/S.dat", ofstream::out);
	plot_C.open("snaps/C.dat", ofstream::out);
	plot_qcells.open("snaps/q_cells.dat", ofstream::out);

	t_dim = model->t_dim;

	n = model->Qcell.size();
	dq.Initialize(n);
	q.Initialize(n);

	dpdq.Initialize(n, n - 1);
	mat.Initialize(n - 1, n - 1);
	b.Initialize(n - 1);

	Tt = model->period[model->period.size() - 1];

	jac = new double*[Variable::size];
	for (int i = 0; i < Variable::size; i++)
		jac[i] = new double[15];
}
VPPSolver::~VPPSolver()
{
	plot_P.close();
	plot_S.close();
	plot_C.close();
	plot_qcells.close();

	for (int i = 0; i < Variable::size; i++)
		delete[] jac[i];
	delete[] jac;
}
void VPPSolver::writeData()
{
	double p = 0.0, s = 0.0, c = 0.0;

	plot_qcells << cur_t * t_dim / 3600.0;

	map<int, double>::iterator it;
	for (it = model->Qcell.begin(); it != model->Qcell.end(); ++it)
	{
		p += model->cells[it->first].u_next.p * model->P_dim;
		s += model->cells[it->first].u_next.s;
		c += model->cells[it->first].u_next.c;
		if (model->leftBoundIsRate)
			plot_qcells << "\t" << it->second * model->Q_dim * 86400.0;
		else
			plot_qcells << "\t" << model->getRate(it->first) * model->Q_dim * 86400.0;
	}

	plot_P << cur_t * t_dim / 3600.0 << "\t" << p / (double)(model->Qcell.size()) << endl;
	plot_S << cur_t * t_dim / 3600.0 << "\t" << s / (double)(model->Qcell.size()) << endl;
	plot_C << cur_t * t_dim / 3600.0 << "\t" << c / (double)(model->Qcell.size()) << endl;

	plot_qcells << endl;
}
void VPPSolver::control()
{
	writeData();

	if (cur_t >= model->period[curTimePeriod])
	{
		curTimePeriod++;
		model->ht = model->ht_min;
		model->setPeriod(curTimePeriod);
	}

	if (model->ht <= model->ht_max && iterations < 6)
		model->ht = model->ht * 1.5;
	else if (iterations > 6 && model->ht > model->ht_min)
		model->ht = model->ht / 1.5;

	if (cur_t + model->ht > model->period[curTimePeriod])
		model->ht = model->period[curTimePeriod] - cur_t;

	cur_t += model->ht;
}
void VPPSolver::doNextStep()
{
	solveStep();

	if (n > 1 && model->Q_sum > EQUALITY_TOLERANCE)
	{
		double H0 = fabs(model->solveH());
		if (H0 > 0.1)
		{
			printWellRates();
			fillq();

			double mult = 0.9;
			double H = H0;

			while (H > H0 / 50.0 || H > 0.05)
			{
				solveDq(mult);

				int i = 0;
				map<int, double>::iterator it = model->Qcell.begin();

				while (it != model->Qcell.end())
				{
					q[i] += mult * dq[i];
					it->second = q[i];
					i++;	++it;
				}

				solveStep();
				printWellRates();

				H = fabs(model->solveH());
			}
		}
	}
}
void VPPSolver::fillq()
{
	int i = 0;
	map<int, double>::iterator it = model->Qcell.begin();
	while (it != model->Qcell.end())
	{
		q[i++] = it->second;
		++it;
	}
}
void VPPSolver::fillDq()
{
	for (int i = 0; i < n; i++)
		dq[i] = 0.0;
}
void VPPSolver::solveDq(double mult)
{
	fillDq();
	filldPdQ(mult);
	solveStep();
	solveSystem();

	int i = 0;
	map<int, double>::iterator it;
	for (it = model->Qcell.begin(); it != model->Qcell.end(); ++it)
	{
		cout << "dq[" << it->first << "] = " << dq[i++] * model->Q_dim * 86400.0 << endl;
	}
	cout << endl;
}
void VPPSolver::solveSystem()
{
	double s = 0.0, p1, p2;
	map<int, double>::iterator it;

	for (int i = 0; i < n - 1; i++)
	{
		for (int j = 0; j < n - 1; j++)
		{
			s = 0.0;
			for (int k = 0; k < n - 1; k++)
				s += (dpdq[k + 1][j] - dpdq[k][j]) * (dpdq[k + 1][i] - dpdq[k][i]);
			mat[i][j] = s;
		}

		s = 0.0;
		it = model->Qcell.begin();
		for (int k = 0; k < n - 1; k++)
		{
			p1 = model->cells[it->first].u_next.p;
			p2 = model->cells[(++it)->first].u_next.p;
			s += (p2 - p1) * (dpdq[k + 1][i] - dpdq[k][i]);
		}
		b[i] = -s;
	}

	MC_LU rateSystem(mat, b);
	rateSystem.LU_Solve();

	s = 0.0;
	for (int i = 0; i < n - 1; i++)
	{
		dq[i + 1] = rateSystem.ptResult[i];
		s += rateSystem.ptResult[i];
	}
	dq[0] = -s;
}
void VPPSolver::filldPdQ(double mult)
{
	double p1, p2, ratio;
	ratio = mult * 0.001 / (double)(n);

	int i = 0, j = 0;
	map<int, double>::iterator it0 = model->Qcell.begin();
	map<int, double>::iterator it1 = model->Qcell.begin();
	map<int, double>::iterator it2 = it1;	++it2;
	while (it1 != model->Qcell.end())
	{
		j = 0;
		it2 = model->Qcell.begin();		++it2;
		while (it2 != model->Qcell.end())
		{
			model->setRateDeviation(it2->first, -ratio);
			model->setRateDeviation(it0->first, ratio);
			solveStep();
			p1 = model->cells[it1->first].u_next.p;

			model->setRateDeviation(it2->first, 2.0 * ratio);
			model->setRateDeviation(it0->first, -2.0 * ratio);
			solveStep();
			p2 = model->cells[it1->first].u_next.p;

			model->setRateDeviation(it2->first, -ratio);
			model->setRateDeviation(it0->first, ratio);

			dpdq[i][j++] = (p2 - p1) / (2.0 * ratio * model->Q_sum);

			++it2;
		}
		i++;
		++it1;
	}
}

void VPPSolver::solveStep()
{
	int cellIdx, varIdx;
	double err_newton = 1.0;
	double averPresPrev = averValue(0);
	double averSatPrev = averValue(1);
	double averСPrev = averValue(2);
	double averPres, averSat, averC;
	double dAverPres = 1.0, dAverSat = 1.0, dAverC = 1.0;

	iterations = 0;
	while (err_newton > 1.e-4 && (dAverSat > 1.e-10 || dAverPres > 1.e-10 || dAverC > 1.e-10) && iterations < 9)
	{
		copyIterLayer();

		Solve(model->cellsNum_r + 1, Variable::size * (model->cellsNum_z + 2), PRES);
		construction_from_fz(model->cellsNum_r + 2, Variable::size * (model->cellsNum_z + 2), PRES);
		model->solveFixVar();

		err_newton = convergance(cellIdx, varIdx);

		averPres = averValue(0);					averSat = averValue(1);					averC = averValue(2);
		dAverPres = fabs(averPres - averPresPrev);	dAverSat = fabs(averSat - averSatPrev);	dAverC = fabs(averC - averСPrev);
		averPresPrev = averPres;					averSatPrev = averSat;					averСPrev = averC;

		iterations++;
	}

	cout << "Newton Iterations = " << iterations << endl;
}
void VPPSolver::construction_from_fz(int N, int n, int key)
{
	vector<Cell>::iterator it;
	if (key == PRES)
	{
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < model->cellsNum_z + 2; j++)
			{
				Variable& var = model->cells[i * (model->cellsNum_z + 2) + j].u_next;
				var.p += fz[i][2 * j + 1];
				var.s += fz[i][2 * j + 2];
				var.c += fz[i][2 * j + 3];
			}
		}
	}
}

void VPPSolver::LeftBoundAppr(int MZ, int key)
{
	for (int i = 0; i < MZ; i++)
	{
		for (int j = 0; j < MZ; j++)
		{
			C[i][j] = 0.0;
			B[i][j] = 0.0;
			A[i][j] = 0.0;
		}
		
		C[i][i] = 1.0;
		B[i][i] = -1.0;
		A[i][i] = 0.0;
		RightSide[i][0] = model->cells[i / Variable::size].u_next.values[i % Variable::size] - 
						model->cells[int(i / Variable::size) + model->cellsNum_z + 2].u_next.values[i % Variable::size];
	}

	map<int, double>::iterator it;
	int idx;
	if (key == PRES)
	{
		for (it = model->Qcell.begin(); it != model->Qcell.end(); ++it)
		{
			idx = Variable::size * it->first;

			model->setVariables(it->first);

			model->solve_eqLeft(it->first);
			for (int i = 0; i < Variable::size; i++)
				RightSide[idx + i][0] = -model->y[i];

			jacobian(left, Variable::size, 6, model->x, jac);
			for (int i = 0; i < Variable::size; i++)
			{
				// i - equation index
				for (int j = 0; j < Variable::size; j++)
				{
					// j - variable index
					C[idx + i][idx + j] = jac[i][j];
					B[idx + i][idx + j] = jac[i][Variable::size + j];
				}
			}
		}
	}

	construction_bz(MZ, 2);
}
void VPPSolver::RightBoundAppr(int MZ, int key)
{
	for (int i = 0; i < MZ; i++)
	{
		for (int j = 0; j < MZ; j++)
		{
			C[i][j] = 0.0;
			B[i][j] = 0.0;
			A[i][j] = 0.0;
		}
		RightSide[i][0] = 0.0;
	}

	if (key == PRES)
	{
		int idx = 0;

		for (int i = (model->cellsNum_r + 1)*(model->cellsNum_z + 2); i < (model->cellsNum_r + 2)*(model->cellsNum_z + 2); i++)
		{
			model->setVariables(i);

			model->solve_eqRight(i);
			for (int i = 0; i < Variable::size; i++)
				RightSide[idx + i][0] = -model->y[i];

			jacobian(right, Variable::size, 9, model->x, jac);
			for (int i = 0; i < Variable::size; i++)
			{
				// i - equation index
				for (int j = 0; j < Variable::size; j++)
				{
					// j - variable index
					A[idx + i][idx + j] = jac[i][j];
					B[idx + i][idx + j] = jac[i][Variable::size + j];
					C[idx + i][idx + j] = jac[i][2 * Variable::size + j];
				}
			}

			idx += Variable::size;
		}
	}

	construction_bz(MZ, 1);
}
void VPPSolver::MiddleAppr(int current, int MZ, int key)
{
	for (int i = 0; i < MZ; i++)
	{
		for (int j = 0; j < MZ; j++)
		{
			A[i][j] = 0.0;
			B[i][j] = 0.0;
			C[i][j] = 0.0;
		}
		RightSide[i][0] = 0.0;
	}

	if (key == PRES)
	{
		int idx = 0;
		int i = current * (model->cellsNum_z + 2);

		// Top cell
		model->setVariables(i);

		model->solve_eqVertical(i);
		for (int i = 0; i < Variable::size; i++)
			RightSide[idx + i][0] = -model->y[i];

		jacobian(vertical, Variable::size, 6, model->x, jac);
		for (int i = 0; i < Variable::size; i++)
		{
			// i - equation index
			for (int j = 0; j < Variable::size; j++)
			{
				// j - variable index
				B[idx + i][idx + j] = jac[i][j];
				B[idx + i][idx + j + Variable::size] = jac[i][Variable::size + j];
			}
		}

		idx += Variable::size;

		// Middle cells
		for (i = current * (model->cellsNum_z + 2) + 1; i < (current + 1) * (model->cellsNum_z + 2) - 1; i++)
		{
			model->setVariables(i);

			model->solve_eqMiddle(i);
			for (int i = 0; i < Variable::size; i++)
				RightSide[idx + i][0] = -model->y[i];

			jacobian(mid, Variable::size, 15, model->x, jac);
			for (int i = 0; i < Variable::size; i++)
			{
				// i - equation index
				for (int j = 0; j < Variable::size; j++)
				{
					// j - variable index
					B[idx + i][idx + j] = jac[i][j];
					C[idx + i][idx + j] = jac[i][Variable::size + j];
					A[idx + i][idx + j] = jac[i][Variable::size * 2 + j];
					B[idx + i][idx + j - Variable::size] = jac[i][Variable::size * 3 + j];
					B[idx + i][idx + j + Variable::size] = jac[i][Variable::size * 4 + j];
				}
			}

			idx += Variable::size;
		}

		// Bottom cell
		model->setVariables(i);

		model->solve_eqVertical(i);
		for (int i = 0; i < Variable::size; i++)
			RightSide[idx + i][0] = -model->y[i];

		jacobian(vertical, Variable::size, 6, model->x, jac);
		for (int i = 0; i < Variable::size; i++)
		{
			// i - equation index
			for (int j = 0; j < Variable::size; j++)
			{
				// j - variable index
				B[idx + i][idx + j] = jac[i][j];
				B[idx + i][idx + j - Variable::size] = jac[i][Variable::size + j];
			}
		}
	}

	construction_bz(MZ, 2);
}
