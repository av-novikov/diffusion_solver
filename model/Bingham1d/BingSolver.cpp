#include "model/Bingham1d/BingSolver.hpp"

using namespace std;
using namespace bing1d;

Bing1dSolver::Bing1dSolver(Bingham1d* _model) : AbstractSolver<Bingham1d>(_model)
{
	Initialize(model->cellsNum_r + 2, Variable::size);

	plot_Pdyn.open("snaps/P_dyn.dat", ofstream::out);
	plot_qcells.open("snaps/q_cells.dat", ofstream::out);

	t_dim = model->t_dim;
	Tt = model->period[model->period.size() - 1];

	jac = new double*[Variable::size - 1];
	for (int i = 0; i < Variable::size - 1; i++)
		jac[i] = new double[stencil * Variable::size];
}
Bing1dSolver::~Bing1dSolver()
{
	plot_Pdyn.close();
	plot_qcells.close();

	for (int i = 0; i < Variable::size - 1; i++)
		delete[] jac[i];
	delete[] jac;
}
void Bing1dSolver::writeData()
{
	double p = 0.0;

	plot_qcells << cur_t * t_dim / 3600.0;

	map<int, double>::iterator it;
	for (it = model->Qcell.begin(); it != model->Qcell.end(); ++it)
	{
		p += model->cells[it->first].u_next.p * model->P_dim;
		if (model->leftBoundIsRate)
			plot_qcells << "\t" << it->second * model->Q_dim * 86400.0;
		else
			plot_qcells << "\t" << model->getRate(it->first) * model->Q_dim * 86400.0;
	}

	plot_Pdyn << cur_t * t_dim / 3600.0 << "\t" << p / (double)(model->Qcell.size()) << endl;
	plot_qcells << endl;
}
void Bing1dSolver::control()
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
void Bing1dSolver::doNextStep()
{
	solveStep();
}
void Bing1dSolver::solveStep()
{
	int cellIdx, varIdx;
	double err_newton = 1.0;
	double averPresPrev = averValue(0);
	double averPres;
	double dAverPres = 1.0;

	iterations = 0;
	while (err_newton > 1.e-4 && (dAverPres > 1.e-10) && iterations < 9)
	{
		copyIterLayer();

		Solve(model->cellsNum_r + 1, Variable::size, PRES);
		construction_from_fz(model->cellsNum_r + 2, Variable::size, PRES);


		err_newton = convergance(cellIdx, varIdx);

		averPres = averValue(0);
		dAverPres = fabs(averPres - averPresPrev);
		averPresPrev = averPres;

		iterations++;
	}

	cout << "Newton Iterations = " << iterations << endl;
}
void Bing1dSolver::construction_from_fz(int N, int n, int key)
{
	vector<Cell>::iterator it;
	if (key == PRES)
	{
		for (int i = 0; i < N; i++)
			model->cells[i].u_next.p += fz[i][1];
	}
}

void Bing1dSolver::LeftBoundAppr(int MZ, int key)
{
/*	for (int i = 0; i < MZ; i++)
	{
		for (int j = 0; j < MZ; j++)
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
		}
		else
		{
			B[i][i] = -1.0;
			RightSide[i][0] = -next.p + nebr.p;
		}
	}

	map<int, double>::iterator it;
	int idx;
	if (key == PRES)
	{
		for (it = model->Qcell.begin(); it != model->Qcell.end(); ++it)
		{
			idx = (Variable::size - 1) * it->first;

			model->setVariables(it->first);

			model->solve_eqLeft(it->first);
			for (int i = 0; i < Variable::size - 1; i++)
				RightSide[idx + i][0] = -model->y[i];

			jacobian(left, Variable::size - 1, Variable::size * Lstencil, model->x, jac);
			for (int i = 0; i < Variable::size - 1; i++)
			{
				// i - equation index
				//if(next.SATUR)
				C[idx + i][idx] = jac[i][0];
				B[idx + i][idx] = jac[i][Variable::size];
				A[idx + i][idx] = jac[i][2 * Variable::size];

				if (model->cells[it->first].u_next.SATUR)
					C[idx + i][idx + 1] = jac[i][1];
				else
					C[idx + i][idx + 1] = jac[i][2];
				if (model->cells[it->first + model->cellsNum_z + 2].u_next.SATUR)
					B[idx + i][idx + 1] = jac[i][Variable::size + 1];
				else
					B[idx + i][idx + 1] = jac[i][Variable::size + 2];
				if (model->cells[it->first + 2 * model->cellsNum_z + 4].u_next.SATUR)
					A[idx + i][idx + 1] = jac[i][2 * Variable::size + 1];
				else
					A[idx + i][idx + 1] = jac[i][2 * Variable::size + 2];
			}
		}
	}

	construction_bz(MZ, 2);*/
}
void Bing1dSolver::RightBoundAppr(int MZ, int key)
{
/*	for (int i = 0; i < MZ; i++)
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

		for (int cell_idx = (model->cellsNum_r + 1)*(model->cellsNum_z + 2); cell_idx < (model->cellsNum_r + 2)*(model->cellsNum_z + 2); cell_idx++)
		{
			model->setVariables(cell_idx);

			model->solve_eqRight(cell_idx);
			for (int i = 0; i < Variable::size - 1; i++)
				RightSide[idx + i][0] = -model->y[i];


			jacobian(right, Variable::size - 1, Variable::size * Rstencil, model->x, jac);
			for (int i = 0; i < Variable::size - 1; i++)
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

			idx += Variable::size - 1;
		}
	}

	construction_bz(MZ, 1);*/
}
void Bing1dSolver::MiddleAppr(int current, int MZ, int key)
{
/*	for (int i = 0; i < MZ; i++)
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
		int cell_idx = current * (model->cellsNum_z + 2);

		// Top cell
		model->setVariables(cell_idx);

		model->solve_eqVertical(cell_idx);
		for (int i = 0; i < Variable::size - 1; i++)
			RightSide[idx + i][0] = -model->y[i];

		jacobian(vertical, Variable::size - 1, Variable::size * Vstencil, model->x, jac);
		for (int i = 0; i < Variable::size - 1; i++)
		{
			// i - equation index
			B[idx + i][idx] = jac[i][0];
			B[idx + i][idx + Variable::size - 1] = jac[i][Variable::size];

			if (model->cells[cell_idx].u_next.SATUR)
				B[idx + i][idx + 1] = jac[i][1];
			else
				B[idx + i][idx + 1] = jac[i][2];
			if (model->cells[cell_idx + 1].u_next.SATUR)
				B[idx + i][idx + 1 + Variable::size - 1] = jac[i][Variable::size + 1];
			else
				B[idx + i][idx + 1 + Variable::size - 1] = jac[i][Variable::size + 2];
		}

		idx += Variable::size - 1;

		// Middle cells
		for (cell_idx = current * (model->cellsNum_z + 2) + 1; cell_idx < (current + 1) * (model->cellsNum_z + 2) - 1; cell_idx++)
		{
			model->setVariables(cell_idx);

			model->solve_eqMiddle(cell_idx);
			for (int i = 0; i < Variable::size - 1; i++)
				RightSide[idx + i][0] = -model->y[i];

			jacobian(mid, Variable::size - 1, Variable::size * stencil, model->x, jac);
			for (int i = 0; i < Variable::size - 1; i++)
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

		jacobian(vertical, Variable::size - 1, Variable::size * Vstencil, model->x, jac);
		for (int i = 0; i < Variable::size - 1; i++)
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

	construction_bz(MZ, 2);*/
}
