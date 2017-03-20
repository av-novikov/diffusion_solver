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

	double averVal[3], averValPrev[3], dAverVal[3];
	averValPrev[0] = averValue(0);	averValPrev[1] = averValue(1);	averValPrev[2] = averValue(2);
	dAverVal[0] = 1.0;	dAverVal[1] = 1.0;	averVal[2] = 1.0;

	iterations = 0;
	while (err_newton > 1.e-4 && (dAverVal[0] > 1.e-10 || dAverVal[1] > 1.e-10 || dAverVal[2] > 1.e-10) && iterations < 9)
	{
		copyIterLayer();

		Solve(model->cellsNum_r + 1, BlackOil_RZ::var_size * (model->cellsNum_z + 2), PRES);
		construction_from_fz(model->cellsNum_r + 2, BlackOil_RZ::var_size * (model->cellsNum_z + 2), PRES);
		model->solveP_bub();

		err_newton = convergance(cellIdx, varIdx);

		averVal[0] = averValue(0);	averVal[1] = averValue(1);	averVal[2] = averValue(2);
		dAverVal[0] = fabs(averVal[0] - averValPrev[0]);
		dAverVal[1] = fabs(averVal[1] - averValPrev[1]);
		dAverVal[2] = fabs(averVal[2] - averValPrev[2]);
		averValPrev[0] = averVal[0];	averValPrev[1] = averVal[1];	averValPrev[2] = averVal[2];

		iterations++;
	}

	cout << "Newton Iterations = " << iterations << endl;
}
void BlackOil2dSolver::construction_from_fz(int N, int n, int key)
{
	vector<Cell>::iterator it;
	if (key == PRES)
	{
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < model->cellsNum_z + 2; j++)
			{
				Variable& next = model->cells[i*(model->cellsNum_z + 2) + j].u_next;
				next.p += fz[i][2 * j + 1];
				next.s_w += fz[i][2 * j + 2];
				if (next.SATUR)
				{
					next.s_o += fz[i][2 * j + 3];
					next.p_bub = next.p;
				}
				else
				{
					next.s_o -= fz[i][2 * j + 2];
					next.p_bub += fz[i][2 * j + 3];
				}
			}
		}
	}
}
void BlackOil2dSolver::MiddleAppr(int current, int MZ, int key)
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
		int cell_idx = current * (model->cellsNum_z + 2);

		// Top cell
		model->setVariables(model->cells[cell_idx]);

		for (int i = 0; i < var_size; i++)
		{
			RightSide[idx + i][0] = -model->y[i];

			// p
			B[idx + i][idx] = model->jac[i][0];
			B[idx + i][idx + var_size] = model->jac[i][size];
			// s_w
			B[idx + i][idx + 1] = model->jac[i][1];
			B[idx + i][idx + 1 + var_size] = model->jac[i][1 + size];
			// s_o / p_bub
			if (model->cells[cell_idx].u_next.SATUR)
				B[idx + i][idx + 2] = model->jac[i][2];
			else
				B[idx + i][idx + 2] = model->jac[i][3];
			if (model->cells[cell_idx + 1].u_next.SATUR)
				B[idx + i][idx + 2 + var_size] = model->jac[i][size + 2];
			else
				B[idx + i][idx + 2 + var_size] = model->jac[i][size + 3];
		}
		idx += var_size;

		// Middle cells
		for (cell_idx = current * (model->cellsNum_z + 2) + 1; cell_idx < (current + 1) * (model->cellsNum_z + 2) - 1; cell_idx++)
		{
			model->setVariables(model->cells[cell_idx]);

			for (int i = 0; i < var_size; i++)
			{
				RightSide[idx + i][0] = -model->y[i];

				// p
				B[idx + i][idx] = model->jac[i][0];
				C[idx + i][idx] = model->jac[i][size];
				A[idx + i][idx] = model->jac[i][size * 2];
				B[idx + i][idx - var_size] = model->jac[i][size * 3];
				B[idx + i][idx + var_size] = model->jac[i][size * 4];
				// s_w
				B[idx + i][idx + 1] = model->jac[i][1];
				C[idx + i][idx + 1] = model->jac[i][size + 1];
				A[idx + i][idx + 1] = model->jac[i][size * 2 + 1];
				B[idx + i][idx + 1 - var_size] = model->jac[i][size * 3 + 1];
				B[idx + i][idx + 1 + var_size] = model->jac[i][size * 4 + 1];
				// s_o / p_bub
				if (model->cells[cell_idx].u_next.SATUR)
					B[idx + i][idx + 2] = model->jac[i][2];
				else
					B[idx + i][idx + 2] = model->jac[i][3];
				if (model->cells[cell_idx - model->cellsNum_z - 2].u_next.SATUR)
					C[idx + i][idx + 2] = model->jac[i][size + 2];
				else
					C[idx + i][idx + 2] = model->jac[i][size + 3];
				if (model->cells[cell_idx + model->cellsNum_z + 2].u_next.SATUR)
					A[idx + i][idx + 2] = model->jac[i][size * 2 + 2];
				else
					A[idx + i][idx + 2] = model->jac[i][size * 2 + 3];
				if (model->cells[cell_idx - 1].u_next.SATUR)
					B[idx + i][idx + 2 - var_size] = model->jac[i][size * 3 + 2];
				else
					B[idx + i][idx + 2 - var_size] = model->jac[i][size * 3 + 3];
				if (model->cells[cell_idx + 1].u_next.SATUR)
					B[idx + i][idx + 2 + var_size] = model->jac[i][size * 4 + 2];
				else
					B[idx + i][idx + 2 + var_size] = model->jac[i][size * 4 + 3];
			}

			idx += var_size;
		}

		// Bottom cell
		model->setVariables(model->cells[cell_idx]);

		for (int i = 0; i < Variable::size - 1; i++)
		{
			RightSide[idx + i][0] = -model->y[i];

			// p
			B[idx + i][idx] = model->jac[i][0];
			B[idx + i][idx - var_size] = model->jac[i][size];
			// s_w
			B[idx + i][idx + 1] = model->jac[i][1];
			B[idx + i][idx + 1 - var_size] = model->jac[i][1 + size];
			// s_o / p_bub
			if (model->cells[cell_idx].u_next.SATUR)
				B[idx + i][idx + 2] = model->jac[i][2];
			else
				B[idx + i][idx + 2] = model->jac[i][3];
			if (model->cells[cell_idx - 1].u_next.SATUR)
				B[idx + i][idx + 2 - var_size] = model->jac[i][Variable::size + 2];
			else
				B[idx + i][idx + 2 - var_size] = model->jac[i][Variable::size + 3];
		}
	}

	construction_bz(MZ, 2);
}
void BlackOil2dSolver::LeftBoundAppr(int MZ, int key)
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

		const Variable next = model->cells[int(i / var_size)].u_next;
		const Variable nebr = model->cells[int(i / var_size) + model->cellsNum_z + 2].u_next;

		if (i % var_size == 2)
		{
			B[i][i] = -1.0;
			RightSide[i][0] = -next.values[2 + !next.SATUR] + nebr.values[2 + !nebr.SATUR];
		}
		else
		{
			B[i][i] = -1.0;
			RightSide[i][0] = -next.values[i % var_size] + nebr.values[i % var_size];
		}
	}

	map<int, double>::iterator it;
	int idx;
	if (key == PRES)
	{
		for (it = model->Qcell.begin(); it != model->Qcell.end(); ++it)
		{
			idx = var_size * it->first;
			model->setVariables(model->cells[it->first]);

			for (int i = 0; i < var_size; i++)
			{
				RightSide[idx + i][0] = -model->y[i];

				// p
				C[idx + i][idx] = model->jac[i][0];
				B[idx + i][idx] = model->jac[i][size];
				A[idx + i][idx] = model->jac[i][2 * size];
				// s_w
				C[idx + i][idx + 1] = model->jac[i][1];
				B[idx + i][idx + 1] = model->jac[i][1 + size];
				A[idx + i][idx + 1] = model->jac[i][1 + 2 * size];
				// s_o / p_bub
				if (model->cells[it->first].u_next.SATUR)
					C[idx + i][idx + 2] = model->jac[i][2];
				else
					C[idx + i][idx + 2] = model->jac[i][3];
				if (model->cells[it->first + model->cellsNum_z + 2].u_next.SATUR)
					B[idx + i][idx + 2] = model->jac[i][size + 2];
				else
					B[idx + i][idx + 2] = model->jac[i][size + 3];
				if (model->cells[it->first + 2 * model->cellsNum_z + 4].u_next.SATUR)
					A[idx + i][idx + 2] = model->jac[i][2 * size + 2];
				else
					A[idx + i][idx + 2] = model->jac[i][2 * size + 3];
			}
		}
	}

	construction_bz(MZ, 2);
}
void BlackOil2dSolver::RightBoundAppr(int MZ, int key)
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

		for (int cell_idx = (model->cellsNum_r + 1)*(model->cellsNum_z + 2); cell_idx < (model->cellsNum_r + 2)*(model->cellsNum_z + 2); cell_idx++)
		{
			model->setVariables(model->cells[cell_idx]);

			for (int i = 0; i < var_size; i++)
			{
				RightSide[idx + i][0] = -model->y[i];

				// p
				A[idx + i][idx] = model->jac[i][0];
				B[idx + i][idx] = model->jac[i][size];
				// s_w
				A[idx + i][idx + 1] = model->jac[i][1];
				B[idx + i][idx + 1] = model->jac[i][1 + size];
				// s_o / p_bub
				if (model->cells[cell_idx].u_next.SATUR)
					A[idx + i][idx + 2] = model->jac[i][2];
				else
					A[idx + i][idx + 2] = model->jac[i][3];
				if (model->cells[cell_idx - model->cellsNum_z - 2].u_next.SATUR)
					B[idx + i][idx + 2] = model->jac[i][size + 2];
				else
					B[idx + i][idx + 2] = model->jac[i][size + 3];
			}

			idx += var_size;
		}
	}

	construction_bz(MZ, 1);
}