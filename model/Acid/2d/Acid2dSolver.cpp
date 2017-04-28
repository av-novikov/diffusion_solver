#include "model/Acid/2d/Acid2dSolver.hpp"

#include <map>

using namespace std;
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
void Acid2dSolver::writeData()
{
/*	double p = 0.0, s_w = 0.0, s_o = 0.0;

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
	*/
}
void Acid2dSolver::solveStep()
{
	int cellIdx, varIdx;
	err_newton = 1.0;
	averValue(averValPrev);
	std::fill(dAverVal.begin(), dAverVal.end(), 1.0);
	iterations = 0;

	auto continueIterations = [this]()
	{
		bool result = false;

		for (const auto& val : dAverVal)
			result += (val > 1.e-10);
		
		return result * (err_newton > 1.e-4) * (iterations < 20);
	};

	while (continueIterations())
	{
		copyIterLayer();

		Solve(model->cellsNum_r + 1, var_size * (model->cellsNum_z + 2), PRES);
		construction_from_fz(model->cellsNum_r + 2, var_size * (model->cellsNum_z + 2), PRES);
		//checkStability();

		err_newton = convergance(cellIdx, varIdx);

		averValue(averVal);
		for (int i = 0; i < var_size; i++)
			dAverVal[i] = fabs(averVal[i] - averValPrev[i]);
		averValPrev = averVal;

		iterations++;
	}

	cout << "Newton Iterations = " << iterations << endl;
}

void Acid2dSolver::construction_from_fz(int N, int n, int key)
{
	vector<Cell>::iterator it;
	if (key == PRES)
	{
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < model->cellsNum_z + 2; j++)
			{
				Variable& next = model->cells[i*(model->cellsNum_z + 2) + j].u_next;

				for (int k = 0; k < var_size; k++)
					next.values[k] += fz[i][var_size * j + k + 1];
			}
		}
	}
}
void Acid2dSolver::MiddleAppr(int current, int MZ, int key)
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

			for (int j = 0; j < var_size; j++)
			{
				B[idx + i][idx + j] = model->jac[i][j];
				B[idx + i][idx + j + var_size] = model->jac[i][j + size];
			}
		}
		idx += var_size;

		// Middle cells
		for (cell_idx = current * (model->cellsNum_z + 2) + 1; cell_idx < (current + 1) * (model->cellsNum_z + 2) - 1; cell_idx++)
		{
			model->setVariables(model->cells[cell_idx]);

			for (int i = 0; i < var_size; i++)
			{
				RightSide[idx + i][0] = -model->y[i];

				for (int j = 0; j < var_size; j++)
				{
					B[idx + i][idx + j] = model->jac[i][j];
					C[idx + i][idx + j] = model->jac[i][j + size];
					A[idx + i][idx + j] = model->jac[i][j + size * 2];
					B[idx + i][idx + j - var_size] = model->jac[i][j + size * 3];
					B[idx + i][idx + j + var_size] = model->jac[i][j + size * 4];
				}
			}

			idx += var_size;
		}

		// Bottom cell
		model->setVariables(model->cells[cell_idx]);

		for (int i = 0; i < Variable::size; i++)
		{
			RightSide[idx + i][0] = -model->y[i];

			for (int j = 0; j < var_size; j++)
			{
				B[idx + i][idx + j] = model->jac[i][j];
				B[idx + i][idx + j - var_size] = model->jac[i][j + size];
			}
		}
	}

	construction_bz(MZ, 2);
}
void Acid2dSolver::LeftBoundAppr(int MZ, int key)
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

		const Variable& next = model->cells[int(i / var_size)].u_next;
		const Variable& nebr = model->cells[int(i / var_size) + model->cellsNum_z + 2].u_next;
		RightSide[i][0] = -next.values[i % var_size] + nebr.values[i % var_size];
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

				for (int j = 0; j < var_size; j++)
				{
					C[idx + i][idx + j] = model->jac[i][j];
					B[idx + i][idx + j] = model->jac[i][j + size];
					A[idx + i][idx + j] = model->jac[i][j + 2 * size];
				}
			}
		}
	}

	construction_bz(MZ, 2);
}
void Acid2dSolver::RightBoundAppr(int MZ, int key)
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

				for (int j = 0; j < var_size; j++)
				{
					A[idx + i][idx + j] = model->jac[i][j];
					B[idx + i][idx + j] = model->jac[i][j + size];
				}
			}

			idx += var_size;
		}
	}

	construction_bz(MZ, 1);
}
void Acid2dSolver::writeMatrixes()
{
	const int MZ = (model->cellsNum_z + 2) * var_size;

	// Left
	mat_a.open("snaps/a_left.dat", ofstream::out);
	mat_b.open("snaps/b_left.dat", ofstream::out);
	mat_c.open("snaps/c_left.dat", ofstream::out);
	rhs.open("snaps/rhs_left.dat", ofstream::out);

	LeftBoundAppr(MZ, PRES);
	for (int i = 0; i < MZ; i++)
	{
		for (int j = 0; j < MZ; j++)
		{
			mat_a << i << "\t" << j << "\t" << A[i][j] << endl;
			mat_b << i << "\t" << j << "\t" << B[i][j] << endl;
			mat_c << i << "\t" << j << "\t" << C[i][j] << endl;
		}
		rhs << i << "\t" << RightSide[i][0] << endl;
	}
	mat_a.close();		mat_b.close();		mat_c.close();		rhs.close();

	// Middle
	for (int k = 0; k < model->cellsNum_r; k++)
	{
		const string afilename = "snaps/a_" + to_string(k + 1) + ".dat";
		const string bfilename = "snaps/b_" + to_string(k + 1) + ".dat";
		const string cfilename = "snaps/c_" + to_string(k + 1) + ".dat";
		const string rhsfilename = "snaps/rhs_" + to_string(k + 1) + ".dat";
		mat_a.open(afilename.c_str(), ofstream::out);
		mat_b.open(bfilename.c_str(), ofstream::out);
		mat_c.open(cfilename.c_str(), ofstream::out);
		rhs.open(rhsfilename.c_str(), ofstream::out);

		MiddleAppr(k + 1, MZ, PRES);
		for (int i = 0; i < MZ; i++)
		{
			for (int j = 0; j < MZ; j++)
			{
				mat_a << i << "\t" << j << "\t" << A[i][j] << endl;
				mat_b << i << "\t" << j << "\t" << B[i][j] << endl;
				mat_c << i << "\t" << j << "\t" << C[i][j] << endl;
			}
			rhs << i << "\t" << RightSide[i][0] << endl;
		}

		mat_a.close();		mat_b.close();		mat_c.close();		rhs.close();
	}

	// Right
	mat_a.open("snaps/a_right.dat", ofstream::out);
	mat_b.open("snaps/b_right.dat", ofstream::out);
	mat_c.open("snaps/c_right.dat", ofstream::out);
	rhs.open("snaps/rhs_right.dat", ofstream::out);

	RightBoundAppr(MZ, PRES);
	for (int i = 0; i < MZ; i++)
	{
		for (int j = 0; j < MZ; j++)
		{
			mat_a << i << "\t" << j << "\t" << A[i][j] << endl;
			mat_b << i << "\t" << j << "\t" << B[i][j] << endl;
			mat_c << i << "\t" << j << "\t" << C[i][j] << endl;
		}
		rhs << i << "\t" << RightSide[i][0] << endl;
	}
	mat_a.close();		mat_b.close();		mat_c.close();		rhs.close();
};