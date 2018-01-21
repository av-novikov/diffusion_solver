#include "model/Acid/1d/Acid1dSolver.hpp"

#include <map>

using namespace std;
using namespace acid1d;

Acid1dSolver::Acid1dSolver(Acid1d* _model) : basic1d::Basic1dSolver<Acid1d>(_model)
{
	//P.open("snaps/P.dat", ofstream::out);
	//S.open("snaps/S.dat", ofstream::out);
	//qcells.open("snaps/q_cells.dat", ofstream::out);
	const int strNum = var_size * model->cellsNum;
	
	CHOP_MULT = 0.1;
	MAX_SAT_CHANGE = 1.0;
	
	CONV_W2 = 1.e-4;		CONV_VAR = 1.e-10;
	MAX_ITER = 20;
}
Acid1dSolver::~Acid1dSolver()
{
	//P.close();
	//S.close();
	//qcells.close();
}
void Acid1dSolver::writeData()
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
void Acid1dSolver::start()
{
	int counter = 0;
	iterations = 8;

	model->setPeriod(curTimePeriod);

	while (cur_t < Tt)
	{
		control();
		if (model->isWriteSnaps)
			model->snapshot_all(counter++);
		doNextStep();
		copyTimeLayer();
		cout << "---------------------NEW TIME STEP---------------------" << endl;
	}
	if (model->isWriteSnaps)
		model->snapshot_all(counter++);
	writeData();
}

void Acid1dSolver::checkStability()
{
	auto barelyMobilLeft = [this](double s_cur, double s_crit) -> double
	{
		return s_crit + fabs(s_cur - s_crit) * CHOP_MULT;
	};
	auto barelyMobilRight = [this](double s_cur, double s_crit) -> double
	{
		return s_crit - fabs(s_crit - s_cur) * CHOP_MULT;
	};
	auto checkCritPoints = [=, this](auto& next, auto& iter, auto& props)
	{
		// Oil
		if ((next.so - props.s_oc) * (iter.so - props.s_oc) < 0.0)
			next.so = barelyMobilLeft(next.so, props.s_oc);
		if ((next.so - (1.0 - props.s_wc - props.s_gc)) * (iter.so - (1.0 - props.s_wc - props.s_gc)) < 0.0)
			next.so = barelyMobilRight(next.so, 1.0 - props.s_wc - props.s_gc);
		// Water
		if ((next.sw - props.s_wc) * (iter.sw - props.s_wc) < 0.0)
			next.sw = barelyMobilLeft(next.sw, props.s_wc);
		if ((next.sw - (1.0 - props.s_oc - props.s_gc)) * (iter.sw - (1.0 - props.s_oc - props.s_gc)) < 0.0)
			next.sw = barelyMobilRight(next.sw, 1.0 - props.s_oc - props.s_gc);
		// Gas
		if ((1.0 - next.sw - next.so - props.s_gc) * (1.0 - iter.sw - iter.so - props.s_gc) < 0.0)
			if (fabs(next.so - iter.so) > fabs(next.sw - iter.sw))
				next.so = 1.0 - next.sw - barelyMobilLeft(1.0 - next.so - next.sw, props.s_gc);
			else
				next.sw = 1.0 - next.so - barelyMobilLeft(1.0 - next.so - next.sw, props.s_gc);
		if ((props.s_wc - next.sw + props.s_oc - next.so) * (props.s_wc - iter.sw + props.s_oc - iter.so) < 0.0)
			if (fabs(next.so - iter.so) > fabs(next.sw - iter.sw))
				next.so = 1.0 - next.sw - barelyMobilRight(1.0 - next.so - next.sw, 1.0 - props.s_oc - props.s_wc);
			else
				next.sw = 1.0 - next.so - barelyMobilRight(1.0 - next.so - next.sw, 1.0 - props.s_oc - props.s_wc);
	};
	auto checkMaxResidual = [=, this](auto& next, auto& iter)
	{
		if (fabs(next.sw - iter.sw) > MAX_SAT_CHANGE)
			next.sw = iter.sw + sign(next.sw - iter.sw) * MAX_SAT_CHANGE;
	};

	for (auto& cell : model->cells)
	{
		const Skeleton_Props& props = *cell.props;
		Variable& next = cell.u_next;
		const Variable& iter = cell.u_iter;

		checkMaxResidual(next, iter);
	}
}
void Acid1dSolver::solveStep()
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
			result += (val > CONV_VAR);
		
		return result * (err_newton > CONV_W2) * (iterations < MAX_ITER);
	};

	while (continueIterations())
	{
		copyIterLayer();

		//fill();
		//solver.Assemble(ind_i, ind_j, a, elemNum, ind_rhs, rhs);
		//solver.Solve(PRECOND::ILU_SIMPLE);
		//copySolution(solver.getSolution());

		//writeMatrixes();
		Solve(model->cellsNum_x + 1, var_size, PRES);
		construction_from_fz(model->cellsNum_x + 2, var_size, PRES);
		checkStability();
		err_newton = convergance(cellIdx, varIdx);

		averValue(averVal);
		for (int i = 0; i < var_size; i++)
			dAverVal[i] = fabs(averVal[i] - averValPrev[i]);
		averValPrev = averVal;

		model->snapshot_all(iterations + 1);
		iterations++;
	}

	cout << "Newton Iterations = " << iterations << endl;
}

void Acid1dSolver::construction_from_fz(int N, int n, int key)
{
	vector<Cell>::iterator it;
	if (key == PRES)
	{
		for (int i = 0; i < N; i++)
		{
				Variable& next = model->cells[i].u_next;
				next.m += fz[i][1];
				next.p += fz[i][2];
				next.sw += fz[i][3];
				next.xa += fz[i][4];
				next.xw += fz[i][5];
				next.xs += fz[i][6];
		}
	}
}
void Acid1dSolver::MiddleAppr(int current, int MZ, int key)
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
		int cell_idx = current;

		model->setVariables(model->cells[cell_idx]);
		for (int i = 0; i < var_size; i++)
		{
			RightSide[idx + i][0] = -model->y[i];

			for (int j = 0; j < var_size; j++)
			{
				B[idx + i][idx + j] = model->jac[i][j];
				C[idx + i][idx + j] = model->jac[i][j + size];
				A[idx + i][idx + j] = model->jac[i][j + size * 2];
			}
		}
	}

	construction_bz(MZ, 2);
}
void Acid1dSolver::LeftBoundAppr(int MZ, int key)
{
	for (int i = 0; i < MZ; i++)
	{
		for (int j = 0; j < MZ; j++)
		{
			C[i][j] = 0.0;
			B[i][j] = 0.0;
			A[i][j] = 0.0;
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
void Acid1dSolver::RightBoundAppr(int MZ, int key)
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
		int cell_idx = model->cellsNum_x + 1;
		
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
	}

	construction_bz(MZ, 1);
}
void Acid1dSolver::writeMatrixes()
{
	const int MZ = var_size;

	// Left
	mat_a.open("snaps/a_left.dat", ofstream::out);
	mat_b.open("snaps/b_left.dat", ofstream::out);
	mat_c.open("snaps/c_left.dat", ofstream::out);
	rhs_os.open("snaps/rhs_left.dat", ofstream::out);

	LeftBoundAppr(MZ, PRES);
	for (int i = 0; i < MZ; i++)
	{
		for (int j = 0; j < MZ; j++)
		{
			mat_a << i << "\t" << j << "\t" << A[i][j] << endl;
			mat_b << i << "\t" << j << "\t" << B[i][j] << endl;
			mat_c << i << "\t" << j << "\t" << C[i][j] << endl;
		}
		rhs_os << i << "\t" << RightSide[i][0] << endl;
	}
	mat_a.close();		mat_b.close();		mat_c.close();		rhs_os.close();

	// Middle
	for (int k = 0; k < model->cellsNum_x; k++)
	{
		const string afilename = "snaps/a_" + to_string(k + 1) + ".dat";
		const string bfilename = "snaps/b_" + to_string(k + 1) + ".dat";
		const string cfilename = "snaps/c_" + to_string(k + 1) + ".dat";
		const string rhsfilename = "snaps/rhs_" + to_string(k + 1) + ".dat";
		mat_a.open(afilename.c_str(), ofstream::out);
		mat_b.open(bfilename.c_str(), ofstream::out);
		mat_c.open(cfilename.c_str(), ofstream::out);
		rhs_os.open(rhsfilename.c_str(), ofstream::out);

		MiddleAppr(k + 1, MZ, PRES);
		for (int i = 0; i < MZ; i++)
		{
			for (int j = 0; j < MZ; j++)
			{
				mat_a << i << "\t" << j << "\t" << A[i][j] << endl;
				mat_b << i << "\t" << j << "\t" << B[i][j] << endl;
				mat_c << i << "\t" << j << "\t" << C[i][j] << endl;
			}
			rhs_os << i << "\t" << RightSide[i][0] << endl;
		}

		mat_a.close();		mat_b.close();		mat_c.close();		rhs_os.close();
	}

	// Right
	mat_a.open("snaps/a_right.dat", ofstream::out);
	mat_b.open("snaps/b_right.dat", ofstream::out);
	mat_c.open("snaps/c_right.dat", ofstream::out);
	rhs_os.open("snaps/rhs_right.dat", ofstream::out);

	RightBoundAppr(MZ, PRES);
	for (int i = 0; i < MZ; i++)
	{
		for (int j = 0; j < MZ; j++)
		{
			mat_a << i << "\t" << j << "\t" << A[i][j] << endl;
			mat_b << i << "\t" << j << "\t" << B[i][j] << endl;
			mat_c << i << "\t" << j << "\t" << C[i][j] << endl;
		}
		rhs_os << i << "\t" << RightSide[i][0] << endl;
	}
	mat_a.close();		mat_b.close();		mat_c.close();		rhs_os.close();
};