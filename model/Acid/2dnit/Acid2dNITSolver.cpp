#include "model/Acid/2dnit/Acid2dNITSolver.hpp"

#include <map>

using namespace std;
using namespace acid2dnit;

Acid2dNITSolver::Acid2dNITSolver(Acid2dNIT* _model) : basic2d::Basic2dSolver<Acid2dNIT>(_model)
{
	//P.open("snaps/P.dat", ofstream::out);
	//S.open("snaps/S.dat", ofstream::out);
	//qcells.open("snaps/q_cells.dat", ofstream::out);
	//temp.open("snaps/temp.dat", ofstream::out);
	const int strNum = var_size * model->cellsNum;
	ind_i = new int[stencil * var_size * strNum];
	ind_j = new int[stencil * var_size * strNum];
	cols = new int[strNum];
	a = new double[stencil * var_size * strNum];
	ind_rhs = new int[strNum];
	rhs = new double[strNum];

	CHOP_MULT = 0.1;
	MAX_SAT_CHANGE = 1.0;
	
	CONV_W2 = 1.e-4;		CONV_VAR = 1.e-10;
	MAX_ITER = 20;

	pvd.open("snaps/Acid2dNIT.pvd", std::ofstream::out);
	pvd << "<VTKFile type = \"Collection\" version = \"1.0\" byte_order = \"LittleEndian\" header_type = \"UInt64\">\n";
	pvd << "\t<Collection>\n";
}
Acid2dNITSolver::~Acid2dNITSolver()
{
	//P.close();
	//S.close();
	//qcells.close();
	//temp.close();

	delete[] ind_i, ind_j, ind_rhs;
	delete[] cols;
	delete[] a, rhs;

	pvd << "\t</Collection>\n";
	pvd << "</VTKFile>\n";
	pvd.close();
}
void Acid2dNITSolver::writeData()
{
	pvd << "\t\t<DataSet part=\"0\" timestep=\"";
	if (cur_t == 0.0)
		pvd << std::to_string(0.0);
	else
		pvd << cur_t * t_dim / 3600.0;
	pvd << "0\" file=\"Acid2dNIT_" + std::to_string(step_idx) + ".vtp\"/>\n";


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
void Acid2dNITSolver::start()
{
	step_idx = 0;
	iterations = 8;

	fillIndices();
	solver.Init(var_size * model->cellsNum, 1.e-15, 1.e-15);

	model->setPeriod(curTimePeriod);

	while (cur_t < Tt)
	{
		control();
		if (model->isWriteSnaps)
			model->snapshot_all(step_idx++);
		doNextStep();
		copyTimeLayer();
		cout << "---------------------NEW TIME STEP---------------------" << endl;
	}
	if (model->isWriteSnaps)
		model->snapshot_all(step_idx);
	writeData();
}
void Acid2dNITSolver::copySolution(const Vector& sol)
{
	for (int i = 0; i < model->cellsNum; i++)
		for(int j = 0; j < var_size; j++)
			model->cells[i].u_next.values[j] += sol[i * var_size + j];
}

void Acid2dNITSolver::checkStability()
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
		if ((1.0 - next.sw - props.s_oc) * (1.0 - iter.sw - props.s_oc) < 0.0)
			next.sw = 1.0 - barelyMobilLeft(1.0 - next.sw, props.s_oc);
		if ((props.s_wc - next.sw) * (props.s_wc - iter.sw) < 0.0)
			next.sw = 1.0 - barelyMobilRight(1.0 - next.sw, 1.0 - props.s_wc);
		// Water
		if ((next.sw - props.s_wc) * (iter.sw - props.s_wc) < 0.0)
			next.sw = barelyMobilLeft(next.sw, props.s_wc);
		if ((next.sw - (1.0 - props.s_oc)) * (iter.sw - (1.0 - props.s_oc)) < 0.0)
			next.sw = barelyMobilRight(next.sw, 1.0 - props.s_oc);
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

		checkCritPoints(next, iter, props);
		checkMaxResidual(next, iter);
	}
}
void Acid2dNITSolver::solveStep()
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
		//writeMatrixes();
		fill();
		solver.Assemble(ind_i, ind_j, a, elemNum, ind_rhs, rhs);
		solver.Solve(PRECOND::ILU_SIMPLE);
		copySolution(solver.getSolution());

		checkStability();
		err_newton = convergance(cellIdx, varIdx);

		averValue(averVal);
		for (int i = 0; i < var_size; i++)
			dAverVal[i] = fabs(averVal[i] - averValPrev[i]);
		averValPrev = averVal;

		iterations++;
	}

	cout << "Newton Iterations = " << iterations << "\t cur_t = " << cur_t << endl;
}

void Acid2dNITSolver::fillIndices()
{
	int pres_counter = 0, temp_counter = 0;
	int col_idx = 0;

	for (const auto& cell : model->cells)
	{
		cols[col_idx] = 0;
		auto& pres_idx = getMatrixStencil(cell);

		for (int i = 0; i < var_size; i++)
		{
			const int str_idx = var_size * cell.num + i;
			for (const int idx : pres_idx)
			{
				for (int j = 0; j < var_size; j++)
				{
					ind_i[pres_counter] = str_idx;			ind_j[pres_counter++] = var_size * idx + j;
				}
			}
		}

		cols[col_idx++] += var_size * pres_idx.size();
		pres_idx.clear();
	}

	elemNum = pres_counter;

	for (int i = 0; i < model->cellsNum * var_size; i++)
		ind_rhs[i] = i;
}
void Acid2dNITSolver::fill()
{
	int counter = 0;
	int nebr_idx;

	for (const auto& cell : model->cells)
	{
		model->setVariables(cell);

		auto& mat_idx = getMatrixStencil(cell);
		for (int i = 0; i < var_size; i++)
		{
			const int str_idx = var_size * cell.num + i;
			nebr_idx = 0;
			for (const int idx : mat_idx)
			{
				for (int j = 0; j < var_size; j++)
					a[counter++] = model->jac[i][var_size * nebr_idx + j];

				nebr_idx++;
			}

			rhs[str_idx] = -model->y[i];
		}
		mat_idx.clear();
	}

}