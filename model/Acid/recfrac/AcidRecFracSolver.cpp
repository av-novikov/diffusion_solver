#include <iomanip>
#include <iostream>

#include "model/Acid/recfrac/AcidRecFracSolver.hpp"

#include "adolc/sparse/sparsedrivers.h"
#include "adolc/drivers/drivers.h"

using std::endl;
using std::setprecision;
using std::map;
using std::ofstream;
using namespace acidrecfrac;

AcidRecFracSolver::AcidRecFracSolver(AcidRecFrac* _model) : AbstractSolver<AcidRecFrac>(_model)
{
	strNum = model->cellsNum_frac * var_frac_size + model->cellsNum_poro * var_poro_size;
	x = new double[strNum];
	y = new double[strNum];

	const int poro_elems = (poro_stencil * var_poro_size) * (model->cellsNum_poro * var_poro_size);
	const int frac_elems = (frac_stencil * var_frac_size) * (model->cellsNum_frac * var_frac_size);
	ind_i = new int[poro_elems + frac_elems];
	ind_j = new int[poro_elems + frac_elems];
	cols = new int[strNum];
	a = new double[poro_elems + frac_elems];
	ind_rhs = new int[strNum];
	rhs = new double[strNum];

	options[0] = 0;          /* sparsity pattern by index domains (default) */
	options[1] = 0;          /*                         safe mode (default) */
	options[2] = 0;          /*              not required if options[0] = 0 */
	options[3] = 0;          /*                column compression (default) */
	repeat = 0;
	//P.open("snaps/P.dat", ofstream::out);
	//S.open("snaps/S.dat", ofstream::out);
	//qcells.open("snaps/Q.dat", ofstream::out);

	CHOP_MULT = 0.1;
	MAX_SAT_CHANGE = 1.0;

	CONV_W2 = 1.e-4;		CONV_VAR = 1.e-10;
	MAX_ITER = 20;

	P.open("snaps/P.dat", ofstream::out);
	qcells.open("snaps/q_cells.dat", ofstream::out);
	pvd_frac.open("snaps/AcidRecFrac_frac.pvd", std::ofstream::out);
	pvd_poro.open("snaps/AcidRecFrac_poro.pvd", std::ofstream::out);
	pvd_frac << "<VTKFile type = \"Collection\" version = \"1.0\" byte_order = \"LittleEndian\" header_type = \"UInt64\">\n";
	pvd_poro << "<VTKFile type = \"Collection\" version = \"1.0\" byte_order = \"LittleEndian\" header_type = \"UInt64\">\n";
	pvd_frac << "\t<Collection>\n";
	pvd_poro << "\t<Collection>\n";
}
AcidRecFracSolver::~AcidRecFracSolver()
{
	delete[] x, y;
	delete[] ind_i, ind_j, ind_rhs;
	delete[] cols;
	delete[] a, rhs;

	P.close();
	qcells.close();

	pvd_frac << "\t</Collection>\n";
	pvd_frac << "</VTKFile>\n";
	pvd_frac.close();

	pvd_poro << "\t</Collection>\n";
	pvd_poro << "</VTKFile>\n";
	pvd_poro.close();
}
void AcidRecFracSolver::writeData()
{
	pvd_frac << "\t\t<DataSet part=\"0\" timestep=\"";
	if (cur_t == 0.0)
		pvd_frac << std::to_string(0.0);
	else
		pvd_frac << cur_t * t_dim / 3600.0;
	pvd_frac << "0\" file=\"AcidRecFrac_frac_" + std::to_string(step_idx) + ".vtu\"/>\n";
	pvd_poro << "\t\t<DataSet part=\"0\" timestep=\"";
	if (cur_t == 0.0)
		pvd_poro << std::to_string(0.0);
	else
		pvd_poro << cur_t * t_dim / 3600.0;
	pvd_poro << "0\" file=\"AcidRecFrac_poro_" + std::to_string(step_idx) + ".vtu\"/>\n";

	double p = 0.0;
	qcells << cur_t * t_dim / 3600.0;

	map<int, double>::iterator it;
	for (it = model->Qcell.begin(); it != model->Qcell.end(); ++it)
	{
		p += model->cells_frac[it->first].u_next.p * model->P_dim;
		if (model->leftBoundIsRate)
			qcells << "\t" << it->second * model->Q_dim * 86400.0;
		else
			qcells << "\t" << model->getRate(it->first) * model->Q_dim * 86400.0;
	}

	P << cur_t * t_dim / 3600.0 <<
		"\t" << p / (double)(model->Qcell.size()) << endl;

	qcells << endl;
}
void AcidRecFracSolver::control()
{
	writeData();

	if (cur_t >= model->period[curTimePeriod])
	{
		curTimePeriod++;
		model->ht = model->ht_min;
		model->setPeriod(curTimePeriod);
	}

	if (model->ht <= model->ht_max && iterations < 5 && (curTimePeriod == 0 || err_newton_first < 1.0))
		model->ht = model->ht * 1.5;
	else if (iterations > 5 && model->ht > model->ht_min)
		model->ht = model->ht / 1.5;

	if (cur_t + model->ht > model->period[curTimePeriod])
		model->ht = model->period[curTimePeriod] - cur_t;

	cur_t += model->ht;
}
void AcidRecFracSolver::start()
{
	iterations = 8;
	step_idx = 0;

	fillIndices();
	solver.Init(strNum, 1.e-15, 1.e-15);

	model->setPeriod(curTimePeriod);

	while (cur_t < Tt)
	{
		control();
		model->snapshot_all(step_idx++);
		doNextStep();
		copyTimeLayer();
		cout << "---------------------NEW TIME STEP---------------------" << endl;
		cout << setprecision(6);
		cout << "time = " << cur_t << endl;
	}
	model->snapshot_all(step_idx);
	writeData();
}
void AcidRecFracSolver::doNextStep()
{
	solveStep();
	model->calculateTrans();
}
void AcidRecFracSolver::copySolution(const paralution::LocalVector<double>& sol)
{
	for (int i = 0; i < model->cellsNum_frac; i++)
	{
		auto& cell = model->cells_frac[i];
		for (int j = 0; j < var_frac_size; j++)
			cell.u_next.values[j] += sol[i * var_frac_size + j];
	}

	const int startIdx = var_frac_size * model->cellsNum_frac;
	for (int i = 0; i < model->cellsNum_poro; i++)
	{
		auto& cell = model->cells_poro[i];
		for (int j = 0; j < var_poro_size; j++)
			cell.u_next.values[j] += sol[startIdx + i * var_poro_size + j];
	}
}
void AcidRecFracSolver::checkStability()
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
		if ((1.0 - next.s - props.s_oc) * (1.0 - iter.s - props.s_oc) < 0.0)
			next.s = 1.0 - barelyMobilLeft(1.0 - next.s, props.s_oc);
		if ((1.0 - next.s - (1.0 - props.s_wc)) * (1.0 - iter.s - (1.0 - props.s_wc)) < 0.0)
			next.s = 1.0 - barelyMobilRight(1.0 - next.s, 1.0 - props.s_wc);
		// Water
		if ((next.s - props.s_wc) * (iter.s - props.s_wc) < 0.0)
			next.s = barelyMobilLeft(next.s, props.s_wc);
		if ((next.s - (1.0 - props.s_oc)) * (iter.s - (1.0 - props.s_oc)) < 0.0)
			next.s = barelyMobilRight(next.s, 1.0 - props.s_oc);
	};
	auto checkMaxResidual = [=, this](auto& next, auto& iter)
	{
		if (fabs(next.s - iter.s) > MAX_SAT_CHANGE)
			next.s = iter.s + sign(next.s - iter.s) * MAX_SAT_CHANGE;
	};
	/*for (size_t i = 0; i < size; i++)
	{
		auto& data = (*model)[i];
		checkCritPoints(data.u_next, data.u_iter, model->props_sk[0]);
		checkMaxResidual(data.u_next, data.u_iter);
	}*/
}

void AcidRecFracSolver::solveStep()
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

		computeJac();
		fill();
		solver.Assemble(ind_i, ind_j, a, elemNum, ind_rhs, rhs);
		solver.Solve(PRECOND::ILU_SIMPLE);
		copySolution(solver.getSolution());

		checkStability();
		auto err = convergance();
		err_newton = fmax(err[0].err, err[1].err);
		if (iterations == 0)
			err_newton_first = err_newton;

		averValue(averVal);
		for (int i = 0; i < AcidRecFrac::var_size; i++)
			dAverVal[i] = fabs(averVal[i] - averValPrev[i]);
		averValPrev = averVal;
		iterations++;
	}

	cout << "Newton Iterations = " << iterations << endl;
}
void AcidRecFracSolver::computeJac()
{
	const Regime reg = (curTimePeriod == 0) ? INJECTION : STOP;
	trace_on(0);
	const int start_idx = model->cellsNum_frac * var_frac_size;

	for (const auto& cell : model->cells_frac)
	{
		auto& cur_x = model->x_frac[cell.num];
		cur_x.p <<= cell.u_next.p;
		cur_x.c <<= cell.u_next.c;

		x[cell.num * var_frac_size] = cell.u_next.p;
		x[cell.num * var_frac_size + 1] = cell.u_next.c;
	}
	for (const auto& cell : model->cells_poro)
	{
		auto& cur_x = model->x_poro[cell.num];
		cur_x.m <<= cell.u_next.m;
		cur_x.p <<= cell.u_next.p;
		cur_x.sw <<= cell.u_next.sw;
		cur_x.xw <<= cell.u_next.xw;
		cur_x.xa <<= cell.u_next.xa;
		cur_x.xs <<= cell.u_next.xs;

		x[start_idx + cell.num * var_poro_size] = cell.u_next.m;
		x[start_idx + cell.num * var_poro_size + 1] = cell.u_next.p;
		x[start_idx + cell.num * var_poro_size + 2] = cell.u_next.sw;
		x[start_idx + cell.num * var_poro_size + 3] = cell.u_next.xw;
		x[start_idx + cell.num * var_poro_size + 4] = cell.u_next.xa;
		x[start_idx + cell.num * var_poro_size + 5] = cell.u_next.xs;
	}
	// Fracture
	for (const auto& cell : model->cells_frac)
	{
		const auto res = model->solveFrac(cell, reg);
		model->h[cell.num * var_frac_size] = res.p;
		model->h[cell.num * var_frac_size + 1] = res.c;
	}
	// Porous medium
	for (const auto& cell : model->cells_poro)
	{
		const auto res = model->solvePoro(cell, reg);
		model->h[start_idx + cell.num * var_poro_size] = res.m;
		model->h[start_idx + cell.num * var_poro_size + 1] = res.p;
		model->h[start_idx + cell.num * var_poro_size + 2] = res.sw;
		model->h[start_idx + cell.num * var_poro_size + 3] = res.xw;
		model->h[start_idx + cell.num * var_poro_size + 4] = res.xa;
		model->h[start_idx + cell.num * var_poro_size + 5] = res.xs;
	}

	for (int i = 0; i < strNum; i++)
		model->h[i] >>= y[i];

	trace_off();
}
void AcidRecFracSolver::fill()
{
	sparse_jac(0, strNum, strNum, repeat,
		x, &elemNum, (unsigned int**)(&ind_i), (unsigned int**)(&ind_j), &a, options);

	for(int i = 0; i < strNum; i++)
		rhs[i] = -y[i];
}