#include <iomanip>
#include <iostream>

#include "model/Acid/frac/AcidFracSolver.hpp"

#include "adolc/sparse/sparsedrivers.h"
#include "adolc/drivers/drivers.h"

using std::endl;
using std::setprecision;
using namespace acidfrac;

AcidFracSolver::AcidFracSolver(AcidFrac* _model) : AbstractSolver<AcidFrac>(_model)
{
	strNum = model->cellsNum * var_frac_size + model->cellsPoroNum * var_poro_size;
	x = new double[strNum];
	y = new double[strNum];

	const int poro_elems = (poro_stencil * var_poro_size) * (model->cellsPoroNum * var_poro_size);
	const int frac_elems = (frac_stencil * var_frac_size) * (model->cellsNum * var_frac_size);
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
}
AcidFracSolver::~AcidFracSolver()
{
	delete[] x, y;
	delete[] ind_i, ind_j, ind_rhs;
	delete[] cols;
	delete[] a, rhs;
}
void AcidFracSolver::writeData()
{
}
void AcidFracSolver::control()
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
void AcidFracSolver::start()
{
	int counter = 0;
	iterations = 8;

	fillIndices();
	solver.Init(strNum, 1.e-15, 1.e-15);

	model->setPeriod(curTimePeriod);

	while (cur_t < Tt)
	{
		control();
		model->snapshot_all(counter++);
		doNextStep();
		copyTimeLayer();
		cout << "---------------------NEW TIME STEP---------------------" << endl;
		cout << setprecision(6);
		cout << "time = " << cur_t << endl;
	}
	model->snapshot_all(counter++);
	writeData();
}
void AcidFracSolver::doNextStep()
{
	solveStep();
}
void AcidFracSolver::copySolution(const paralution::LocalVector<double>& sol)
{
	for (int i = 0; i < model->cellsNum; i++)
	{
		auto& cell = model->cells_frac[i];
		for (int j = 0; j < var_frac_size; j++)
			cell.u_next.values[j] += sol[i * var_frac_size + j];
	}

	for (auto& grid : model->poro_grids)
	{
		const int start_idx = var_frac_size * model->cellsNum + var_poro_size * grid.start_idx;
		for (int i = 0; i < grid.cellsNum + 2; i++)
		{
			auto& cell = grid.cells[i];
			for (int j = 0; j < var_poro_size; j++)
				cell.u_next.values[j] += sol[start_idx + i * var_poro_size + j];
		}
	}
}
void AcidFracSolver::checkStability()
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

void AcidFracSolver::solveStep()
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
		solver.Solve(PRECOND::ILU_SERIOUS);
		copySolution(solver.getSolution());

		checkStability();
		err_newton = convergance(cellIdx, varIdx);

		averValue(averVal);
		for (int i = 0; i < AcidFrac::var_size; i++)
			dAverVal[i] = fabs(averVal[i] - averValPrev[i]);
		averValPrev = averVal;

		model->snapshot_all(iterations + 1);
		iterations++;
	}

	cout << "Newton Iterations = " << iterations << endl;
}
void AcidFracSolver::computeJac()
{
	trace_on(0);

	for (int i = 0; i < model->cellsNum; i++)
	{
		model->x_frac[i].p <<= model->cells_frac[i].u_next.p;
		model->x_frac[i].c <<= model->cells_frac[i].u_next.c;

		x[i * var_frac_size] = model->cells_frac[i].u_next.p;
		x[i * var_frac_size + 1] = model->cells_frac[i].u_next.c;
	}
	for (const auto& grid : model->poro_grids)
	{
		for (const auto& cell : grid.cells)
		{
			auto& cur_x = model->x_poro[grid.start_idx + cell.num];
			cur_x.m <<= cell.u_next.m;
			cur_x.p <<= cell.u_next.p;
			cur_x.sw <<= cell.u_next.sw;
			cur_x.xw <<= cell.u_next.xw;
			cur_x.xa <<= cell.u_next.xa;
			cur_x.xs <<= cell.u_next.xs;

			const int start_idx = model->cellsNum * var_frac_size + grid.start_idx * var_poro_size;
			x[start_idx + cell.num * var_poro_size] = cell.u_next.m;
			x[start_idx + cell.num * var_poro_size + 1] = cell.u_next.p;
			x[start_idx + cell.num * var_poro_size + 2] = cell.u_next.sw;
			x[start_idx + cell.num * var_poro_size + 3] = cell.u_next.xw;
			x[start_idx + cell.num * var_poro_size + 4] = cell.u_next.xa;
			x[start_idx + cell.num * var_poro_size + 5] = cell.u_next.xs;
		}
	}

	// Fracture
	for (const auto cell : model->cells_frac)
	{
		const auto res = model->solveFrac(cell);
		model->h[cell.num * var_frac_size] = res.p;
		model->h[cell.num * var_frac_size + 1] = res.c;
	}
	// Porous medium
	for (const auto& grid : model->poro_grids)
	{
		const int start_idx = model->cellsNum * var_frac_size + grid.start_idx * var_poro_size;
		for (const auto cell : grid.cells)
		{
			const auto res = model->solvePoro(cell);
			model->h[start_idx + cell.num * var_poro_size] = res.m;
			model->h[start_idx + cell.num * var_poro_size + 1] = res.p;
			model->h[start_idx + cell.num * var_poro_size + 2] = res.sw;
			model->h[start_idx + cell.num * var_poro_size + 3] = res.xw;
			model->h[start_idx + cell.num * var_poro_size + 4] = res.xa;
			model->h[start_idx + cell.num * var_poro_size + 5] = res.xs;
		}
	}

	for (int i = 0; i < strNum; i++)
		model->h[i] >>= y[i];

	trace_off();
}
void AcidFracSolver::fill()
{
	sparse_jac(0, strNum, strNum, repeat,
		x, &elemNum, (unsigned int**)(&ind_i), (unsigned int**)(&ind_j), &a, options);

	for(int i = 0; i < strNum; i++)
		rhs[i] = -y[i];
}