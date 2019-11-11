#include <iomanip>
#include <iostream>

#include "model/Acid/2drec/Acid2dRecSolver.hpp"

#include "adolc/sparse/sparsedrivers.h"
#include "adolc/drivers/drivers.h"

using std::endl;
using std::setprecision;
using std::map;
using std::ofstream;
using namespace acid2drec;

template<class SolType>
Acid2dRecSolver<SolType>::Acid2dRecSolver(Acid2dRecModel* _model) : AbstractSolver<Acid2dRecModel>(_model)
{
	strNum = model->cellsNum * var_size;
	x = new double[strNum];
	y = new double[strNum];

	const int elems = (stencil * var_size) * (model->cellsNum * var_size);
	ind_i = new int[elems];
	ind_j = new int[elems];
	cols = new int[strNum];
	a = new double[elems];
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

	CONV_W2 = 1.e-5;		CONV_VAR = 1.e-10;
	MAX_ITER = 20;

    MAX_INIT_RES[0].first = 5.E-9;	MAX_INIT_RES[0].second = 2.E-8;         
    MAX_INIT_RES[1].first = 5.E-8;	MAX_INIT_RES[1].second = 2.E-7;

	MULT_UP = MULT_DOWN = 1.5;

	P.open(model->prefix + "P.dat", ofstream::out);
	qcells.open(model->prefix + "q_cells.dat", ofstream::out);
	pvd.open(model->prefix + "Acid2dRec.pvd", std::ofstream::out);
	pvd << "<VTKFile type = \"Collection\" version = \"1.0\" byte_order = \"LittleEndian\" header_type = \"UInt64\">\n";
	pvd << "\t<Collection>\n";
}
template<class SolType>
Acid2dRecSolver<SolType>::~Acid2dRecSolver()
{
	delete[] x, y;
	delete[] ind_i, ind_j, ind_rhs;
	delete[] cols;
	delete[] a, rhs;

	P.close();
	qcells.close();

	pvd << "\t</Collection>\n";
	pvd << "</VTKFile>\n";
	pvd.close();
}
template<class SolType>
void Acid2dRecSolver<SolType>::writeData()
{
	pvd << "\t\t<DataSet part=\"0\" timestep=\"";
	if (cur_t == 0.0)
		pvd << std::to_string(0.0);
	else
		pvd << cur_t * t_dim / 3600.0;
	pvd << "0\" file=\"Acid2dRec_" + std::to_string(step_idx) + ".vtu\"/>\n";

	double rate = 0.0, acid_rate = 0.0, cur_rate, pres = 0.0;
	for (int i = model->cellsNum_y + 2; i < (model->cellsNum_y + 2) * (model->cellsNum_x + 1); i += model->cellsNum_y + 2)
	{
		const auto& cell = model->cells[i];
		assert(cell.type == Type::BOTTOM);
		cur_rate = model->getRate(cell.num);
		rate += cur_rate;
		acid_rate += cur_rate * cell.u_next.xa;
		pres += cell.u_next.p;
	}
	model->injected_sol_volume += rate * model->ht;
	model->injected_acid_volume += acid_rate * model->ht;

	qcells << cur_t * t_dim / 3600.0;
	qcells << "\t" << rate * model->Q_dim * 86400.0;
	qcells << "\t" << model->injected_sol_volume * model->Q_dim * model->t_dim;
	qcells << "\t" << model->injected_acid_volume * model->Q_dim * model->t_dim;

	P << cur_t * t_dim / 3600.0 << "\t" << pres / model->cellsNum_x * model->P_dim / BAR_TO_PA << endl;

	double perm = 0.0;
	for (int i = 1; i < model->cellsNum_y + 1; i++)
	{
		const auto& cell = model->cells[model->cellsNum_y + 2 + i];
		perm += cell.hy / cell.props->getPermCoseni(cell.u_next.m, cell.u_next.p).value();
	}

	double cum_vol = 0.0;
	for (const auto& cell : model->cells)
	{
		cum_vol += (cell.props->m_init - cell.u_next.m) * cell.V * model->R_dim * model->R_dim * model->R_dim;
	}
	qcells << "\t" << cum_vol;
	qcells << "\t" << model->hy / model->props_sk[0].perm / perm;

	qcells << endl;
}
template<class SolType>
void Acid2dRecSolver<SolType>::analyzeNewtonConvergence()
{
    int control_period = 0;
    const double control_time = 1.E-5 * 3600.0 / t_dim;
    DECREASE_STEP = false;
    INCREASE_STEP = true;
    for (int i = 0; i < int(init_step_res.size()) - 1; i++)
    {
        if (cur_t < control_time)
            control_period = 0;
        else
            control_period = 1;

        if (init_step_res[i] > MAX_INIT_RES[control_period].second)
        {
            DECREASE_STEP = true;
            INCREASE_STEP = false;
            break;
        }
        if (init_step_res[i] > MAX_INIT_RES[control_period].first)
        {
            DECREASE_STEP = false;
            INCREASE_STEP = false;
            break;
        }

        if (init_step_res[i + 1] >= init_step_res[i])
        {
            INCREASE_STEP = false;
            if (iterations > 4)
                DECREASE_STEP = true;
            break;
        }
    }
}
template<class SolType>
void Acid2dRecSolver<SolType>::control()
{
	writeData();

	if (/*cur_t >= model->period[curTimePeriod]*/ model->injected_acid_volume >= model->max_acid_volume && curTimePeriod == 0)
	{
		curTimePeriod++;
		model->ht = model->ht_min;
		model->setPeriod(curTimePeriod);
	}

    analyzeNewtonConvergence();

	//if (INCREASE_STEP && model->ht > 0.05 * cur_t)
	//	INCREASE_STEP = false;

    MULT_UP = MULT_DOWN = 1.5;
    if(INCREASE_STEP && ((model->ht <= 5.0 * model->ht_max && curTimePeriod == 0) || 
                        (model->ht <= 0.5 * model->ht_max && curTimePeriod == 1)))
        model->ht = model->ht * MULT_UP;
    else if (DECREASE_STEP)
        model->ht = model->ht / MULT_DOWN;

	cur_t += model->ht;
}
template<class SolType>
void Acid2dRecSolver<SolType>::start()
{
	step_idx = 0;

	fillIndices();
	solver.Init(strNum, 1.e-30, 1.e-14, 1.e+4);

	model->setPeriod(curTimePeriod);

	while (cur_t < Tt && model->LeftBoundIsRate.size() > curTimePeriod)
	{
        model->snapshot_all(step_idx++);
        cout << "---------------------NEW TIME STEP---------------------" << endl;
		control();
        cout << setprecision(6);
        cout << "time = " << cur_t * t_dim / 3600.0 << "\tht = " << model->ht * t_dim / 3600.0 << endl;
        while (!doNextSmartStep())
        {
            cout << "------------------REPEATED TIME STEP------------------" << endl;
            cur_t -= model->ht;
            model->ht /= 1.5;
            cur_t += model->ht;
            cout << setprecision(6);
            cout << "time = " << cur_t * t_dim / 3600.0 << endl;
        }
		copyTimeLayer();
	}
	model->snapshot_all(step_idx);
	writeData();
}
template<class SolType>
bool Acid2dRecSolver<SolType>::doNextSmartStep()
{
    init_step_res.clear();
    final_step_res.clear();
    iter_num.clear();
    if (!solveSmartStep())
        return false;
    return true;
}
template<class SolType>
template<class VecType>
void Acid2dRecSolver<SolType>::copySolution(const VecType& sol)
{
}
template<>
template<>
void Acid2dRecSolver<ParSolver>::copySolution(const paralution::LocalVector<double>& sol)
{
	for (int i = 0; i < model->cellsNum; i++)
	{
		auto& cell = model->cells[i];
		for (int j = 0; j < var_size; j++)
			cell.u_next.values[j] += sol[i * var_size + j];
	}
}
template<>
template<>
void Acid2dRecSolver<HypreSolver>::copySolution(const HypreSolver::Vector& sol)
{
	for (int i = 0; i < model->cellsNum; i++)
	{
		auto& cell = model->cells[i];
		for (int j = 0; j < var_size; j++)
			cell.u_next.values[j] += sol[i * var_size + j];
	}
}
template<class SolType>
void Acid2dRecSolver<SolType>::checkStability()
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
template<class SolType>
void Acid2dRecSolver<SolType>::checkVariables()
{
	for (auto cell : model->cells)
	{
		if (cell.u_next.m > cell.props->m_max)
			cell.u_next.m = cell.props->m_max;
		if (cell.u_next.sw > 1.0)
			cell.u_next.sw = 1.0;
		if (cell.u_next.sw < 0.0)
			cell.u_next.sw = 0.0;
		if (cell.u_next.xa < 0.0)
			cell.u_next.xa = 0.0;
		if (cell.u_next.xs < 0.0)
			cell.u_next.xs = 0.0;
		if (cell.u_next.xw < 0.0)
			cell.u_next.xw = 0.0;
	}
}
template<class SolType>
bool Acid2dRecSolver<SolType>::solveSmartStep()
{
	return true;
}
template<>
bool Acid2dRecSolver<HypreSolver>::solveSmartStep()
{
	int cellIdx, varIdx;
	err_newton = 1.0;
	averValue(averValPrev);
	std::fill(dAverVal.begin(), dAverVal.end(), 1.0);
	iterations = 0;
	bool isHarder = (curTimePeriod == 0) ? false : true;
	bool isInit = false;// (cur_t < 0.001 / t_dim) ? true : false;

	auto continueIterations = [this]()
	{
		bool result = false;

		for (const auto& val : dAverVal)
			if (val > CONV_VAR)
				result = true;

		return result * (err_newton > CONV_W2) * (iterations < MAX_ITER);
	};
	while (continueIterations())
	{
		copyIterLayer();

		computeJac();
		fill();
		solver.Assemble(cols, ind_i, ind_j, a, elemNum, ind_rhs, rhs);
		solver.Solve();
		init_step_res.push_back(solver.init_res);
		final_step_res.push_back(solver.final_res);
		//iter_num.push_back(solver.iter_num);

		//if (init_step_res.back() >= 2.0 * final_step_res.back())
		//{
			copySolution(solver.getSolution());

			checkStability();
			err_newton = convergance(cellIdx, varIdx);
			if (iterations == 0)
				err_newton_first = err_newton;

			averValue(averVal);
			for (int i = 0; i < Acid2dRecModel::var_size; i++)
				dAverVal[i] = fabs(averVal[i] - averValPrev[i]);
			averValPrev = averVal;
		//}
		//else
		//	return false;
		iterations++;
		//model->snapshot_all(iterations+1);
	}

	checkVariables();

	cout << "Newton Iterations = " << iterations << endl;
	return true;
}
template<>
bool Acid2dRecSolver<ParSolver>::solveSmartStep()
{
	int cellIdx, varIdx;
	err_newton = 1.0;
	averValue(averValPrev);
	std::fill(dAverVal.begin(), dAverVal.end(), 1.0);
	iterations = 0;
	bool isHarder = (curTimePeriod == 0) ? false : true;
	bool isInit = true;// (cur_t < 0.001 / t_dim) ? true : false;

	auto continueIterations = [this]()
	{
		bool result = false;

		for (const auto& val : dAverVal)
            if (val > CONV_VAR)
                result = true;

		return result * (err_newton > CONV_W2) * (iterations < MAX_ITER);
	};
	while (continueIterations())
	{
		copyIterLayer();

		computeJac();
		fill();
		solver.Assemble(ind_i, ind_j, a, elemNum, ind_rhs, rhs);
		solver.Solve(PRECOND::ILU_SIMPLE, isInit, 1);
        init_step_res.push_back(solver.init_res);
        final_step_res.push_back(solver.final_res);
        iter_num.push_back(solver.iter_num);
		
        if (init_step_res.back() >= 2.0 * final_step_res.back())
        {
            copySolution(solver.getSolution());

            checkStability();
			err_newton = convergance(cellIdx, varIdx);
            if (iterations == 0)
                err_newton_first = err_newton;

            averValue(averVal);
            for (int i = 0; i < Acid2dRecModel::var_size; i++)
                dAverVal[i] = fabs(averVal[i] - averValPrev[i]);
            averValPrev = averVal;
        }
        else
            return false;
		iterations++;
		//model->snapshot_all(iterations);
	}

	checkVariables();

	cout << "Newton Iterations = " << iterations << endl;
    return true;
}
template<class SolType>
void Acid2dRecSolver<SolType>::computeJac()
{
	const Regime reg = (curTimePeriod == 0) ? INJECTION : STOP;
	trace_on(0);

	for (const auto& cell : model->cells)
	{
		auto& cur_x = model->x[cell.num];
		cur_x.m <<= cell.u_next.m;
		cur_x.p <<= cell.u_next.p;
		cur_x.sw <<= cell.u_next.sw;
		cur_x.xw <<= cell.u_next.xw;
		cur_x.xa <<= cell.u_next.xa;
		cur_x.xs <<= cell.u_next.xs;

		x[cell.num * var_size] = cell.u_next.m;
		x[cell.num * var_size + 1] = cell.u_next.p;
		x[cell.num * var_size + 2] = cell.u_next.sw;
		x[cell.num * var_size + 3] = cell.u_next.xw;
		x[cell.num * var_size + 4] = cell.u_next.xa;
		x[cell.num * var_size + 5] = cell.u_next.xs;
	}
	for (const auto& cell : model->cells)
	{
		const auto res = model->solve(cell, reg);
		model->h[cell.num * var_size] = res.m;
		model->h[cell.num * var_size + 1] = res.p;
		model->h[cell.num * var_size + 2] = res.sw;
		model->h[cell.num * var_size + 3] = res.xw;
		model->h[cell.num * var_size + 4] = res.xa;
		model->h[cell.num * var_size + 5] = res.xs;
	}

	for (int i = 0; i < strNum; i++)
		model->h[i] >>= y[i];

	trace_off();
}
template<class SolType>
void Acid2dRecSolver<SolType>::fill()
{
	sparse_jac(0, strNum, strNum, repeat,
		x, &elemNum, (unsigned int**)(&ind_i), (unsigned int**)(&ind_j), &a, options);

	for (int i = 0; i < strNum; i++)
	{
		cols[i] = 0;
		rhs[i] = -y[i];
	}

	for (int i = 0; i < elemNum; i++)
		cols[ind_i[i]]++;

}

template class Acid2dRecSolver<ParSolver>;
template class Acid2dRecSolver<HypreSolver>;