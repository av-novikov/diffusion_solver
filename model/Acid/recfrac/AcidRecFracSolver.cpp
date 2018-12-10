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

	CONV_W2 = 1.e-5;		CONV_VAR = 1.e-12;
	MAX_ITER = 20;

    MAX_INIT_RES1 = 3.E-9;
    MAX_INIT_RES2 = 1.E-9;

	P.open("snaps/P.dat", ofstream::out);
	qcells.open("snaps/q_cells.dat", ofstream::out);
	trans.open("snaps/Trans.dat", ofstream::out);
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
	trans.close();

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
	pvd_poro << "\t\t<DataSet part=\"0\" timestep=\"";
	if (cur_t == 0.0)
	{
		pvd_frac << std::to_string(0.0);
		pvd_poro << std::to_string(0.0);
	}
	else
	{
		pvd_frac << cur_t * t_dim / 3600.0;
		pvd_poro << cur_t * t_dim / 3600.0;
	}
	pvd_frac << "0\" file=\"AcidRecFrac_frac_" + std::to_string(step_idx) + ".vtu\"/>\n";
	pvd_poro << "0\" file=\"AcidRecFrac_poro_" + std::to_string(step_idx) + ".vtu\"/>\n";

	double q = 0.0;
	for (int i = 0; i < (model->cellsNum_y_frac + 1) * (model->cellsNum_z + 2); i++)
	{
		const auto& cell = model->cells_frac[i];
		if (cell.type == FracType::FRAC_IN && cell.hy != 0.0 && cell.hz != 0.0)
			q += model->getRate(cell.num);
	}
	model->injected_sol_volume -= q * model->ht;
	model->injected_acid_volume -= q * model->c * model->ht;
	qcells << cur_t * t_dim / 3600.0;
	qcells << "\t" << -q * model->Q_dim * 86400.0;
	qcells << "\t" << model->injected_sol_volume * model->Q_dim * model->t_dim;
	qcells << "\t" << model->injected_acid_volume * model->Q_dim * model->t_dim;

	double p = model->cells_frac[model->cellsNum_y_frac + 2].u_next.p;
	P << cur_t * t_dim / 3600.0 << "\t" << p * model->P_dim / BAR_TO_PA << endl;

	const double k0 = M2toMilliDarcy(model->props_sk[0].perm * model->R_dim * model->R_dim);
	double mean_trans = 0.0, harm_mean_trans = 0.0, cur_trans;
	for (int i = 0; i < model->trans.size(); i++)
	{
		cur_trans = model->trans[i] * model->widths[i] * model->R_dim * k0;
		mean_trans += cur_trans / model->trans.size();
		harm_mean_trans += 1.0 / cur_trans;
	}
	harm_mean_trans = model->trans.size() / harm_mean_trans;

	trans << cur_t * t_dim / 3600.0 << "\t" << mean_trans << "\t" << harm_mean_trans << endl;

	qcells << endl;
}
void AcidRecFracSolver::analyzeNewtonConvergence()
{
    DECREASE_STEP = false;
    INCREASE_STEP = true;
    for (int i = 0; i < int(init_step_res.size()) - 1; i++)
    {
        if (init_step_res[i] > MAX_INIT_RES1)
        {
            DECREASE_STEP = true;
            INCREASE_STEP = false;
            break;
        }
        if (init_step_res[i] > MAX_INIT_RES2)
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
void AcidRecFracSolver::control()
{
	writeData();

	if (/*cur_t >= model->period[curTimePeriod]*/ model->injected_sol_volume >= model->max_sol_volume && curTimePeriod == 0)
	{
		curTimePeriod++;
		model->ht = model->ht_min;
		model->setPeriod(curTimePeriod);
	}

    analyzeNewtonConvergence();
    if(INCREASE_STEP && ((model->ht <= 5.0 * model->ht_max && curTimePeriod == 0) || 
                        (model->ht <= 0.5 * model->ht_max && curTimePeriod == 1)))
        model->ht = model->ht * 1.5;
    else if (DECREASE_STEP)
        model->ht = model->ht / 1.5;

	/*if (model->ht <= model->ht_max * 5.0 && iterations < 4 && err_newton_first < 1.0 && curTimePeriod == 0)
		model->ht = model->ht * 1.5;
	else if (model->ht <= model->ht_max && iterations < 5 && err_newton_first < 1.0 && curTimePeriod == 0)
		model->ht = model->ht * 1.5;
    else if (model->ht <= model->ht_max / 2.0 && iterations < 4 && err_newton_first < 1.0 && curTimePeriod == 1)
        model->ht = model->ht * 1.5;
	else if (iterations > 5 && model->ht > model->ht_min)
		model->ht = model->ht / 1.5;*/

	//if (cur_t + model->ht > model->period[curTimePeriod])
	//	model->ht = model->period[curTimePeriod] - cur_t;

	cur_t += model->ht;
}
void AcidRecFracSolver::start()
{
	step_idx = 0;

	fillIndices();
	solver.Init(strNum, 1.e-15, 1.e-15);

	model->setPeriod(curTimePeriod);

	while (cur_t < Tt)
	{
        model->snapshot_all(step_idx++);
        cout << "---------------------NEW TIME STEP---------------------" << endl;
		control();
        cout << setprecision(6);
        cout << "time = " << cur_t * t_dim / 3600.0 << endl;
		cfl_x = model->max_vel_x * model->ht / (model->props_frac.l2 / model->cellsNum_x);
		cfl_y = model->max_vel_y * model->ht / (model->props_frac.w2 / model->cellsNum_y_frac);
		model->max_vel_x = model->max_vel_y = 0.0;
		cfl = cfl_x + cfl_y;
		std::cout << "cfl_x = " << cfl_x << "\tcfl_y = " << cfl_y << "\tcfl = " << cfl << std::endl;
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
bool AcidRecFracSolver::doNextSmartStep()
{
    init_step_res.clear();
    final_step_res.clear();
    iter_num.clear();
    if (!solveSmartStep())
        return false;
	model->calculateTrans();
    return true;
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
void AcidRecFracSolver::checkVariables()
{
	for (auto cell : model->cells_poro)
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
bool AcidRecFracSolver::solveSmartStep()
{
	int cellIdx, varIdx;
	err_newton = 1.0;
	averValue(averValPrev);
	std::fill(dAverVal.begin(), dAverVal.end(), 1.0);
	iterations = 0;
	bool isHarder = (curTimePeriod == 0) ? false : true;
    bool isInit = false;// (cur_t < 0.001 / t_dim) ? false : true;

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
		solver.Solve(PRECOND::ILU_SIMPLE, isInit);
        init_step_res.push_back(solver.init_res);
        final_step_res.push_back(solver.final_res);
        iter_num.push_back(solver.iter_num);
		
        if (init_step_res.back() >= 2.0 * final_step_res.back())
        {
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
        }
        else
            return false;
		iterations++;
	}

	checkVariables();

	cout << "Newton Iterations = " << iterations << endl;
    return true;
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