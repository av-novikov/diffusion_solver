#include <iomanip>
#include <iostream>

#include "model/Acid/recfrac/RecFracProdSolver.hpp"

#include "adolc/sparse/sparsedrivers.h"
#include "adolc/drivers/drivers.h"

using std::endl;
using std::setprecision;
using std::ofstream;
using namespace acidrecfrac_prod;

RecFracProdSolver::RecFracProdSolver(RecFracProd* _model) : AbstractSolver<RecFracProd>(_model)
{
	strNum = model->cellsNum * var_size;
	x = new double[strNum];
	y = new double[strNum];

	ind_i = new int[stencil * var_size * strNum];
	ind_j = new int[stencil * var_size * strNum];
	cols = new int[strNum];
	a = new double[stencil * var_size * strNum];
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

	CONV_W2 = 1.e-5;		CONV_VAR = 1.e-8;
	MAX_ITER = 20;

    MAX_INIT_RES1 = 5.E-7;
	MAX_INIT_RES2 = 1.E-8;
    //MAX_INIT_RES1 = 4.E-8; 
    //MAX_INIT_RES2 = 1.E-9;

	MULT_UP = MULT_DOWN = 2;

	P.open("snaps/P.dat", ofstream::out);
	qcells.open("snaps/q_prod.dat", ofstream::out);
	pvd.open("snaps/AcidRecFrac_prod.pvd", std::ofstream::out);
	pvd << "<VTKFile type = \"Collection\" version = \"1.0\" byte_order = \"LittleEndian\" header_type = \"UInt64\">\n";
	pvd << "\t<Collection>\n";
}
RecFracProdSolver::~RecFracProdSolver()
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
void RecFracProdSolver::writeData()
{
	pvd << "\t\t<DataSet part=\"0\" timestep=\"";
	if (cur_t == 0.0)
		pvd << std::to_string(0.0);
	else
		pvd << cur_t * t_dim / 3600.0;
	pvd << "0\" file=\"AcidRecFrac_prod_" + std::to_string(step_idx) + ".vtu\"/>\n";

	double q = 0.0, s = 0.0;
	for (int i = 0; i < (model->cellsNum_y + 2) * (model->cellsNum_z + 2); i++)
	{
		const auto& cell = model->cells[i];
		if (cell.type == Type::FRAC_IN && cell.hy != 0.0 && cell.hz != 0.0)
		{
			q += model->getRate(cell);
			s += cell.hy;
		}
	}

	qcells << cur_t * t_dim / 3600.0;
	qcells << "\t" << (model->props_sk.back().p_out - model->Pwf) * model->P_dim / BAR_TO_PA;
	qcells << "\t" << s * model->R_dim;
	qcells << "\t" << 4.0 * q * model->Q_dim * 86400.0;
	qcells << "\t" << 4.0 * q * model->Q_dim * 86400.0 / ((model->props_sk.back().p_out - model->Pwf) * model->P_dim / BAR_TO_PA);
	qcells << endl;
}
void RecFracProdSolver::analyzeNewtonConvergence()
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
        //INCREASE_STEP = true;
        //DECREASE_STEP = false;
        if (init_step_res[i + 1] >= init_step_res[i])
        {
            INCREASE_STEP = false;
            if (iterations > 4)
                DECREASE_STEP = true;
            break;
        }
    }
}
void RecFracProdSolver::control()
{
	writeData();

	if (cur_t >= model->period[curTimePeriod])
	{
		curTimePeriod++;
		model->ht = model->ht_min;
		model->setPeriod(curTimePeriod);
	}

    analyzeNewtonConvergence();
    
    if(INCREASE_STEP && model->ht <= model->ht_max)
        model->ht = model->ht * MULT_UP;
    else if (DECREASE_STEP)
        model->ht = model->ht / MULT_DOWN;

	if (cur_t + model->ht > model->period[curTimePeriod])
		model->ht = model->period[curTimePeriod] - cur_t;

	cur_t += model->ht;
}
void RecFracProdSolver::start()
{
	step_idx = 0;

	fillIndices();
	solver.Init(strNum, 1.e-15, 1.e-8, 1.e+4);

	model->setPeriod(curTimePeriod);

	while (cur_t < Tt)
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
bool RecFracProdSolver::doNextSmartStep()
{
    init_step_res.clear();
    final_step_res.clear();
    iter_num.clear();
    if (!solveSmartStep())
        return false;
    return true;
}
void RecFracProdSolver::copySolution(const paralution::LocalVector<double>& sol)
{
	for (int i = 0; i < model->cellsNum; i++)
	{
		auto& cell = model->cells[i];
		for (int j = 0; j < var_size; j++)
			cell.u_next.values[j] += sol[i * var_size + j];
	}
}
void RecFracProdSolver::checkStability()
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
void RecFracProdSolver::checkVariables()
{
	/*for (auto cell : model->cells)
	{
		if (cell.u_next.s > 1.0)
			cell.u_next.s = 1.0;
		if (cell.u_next.s < 0.0)
			cell.u_next.s = 0.0;
	}*/
}
bool RecFracProdSolver::solveSmartStep()
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
		solver.Solve(PRECOND::ILU_SERIOUS, isInit);
        init_step_res.push_back(solver.init_res);
        final_step_res.push_back(solver.final_res);
        iter_num.push_back(solver.iter_num);
		
        if (init_step_res.back() >= 2.0 * final_step_res.back())
        {
            copySolution(solver.getSolution());

            checkStability();
            double err = convergance(cellIdx, varIdx);
            if (iterations == 0)
                err_newton_first = err_newton;

            averValue(averVal);
            for (int i = 0; i < var_size; i++)
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
void RecFracProdSolver::computeJac()
{
	trace_on(0);

	for (const auto& cell : model->cells)
	{
		auto& cur_x = model->x[cell.num];
		cur_x.p <<= cell.u_next.p;
		//cur_x.s <<= cell.u_next.s;

		x[cell.num * var_size] = cell.u_next.p;
		//x[cell.num * var_size + 1] = cell.u_next.s;
	}
	for (const auto& cell : model->cells)
	{
		const auto res = model->solveFrac(cell);
		model->h[cell.num * var_size] = res.p;
		//model->h[cell.num * var_size + 1] = res.s;
	}

	for (int i = 0; i < strNum; i++)
		model->h[i] >>= y[i];

	trace_off();
}
void RecFracProdSolver::fill()
{
	sparse_jac(0, strNum, strNum, repeat,
		x, &elemNum, (unsigned int**)(&ind_i), (unsigned int**)(&ind_j), &a, options);

	for(int i = 0; i < strNum; i++)
		rhs[i] = -y[i];
}