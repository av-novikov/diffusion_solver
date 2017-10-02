#include "model/Bingham1d/Bingham1d.hpp"
#include "util/utils.h"

#include "adolc/drivers/drivers.h"
#include "adolc/adolc.h"

using namespace std;
using namespace bing1d;

Bingham1d::Bingham1d()
{
	x = new double[stencil * Variable::size];
	y = new double[Variable::size];

	jac = new double*[var_size];
	for (int i = 0; i < var_size; i++)
		jac[i] = new double[stencil * Variable::size];
}
Bingham1d::~Bingham1d()
{
	delete x;
	delete y;

	for (int i = 0; i < var_size; i++)
		delete[] jac[i];
	delete[] jac;
}
void Bingham1d::setProps(Properties& props)
{
	leftBoundIsRate = props.leftBoundIsRate;
	rightBoundIsPres = props.rightBoundIsPres;

	// Setting grid properties
	r_w = props.r_w;
	r_e = props.r_e;
	cellsNum_r = props.cellsNum_r;
	cellsNum = cellsNum_r + 2;

	props_sk = props.props_sk;
	props_sk.perm_r = MilliDarcyToM2(props_sk.perm_r);
	props_sk.perm_z = MilliDarcyToM2(props_sk.perm_z);

	perfIntervals = props.perfIntervals;

	periodsNum = props.timePeriods.size();
	for (int i = 0; i < periodsNum; i++)
	{
		period.push_back(props.timePeriods[i]);
		if (leftBoundIsRate)
			rate.push_back(props.rates[i] / 86400.0);
		else
			pwf.push_back(props.pwf[i]);

		if (props_sk.radiuses_eff[i] > props.r_w)
			props_sk.perms_eff.push_back(MilliDarcyToM2(props.props_sk.perm_r * log(props.props_sk.radiuses_eff[i] / props.r_w) / (log(props.props_sk.radiuses_eff[i] / props.r_w) + props.props_sk.skins[i])));
		else
			props_sk.perms_eff.push_back(MilliDarcyToM2(props.props_sk.perm_r));
	}

	// Temporal properties
	ht = props.ht;
	ht_min = props.ht_min;
	ht_max = props.ht_max;

	// Oil properties
	props_oil = props.props_oil;
	props_oil.visc = cPToPaSec(props_oil.visc);

	makeDimLess();

	props_oil.u_dp = setDataset(props.u_dp_dimless, props_oil.d / props_oil.tau0, 1);
}
void Bingham1d::makeDimLess()
{
	// Main units
	R_dim = r_w;
	t_dim = 3600.0;
	P_dim = props_sk.p_init;

	// Temporal properties
	ht /= t_dim;
	ht_min /= t_dim;
	ht_max /= t_dim;

	// Grid properties
	r_w /= R_dim;
	r_e /= R_dim;

	// Skeleton properties
	props_sk.perm_r /= (R_dim * R_dim);
	props_sk.perm_z /= (R_dim * R_dim);
	props_sk.beta /= (1.0 / P_dim);
	props_sk.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_sk.height /= R_dim;
	props_sk.p_init /= P_dim;
	props_sk.p_out /= P_dim;

	for (int j = 0; j < periodsNum; j++)
	{
		props_sk.perms_eff[j] /= (R_dim * R_dim);
		props_sk.radiuses_eff[j] /= R_dim;
	}

	Q_dim = R_dim * R_dim * R_dim / t_dim;
	for (int i = 0; i < periodsNum; i++)
	{
		period[i] /= t_dim;
		if (leftBoundIsRate)
			rate[i] /= Q_dim;
		else
			pwf[i] /= P_dim;
	}

	// Oil properties
	props_oil.visc /= (P_dim * t_dim);
	props_oil.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_oil.beta /= (1.0 / P_dim);
	props_oil.p_ref /= P_dim;

	props_oil.tau0 /= (P_dim);
	props_oil.d /= (R_dim);
}
void Bingham1d::buildGridLog()
{
	cells.reserve(cellsNum);

	Volume = 0.0;
	int counter = 0;

	double r_prev = r_w;
	double logMax = log(r_e / r_w);
	double logStep = logMax / (double)cellsNum_r;

	double hr = r_prev * (exp(logStep) - 1.0);
	double cm = r_w;

	cells.push_back(Cell(counter++, cm, 0.0, props_sk.height));
	cm += hr / 2.0;
	for (int i = 0; i < cellsNum_r; i++)
	{
		cm = r_prev * (1.0 + exp(logStep)) / 2.0;
		hr = r_prev * (exp(logStep) - 1.0);
		cells.push_back(Cell(counter++, cm, hr, props_sk.height));

		Volume += cells[cells.size() - 1].V;
		r_prev = r_prev * exp(logStep);
	}
	cells.push_back(Cell(counter++, r_e, 0.0, props_sk.height));
}
void Bingham1d::setPeriod(int period)
{
	if (leftBoundIsRate)
	{
		Q_sum = rate[period];

		if (period == 0 || rate[period - 1] < EQUALITY_TOLERANCE)
				Qcell[0] = Q_sum;
		else
			Qcell[0] *= (Q_sum / rate[period - 1]);
	}
	else
	{
		Pwf = pwf[period];
		Q_sum = 0.0;
	}

	props_sk.radius_eff = props_sk.radiuses_eff[period];
	props_sk.perm_eff = props_sk.perms_eff[period];
	props_sk.skin = props_sk.skins[period];
}
void Bingham1d::setPerforated()
{
	Qcell[0] = 0.0;
	height_perf = props_sk.height;
}
void Bingham1d::setInitialState()
{
	vector<Cell>::iterator it;
	for (it = cells.begin(); it != cells.end(); ++it)
		it->u_prev.p = it->u_iter.p = it->u_next.p = props_sk.p_init;
}
double Bingham1d::getRate(int cur)
{
	int neighbor = cur + 1;
	Variable& upwd = cells[getUpwindIdx(cur, neighbor)].u_next;
	Variable& next = cells[cur].u_next;
	double dp = cells[neighbor].u_next.p - next.p;
	return getTrans(cells[cur], cells[neighbor]) / props_oil.visc / props_oil.getB(next.p).value() *
		dp * props_oil.getU(fabs(dp / (cells[neighbor].r - cells[cur].r))).value();
}

void Bingham1d::solve_eqMiddle(const Cell& cell)
{
	trace_on(mid);
	adouble h[Variable::size];
	TapeVariable var[stencil];
	const Variable& prev = cell.u_prev;

	for (int i = 0; i < stencil; i++)
		var[i].p <<= x[i * Variable::size];

	const TapeVariable& next = var[0];
	h[0] = props_sk.getPoro(next.p) * props_oil.getDensity(next.p) -
		props_sk.getPoro(prev.p) * props_oil.getDensity(prev.p);

	int neighbor[2];
	getNeighborIdx(cell.num, neighbor);
	for (int i = 0; i < 2; i++)
	{
		Cell& beta = cells[neighbor[i]];
		const TapeVariable& nebr = var[i + 1];

		h[0] += ht / cell.V * getTrans(cell, beta) * (next.p - nebr.p) *
			linearInterp(props_oil.getDensity(next.p), cell, props_oil.getDensity(nebr.p), beta) / 
			props_oil.visc * props_oil.getU(fabs((next.p - nebr.p) / (cell.r - beta.r)));
	}
	h[0] >>= y[0];

	trace_off();
}
void Bingham1d::solve_eqLeft(const Cell& cell)
{
	const Cell& beta = cells[1];

	trace_on(left);
	adouble h[var_size];
	TapeVariable var[Lstencil];
	for (int i = 0; i < Lstencil; i++)
		var[i].p <<= x[i * Variable::size];

	const TapeVariable& next = var[0];
	const TapeVariable& nebr = var[1];

	adouble leftIsRate = leftBoundIsRate;
	condassign(h[0], leftIsRate,
		getTrans(cell, beta) / props_oil.visc / props_oil.getB(next.p).value() *
		(nebr.p - next.p) * props_oil.getU(fabs((next.p - nebr.p) / (cell.r - beta.r))) - Qcell[cell.num],
		next.p - Pwf);

	h[0] >>= y[0];

	trace_off();
}
void Bingham1d::solve_eqRight(const Cell& cell)
{
	const Cell& beta = cells[cellsNum - 2];

	trace_on(right);
	adouble h[var_size];
	TapeVariable var[Rstencil];
	for (int i = 0; i < Rstencil; i++)
		var[i].p <<= x[i * Variable::size];

	const TapeVariable& next = var[0];
	const TapeVariable& nebr = var[1];

	adouble rightIsPres = rightBoundIsPres;
	condassign(h[0], rightIsPres, next.p - (adouble)(props_sk.p_out), next.p - (adouble)(nebr.p));

	h[0] >>= y[0];

	trace_off();
}
void Bingham1d::setVariables(const Cell& cell)
{
	if (cell.num == 0)
	{
		// Left
		const Variable& next = cells[cell.num].u_next;
		const Variable& nebr1 = cells[cell.num + 1].u_next;

		for (int i = 0; i < Variable::size; i++)
		{
			x[i] = next.values[i];
			x[Variable::size + i] = nebr1.values[i];
		}

		solve_eqLeft(cell);
		jacobian(left, var_size, Variable::size * Lstencil, x, jac);
	}
	else if (cell.num == cellsNum_r + 1)
	{
		// Right
		const Variable& next = cells[cell.num].u_next;
		const Variable& nebr = cells[cell.num - 1].u_next;

		for (int i = 0; i < Variable::size; i++)
		{
			x[i] = next.values[i];
			x[Variable::size + i] = nebr.values[i];
		}

		solve_eqRight(cell);
		jacobian(right, var_size, Variable::size * Rstencil, x, jac);
	}
	else
	{
		// Middle
		const Variable& next = cells[cell.num].u_next;
		int neighbor[2];
		getNeighborIdx(cell.num, neighbor);

		for (int i = 0; i < Variable::size; i++)
		{
			x[i] = next.values[i];

			for (int j = 0; j < stencil - 1; j++)
			{
				const Variable& nebr = cells[neighbor[j]].u_next;
				x[(j + 1) * Variable::size + i] = nebr.values[i];
			}
		}

		solve_eqMiddle(cell);
		jacobian(mid, var_size, Variable::size * stencil, x, jac);
	}
}