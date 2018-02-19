#include "model/Acid/1d/Acid1d.hpp"

#include <cassert>

#include "adolc/drivers/drivers.h"
#include "adolc/adolc.h"

using namespace acid1d;

double acid1d::Component::R = 8.3144598;
double acid1d::Component::p_std = 101325.0;

Acid1d::Acid1d()
{
	x = new double[stencil * Variable::size];
	y = new double[var_size];

	jac = new double*[var_size];
	for (int i = 0; i < var_size; i++)
		jac[i] = new double[stencil * Variable::size];
}
Acid1d::~Acid1d()
{
	delete x;
	delete y;

	for (int i = 0; i < var_size; i++)
		delete[] jac[i];
	delete[] jac;
}
void Acid1d::setProps(Properties& props)
{
	props_sk = props.props_sk;
	setBasicProps(props);
	for (int i = 0; i < periodsNum; i++)
		xas.push_back(props.xa[i]);

	props_w = props.props_w;
	props_w.visc = cPToPaSec(props_w.visc);

	props_o = props.props_o;
	props_o.visc = cPToPaSec(props_o.visc);

	props_g = props.props_g;
	props_g.visc = cPToPaSec(props_g.visc);
	props_g.co2.mol_weight = gramToKg(props_g.co2.mol_weight);

	props_o.gas_dens_stc = props_g.co2.rho_stc;

	for (auto& comp : reac.comps)
		comp.mol_weight = gramToKg(comp.mol_weight);

	makeBasicDimLess();
	makeDimLess();

	// Data sets
	props_o.b = setDataset(props.B_oil, P_dim / BAR_TO_PA, 1.0);
	props_o.Rs = setDataset(props.Rs, P_dim / BAR_TO_PA, 1.0);
	props_g.rho = setDataset(props.rho_co2, P_dim / BAR_TO_PA, (P_dim * t_dim * t_dim / R_dim / R_dim));
}
void Acid1d::makeDimLess()
{
	T_dim = Component::T;

	Component::p_std /= P_dim;
	Component::R /= (P_dim * R_dim * R_dim * R_dim / T_dim);
	Component::T /= T_dim;

	props_sk.p_sat /= P_dim;

	props_w.visc /= (P_dim * t_dim);
	props_w.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_w.beta /= (1.0 / P_dim);
	props_o.visc /= (P_dim * t_dim);
	props_o.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_o.gas_dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_o.beta /= (1.0 / P_dim);
	props_o.p_ref /= P_dim;
	props_g.visc /= (P_dim * t_dim);
	props_g.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_g.co2.mol_weight /= (P_dim * t_dim * t_dim * R_dim);
	props_g.co2.rho_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);

	reac.activation_energy /= (P_dim * R_dim * R_dim * R_dim);
	reac.surf_init /= (1.0 / R_dim);
	reac.reaction_const /= (R_dim / t_dim);
	for (auto& comp : reac.comps)
	{
		comp.rho_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
		comp.mol_weight /= (P_dim * t_dim * t_dim * R_dim);
	}
}
void Acid1d::setPeriod(int period)
{
	if (leftBoundIsRate)
	{
		Q_sum = rate[period];

		if (period == 0 || rate[period - 1] < EQUALITY_TOLERANCE) {
			std::map<int, double>::iterator it;
			for (it = Qcell.begin(); it != Qcell.end(); ++it)
				it->second = Q_sum;
		}
		else {
			std::map<int, double>::iterator it;
			for (it = Qcell.begin(); it != Qcell.end(); ++it)
				it->second = it->second * Q_sum / rate[period - 1];
		}
	}
	else
	{
		Pwf = pwf[period];
		Q_sum = 0.0;
	}

	xa = xas[period];

	props_sk.radius_eff = props_sk.radiuses_eff[period];
	props_sk.perm_eff = props_sk.perms_eff[period];
	props_sk.skin = props_sk.skins[period];
}
void Acid1d::setInitialState()
{
	vector<Cell>::iterator it;
	for (it = cells.begin(); it != cells.end(); ++it)
	{
		it->u_prev.m = it->u_iter.m = it->u_next.m = props_sk.m_init;
		it->u_prev.p = it->u_iter.p = it->u_next.p = props_sk.p_init;
		it->u_prev.sw = it->u_iter.sw = it->u_next.sw = props_sk.sw_init;
		it->u_prev.xa = it->u_iter.xa = it->u_next.xa = props_sk.xa_init;
		it->u_prev.xw = it->u_iter.xw = it->u_next.xw = props_sk.xw_init;
		it->u_prev.xs = it->u_iter.xs = it->u_next.xs = 0.0;
	}
}
double Acid1d::getRate(int cur)
{
	const Cell& cell = cells[cur];
	const Cell& beta = cells[cur + 1];
	const Variable& next = cell.u_next;
	const Variable& nebr = beta.u_next;

	return props_w.getDensity(next.p, next.xa, next.xw, next.xs).value() / props_w.getDensity(Component::p_std, next.xa, next.xw, next.xs).value() *
		getTrans(cell, next.m, beta, nebr.m).value() / props_w.getViscosity(next.p, next.xa, next.xw, next.xs).value() * 
		(nebr.p - next.p);
};

void Acid1d::solve_eqMiddle(const Cell& cell)
{
	trace_on(mid);
	adouble h[var_size], tmp;
	TapeVariable var[stencil];
	const Variable& prev = cell.u_prev;

	for (int i = 0; i < stencil; i++)
	{
		var[i].m <<= x[i * Variable::size];
		var[i].p <<= x[i * Variable::size + 1];
		var[i].sw <<= x[i * Variable::size + 2];
		var[i].xw <<= x[i * Variable::size + 3];
		var[i].xa <<= x[i * Variable::size + 4];
		var[i].xs <<= x[i * Variable::size + 5];
	}

	TapeVariable& next = var[0];
	adouble rate = getReactionRate(next);

	h[0] = (1.0 - next.m) * props_sk.getDensity(next.p) -
			(1.0 - prev.m) * props_sk.getDensity(prev.p) - 
			ht * reac.indices[REACTS::CALCITE] * reac.comps[REACTS::CALCITE].mol_weight * rate;
	/*h[1] = next.m * (1.0 - next.sw - next.so) * props_g.getDensity(next.p) -
		prev.m * (1.0 - prev.sw - prev.so) * props_g.getDensity(prev.p) -
		ht * reac.indices[REACTS::CO2] * reac.comps[REACTS::CO2].mol_weight * rate;

	condassign(h[1], satur,
				next.m * ((1.0 - next.sw - next.so) * props_g.getDensity(next.p) +
			next.so * props_o.gas_dens_stc * props_o.getRs(next.p, next.p_bub, satur) / props_o.getB(next.p, next.p_bub, satur)) -
				prev.m * ((1.0 - prev.sw - prev.so) * props_g.getDensity(prev.p) + 
			prev.so * props_o.gas_dens_stc * props_o.getRs(prev.p, prev.p_bub, prev.SATUR) / props_o.getB(prev.p, prev.p_bub, prev.SATUR)) -
				ht * reac.indices[REACTS::CO2] * reac.comps[REACTS::CO2].mol_weight * rate,
				next.m * next.so * props_o.gas_dens_stc * props_o.getRs(next.p, next.p_bub, satur) / props_o.getB(next.p, next.p_bub, satur) -
			prev.m * prev.so * props_o.gas_dens_stc * props_o.getRs(prev.p, prev.p_bub, prev.SATUR) / props_o.getB(prev.p, prev.p_bub, prev.SATUR) - 
				ht * reac.indices[REACTS::CO2] * reac.comps[REACTS::CO2].mol_weight * rate);*/
	h[1] = next.m * next.sw * props_w.getDensity(next.p, next.xa, next.xw, next.xs) -
		prev.m * prev.sw * props_w.getDensity(prev.p, prev.xa, prev.xw, prev.xs) -
		ht * (reac.indices[REACTS::ACID] * reac.comps[REACTS::ACID].mol_weight +
			reac.indices[REACTS::WATER] * reac.comps[REACTS::WATER].mol_weight +
			reac.indices[REACTS::SALT] * reac.comps[REACTS::SALT].mol_weight) * rate;
	h[2] = next.m * (1.0 - next.sw) * props_o.getDensity(next.p) -
		prev.m * (1.0 - prev.sw) * props_o.getDensity(prev.p);
	h[3] = next.m * next.sw * props_w.getDensity(next.p, next.xa, next.xw, next.xs) * next.xw -
		prev.m * prev.sw * props_w.getDensity(prev.p, prev.xa, prev.xw, prev.xw) * prev.xw -
		ht * reac.indices[REACTS::WATER] * reac.comps[REACTS::WATER].mol_weight * rate;
	h[4] = next.m * next.sw * props_w.getDensity(next.p, next.xa, next.xw, next.xs) * next.xa -
		prev.m * prev.sw * props_w.getDensity(prev.p, prev.xa, prev.xw, prev.xs) * prev.xa -
		ht * reac.indices[REACTS::ACID] * reac.comps[REACTS::ACID].mol_weight * rate;
	h[5] = next.m * next.sw * props_w.getDensity(next.p, next.xa, next.xw, next.xs) * next.xs -
		prev.m * prev.sw * props_w.getDensity(prev.p, prev.xa, prev.xw, prev.xs) * prev.xs -
		ht * reac.indices[REACTS::SALT] * reac.comps[REACTS::SALT].mol_weight * rate;

	int neighbor[4];
	getNeighborIdx(cell.num, neighbor);
	for (int i = 0; i < 2; i++)
	{
		const Cell& beta = cells[neighbor[i]];
		const int upwd_idx = (getUpwindIdx(cell.num, neighbor[i]) == cell.num) ? 0 : i + 1;
		const TapeVariable& nebr = var[i + 1];
		TapeVariable& upwd = var[upwd_idx];
		adouble dens_w = getAverage(props_w.getDensity(next.p, next.xa, next.xw, next.xs), cell, 
									props_w.getDensity(nebr.p, nebr.xa, nebr.xw, nebr.xs), beta);
		adouble dens_o = getAverage(props_o.getDensity(next.p), cell, 
									props_o.getDensity(nebr.p), beta);
		adouble buf_w = ht / cell.V * getTrans(cell, next.m, beta, nebr.m) * (next.p - nebr.p) *
						dens_w * props_w.getKr(upwd.sw, &props_sk) / props_w.getViscosity(upwd.p, upwd.xa, upwd.xw, upwd.xs);
		adouble buf_o = ht / cell.V * getTrans(cell, next.m, beta, nebr.m) * (next.p - nebr.p) *
						dens_o * props_o.getKr(upwd.sw, &props_sk) / props_o.getViscosity(upwd.p);

		h[1] += buf_w;
		h[2] += buf_o;
		h[3] += buf_w * upwd.xw;
		h[4] += buf_w * upwd.xa;
		h[5] += buf_w * upwd.xs;
	}

	//for (int i = 0; i < var_size; i++)
	//	h[i] /= P_dim;

	for (int i = 0; i < var_size; i++)
		h[i] >>= y[i];

	trace_off();
}
void Acid1d::solve_eqLeft(const Cell& cell)
{
	const Cell& beta1 = cells[cell.num + 1];
	const Cell& beta2 = cells[cell.num + 2];

	trace_on(left);
	adouble h[var_size];
	TapeVariable var[Lstencil];
	for (int i = 0; i < Lstencil; i++)
	{
		var[i].m <<= x[i * Variable::size];
		var[i].p <<= x[i * Variable::size + 1];
		var[i].sw <<= x[i * Variable::size + 2];
		var[i].xw <<= x[i * Variable::size + 3];
		var[i].xa <<= x[i * Variable::size + 4];
		var[i].xs <<= x[i * Variable::size + 5];
	}

	const adouble leftIsRate = leftBoundIsRate;
	TapeVariable& next = var[0];
	TapeVariable& nebr1 = var[1];
	TapeVariable& nebr2 = var[2];
	const Variable& prev = cell.u_prev;

	adouble isPerforated, leftPresPerforated, leftPresNonPerforated;
	auto it = Qcell.find(cell.num);
	double rate;
	if (it != Qcell.end())
	{
		leftPresPerforated = !leftBoundIsRate;
		leftPresNonPerforated = false;
		isPerforated = true;
		rate = Qcell[cell.num];
	}
	else
	{
		isPerforated = false;
		leftPresPerforated = false;
		leftPresNonPerforated = !leftBoundIsRate;
		rate = 0.0;
	}

	h[0] = (1.0 - next.m) * props_sk.getDensity(next.p) - (1.0 - prev.m) * props_sk.getDensity(prev.p) -
		ht * reac.indices[REACTS::CALCITE] * reac.comps[REACTS::CALCITE].mol_weight * getReactionRate(next);
	condassign(h[1], leftIsRate, props_w.getDensity(next.p, next.xa, next.xw, next.xs) * getTrans(cell, next.m, beta1, nebr1.m) /
		props_w.getViscosity(next.p, next.xa, next.xw, next.xs) * (nebr1.p - next.p) +
		props_w.getDensity(Component::p_std, next.xa, next.xw, next.xs) * rate,
									next.p - Pwf);
	h[2] = next.sw - (1.0 - props_sk.s_oc);
	h[3] = next.xw - (1.0 - xa);
	h[4] = next.xa - xa;
	h[5] = next.xs;

	//for (int i = 0; i < var_size; i++)
	//	h[i] /= P_dim;

	for (int i = 0; i < var_size; i++)
		h[i] >>= y[i];

	trace_off();
}
void Acid1d::solve_eqRight(const Cell& cell)
{
	trace_on(right);
	adouble h[var_size];
	TapeVariable var[Rstencil];
	adouble rightIsPres = rightBoundIsPres;
	for (int i = 0; i < Lstencil; i++)
	{
		var[i].m <<= x[i * Variable::size];
		var[i].p <<= x[i * Variable::size + 1];
		var[i].sw <<= x[i * Variable::size + 2];
		var[i].xw <<= x[i * Variable::size + 3];
		var[i].xa <<= x[i * Variable::size + 4];
		var[i].xs <<= x[i * Variable::size + 5];
	}

	TapeVariable& next = var[0];
	TapeVariable& nebr = var[1];
	const Variable& prev = cell.u_prev;
	const Skeleton_Props& props = *cell.props;

	h[0] = next.m - nebr.m;
	condassign(h[1], rightIsPres, next.p - (adouble)(props_sk.p_out), next.p - nebr.p);
	h[2] = next.sw - nebr.sw;
	h[3] = next.xw - nebr.xw;
	h[4] = next.xa - nebr.xa;
	h[5] = next.xs - nebr.xs;
	
	//for (int i = 0; i < var_size; i++)
	//	h[i] /= P_dim;

	for (int i = 0; i < var_size; i++)
		h[i] >>= y[i];

	trace_off();
}
void Acid1d::setVariables(const Cell& cell)
{
	if (cell.type == Type::WELL_LAT)
	{
		// Left
		const Variable& next = cell.u_next;
		const Variable& nebr1 = cells[cell.num + 1].u_next;
		const Variable& nebr2 = cells[cell.num + 2].u_next;
		for (int i = 0; i < Variable::size; i++)
		{
			x[i] = next.values[i];
			x[Variable::size + i] = nebr1.values[i];
			x[2 * Variable::size + i] = nebr2.values[i];
		}

		solve_eqLeft(cell);
		jacobian(left, var_size, Variable::size * Lstencil, x, jac);
	}
	else if (cell.type == Type::RIGHT)
	{
		// Right
		const Variable& next = cells[cell.num].u_next;
		const Variable& nebr1 = cells[cell.num - 1].u_next;
		const Variable& nebr2 = cells[cell.num - 2].u_next;
		for (int i = 0; i < Variable::size; i++)
		{
			x[i] = next.values[i];
			x[Variable::size + i] = nebr1.values[i];
			x[2 * Variable::size + i] = nebr2.values[i];
		}

		solve_eqRight(cell);
		jacobian(right, var_size, Variable::size * Rstencil, x, jac);
	}
	else if (cell.type == Type::MIDDLE)
	{
		// Middle
		const Variable& next = cells[cell.num].u_next;
		int neighbor[stencil - 1];
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