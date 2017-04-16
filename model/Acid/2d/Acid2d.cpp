#include "model/Acid/2d/Acid2d.hpp"

#include <cassert>

#include "adolc/drivers/drivers.h"
#include "adolc/adolc.h"

using namespace acid2d;

double acid2d::Component::R = 8.3144598;
double acid2d::Component::p_std = 101325.0;

Acid2d::Acid2d()
{
	x = new double[stencil * Variable::size];
	y = new double[Variable::size];

	jac = new double*[Variable::size];
	for (int i = 0; i < Variable::size; i++)
		jac[i] = new double[stencil * Variable::size];
}
Acid2d::~Acid2d()
{
	delete x;
	delete y;

	for (int i = 0; i < Variable::size; i++)
		delete[] jac[i];
	delete[] jac;
}
void Acid2d::setProps(Properties& props)
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

	for (auto& comp : reac.comps)
		comp.mol_weight = gramToKg(comp.mol_weight);

	makeBasicDimLess();
	makeDimLess();

	// Data sets
	//props_w.kr = setDataset(props.kr_l, 1.0, 1.0);
	//props_g.kr = setDataset(props.kr_g, 1.0, 1.0);
	//props_oil.b = setDataset(props.B_oil, P_dim / BAR_TO_PA, 1.0);
	//props_gas.b = setDataset(props.B_gas, P_dim / BAR_TO_PA, 1.0);
	//props_oil.Rs = setDataset(props.Rs, P_dim / BAR_TO_PA, 1.0);
}
void Acid2d::makeDimLess()
{
	T_dim = Component::T;

	Component::p_std /= P_dim;
	Component::R /= (P_dim * R_dim * R_dim * R_dim / T_dim);
	Component::T /= T_dim;

	props_w.visc /= (P_dim * t_dim);
	props_w.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_w.beta /= (1.0 / P_dim);
	props_o.visc /= (P_dim * t_dim);
	props_o.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_o.beta /= (1.0 / P_dim);
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
void Acid2d::setPeriod(int period)
{
	if (leftBoundIsRate)
	{
		Q_sum = rate[period];

		if (period == 0 || rate[period - 1] < EQUALITY_TOLERANCE) {
			std::map<int, double>::iterator it;
			for (it = Qcell.begin(); it != Qcell.end(); ++it)
				it->second = Q_sum * cells[it->first].hz / height_perf;
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

	for (int i = 0; i < skeletonsNum; i++)
	{
		props_sk[i].radius_eff = props_sk[i].radiuses_eff[period];
		props_sk[i].perm_eff = props_sk[i].perms_eff[period];
		props_sk[i].skin = props_sk[i].skins[period];
	}
}
void Acid2d::setInitialState()
{
	vector<Cell>::iterator it;
	for (it = cells.begin(); it != cells.end(); ++it)
	{
		Skeleton_Props* props = &props_sk[getSkeletonIdx(*it)];
		it->u_prev.m = it->u_iter.m = it->u_next.m = props->m_init;
		it->u_prev.p = it->u_iter.p = it->u_next.p = props->p_init;
		it->u_prev.sw = it->u_iter.sw = it->u_next.sw = props->sw_init;
		it->u_prev.so = it->u_iter.so = it->u_next.so = props->so_init;
		it->u_prev.xa = it->u_iter.xa = it->u_next.xa = props->xa_init;
		it->u_prev.xw = it->u_iter.xw = it->u_next.xw = props->xw_init;

		it->props = props;
	}
}
double Acid2d::getRate(int cur)
{
	const Cell& cell = cells[cur];
	const Cell& beta = cells[cur + cellsNum_z + 2];
	const Variable& next = cell.u_next;
	const Variable& nebr = beta.u_next;

	return props_w.getDensity(next.p, next.xa, next.xw).value() / props_w.getDensity(Component::p_std, next.xa, next.xw).value() *
		getTrans(cell, next.m, beta, nebr.m).value() / props_w.getViscosity(next.p, next.xa, next.xw).value() * 
		(nebr.p - next.p);
};

void Acid2d::solve_eqMiddle(const Cell& cell)
{
	const Skeleton_Props& props = *cell.props;

	trace_on(mid);
	adouble h[var_size], tmp;
	TapeVariable var[stencil];
	const Variable& prev = cell.u_prev;

	for (int i = 0; i < stencil; i++)
	{
		var[i].m <<= x[i * Variable::size];
		var[i].p <<= x[i * Variable::size + 1];
		var[i].sw <<= x[i * Variable::size + 2];
		var[i].so <<= x[i * Variable::size + 3];
		var[i].xa <<= x[i * Variable::size + 4];
		var[i].xw <<= x[i * Variable::size + 5];
	}

	TapeVariable& next = var[0];
	adouble rate = getReactionRate(next, props);
	
	h[0] = (1.0 - next.m) * props.getDensity(next.p) -
			(1.0 - prev.m) * props.getDensity(prev.p) - 
			ht * reac.indices[REACTS::CALCITE] * reac.comps[REACTS::CALCITE].mol_weight * rate;
	h[1] = next.m * (1.0 - next.sw - next.so) * props_g.getDensity(next.p) -
		prev.m * (1.0 - prev.sw - prev.so) * props_g.getDensity(prev.p) -
		ht * reac.indices[REACTS::CO2] * reac.comps[REACTS::CO2].mol_weight * rate;
	h[2] = next.m * next.sw * props_w.getDensity(next.p, next.xa, next.xw) -
		prev.m * prev.sw * props_w.getDensity(prev.p, prev.xa, prev.xw) -
		ht * (reac.indices[REACTS::ACID] * reac.comps[REACTS::ACID].mol_weight +
			reac.indices[REACTS::WATER] * reac.comps[REACTS::WATER].mol_weight +
			reac.indices[REACTS::SALT] * reac.comps[REACTS::SALT].mol_weight) * rate;
	h[3] = next.m * next.so * props_o.getDensity(next.p) - prev.m * prev.so * props_o.getDensity(prev.p);
	h[4] = next.m * next.sw * props_w.getDensity(next.p, next.xa, next.xw) * next.xa -
		prev.m * prev.sw * props_w.getDensity(prev.p, prev.xa, prev.xw) * prev.xa -
		ht * reac.indices[REACTS::ACID] * reac.comps[REACTS::ACID].mol_weight * rate;
	h[5] = next.m * next.sw * props_w.getDensity(next.p, next.xa, next.xw) * next.xw -
		prev.m * prev.sw * props_w.getDensity(prev.p, prev.xa, prev.xw) * prev.xw -
		ht * reac.indices[REACTS::WATER] * reac.comps[REACTS::WATER].mol_weight * rate;

	int neighbor[4];
	getNeighborIdx(cell.num, neighbor);
	for (int i = 0; i < 4; i++)
	{
		const Cell& beta = cells[neighbor[i]];
		const int upwd_idx = (getUpwindIdx(cell.num, neighbor[i]) == cell.num) ? 0 : i + 1;
		const TapeVariable& nebr = var[i + 1];
		TapeVariable& upwd = var[upwd_idx];
		
		adouble dens_w = getAvarage(props_w.getDensity(next.p, next.xa, next.xw), cell, props_w.getDensity(nebr.p, nebr.xa, nebr.xw), beta);
		adouble dens_o = getAvarage(props_o.getDensity(next.p), cell, props_o.getDensity(nebr.p), beta);
		adouble dens_g = getAvarage(props_g.getDensity(next.p), cell, props_g.getDensity(nebr.p), beta);
		adouble buf = ht / cell.V * getTrans(cell, next.m, beta, nebr.m) * (next.p - nebr.p);
		adouble buf_w = buf * dens_w * props_w.getKr(upwd.sw, upwd.so, cells[upwd_idx].props) / props_w.getViscosity(upwd.p, upwd.xa, upwd.xw);

		h[1] += buf * dens_g * props_g.getKr(upwd.sw, upwd.so, cells[upwd_idx].props) / props_g.getViscosity(upwd.p);
		h[2] += buf_w;
		h[3] += buf * dens_o * props_o.getKr(upwd.sw, upwd.so, cells[upwd_idx].props) / props_o.getViscosity(upwd.p);
		h[4] += buf_w * upwd.xa;
		h[5] += buf_w * upwd.xw;
	}

	for (int i = 0; i < var_size; i++)
		h[i] >>= y[i];

	trace_off();
}
void Acid2d::solve_eqLeft(const Cell& cell)
{
	const Cell& beta = cells[cell.num + cellsNum_z + 2];

	trace_on(left);
	adouble h[var_size];
	TapeVariable var[Lstencil];
	for (int i = 0; i < Lstencil; i++)
	{
		var[i].m <<= x[i * Variable::size];
		var[i].p <<= x[i * Variable::size + 1];
		var[i].sw <<= x[i * Variable::size + 2];
		var[i].so <<= x[i * Variable::size + 3];
		var[i].xa <<= x[i * Variable::size + 4];
		var[i].xw <<= x[i * Variable::size + 5];
	}

	const adouble leftIsRate = leftBoundIsRate;
	TapeVariable& next = var[0];
	TapeVariable& nebr = var[1];
	const Variable& prev = cell.u_prev;
	const Skeleton_Props& props = *cell.props;

	h[0] = (1.0 - next.m) * props.getDensity(next.p) - (1.0 - prev.m) * props.getDensity(prev.p) -
		ht * reac.indices[REACTS::CALCITE] * getReactionRate(next, props);
	condassign(h[1], leftIsRate, props_w.getDensity(next.p, next.xa, next.xw) * getTrans(cell, next.m, beta, nebr.m) /
		props_w.getViscosity(next.p, next.xa, next.xw) * (nebr.p - next.p) -
		props_w.getDensity(Component::p_std, next.xa, next.xw) * Qcell[cell.num],
		next.p - Pwf);
	h[2] = next.sw - (1.0 - props.s_oc);
	h[3] = next.so - props.s_oc;
	h[4] = next.xa - xa;
	h[5] = next.xw - (1.0 - xa);

	for (int i = 0; i < var_size; i++)
		h[i] >>= y[i];

	trace_off();
}
void Acid2d::solve_eqRight(const Cell& cell)
{
	trace_on(right);
	adouble h[var_size];
	TapeVariable var[Rstencil];
	adouble rightIsPres = rightBoundIsPres;
	for (int i = 0; i < Rstencil; i++)
	{
		var[i].m <<= x[i * Variable::size];
		var[i].p <<= x[i * Variable::size + 1];
		var[i].sw <<= x[i * Variable::size + 2];
		var[i].so <<= x[i * Variable::size + 3];
		var[i].xa <<= x[i * Variable::size + 4];
		var[i].xw <<= x[i * Variable::size + 5];
	}

	TapeVariable& next = var[0];
	TapeVariable& nebr = var[1];
	const Variable& prev = cell.u_prev;
	const Skeleton_Props& props = *cell.props;

	h[0] = next.m - nebr.m;
	condassign(h[1], rightIsPres, next.p - (adouble)(cell.props->p_out), next.p - (adouble)(nebr.p));
	h[2] = next.sw - nebr.sw;
	h[3] = next.so - nebr.so;
	h[4] = next.xa - nebr.xa;
	h[5] = next.xw - nebr.xw;
	
	for (int i = 0; i < var_size; i++)
		h[i] >>= y[i];

	trace_off();
}
void Acid2d::solve_eqVertical(const Cell& cell)
{
	trace_on(vertical);
	adouble h[var_size];
	TapeVariable var[Vstencil];

	for (int i = 0; i < Vstencil; i++)
	{
		var[i].m <<= x[i * Variable::size];
		var[i].p <<= x[i * Variable::size + 1];
		var[i].sw <<= x[i * Variable::size + 2];
		var[i].so <<= x[i * Variable::size + 3];
		var[i].xa <<= x[i * Variable::size + 4];
		var[i].xw <<= x[i * Variable::size + 5];
	}

	const TapeVariable& next = var[0];
	const TapeVariable& nebr = var[1];

	h[0] = next.m - nebr.m;
	h[1] = next.p - nebr.p;
	h[2] = next.sw - nebr.sw;
	h[3] = next.so - nebr.so;
	h[4] = nebr.xa - nebr.xa;
	h[5] = nebr.xw - nebr.xw;

	for (int i = 0; i < var_size; i++)
		h[i] >>= y[i];

	trace_off();
}
void Acid2d::setVariables(const Cell& cell)
{
	if (cell.type == Type::WELL_LAT)
	{
		// Left
		const Variable& next = cell.u_next;
		const Variable& nebr1 = cells[cell.num + cellsNum_z + 2].u_next;
		const Variable& nebr2 = cells[cell.num + 2 * cellsNum_z + 4].u_next;
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
		const Variable& nebr1 = cells[cell.num - cellsNum_z - 2].u_next;
		const Variable& nebr2 = cells[cell.num - 2 * cellsNum_z - 4].u_next;
		for (int i = 0; i < Variable::size; i++)
		{
			x[i] = next.values[i];
			x[Variable::size + i] = nebr1.values[i];
			x[2 * Variable::size + i] = nebr2.values[i];
		}

		solve_eqRight(cell);
		jacobian(right, var_size, Variable::size * Rstencil, x, jac);
	}
	else if (cell.type == Type::TOP)
	{
		// Top
		const Variable& next = cells[cell.num].u_next;
		const Variable& nebr = cells[cell.num + 1].u_next;
		for (int i = 0; i < Variable::size; i++)
		{
			x[i] = next.values[i];
			x[Variable::size + i] = nebr.values[i];
		}

		solve_eqVertical(cell);
		jacobian(vertical, var_size, Variable::size * Vstencil, x, jac);
	}
	else if (cell.type == Type::BOTTOM)
	{
		// Bottom
		const Variable& next = cells[cell.num].u_next;
		const Variable& nebr = cells[cell.num - 1].u_next;
		for (int i = 0; i < Variable::size; i++)
		{
			x[i] = next.values[i];
			x[Variable::size + i] = nebr.values[i];
		}

		solve_eqVertical(cell);
		jacobian(vertical, var_size, Variable::size * Vstencil, x, jac);
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