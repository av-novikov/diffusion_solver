#include "model/Acid/2dnit/Acid2dNIT.hpp"

#include <cassert>

#include "adolc/drivers/drivers.h"
#include "adolc/adolc.h"

using namespace acid2dnit;

double acid2dnit::Component::R = 8.3144598;
double acid2dnit::Component::p_std = 101325.0;

Acid2dNIT::Acid2dNIT()
{
	x = new double[stencil * Variable::size];
	y = new double[var_size];

	jac = new double*[var_size];
	for (int i = 0; i < var_size; i++)
		jac[i] = new double[stencil * Variable::size];
}
Acid2dNIT::~Acid2dNIT()
{
	delete[] x, y, h;
	delete[] var;

	for (int i = 0; i < var_size; i++)
		delete[] jac[i];
	delete[] jac;
}
void Acid2dNIT::setProps(Properties& props)
{
	props_sk = props.props_sk;
	

	rightBoundIsPres = props.rightBoundIsPres;

	// Setting grid properties
	r_w = props.r_w;
	r_e = props.r_e;
	cellsNum_r = props.cellsNum_r;
	cellsNum_z = props.cellsNum_z;
	cellsNum = (cellsNum_r + 2) * (cellsNum_z + 2);

	// Setting skeleton properties
	perfIntervals = props.perfIntervals;

	skeletonsNum = props.props_sk.size();
	checkSkeletons(props.props_sk);
	for (int j = 0; j < skeletonsNum; j++)
	{
		props_sk[j].perm_r = MilliDarcyToM2(props_sk[j].perm_r);
		props_sk[j].perm_z = MilliDarcyToM2(props_sk[j].perm_z);
	}

	periodsNum = props.timePeriods.size();
	rate.resize(periodsNum);
	pwf.resize(periodsNum);
	int rate_idx = 0, pres_idx = 0;
	for (int i = 0; i < periodsNum; i++)
	{
		period.push_back(props.timePeriods[i]);
		LeftBoundIsRate.push_back(props.LeftBoundIsRate[i]);
		xas.push_back(props.xa[i]);
		temps.push_back(props.temps[i]);

		if (LeftBoundIsRate.back())
		{
			rate[i] = props.rates[rate_idx++] / 86400.0;
			pwf[i] = 0.0;
		}
		else
		{
			pwf[i] = props.pwf[pres_idx++];
			rate[i] = 0.0;
		}
		for (int j = 0; j < skeletonsNum; j++)
		{
			if (props_sk[j].radiuses_eff[i] > props.r_w)
				props_sk[j].perms_eff.push_back(MilliDarcyToM2(props.props_sk[j].perm_r * log(props.props_sk[j].radiuses_eff[i] / props.r_w) / (log(props.props_sk[j].radiuses_eff[i] / props.r_w) + props.props_sk[j].skins[i])));
			else
				props_sk[j].perms_eff.push_back(MilliDarcyToM2(props.props_sk[j].perm_r));
		}
	}

	// Temporal properties
	ht = props.ht;
	ht_min = props.ht_min;
	ht_max = props.ht_max;

	alpha = props.alpha;
	depth_point = props.depth_point;

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

	// Main units
	R_dim = r_e / 100.0;
	t_dim = 3600.0;
	P_dim = props_sk[0].p_init;

	// Temporal properties
	ht /= t_dim;
	ht_min /= t_dim;
	ht_max /= t_dim;

	// Grid properties
	r_w /= R_dim;
	r_e /= R_dim;

	// Skeleton properties
	for (int i = 0; i < skeletonsNum; i++)
	{
		props_sk[i].perm_r /= (R_dim * R_dim);
		props_sk[i].perm_z /= (R_dim * R_dim);

		props_sk[i].beta /= (1.0 / P_dim);
		props_sk[i].dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
		props_sk[i].h1 = (props_sk[i].h1 - depth_point) / R_dim;
		props_sk[i].h2 = (props_sk[i].h2 - depth_point) / R_dim;
		props_sk[i].height /= R_dim;
		props_sk[i].p_init /= P_dim;
		props_sk[i].p_out /= P_dim;
		props_sk[i].p_ref /= P_dim;

		for (int j = 0; j < periodsNum; j++)
		{
			props_sk[i].perms_eff[j] /= (R_dim * R_dim);
			props_sk[i].radiuses_eff[j] /= R_dim;
		}
	}

	Q_dim = R_dim * R_dim * R_dim / t_dim;
	for (int i = 0; i < periodsNum; i++)
	{
		period[i] /= t_dim;
		rate[i] /= Q_dim;
		pwf[i] /= P_dim;
	}

	grav /= (R_dim / t_dim / t_dim);
	// Rest properties
	alpha /= t_dim;
	//depth_point = 0.0;

	makeDimLess();

	// Data sets
	props_o.b = setDataset(props.B_oil, P_dim / BAR_TO_PA, 1.0);
	props_o.Rs = setDataset(props.Rs, P_dim / BAR_TO_PA, 1.0);
	props_g.rho = setDataset(props.rho_co2, P_dim / BAR_TO_PA, (P_dim * t_dim * t_dim / R_dim / R_dim));
}
void Acid2dNIT::makeDimLess()
{
	//T_dim = Component::T;
	T_dim = props_sk[0].t_init;

	Component::p_std /= P_dim;
	Component::R /= (P_dim * R_dim * R_dim * R_dim / T_dim);
	Component::T /= T_dim;

	for (auto& sk : props_sk)
	{
		sk.p_sat /= P_dim;
		sk.t_init /= T_dim;
		sk.c = sk.c / R_dim / R_dim * T_dim * t_dim * t_dim;
		sk.lambda_r = sk.lambda_r * T_dim * t_dim / P_dim / R_dim / R_dim;
		sk.lambda_z = sk.lambda_z * T_dim * t_dim / P_dim / R_dim / R_dim;
	}

	for (int i = 0; i < periodsNum; i++)
		temps[i] /= T_dim;

	props_w.visc /= (P_dim * t_dim);
	props_w.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_w.beta /= (1.0 / P_dim);
	props_w.c = props_w.c / R_dim / R_dim * T_dim * t_dim * t_dim;
	props_w.lambda = props_w.lambda * T_dim * t_dim / P_dim / R_dim / R_dim;
	props_w.jt = props_w.jt * P_dim / T_dim;
	props_w.ad = props_w.ad * P_dim / T_dim;
	props_o.visc /= (P_dim * t_dim);
	props_o.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_o.gas_dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_o.beta /= (1.0 / P_dim);
	props_o.p_ref /= P_dim;
	props_o.c = props_o.c / R_dim / R_dim * T_dim * t_dim * t_dim;
	props_o.lambda = props_o.lambda * T_dim * t_dim / P_dim / R_dim / R_dim;
	props_o.jt = props_o.jt * P_dim / T_dim;
	props_o.ad = props_o.ad * P_dim / T_dim;
	props_g.visc /= (P_dim * t_dim);
	props_g.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_g.co2.mol_weight /= (P_dim * t_dim * t_dim * R_dim);
	props_g.co2.rho_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_g.c = props_g.c / R_dim / R_dim * T_dim * t_dim * t_dim;
	props_g.lambda = props_g.lambda * T_dim * t_dim / P_dim / R_dim / R_dim;
	props_g.jt = props_g.jt * P_dim / T_dim;
	props_g.ad = props_g.ad * P_dim / T_dim;

	reac.activation_energy /= (P_dim * R_dim * R_dim * R_dim);
	reac.surf_init /= (1.0 / R_dim);
	reac.reaction_const /= (R_dim / t_dim);
	for (auto& comp : reac.comps)
	{
		comp.rho_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
		comp.mol_weight /= (P_dim * t_dim * t_dim * R_dim);
	}
}
void Acid2dNIT::setPeriod(int period)
{
	leftBoundIsRate = LeftBoundIsRate[period];
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
	temp = temps[period];

	for (int i = 0; i < skeletonsNum; i++)
	{
		props_sk[i].radius_eff = props_sk[i].radiuses_eff[period];
		props_sk[i].perm_eff = props_sk[i].perms_eff[period];
		props_sk[i].skin = props_sk[i].skins[period];
	}
}
void Acid2dNIT::setInitialState()
{
	vector<Cell>::iterator it;
	for (it = cells.begin(); it != cells.end(); ++it)
	{
		Skeleton_Props* props = &props_sk[getSkeletonIdx(*it)];
		it->u_prev.m = it->u_iter.m = it->u_next.m = props->m_init;
		it->u_prev.p = it->u_iter.p = it->u_next.p = props->p_init;
		it->u_prev.sw = it->u_iter.sw = it->u_next.sw = props->sw_init;
		it->u_prev.xw = it->u_iter.xw = it->u_next.xw = props->xw_init;
		it->u_prev.xa = it->u_iter.xa = it->u_next.xa = props->xa_init;
		it->u_prev.xs = it->u_iter.xs = it->u_next.xs = 0.0;
		it->u_prev.t = it->u_iter.t = it->u_next.t = props->t_init;
		it->props = props;
	}

	var = new TapeVariable[stencil];
	h = new adouble[var_size];
}
double Acid2dNIT::getRate(int cur)
{
	const Cell& cell = cells[cur];
	const Cell& beta = cells[cur + cellsNum_z + 2];
	const Variable& next = cell.u_next;
	const Variable& nebr = beta.u_next;

	return props_w.getDensity(next.p, next.xa, next.xw, next.xs).value() / props_w.getDensity(Component::p_std, next.xa, next.xw, next.xs).value() *
		getTrans(cell, next.m, beta, nebr.m).value() / props_w.getViscosity(next.p, next.xa, next.xw, next.xs).value() * 
		(nebr.p - next.p);
};

void Acid2dNIT::solve_eqMiddle(const Cell& cell)
{
	const Skeleton_Props& props = *cell.props;

	trace_on(mid);
	adouble tmp;
	const Variable& prev = cell.u_prev;

	for (int i = 0; i < stencil; i++)
	{
		var[i].m <<= x[i * Variable::size];
		var[i].p <<= x[i * Variable::size + 1];
		var[i].sw <<= x[i * Variable::size + 2];
		var[i].xw <<= x[i * Variable::size + 3];
		var[i].xa <<= x[i * Variable::size + 4];
		var[i].xs <<= x[i * Variable::size + 5];
		var[i].t <<= x[i * Variable::size + 6];
	}

	TapeVariable& next = var[0];
	adouble rate = getReactionRate(next, props);
	//double rate0 = getReactionRate(next, props).value();

	h[0] = (1.0 - next.m) * props.getDensity(next.p) - (1.0 - prev.m) * props.getDensity(prev.p) - 
			ht * reac.indices[REACTS::CALCITE] * reac.comps[REACTS::CALCITE].mol_weight * rate;
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
	h[6] = getCn(cell) * (next.t - prev.t) - getAd(cell) * (next.p - prev.p);

	int neighbor[4];
	getNeighborIdx(cell.num, neighbor);
	for (int i = 0; i < 4; i++)
	{
		const Cell& beta = cells[neighbor[i]];
		const int upwd_idx = (getUpwindIdx(cell.num, neighbor[i]) == cell.num) ? 0 : i + 1;
		const TapeVariable& nebr = var[i + 1];
		TapeVariable& upwd = var[upwd_idx];
		
		adouble dens_w = getAverage(props_w.getDensity(next.p, next.xa, next.xw, next.xs), cell,
			props_w.getDensity(nebr.p, nebr.xa, nebr.xw, nebr.xs), beta);
		adouble dens_o = getAverage(props_o.getDensity(next.p), cell,
			props_o.getDensity(nebr.p), beta);
		adouble buf_w = ht / cell.V * getTrans(cell, next.m, beta, nebr.m) * ((next.p - nebr.p) - dens_w * grav * (cell.z - beta.z)) *
			dens_w * props_w.getKr(upwd.sw, cells[upwd_idx].props) / props_w.getViscosity(upwd.p, upwd.xa, upwd.xw, upwd.xs);
		adouble buf_o = ht / cell.V * getTrans(cell, next.m, beta, nebr.m) * ((next.p - nebr.p) - dens_o * grav * (cell.z - beta.z)) *
			dens_o * props_o.getKr(upwd.sw, cells[upwd_idx].props) / props_o.getViscosity(upwd.p);
		const auto mult = getDivCoeff(const_cast<Cell&>(cell), const_cast<Cell&>(beta));

		h[1] += buf_w;
		h[2] += buf_o;
		h[3] += buf_w * upwd.xw;
		h[4] += buf_w * upwd.xa;
		h[5] += buf_w * upwd.xs;
		h[6] += ht * (mult.ther * (next.t - nebr.t) + mult.pres * (next.p - nebr.p)) / fabs(getDistance(beta, cell));
	}

	//for (int i = 0; i < var_size; i++)
	//	h[i] /= P_dim;

	for (int i = 0; i < var_size; i++)
		h[i] >>= y[i];

	trace_off();
}
void Acid2dNIT::solve_eqLeft(const Cell& cell)
{
	const Cell& beta1 = cells[cell.num + cellsNum_z + 2];
	const Cell& beta2 = cells[cell.num + 2 * cellsNum_z + 4];

	trace_on(left);
	for (int i = 0; i < Lstencil; i++)
	{
		var[i].m <<= x[i * Variable::size];
		var[i].p <<= x[i * Variable::size + 1];
		var[i].sw <<= x[i * Variable::size + 2];
		var[i].xw <<= x[i * Variable::size + 3];
		var[i].xa <<= x[i * Variable::size + 4];
		var[i].xs <<= x[i * Variable::size + 5];
		var[i].t <<= x[i * Variable::size + 6];
	}

	const adouble leftIsRate = leftBoundIsRate;
	TapeVariable& next = var[0];
	TapeVariable& nebr1 = var[1];
	TapeVariable& nebr2 = var[2];
	const Variable& prev = cell.u_prev;
	const Skeleton_Props& props = *cell.props;

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

	h[0] = (1.0 - next.m) * props.getDensity(next.p) - (1.0 - prev.m) * props.getDensity(prev.p) -
		ht * reac.indices[REACTS::CALCITE] * reac.comps[REACTS::CALCITE].mol_weight * getReactionRate(next, props);
	if (leftBoundIsRate)
	{
		double tmp = props_w.getDensity(next.p, next.xa, next.xw, next.xs).value() * getTrans(cell, next.m, beta1, nebr1.m).value() /
			props_w.getViscosity(next.p, next.xa, next.xw, next.xs).value() * (nebr1.p - next.p).value() +
			props_w.getDensity(Component::p_std, next.xa, next.xw, next.xs).value() * rate;
		h[1] = /*props_w.getDensity(next.p, next.xa, next.xw, next.xs) * getTrans(cell, next.m, beta1, nebr1.m) /
			props_w.getViscosity(next.p, next.xa, next.xw, next.xs) */ (nebr1.p - next.p) /*+
			props_w.getDensity(Component::p_std, next.xa, next.xw, next.xs) * rate*/;
	}
	else
		h[1] = next.p - Pwf;
	condassign(h[1], leftPresNonPerforated, next.p - nebr1.p);
	condassign(h[2], isPerforated, next.sw - (1.0 - props.s_oc), next.sw - nebr1.sw);
	condassign(h[3], isPerforated, next.xw - (1.0 - xa), next.xw - nebr1.xw);
	condassign(h[4], isPerforated, next.xa - xa, next.xa - nebr1.xa);
	condassign(h[5], isPerforated, next.xs, next.xs - nebr1.xs);
	condassign(h[6], isPerforated, next.t - temp, next.t - nebr1.t);
	//for (int i = 0; i < var_size; i++)
	//	h[i] /= P_dim;

	for (int i = 0; i < var_size; i++)
		h[i] >>= y[i];

	trace_off();
}
void Acid2dNIT::solve_eqRight(const Cell& cell)
{
	trace_on(right);
	adouble rightIsPres = rightBoundIsPres;
	for (int i = 0; i < Rstencil; i++)
	{
		var[i].m <<= x[i * Variable::size];
		var[i].p <<= x[i * Variable::size + 1];
		var[i].sw <<= x[i * Variable::size + 2];
		var[i].xw <<= x[i * Variable::size + 3];
		var[i].xa <<= x[i * Variable::size + 4];
		var[i].xs <<= x[i * Variable::size + 5];
		var[i].t <<= x[i * Variable::size + 6];
	}

	TapeVariable& next = var[0];
	TapeVariable& nebr = var[1];
	const Variable& prev = cell.u_prev;
	const Skeleton_Props& props = *cell.props;

	h[0] = next.m - nebr.m;
	condassign(h[1], rightIsPres, next.p - (adouble)(props.p_out), next.p - (adouble)(nebr.p));
	h[2] = next.sw - nebr.sw;
	h[3] = next.xw - nebr.xw;
	h[4] = next.xa - nebr.xa;
	h[5] = next.xs - nebr.xs;
	h[6] = next.t - props.t_init;
	
	//for (int i = 0; i < var_size; i++)
	//	h[i] /= P_dim;

	for (int i = 0; i < var_size; i++)
		h[i] >>= y[i];

	trace_off();
}
void Acid2dNIT::solve_eqVertical(const Cell& cell)
{
	trace_on(vertical);
	for (int i = 0; i < Vstencil; i++)
	{
		var[i].m <<= x[i * Variable::size];
		var[i].p <<= x[i * Variable::size + 1];
		var[i].sw <<= x[i * Variable::size + 2];
		var[i].xw <<= x[i * Variable::size + 3];
		var[i].xa <<= x[i * Variable::size + 4];
		var[i].xs <<= x[i * Variable::size + 5];
		var[i].t <<= x[i * Variable::size + 6];
	}

	const TapeVariable& next = var[0];
	const TapeVariable& nebr = var[1];

	h[0] = next.m - nebr.m;
	h[1] = next.p - nebr.p;
	h[2] = next.sw - nebr.sw;
	h[3] = next.xw - nebr.xw;
	h[4] = next.xa - nebr.xa;
	h[5] = next.xs - nebr.xs;
	h[6] = next.t - nebr.t;

	//for (int i = 0; i < var_size; i++)
	//	h[i] /= P_dim;

	for (int i = 0; i < var_size; i++)
		h[i] >>= y[i];

	trace_off();
}
void Acid2dNIT::setVariables(const Cell& cell)
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