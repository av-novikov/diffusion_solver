#include "model/WaxNIT/WaxNIT.hpp"

#include "adolc/drivers/drivers.h"
#include "adolc/adolc.h"

using namespace wax_nit;

WaxNIT::WaxNIT()
{
	x = new double[stencil * Variable::size];
	y = new double[var_size];

	jac = new double*[var_size];
	for (int i = 0; i < var_size; i++)
		jac[i] = new double[stencil * Variable::size];
}
WaxNIT::~WaxNIT()
{
	delete[] y;

	for (int i = 0; i < var_size; i++)
		delete[] jac[i];
	delete[] jac;

	delete[] x, h;
	delete[] var;
}
void WaxNIT::setProps(Properties& props)
{
	props_sk = props.props_sk;
	setBasicProps(props);

	props_wat = props.props_wat;
	props_wat.visc = cPToPaSec(props_wat.visc);
	props_oil = props.props_oil;
	props_oil.visc = cPToPaSec(props_oil.visc);
	props_gas = props.props_gas;
	props_gas.visc = cPToPaSec(props_gas.visc);
	props_wax = props.props_wax;
	L = props.L;

	makeBasicDimLess();
	makeDimLess();

	// Data sets
	//props_wat.kr = setDataset(props.kr_wat, 1.0, 1.0);
	//props_oil.kr = setDataset(props.kr_oil, 1.0, 1.0);
	//props_gas.kr = setDataset(props.kr_gas, 1.0, 1.0);
	//props_wat.b = setDataset(props.B_wat, P_dim / BAR_TO_PA, 1.0);
	props_oil.b = setDataset(props.B_oil, P_dim / BAR_TO_PA, 1.0);
	props_gas.b = setDataset(props.B_gas, P_dim / BAR_TO_PA, 1.0);
	props_oil.Rs = setDataset(props.Rs, P_dim / BAR_TO_PA, 1.0);
	props_oil.lp = setDataset(props.lp, T_dim, 1.0);
}
void WaxNIT::makeDimLess()
{
	if (props_sk[0].t_init != 0.0)
		T_dim = fabs(props_sk[0].t_init);
	else
		T_dim = 1.0;

	for (auto& sk : props_sk)
	{
		sk.p_sat /= P_dim;
		sk.t_init /= T_dim;
		sk.t_sat /= T_dim;

		sk.c = sk.c / R_dim / R_dim * T_dim * t_dim * t_dim;
		sk.lambda_r = sk.lambda_r * T_dim * t_dim / P_dim / R_dim / R_dim;
		sk.lambda_z = sk.lambda_z * T_dim * t_dim / P_dim / R_dim / R_dim;
	}

	props_wat.visc /= (P_dim * t_dim);
	props_wat.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_wat.beta /= (1.0 / P_dim);
	props_wat.p_ref /= P_dim;
	props_wat.c = props_wat.c / R_dim / R_dim * T_dim * t_dim * t_dim;
	props_wat.lambda = props_wat.lambda * T_dim * t_dim / P_dim / R_dim / R_dim;
	props_wat.jt = props_wat.jt * P_dim / T_dim;
	props_wat.ad = props_wat.ad * P_dim / T_dim;

	props_oil.visc /= (P_dim * t_dim);
	props_oil.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_oil.dens_gas_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_oil.dens_wax_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_oil.gamma /= (1.0 / R_dim);
	props_oil.beta /= (1.0 / P_dim);
	props_oil.p_ref /= P_dim;
	props_oil.c = props_oil.c / R_dim / R_dim * T_dim * t_dim * t_dim;
	props_oil.lambda = props_oil.lambda * T_dim * t_dim / P_dim / R_dim / R_dim;
	props_oil.jt = props_oil.jt * P_dim / T_dim;
	props_oil.ad = props_oil.ad * P_dim / T_dim;

	props_gas.visc /= (P_dim * t_dim);
	props_gas.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_gas.c = props_gas.c / R_dim / R_dim * T_dim * t_dim * t_dim;
	props_gas.lambda = props_gas.lambda * T_dim * t_dim / P_dim / R_dim / R_dim;
	props_gas.jt = props_gas.jt * P_dim / T_dim;
	props_gas.ad = props_gas.ad * P_dim / T_dim;

	props_wax.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);

	L = L / R_dim / R_dim * t_dim * t_dim;	
}
void WaxNIT::setInitialState()
{
	vector<Cell>::iterator it;
	for (it = cells.begin(); it != cells.end(); ++it)
	{
		Skeleton_Props* props = &props_sk[getSkeletonIdx(*it)];
		it->u_prev.m = it->u_iter.m = it->u_next.m = props->m_init;
		it->u_prev.t = it->u_iter.t = it->u_next.t = props->t_init;
		it->u_prev.p = it->u_iter.p = it->u_next.p = props->p_init;
		it->u_prev.s_w = it->u_iter.s_w = it->u_next.s_w = props->sw_init;
		it->u_prev.s_o = it->u_iter.s_o = it->u_next.s_o = props->so_init;
		it->u_prev.s_g = it->u_iter.s_g = it->u_next.s_g = props->sg_init;

		it->u_prev.p_bub = it->u_iter.p_bub = it->u_next.p_bub = props->p_sat;
		it->u_prev.t_bub = it->u_iter.t_bub = it->u_next.t_bub = props->t_sat;

		if (props->p_init > props->p_sat)
			it->u_prev.satur_gas = it->u_iter.satur_gas = it->u_next.satur_gas = false;
		else
			it->u_prev.satur_gas = it->u_iter.satur_gas = it->u_next.satur_gas = true;
		if (props->t_init > props->t_sat)
			it->u_prev.satur_wax = it->u_iter.satur_wax = it->u_next.satur_wax = false;
		else
			it->u_prev.satur_wax = it->u_iter.satur_wax = it->u_next.satur_wax = true;

		it->props = props;
	}

	var = new TapeVariable[stencil];
	h = new adouble[var_size];
}
double WaxNIT::getRate(int cur) const
{
	const Cell& cell = cells[cur];
	const Cell& beta = cells[cur + cellsNum_z + 2];
	const Variable& next = cell.u_next;
	const Variable& nebr = beta.u_next;
	const Variable& upwd = cells[getUpwindIdx(cur, cur + cellsNum_z + 2)].u_next;

	return getTrans(cell, next.m, beta, nebr.m).value() * props_oil.getKr(upwd.s_w, upwd.s_o, cell.props).value() * 
		(1.0 - props_oil.getfp(next.p, next.p_bub, next.satur_gas, next.t, next.t_bub, next.satur_wax).value()) /	
		props_oil.getB(next.p, next.p_bub, next.satur_gas).value() / props_oil.getViscosity(next.p).value() * (nebr.p - next.p);
}

adouble WaxNIT::phaseTrans(const Cell& cell)
{
	const Skeleton_Props& props = *cell.props;
	const TapeVariable& next = var[0];
	const Variable& prev = cell.u_prev;

	adouble fp = props_oil.getfp(next.p, next.p_bub, next.satur_gas, next.t, next.t_bub, next.satur_wax);
	adouble H = (next.m * next.s_o * (1.0 - fp) * props_oil.getRhoTilde(next.p, next.p_bub, next.satur_gas) -
					prev.m * prev.s_o * (1.0 - props_oil.getfp(prev.p, prev.p_bub, prev.satur_gas, prev.t, prev.t_bub, prev.satur_wax)) * 
											props_oil.getRhoTilde(prev.p, prev.p_bub, prev.satur_gas)) / ht;
	int neighbor[4];
	getNeighborIdx(cell.num, neighbor);
	for (int i = 0; i < 4; i++)
	{
		const Cell& beta = cells[neighbor[i]];
		const Cell& upwd_cell = cells[getUpwindIdx(cell.num, neighbor[i])];
		const int upwd_idx = (upwd_cell.num == cell.num) ? 0 : i + 1;
		const TapeVariable& nebr = var[i + 1];
		TapeVariable& upwd = var[upwd_idx];

		adouble fp_nebr = props_oil.getfp(nebr.p, nebr.p_bub, nebr.satur_gas, nebr.t, nebr.t_bub, nebr.satur_wax);
		H += 1.0 / cell.V * getTrans(cell, next.m, beta, nebr.m) * (next.p - nebr.p) *
			props_gas.getKr(upwd.s_w, upwd.s_o, upwd_cell.props) * 
			getAverage((1.0 - fp) * props_oil.getRhoTilde(next.p, next.p_bub, next.satur_gas) / props_oil.getViscosity(next.p), cell,
						(1.0 - fp_nebr) * props_oil.getRhoTilde(nebr.p, nebr.p_bub, nebr.satur_gas) / props_oil.getViscosity(nebr.p), beta);
	}
	return -H;
}
void WaxNIT::solve_eqMiddle(const Cell& cell)
{
	const Skeleton_Props& props = *cell.props;

	trace_on(mid);
	adouble tmp, tmp_o;
	const Variable& prev = cell.u_prev;

	for (int i = 0; i < stencil; i++)
	{
		var[i].m <<= x[i * Variable::size];
		var[i].t <<= x[i * Variable::size + 1];
		var[i].p <<= x[i * Variable::size + 2];
		var[i].s_w <<= x[i * Variable::size + 3];
		var[i].s_o <<= x[i * Variable::size + 4];
		var[i].s_g <<= x[i * Variable::size + 5];
		var[i].p_bub <<= x[i * Variable::size + 6];
		var[i].t_bub <<= x[i * Variable::size + 7];
	}
	const TapeVariable& next = var[0];
	adouble satur_gas = cell.u_next.satur_gas;
	adouble satur_wax = cell.u_next.satur_wax;

	adouble dadt = props_oil.gamma * (1.0 - next.s_w - next.s_o - next.s_g) * props_wax.getRho() * getOilVelocityAbs(cell);
	adouble fp = props_oil.getfp(next.p, next.p_bub, satur_gas, next.t, next.t_bub, satur_wax);
	adouble fp_prev = props_oil.getfp(prev.p, prev.p_bub, prev.satur_gas, prev.t, prev.t_bub, prev.satur_wax);
	h[0] = props.dens_stc * ((1.0 - next.m) - (1.0 - prev.m)) - dadt;
	h[1] = getCn(cell) * (next.t - prev.t) - getAd(cell) * (next.p - prev.p) - ht * L * phaseTrans(cell);
	condassign(h[2], satur_wax,
		next.m * (next.s_o * fp * props_oil.dens_wax_stc + (1.0 - next.s_w - next.s_o - next.s_g) * props_wax.getRho()) -
		prev.m * (prev.s_o * fp_prev * props_oil.dens_wax_stc + (1.0 - prev.s_w - prev.s_o - prev.s_g) * props_wax.getRho()) + dadt,
		next.m * next.s_o * fp * props_oil.dens_wax_stc - prev.m * prev.s_o * fp_prev * props_oil.dens_wax_stc + dadt);
	h[3] = next.m * next.s_w * props_wat.getRho(next.p) - prev.m * prev.s_w * props_wat.getRho(prev.p);
	h[4] = next.m * next.s_o * (1.0 - fp) * props_oil.dens_stc / props_oil.getB(next.p, next.p_bub, satur_gas) -
		prev.m * prev.s_o * (1.0 - fp_prev) * props_oil.dens_stc / props_oil.getB(prev.p, prev.p_bub, prev.satur_gas);
	condassign(h[5], satur_gas, next.m * (next.s_o * (1.0 - fp) * props_oil.getRs(next.p, next.p_bub, satur_gas) * props_oil.dens_gas_stc / props_oil.getB(next.p, next.p_bub, satur_gas) +
									next.s_g * props_gas.getRho(next.p)) -
								prev.m * (prev.s_o * (1.0 - fp_prev) * props_oil.getRs(prev.p, prev.p_bub, prev.satur_gas) * props_oil.dens_gas_stc / props_oil.getB(prev.p, prev.p_bub, prev.satur_gas) +
									prev.s_g * props_gas.getRho(prev.p)),
								next.m * next.s_o * (1.0 - fp) * props_oil.getRs(next.p, next.p_bub, satur_gas) * props_oil.dens_gas_stc / props_oil.getB(next.p, next.p_bub, satur_gas) - 
								prev.m * prev.s_o * (1.0 - fp_prev) * props_oil.getRs(prev.p, prev.p_bub, prev.satur_gas) * props_oil.dens_gas_stc / props_oil.getB(prev.p, prev.p_bub, prev.satur_gas));
	int neighbor[4];
	getNeighborIdx(cell.num, neighbor);
	for (int i = 0; i < 4; i++)
	{
		const Cell& beta = cells[neighbor[i]];
		const Cell& upwd_cell = cells[getUpwindIdx(cell.num, neighbor[i])];
		const int upwd_idx = (upwd_cell.num == cell.num) ? 0 : i + 1;
		const TapeVariable& nebr = var[i + 1];
		TapeVariable& upwd = var[upwd_idx];

		adouble fp_nebr = props_oil.getfp(nebr.p, nebr.p_bub, beta.u_next.satur_gas, nebr.t, nebr.t_bub, beta.u_next.satur_wax);
		adouble fp_avg = getAverage(fp, cell, fp_nebr, beta);
		const auto mult = getDivCoeff(const_cast<Cell&>(cell), const_cast<Cell&>(beta));
		tmp_o = ht / cell.V * getTrans(cell, next.m, beta, nebr.m) * (next.p - nebr.p) *
			props_oil.getKr(upwd.s_w, upwd.s_o, cells[upwd_idx].props) / props_oil.getViscosity(upwd.p);
		h[1] += ht * (mult.ther * (next.t - nebr.t) + mult.pres * (next.p - nebr.p)) / fabs(getDistance(beta, cell));
		condassign(tmp, satur_wax, (props_oil.dens_wax_stc * fp_avg + 
									props_wax.getRho() * (1.0 - upwd.s_w - upwd.s_o - upwd.s_g) / upwd.s_o) * tmp_o, 
									props_oil.dens_wax_stc * fp_avg * tmp_o);
		h[2] += tmp;
		h[3] += ht / cell.V * getTrans(cell, next.m, beta, nebr.m) * (next.p - nebr.p) *
					props_wat.getKr(upwd.s_w, upwd.s_o, cells[upwd_idx].props) / props_wat.getViscosity(upwd.p) *
			getAverage(props_wat.getRho(next.p), cell, props_wat.getRho(nebr.p), beta);
		h[4] += (1.0 - fp_avg) * tmp_o * props_oil.dens_stc /
			getAverage(props_oil.getB(next.p, next.p_bub, satur_gas), cell, props_oil.getB(nebr.p, nebr.p_bub, beta.u_next.satur_gas), beta);
		condassign(tmp, satur_gas, (1.0 - fp_avg) * tmp_o * props_oil.dens_gas_stc * 
			getAverage(props_oil.getRs(next.p, next.p_bub, satur_gas) / props_oil.getB(next.p, next.p_bub, satur_gas), cell, 
						props_oil.getRs(nebr.p, nebr.p_bub, beta.u_next.satur_gas) / props_oil.getB(nebr.p, nebr.p_bub, beta.u_next.satur_gas), beta) + 
				ht / cell.V * getTrans(cell, next.m, beta, nebr.m) * (next.p - nebr.p) * getAverage(props_gas.getRho(next.p), cell, props_gas.getRho(nebr.p), beta) * 
			props_gas.getKr(upwd.s_w, upwd.s_o, cells[upwd_idx].props) / props_gas.getViscosity(upwd.p),
				(1.0 - fp_avg) * tmp_o * props_oil.dens_gas_stc *
			getAverage(props_oil.getRs(next.p, next.p_bub, satur_gas) / props_oil.getB(next.p, next.p_bub, next.satur_gas), cell,
				props_oil.getRs(nebr.p, nebr.p_bub, beta.u_next.satur_gas) / props_oil.getB(nebr.p, nebr.p_bub, beta.u_next.satur_gas), beta));
		h[5] += tmp;
	}

	for (int i = 0; i < var_size; i++)
		h[i] >>= y[i];

	trace_off();
}
void WaxNIT::solve_eqLeft(const Cell& cell)
{
	const Cell& beta1 = cells[cell.num + cellsNum_z + 2];
	const Cell& beta2 = cells[cell.num + 2 * cellsNum_z + 4];
	
	trace_on(left);
	for (int i = 0; i < Lstencil; i++)
	{
		var[i].m <<= x[i * Variable::size];
		var[i].t <<= x[i * Variable::size + 1];
		var[i].p <<= x[i * Variable::size + 2];
		var[i].s_w <<= x[i * Variable::size + 3];
		var[i].s_o <<= x[i * Variable::size + 4];
		var[i].s_g <<= x[i * Variable::size + 5];
		var[i].p_bub <<= x[i * Variable::size + 6];
		var[i].t_bub <<= x[i * Variable::size + 7];
	}

	adouble leftIsRate = leftBoundIsRate;
	const TapeVariable& next = var[0];
	const TapeVariable& nebr1 = var[1];
	const TapeVariable& nebr2 = var[2];
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

	/*condassign(h[0], isPerforated, props.dens_stc * ((1.0 - next.m) - (1.0 - prev.m)) -
		props_oil.gamma * props_oil.getRho(next.p, next.p_bub, next.SATUR) * props_oil.getlp(next.t) * getOilVelocityAbs(cell),
		next.m - nebr1.m);*/
	h[0] = (next.m - nebr1.m) / (adouble)(cell.r - beta1.r) - (nebr1.m - nebr2.m) / (adouble)(beta1.r - beta2.r);
	h[1] = next.t - nebr1.t;
	condassign(h[2], leftIsRate, props_oil.getRho(next.p, next.p_bub, next.satur_gas, next.t, next.t_bub, next.satur_wax) * getTrans(cell, next.m, beta1, nebr1.m) * 
		props_oil.getKr(next.s_w, next.s_o, cell.props) / props_oil.getViscosity(next.p) * (nebr1.p - next.p) + props_oil.dens_stc * rate,
		next.p - Pwf);
	condassign(h[2], leftPresNonPerforated, next.p - nebr1.p);
	h[3] = (next.s_w - nebr1.s_w) / (adouble)(cell.r - beta1.r) - (nebr1.s_w - nebr2.s_w) / (adouble)(beta1.r - beta2.r);
	adouble satur_gas = cell.u_next.satur_gas;
	adouble satur_wax = cell.u_next.satur_wax;
	condassign(h[4], satur_gas,
		(next.s_o - nebr1.s_o) / (adouble)(cell.r - beta1.r) - (nebr1.s_o - nebr2.s_o) / (adouble)(beta1.r - beta2.r),
		(next.p_bub - nebr1.p_bub) / (adouble)(cell.r - beta1.r) - (nebr1.p_bub - nebr2.p_bub) / (adouble)(beta1.r - beta2.r));
	condassign(h[5], satur_wax,
		(next.s_g - nebr1.s_g) / (adouble)(cell.r - beta1.r) - (nebr1.s_g - nebr2.s_g) / (adouble)(beta1.r - beta2.r),
		(next.t_bub - nebr1.t_bub) / (adouble)(cell.r - beta1.r) - (nebr1.t_bub - nebr2.t_bub) / (adouble)(beta1.r - beta2.r));

	for (int i = 0; i < var_size; i++)
		h[i] /= P_dim;

	for (int i = 0; i < var_size; i++)
		h[i] >>= y[i];

	trace_off();
}
void WaxNIT::solve_eqRight(const Cell& cell)
{
	trace_on(right);
	adouble rightIsPres = rightBoundIsPres;
	for (int i = 0; i < Rstencil; i++)
	{
		var[i].m <<= x[i * Variable::size];
		var[i].t <<= x[i * Variable::size + 1];
		var[i].p <<= x[i * Variable::size + 2];
		var[i].s_w <<= x[i * Variable::size + 3];
		var[i].s_o <<= x[i * Variable::size + 4];
		var[i].s_g <<= x[i * Variable::size + 5];
		var[i].p_bub <<= x[i * Variable::size + 6];
		var[i].t_bub <<= x[i * Variable::size + 7];
	}

	const TapeVariable& next = var[0];
	const TapeVariable& nebr = var[1];

	h[0] = next.m - cell.props->m_init;
	h[1] = next.t - cell.props->t_init;
	condassign(h[2], rightIsPres, next.p - (adouble)(cell.props->p_out), next.p - (adouble)(nebr.p));
	h[3] = next.s_w - nebr.s_w;
	adouble satur_gas = next.satur_gas;
	adouble satur_wax = next.satur_gas;
	condassign(h[4], satur_gas, next.s_o - nebr.s_o, next.p_bub - nebr.p_bub);
	condassign(h[5], satur_wax, next.s_g - nebr.s_g, next.t_bub - nebr.t_bub);

	for (int i = 0; i < var_size; i++)
		h[i] /= P_dim;

	for (int i = 0; i < var_size; i++)
		h[i] >>= y[i];

	trace_off();
}
void WaxNIT::solve_eqVertical(const Cell& cell)
{
	trace_on(vertical);
	for (int i = 0; i < Vstencil; i++)
	{
		var[i].m <<= x[i * Variable::size];
		var[i].t <<= x[i * Variable::size + 1];
		var[i].p <<= x[i * Variable::size + 2];
		var[i].s_w <<= x[i * Variable::size + 3];
		var[i].s_o <<= x[i * Variable::size + 4];
		var[i].s_g <<= x[i * Variable::size + 5];
		var[i].p_bub <<= x[i * Variable::size + 6];
		var[i].t_bub <<= x[i * Variable::size + 7];
	}
	const TapeVariable& next = var[0];
	const TapeVariable& nebr = var[1];

	adouble satur_gas = cell.u_next.satur_gas;
	adouble satur_wax = cell.u_next.satur_wax;
	h[0] = next.m - nebr.m;
	h[1] = next.t - nebr.t;
	h[2] = next.p - nebr.p;
	h[3] = next.s_w - nebr.s_w;
	condassign(h[4], satur_gas, next.s_o - nebr.s_o, next.p_bub - nebr.p_bub);
	condassign(h[5], satur_wax, next.s_g - nebr.s_g, next.t_bub - nebr.t_bub);

	for (int i = 0; i < var_size; i++)
		h[i] /= P_dim;

	for (int i = 0; i < var_size; i++)
		h[i] >>= y[i];

	trace_off();
}
void WaxNIT::setVariables(const Cell& cell)
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