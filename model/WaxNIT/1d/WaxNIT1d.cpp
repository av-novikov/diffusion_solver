#include "model/WaxNIT/1d/WaxNIT1d.hpp"

#include "adolc/drivers/drivers.h"
#include "adolc/adolc.h"

using namespace wax_nit1d;

WaxNIT1d::WaxNIT1d()
{
	x = new double[stencil * Variable::size];
	y = new double[var_size];

	jac = new double*[var_size];
	for (int i = 0; i < var_size; i++)
		jac[i] = new double[stencil * Variable::size];
}
WaxNIT1d::~WaxNIT1d()
{
	delete[] y;

	for (int i = 0; i < var_size; i++)
		delete[] jac[i];
	delete[] jac;

	delete[] x, h;
	delete[] var;
}
void WaxNIT1d::setProps(Properties& props)
{
	props_sk = props.props_sk;
	setBasicProps(props);

	props_oil = props.props_oil;
	props_oil.visc = cPToPaSec(props_oil.visc);
	props_wax = props.props_wax;

	makeBasicDimLess();
	makeDimLess();

	// Data sets
	//props_wat.kr = setDataset(props.kr_wat, 1.0, 1.0);
	//props_oil.kr = setDataset(props.kr_oil, 1.0, 1.0);
	//props_gas.kr = setDataset(props.kr_gas, 1.0, 1.0);
	//props_wat.b = setDataset(props.B_wat, P_dim / BAR_TO_PA, 1.0);
	props_oil.b = setDataset(props.B_oil, P_dim / BAR_TO_PA, 1.0);
	props_oil.Rs = setDataset(props.Rs, P_dim / BAR_TO_PA, 1.0);
	props_oil.lp = setDataset(props.lp, T_dim, 1.0);
}
void WaxNIT1d::makeDimLess()
{
	props_sk.p_sat /= P_dim;

	props_oil.visc /= (P_dim * t_dim);
	props_oil.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_oil.dens_gas_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_oil.dens_wax_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_oil.gamma /= (1.0 / R_dim);
	props_oil.beta /= (1.0 / P_dim);
	props_oil.p_ref /= P_dim;

	props_wax.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
}
void WaxNIT1d::setInitialState()
{
	vector<Cell>::iterator it;
	for (it = cells.begin(); it != cells.end(); ++it)
	{
		it->u_prev.m = it->u_iter.m = it->u_next.m = props_sk.m_init;
		it->u_prev.p = it->u_iter.p = it->u_next.p = props_sk.p_init;
		it->u_prev.s_p = it->u_iter.s_p = it->u_next.s_p = 1.0 - props_sk.so_init;
		
		it->props = &props_sk;
	}

	var = new TapeVariable[stencil];
	h = new adouble[var_size];
}
double WaxNIT1d::getRate(int cur_idx, int beta_idx) const
{
	const Cell& cell = cells[cur_idx];
	const Cell& beta = cells[beta_idx];
	const Variable& next = cell.u_next;
	const Variable& nebr = beta.u_next;
	const Variable& upwd = cells[getUpwindIdx(cur_idx, beta_idx)].u_next;

	return getTrans(cell, next.m, beta, nebr.m).value() * props_oil.getKr(1.0 - upwd.s_p, cell.props).value() /	
		props_oil.getB(next.p).value() / props_oil.getViscosity(next.p).value() * (nebr.p - next.p);
}

void WaxNIT1d::solve_eqMiddle(const Cell& cell)
{
	const Skeleton_Props& props = *cell.props;

	trace_on(mid);
	adouble tmp, tmp_o;
	const Variable& prev = cell.u_prev;

	for (int i = 0; i < stencil; i++)
	{
		var[i].m <<= x[i * Variable::size];
		var[i].p <<= x[i * Variable::size + 1];
		var[i].s_p <<= x[i * Variable::size + 2];
	}
	const TapeVariable& next = var[0];

	adouble dadt = props_oil.gamma * next.s_p * props_wax.getRho() * getOilVelocityAbs(cell);
	h[0] = props.dens_stc * ((1.0 - next.m) - (1.0 - prev.m)) - dadt;
	h[1] = next.m * (next.s_p * props_wax.getRho()) -
		prev.m * (prev.s_p * props_wax.getRho()) + dadt;
	h[2] = next.m * (1.0 - next.s_p) * props_oil.getRho(next.p) - prev.m * (1.0 - next.s_p) * props_oil.getRho(prev.p);

	int neighbor[2];
	getNeighborIdx(cell.num, neighbor);
	for (int i = 0; i < 2; i++)
	{
		const Cell& beta = cells[neighbor[i]];
		const Cell& upwd_cell = cells[getUpwindIdx(cell.num, neighbor[i])];
		const int upwd_idx = (upwd_cell.num == cell.num) ? 0 : i + 1;
		const TapeVariable& nebr = var[i + 1];
		TapeVariable& upwd = var[upwd_idx];

		tmp_o = ht / cell.V * getTrans(cell, next.m, beta, nebr.m) * (next.p - nebr.p) *
			props_oil.getKr(1.0 - upwd.s_p, cells[upwd_idx].props) / props_oil.getViscosity(upwd.p);
		h[1] += (props_wax.getRho() *  upwd.s_p / (1.0 - upwd.s_p)) * tmp_o;
		h[2] += tmp_o * getAverage(props_oil.getRho(next.p), cell, props_oil.getRho(nebr.p), beta);
	}

	for (int i = 0; i < var_size; i++)
		h[i] >>= y[i];

	trace_off();
}
void WaxNIT1d::solve_eqLeft(const Cell& cell)
{
	const Cell& beta1 = cells[cell.num + 1];
	const Cell& beta2 = cells[cell.num + 2];
	
	trace_on(left);
	for (int i = 0; i < Lstencil; i++)
	{
		var[i].m <<= x[i * Variable::size];
		var[i].p <<= x[i * Variable::size + 1];
		var[i].s_p <<= x[i * Variable::size + 2];
	}

	adouble leftIsRate = leftBoundIsRate;
	const TapeVariable& next = var[0];
	const TapeVariable& nebr1 = var[1];
	const TapeVariable& nebr2 = var[2];
	const Variable& prev = cell.u_prev;
	const Skeleton_Props& props = *cell.props;

	double rate = Qcell[0];

	h[0] = (next.m - nebr1.m) / (adouble)(cell.x - beta1.x) - (nebr1.m - nebr2.m) / (adouble)(beta1.x - beta2.x);
	condassign(h[1], leftIsRate, props_oil.getRho(next.p) * getTrans(cell, next.m, beta1, nebr1.m) * 
		props_oil.getKr(1.0 - next.s_p, cell.props) / props_oil.getViscosity(next.p) * (nebr1.p - next.p) + props_oil.dens_stc * rate,
		next.p - Pwf);
	h[2] = next.s_p - 0.03;
	
	for (int i = 0; i < var_size; i++)
		h[i] /= P_dim;

	for (int i = 0; i < var_size; i++)
		h[i] >>= y[i];

	trace_off();
}
void WaxNIT1d::solve_eqRight(const Cell& cell)
{
	trace_on(right);
	adouble rightIsPres = rightBoundIsPres;
	for (int i = 0; i < Rstencil; i++)
	{
		var[i].m <<= x[i * Variable::size];
		var[i].p <<= x[i * Variable::size + 1];
		var[i].s_p <<= x[i * Variable::size + 2];
	}

	const Cell& beta1 = cells[cell.num - 1];
	const Cell& beta2 = cells[cell.num - 2];
	const TapeVariable& next = var[0];
	const TapeVariable& nebr1 = var[1];
	const TapeVariable& nebr2 = var[2];

	h[0] = next.m - cell.props->m_init;
	condassign(h[1], rightIsPres, next.p - (adouble)(cell.props->p_out), next.p - (adouble)(nebr1.p));
	h[2] = (next.s_p - nebr1.s_p) / (adouble)(cell.x - beta1.x) - (nebr1.s_p - nebr2.s_p) / (adouble)(beta1.x - beta2.x);
	
	for (int i = 0; i < var_size; i++)
		h[i] /= P_dim;

	for (int i = 0; i < var_size; i++)
		h[i] >>= y[i];

	trace_off();
}
void WaxNIT1d::setVariables(const Cell& cell)
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