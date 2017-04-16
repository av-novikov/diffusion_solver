#include "model/BlackOil_RZ/BlackOil_RZ.hpp"

#include "adolc/drivers/drivers.h"
#include "adolc/adolc.h"

using namespace blackoil_rz;

BlackOil_RZ::BlackOil_RZ()
{
	x = new double[stencil * Variable::size];
	y = new double[var_size];

	jac = new double*[var_size];
	for (int i = 0; i < var_size; i++)
		jac[i] = new double[stencil * Variable::size];
};
BlackOil_RZ::~BlackOil_RZ()
{
	delete x;
	delete y;

	for (int i = 0; i < var_size; i++)
		delete[] jac[i];
	delete[] jac;
};
void BlackOil_RZ::setProps(Properties& props)
{
	props_sk = props.props_sk;
	setBasicProps(props);

	props_wat = props.props_wat;
	props_wat.visc = cPToPaSec(props_wat.visc);
	props_oil = props.props_oil;
	props_oil.visc = cPToPaSec(props_oil.visc);
	props_gas = props.props_gas;
	props_gas.visc = cPToPaSec(props_gas.visc);

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
}
void BlackOil_RZ::makeDimLess()
{
	for (int j = 0; j < skeletonsNum; j++)
	{
		props_sk[j].p_sat /= P_dim;
	}

	props_wat.visc /= (P_dim * t_dim);
	props_wat.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_wat.beta /= (1.0 / P_dim);
	props_wat.p_ref /= P_dim;
	props_oil.visc /= (P_dim * t_dim);
	props_oil.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_oil.beta /= (1.0 / P_dim);
	props_oil.p_ref /= P_dim;
	props_gas.visc /= (P_dim * t_dim);
	props_gas.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
}
void BlackOil_RZ::setInitialState()
{
	vector<Cell>::iterator it;
	for (it = cells.begin(); it != cells.end(); ++it)
	{
		Skeleton_Props* props = &props_sk[getSkeletonIdx(*it)];
		it->u_prev.p = it->u_iter.p = it->u_next.p = props->p_init;
		it->u_prev.p_bub = it->u_iter.p_bub = it->u_next.p_bub = props->p_sat;
		it->u_prev.s_w = it->u_iter.s_w = it->u_next.s_w = props->sw_init;
		it->u_prev.s_o = it->u_iter.s_o = it->u_next.s_o = props->so_init;
		if (props->p_init > props->p_sat)
			it->u_prev.SATUR = it->u_iter.SATUR = it->u_next.SATUR = false;
		else
			it->u_prev.SATUR = it->u_iter.SATUR = it->u_next.SATUR = true;

		it->props = props;
	}
}
double BlackOil_RZ::getRate(int cur) const
{
	const Cell& cell = cells[cur];
	const Cell& beta = cells[cur + cellsNum_z + 2];
	const Variable& next = cell.u_next;
	const Variable& upwd = cells[getUpwindIdx(cur, cur + cellsNum_z + 2)].u_next;

	return getTrans(cell, beta) * props_oil.getKr(upwd.s_w, upwd.s_o, cell.props).value() /
		props_oil.getViscosity(next.p).value() /
		props_oil.getB(next.p, next.p_bub, next.SATUR).value() *
		(beta.u_next.p - next.p);
}

void BlackOil_RZ::solve_eqMiddle(const Cell& cell)
{
	const Skeleton_Props& props = *cell.props;

	trace_on(mid);
	adouble h[var_size], tmp;
	TapeVariable var[stencil];
	const Variable& prev = cell.u_prev;

	for (int i = 0; i < stencil; i++)
	{
		var[i].p <<= x[i * Variable::size];
		var[i].s_w <<= x[i * Variable::size + 1];
		var[i].s_o <<= x[i * Variable::size + 2];
		var[i].p_bub <<= x[i * Variable::size + 3];
	}

	const TapeVariable& next = var[0];
	adouble satur = cell.u_next.SATUR;

	h[0] = props.getPoro(next.p) * next.s_w / props_wat.getB(next.p, next.p_bub, next.SATUR) -
		props.getPoro(prev.p) * prev.s_w / props_wat.getB(prev.p, prev.p_bub, prev.SATUR);
	h[1] = props.getPoro(next.p) * next.s_o / props_oil.getB(next.p, next.p_bub, next.SATUR) -
		props.getPoro(prev.p) * prev.s_o / props_oil.getB(prev.p, prev.p_bub, prev.SATUR);
	condassign(h[2], satur,
		props.getPoro(next.p) * ((1.0 - next.s_o - next.s_w) / props_gas.getB(next.p) +
			next.s_o * props_oil.getRs(next.p, next.p_bub, next.SATUR) / props_oil.getB(next.p, next.p_bub, next.SATUR)) -
		props.getPoro(prev.p) * ((1.0 - prev.s_o - prev.s_w) / props_gas.getB(prev.p) +
			prev.s_o * props_oil.getRs(prev.p, prev.p_bub, prev.SATUR) / props_oil.getB(prev.p, prev.p_bub, prev.SATUR)),
		props.getPoro(next.p) *	next.s_o * props_oil.getRs(next.p, next.p_bub, next.SATUR) / props_oil.getB(next.p, next.p_bub, next.SATUR) - 
			props.getPoro(prev.p) *	prev.s_o * props_oil.getRs(prev.p, prev.p_bub, prev.SATUR) / props_oil.getB(prev.p, prev.p_bub, prev.SATUR));

	int neighbor[4];
	getNeighborIdx(cell.num, neighbor);
	for (int i = 0; i < 4; i++)
	{
		const Cell& beta = cells[neighbor[i]];
		const int upwd_idx = (getUpwindIdx(cell.num, neighbor[i]) == cell.num) ? 0 : i + 1;
		const TapeVariable& nebr = var[i + 1];
		TapeVariable& upwd = var[upwd_idx];

		h[0] += ht / cell.V * getTrans(cell, beta) * (next.p - nebr.p) *
			props_wat.getKr(upwd.s_w, upwd.s_o, cells[upwd_idx].props) / props_wat.getViscosity(upwd.p) /
			props_wat.getB(upwd.p, upwd.p_bub, upwd.SATUR);

		h[1] += ht / cell.V * getTrans(cell, beta) * (next.p - nebr.p) *
			props_oil.getKr(upwd.s_w, upwd.s_o, cells[upwd_idx].props) / props_oil.getViscosity(upwd.p) /
			props_oil.getB(upwd.p, upwd.p_bub, upwd.SATUR);

		condassign(tmp, satur, 
				ht / cell.V * getTrans(cell, beta) * (next.p - nebr.p) *
			(props_oil.getKr(upwd.s_w, upwd.s_o, cells[upwd_idx].props) * props_oil.getRs(upwd.p, upwd.p_bub, upwd.SATUR) /
			props_oil.getViscosity(upwd.p) / props_oil.getB(upwd.p, upwd.p_bub, upwd.SATUR) +
			props_gas.getKr(upwd.s_w, upwd.s_o, cells[upwd_idx].props) / props_gas.getViscosity(upwd.p) / props_gas.getB(upwd.p)),
				ht / cell.V * getTrans(cell, beta) * (next.p - nebr.p) *
			props_oil.getKr(upwd.s_w, upwd.s_o, cells[upwd_idx].props) * props_oil.getRs(upwd.p, upwd.p_bub, upwd.SATUR) /
			props_oil.getViscosity(upwd.p) / props_oil.getB(upwd.p, upwd.p_bub, upwd.SATUR));
		h[2] += tmp;
	}

	for (int i = 0; i < var_size; i++)
		h[i] >>= y[i];

	trace_off();
}
void BlackOil_RZ::solve_eqLeft(const Cell& cell)
{
	const Cell& beta1 = cells[cell.num + cellsNum_z + 2];
	const Cell& beta2 = cells[cell.num + 2 * cellsNum_z + 4];

	trace_on(left);
	adouble h[var_size];
	TapeVariable var[Lstencil];
	for (int i = 0; i < Lstencil; i++)
	{
		var[i].p <<= x[i * Variable::size];
		var[i].s_w <<= x[i * Variable::size + 1];
		var[i].s_o <<= x[i * Variable::size + 2];
		var[i].p_bub <<= x[i * Variable::size + 3];
	}

	adouble leftIsRate = leftBoundIsRate;
	const TapeVariable& next = var[0];
	const TapeVariable& nebr1 = var[1];
	const TapeVariable& nebr2 = var[2];
	const int upwd_idx = (getUpwindIdx(cell.num, cell.num + cellsNum_z + 2) == cell.num) ? 0 : 1;
	TapeVariable& upwd = var[upwd_idx];

	condassign(h[0], leftIsRate,
		getTrans(cell, beta1) * props_oil.getKr(upwd.s_w, upwd.s_o, cells[upwd_idx].props) /
		props_oil.getViscosity(next.p) / props_oil.getB(next.p, next.p_bub, next.SATUR) *
		(nebr1.p - next.p) - Qcell[cell.num],
		next.p - Pwf);

	h[1] = (next.s_w - nebr1.s_w) / (adouble)(cell.r - beta1.r) - (nebr1.s_w - nebr2.s_w) / (adouble)(beta1.r - beta2.r);

	adouble satur = cell.u_next.SATUR;
	condassign(h[2], satur,
		(next.s_o - nebr1.s_o) / (adouble)(cell.r - beta1.r) - (nebr1.s_o - nebr2.s_o) / (adouble)(beta1.r - beta2.r),
		(next.p_bub - nebr1.p_bub) / (adouble)(cell.r - beta1.r) - (nebr1.p_bub - nebr2.p_bub) / (adouble)(beta1.r - beta2.r));

	for (int i = 0; i < var_size; i++)
		h[i] >>= y[i];

	trace_off();
}
void BlackOil_RZ::solve_eqRight(const Cell& cell)
{
	trace_on(right);
	adouble h[var_size];
	TapeVariable var[Rstencil];
	adouble rightIsPres = rightBoundIsPres;
	for (int i = 0; i < Rstencil; i++)
	{
		var[i].p <<= x[i * Variable::size];
		var[i].s_w <<= x[i * Variable::size + 1];
		var[i].s_o <<= x[i * Variable::size + 2];
		var[i].p_bub <<= x[i * Variable::size + 3];
	}

	const TapeVariable& next = var[0];
	const TapeVariable& nebr = var[1];

	condassign(h[0], rightIsPres, next.p - (adouble)(cell.props->p_out), next.p - (adouble)(nebr.p));
	h[1] = next.s_w - nebr.s_w;
	adouble satur = cell.u_next.SATUR;
	condassign(h[2], satur, next.s_o - nebr.s_o, next.p_bub - nebr.p_bub);

	for (int i = 0; i < var_size; i++)
		h[i] >>= y[i];

	trace_off();
}
void BlackOil_RZ::solve_eqVertical(const Cell& cell)
{
	trace_on(vertical);
	adouble h[var_size];
	TapeVariable var[Vstencil];

	for (int i = 0; i < Vstencil; i++)
	{
		var[i].p <<= x[i * Variable::size];
		var[i].s_w <<= x[i * Variable::size + 1];
		var[i].s_o <<= x[i * Variable::size + 2];
		var[i].p_bub <<= x[i * Variable::size + 3];
	}

	const TapeVariable& next = var[0];
	const TapeVariable& nebr = var[1];

	h[0] = next.p - nebr.p;
	h[1] = next.s_w - nebr.s_w;
	adouble satur = cell.u_next.SATUR;
	condassign(h[2], satur, next.s_o - nebr.s_o, next.p_bub - nebr.p_bub);

	for (int i = 0; i < var_size; i++)
		h[i] >>= y[i];

	trace_off();
}
void BlackOil_RZ::setVariables(const Cell& cell)
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