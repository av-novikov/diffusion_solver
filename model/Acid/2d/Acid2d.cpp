#include "model/Acid/2d/Acid2d.hpp"

#include <cassert>

#include "adolc/drivers/drivers.h"
#include "adolc/adolc.h"

using namespace acid2d;

template <>
const int basic2d::Basic2d<Variable, Properties, Skeleton_Props, TCell, Acid2d>::var_size = Variable::size;
const double acid2d::Component::R = 8.3144598;

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
	setBasicProps(props);

	// Water properties
	props_w = props.props_w;
	props_w.visc = cPToPaSec(props_w.visc);
	// Oil properties
	props_o = props.props_o;
	props_o.visc = cPToPaSec(props_o.visc);
	// Gas properties
	props_g = props.props_g;
	props_g.visc = cPToPaSec(props_g.visc);

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
	// Skeleton properties
	for (int i = 0; i < skeletonsNum; i++)
	{
		props_sk[i].d_pore_r /= R_dim;
		props_sk[i].d_pore_z /= R_dim;
	}

	// Liquid properties
	props_w.visc /= (P_dim * t_dim);
	props_w.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_w.beta /= (1.0 / P_dim);
	// Oil properties
	props_o.visc /= (P_dim * t_dim);
	props_o.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_o.beta /= (1.0 / P_dim);
	// Gas properties
	props_g.visc /= (P_dim * t_dim);
	props_g.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
}

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

	h[0] = next.m * next.sw * props_w.getMolarDensity(next.p, next.xa, next.xw) -
			prev.m * prev.sw * props_w.getMolarDensity(prev.p, prev.xa, prev.xw) - 
			ht * (reac.indices[REACTS::ACID] + reac.indices[REACTS::WATER] + reac.indices[REACTS::SALT]) * rate;
	h[1] = next.m * next.sw * props_w.getMolarDensity(next.p, next.xa, next.xw) * next.xa -
			prev.m * prev.sw * props_w.getMolarDensity(prev.p, prev.xa, prev.xw) * prev.xa - 
			ht * reac.indices[REACTS::ACID] * rate;
	h[2] = next.m * next.sw * props_w.getMolarDensity(next.p, next.xa, next.xw) * next.xw -
			prev.m * prev.sw * props_w.getMolarDensity(prev.p, prev.xa, prev.xw) * prev.xw - 
			ht * reac.indices[REACTS::WATER] * rate;
	h[3] = next.m * next.so * props_o.getMolarDensity(next.p) - prev.m * prev.so * props_o.getMolarDensity(prev.p);
	h[4] = next.m * (1.0 - next.sw - next.so) * props_g.getMolarDensity(next.p) - 
			prev.m * (1.0 - prev.sw - prev.so) * props_g.getMolarDensity(prev.p) - 
			ht * reac.indices[REACTS::CO2] * rate;
	h[5] = (1.0 - next.m) * props.getMolarDensity() -
			(1.0 - prev.m) * props.getMolarDensity() - 
			ht * reac.indices[REACTS::CALCITE] * rate;

	int neighbor[4];
	getNeighborIdx(cell.num, neighbor);
	for (int i = 0; i < 4; i++)
	{
		const Cell& beta = cells[neighbor[i]];
		const int upwd_idx = (getUpwindIdx(cell.num, neighbor[i]) == cell.num) ? 0 : i + 1;
		const TapeVariable& nebr = var[i + 1];
		TapeVariable& upwd = var[upwd_idx];
		
		adouble dens_w = getAvarage(props_w.getMolarDensity(next.p, next.xa, next.xw), cell, props_w.getMolarDensity(nebr.p, nebr.xa, nebr.xw), beta);
		adouble dens_o = getAvarage(props_o.getMolarDensity(next.p), cell, props_o.getMolarDensity(nebr.p), beta);
		adouble dens_g = getAvarage(props_g.getMolarDensity(next.p), cell, props_g.getMolarDensity(nebr.p), beta);
		adouble buf = ht / cell.V * getTrans(cell, next.m, beta, nebr.m) * (next.p - nebr.p);
		adouble buf_w = buf * dens_w * props_w.getKr(upwd.sw, upwd.so, cells[upwd_idx].props) / props_w.getViscosity(upwd.p, upwd.xa, upwd.xw);

		h[0] += buf_w;
		h[1] += buf_w * upwd.xa;
		h[2] += buf_w * upwd.xw;
		h[3] += buf * dens_o * props_o.getKr(upwd.sw, upwd.so, cells[upwd_idx].props) / props_o.getViscosity(upwd.p);
		h[4] += buf * dens_g * props_g.getKr(upwd.sw, upwd.so, cells[upwd_idx].props) / props_g.getViscosity(upwd.p);
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

	for (int i = 0; i < var_size; i++)
		h[i] >>= y[i];

	trace_off();
}
void Acid2d::solve_eqRight(const Cell& cell)
{
	trace_on(right);
	adouble h[var_size];
	TapeVariable var[Rstencil];
	for (int i = 0; i < Rstencil; i++)
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
	h[2] = next.so - nebr.so;
	h[3] = next.sw - nebr.sw;
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