#include "model/Acid/2d/Acid2d.hpp"

#include <cassert>

#include "adolc/drivers/drivers.h"
#include "adolc/adolc.h"

using namespace acid2d;

Acid2d::Acid2d()
{
	x = new double[stencil * Variable::size];
	y = new double[Variable::size];
}
Acid2d::~Acid2d()
{
	delete x;
	delete y;
}
void Acid2d::setProps(Properties& props)
{
	setBasicProps(props);

	// Oil properties
	props_l = props.props_l;
	//props_l.visc = cPToPaSec(props_l.visc);

	// Gas properties
	props_g = props.props_g;
	//props_gas.visc = cPToPaSec(props_gas.visc);

	makeBasicDimLess();
	makeDimLess();

	// Data sets
	props_l.kr = setDataset(props.kr_l, 1.0, 1.0);
	props_g.kr = setDataset(props.kr_g, 1.0, 1.0);

	//props_oil.b = setDataset(props.B_oil, P_dim / BAR_TO_PA, 1.0);
	//props_gas.b = setDataset(props.B_gas, P_dim / BAR_TO_PA, 1.0);
	//props_oil.Rs = setDataset(props.Rs, P_dim / BAR_TO_PA, 1.0);
}
void Acid2d::makeDimLess()
{
	// Liquid properties
	props_l.visc /= (P_dim * t_dim);
	props_l.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_l.beta /= (1.0 / P_dim);
	props_l.p_sat /= P_dim;

	// Gas properties
	props_g.visc /= (P_dim * t_dim);
	props_g.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
}

void Acid2d::solve_eqMiddle(int cur)
{

}
void Acid2d::solve_eqLeft(int cur)
{
}
void Acid2d::solve_eqRight(int cur)
{
}
void Acid2d::solve_eqVertical(int cur)
{
}
void Acid2d::setVariables(int cur)
{
	if (cur < cellsNum_z + 2)
	{
		// Left
		const Variable& next = cells[cur].u_next;
		const Variable& nebr = cells[cur + cellsNum_z + 2].u_next;

		for (int i = 0; i < Variable::size; i++)
		{
			x[i] = next.values[i];
			x[Variable::size + i] = nebr.values[i];
		}
	}
	else if (cur >= (cellsNum_z + 2) * (cellsNum_r + 1))
	{
		// Right
		const Variable& next = cells[cur].u_next;
		const Variable& nebr1 = cells[cur - cellsNum_z - 2].u_next;
		const Variable& nebr2 = cells[cur - 2 * cellsNum_z - 4].u_next;

		for (int i = 0; i < Variable::size; i++)
		{
			x[i] = next.values[i];
			x[Variable::size + i] = nebr1.values[i];
			x[2 * Variable::size + i] = nebr2.values[i];
		}
	}
	else if (cur % (cellsNum_z + 2) == 0)
	{
		// Top
		const Variable& next = cells[cur].u_next;
		const Variable& nebr = cells[cur + 1].u_next;

		for (int i = 0; i < Variable::size; i++)
		{
			x[i] = next.values[i];
			x[Variable::size + i] = nebr.values[i];
		}
	}
	else if ((cur + 1) % (cellsNum_z + 2) == 0)
	{
		// Bottom
		const Variable& next = cells[cur].u_next;
		const Variable& nebr = cells[cur - 1].u_next;

		for (int i = 0; i < Variable::size; i++)
		{
			x[i] = next.values[i];
			x[Variable::size + i] = nebr.values[i];
		}
	}
	else
	{
		// Middle
		const Variable& next = cells[cur].u_next;
		int neighbor[4];
		getNeighborIdx(cur, neighbor);

		for (int i = 0; i < Variable::size; i++)
		{
			x[i] = next.values[i];

			for (int j = 0; j < 4; j++)
			{
				const Variable& nebr = cells[neighbor[j]].u_next;
				x[(j + 1) * Variable::size + i] = nebr.values[i];
			}
		}
	}
}