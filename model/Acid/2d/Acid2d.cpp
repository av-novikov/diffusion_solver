#include "model/Acid/2d/Acid2d.hpp"

#include <cassert>

#include "adolc/drivers/drivers.h"
#include "adolc/adolc.h"

using namespace acid2d;

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

	// Oil properties
	props_l = props.props_l;
	props_l.visc = cPToPaSec(props_l.visc);
	
	// Gas properties
	props_g = props.props_g;
	props_g.visc = cPToPaSec(props_g.visc);

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

	// Gas properties
	props_g.visc /= (P_dim * t_dim);
	props_g.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
}

void Acid2d::solve_eqMiddle(const Cell& cell)
{

}
void Acid2d::solve_eqLeft(const Cell& cell)
{
}
void Acid2d::solve_eqRight(const Cell& cell)
{
}
void Acid2d::solve_eqVertical(const Cell& cell)
{
}