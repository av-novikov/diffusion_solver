#include "model/Acid/2d/Acid2d.hpp"

#include <cassert>

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