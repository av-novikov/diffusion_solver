#include "model/BlackOil_RZ/BlackOil_RZ.hpp"

using namespace blackoil_rz;

BlackOil_RZ::BlackOil_RZ()
{
	x = new double[stencil * (Variable::size - 1)];
	y = new double[Variable::size - 1];

	jac = new double*[Variable::size - 1];
	for (int i = 0; i < Variable::size - 1; i++)
		jac[i] = new double[stencil * Variable::size];
};
BlackOil_RZ::~BlackOil_RZ()
{
	delete x;
	delete y;

	for (int i = 0; i < Variable::size - 1; i++)
		delete[] jac[i];
	delete[] jac;
};
void BlackOil_RZ::setProps(Properties& props)
{
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
	props_wat.kr = setDataset(props.kr_wat, 1.0, 1.0);
	props_oil.kr = setDataset(props.kr_oil, 1.0, 1.0);
	props_gas.kr = setDataset(props.kr_gas, 1.0, 1.0);
	props_wat.b = setDataset(props.B_wat, P_dim / BAR_TO_PA, 1.0);
	props_oil.b = setDataset(props.B_oil, P_dim / BAR_TO_PA, 1.0);
	props_gas.b = setDataset(props.B_gas, P_dim / BAR_TO_PA, 1.0);
	props_oil.Rs = setDataset(props.Rs, P_dim / BAR_TO_PA, 1.0);
}
void BlackOil_RZ::makeDimLess()
{
	props_wat.visc /= (P_dim * t_dim);
	props_wat.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_wat.beta /= (1.0 / P_dim);
	props_wat.p_sat /= P_dim;
	props_oil.visc /= (P_dim * t_dim);
	props_oil.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_oil.beta /= (1.0 / P_dim);
	props_oil.p_sat /= P_dim;
	props_gas.visc /= (P_dim * t_dim);
	props_gas.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
}

void BlackOil_RZ::solve_eqMiddle(const Cell& cell)
{

}
void BlackOil_RZ::solve_eqLeft(const Cell& cell)
{
}
void BlackOil_RZ::solve_eqRight(const Cell& cell)
{
}
void BlackOil_RZ::solve_eqVertical(const Cell& cell)
{
}