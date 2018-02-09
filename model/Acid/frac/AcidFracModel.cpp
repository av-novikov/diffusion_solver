#include "model/Acid/frac/AcidFracModel.hpp"

using namespace acidfrac;

AcidFrac::AcidFrac()
{
	isWriteSnaps = true;
	grav = 9.8;
	snapshotter = new VTKSnapshotter<AcidFrac>();
}
AcidFrac::~AcidFrac()
{
	delete snapshotter;
}
void AcidFrac::setProps(Properties& props)
{
	props_sk = props.props_sk;

	leftBoundIsRate = props.leftBoundIsRate;
	rightBoundIsPres = props.rightBoundIsPres;

	// Setting grid properties
	cellsNum_x = props.cellsNum_x;
	cellsNum_y = props.cellsNum_y;
	cellsNum_z = props.cellsNum_z;
	cellsNum = (cellsNum_x + 2) * (cellsNum_z + 2) * (cellsNum_y + 1);

	skeletonsNum = props.props_sk.size();
	for (int j = 0; j < skeletonsNum; j++)
	{
		props_sk[j].perm = MilliDarcyToM2(props_sk[j].perm);
	}

	periodsNum = props.timePeriods.size();
	for (int i = 0; i < periodsNum; i++)
	{
		cs.push_back(props.cs[i]);
		period.push_back(props.timePeriods[i]);
		if (leftBoundIsRate)
			rate.push_back(props.rates[i] / 86400.0);
		else
			pwf.push_back(props.pwf[i]);
	}

	// Temporal properties
	ht = props.ht;
	ht_min = props.ht_min;
	ht_max = props.ht_max;

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
	makeDimLess();

	props_o.b = setDataset(props.B_oil, P_dim / BAR_TO_PA, 1.0);
	props_o.Rs = setDataset(props.Rs, P_dim / BAR_TO_PA, 1.0);
	props_g.rho = setDataset(props.rho_co2, P_dim / BAR_TO_PA, (P_dim * t_dim * t_dim / R_dim / R_dim));
}