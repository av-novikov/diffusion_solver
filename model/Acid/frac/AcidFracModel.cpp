#include "model/Acid/frac/AcidFracModel.hpp"

using namespace acidfrac;

double acidfrac::Component::R = 8.3144598;
double acidfrac::Component::p_std = 101325.0;

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
	props_frac = props.props_frac;
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
void AcidFrac::makeDimLess()
{
	R_dim = props_frac.l2 / 10.0;
	t_dim = 3600.0;
	P_dim = props_sk[0].p_init;
	// Temporal properties
	ht /= t_dim;
	ht_min /= t_dim;
	ht_max /= t_dim;
	// Skeleton properties
	for (auto& sk : props_sk)
	{
		sk.perm /= (R_dim * R_dim);

		sk.beta /= (1.0 / P_dim);
		sk.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
		sk.h1 /= R_dim;
		sk.h2 /= R_dim;
		sk.height /= R_dim;
		sk.p_init /= P_dim;
		sk.p_out /= P_dim;
		sk.p_ref /= P_dim;
		sk.hx /= R_dim;
		sk.hz /= R_dim;
		for (int j = 0; j < periodsNum; j++)
		{
			sk.perms_eff[j] /= (R_dim * R_dim);
			sk.radiuses_eff[j] /= R_dim;
		}
	}

	Q_dim = R_dim * R_dim * R_dim / t_dim;
	for (int i = 0; i < periodsNum; i++)
	{
		period[i] /= t_dim;
		if (leftBoundIsRate)
			rate[i] /= Q_dim;
		else
			pwf[i] /= P_dim;
	}

	grav /= (R_dim / t_dim / t_dim);
	Component::p_std /= P_dim;
	Component::R /= (P_dim * R_dim * R_dim * R_dim / T_dim);
	Component::T /= T_dim;

	props_w.visc /= (P_dim * t_dim);
	props_w.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_w.beta /= (1.0 / P_dim);
	props_o.visc /= (P_dim * t_dim);
	props_o.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_o.gas_dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_o.beta /= (1.0 / P_dim);
	props_o.p_ref /= P_dim;
	props_g.visc /= (P_dim * t_dim);
	props_g.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_g.co2.mol_weight /= (P_dim * t_dim * t_dim * R_dim);
	props_g.co2.rho_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);

	reac.activation_energy /= (P_dim * R_dim * R_dim * R_dim);
	reac.surf_init /= (1.0 / R_dim);
	reac.reaction_const /= (R_dim / t_dim);
	for (auto& comp : reac.comps)
	{
		comp.rho_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
		comp.mol_weight /= (P_dim * t_dim * t_dim * R_dim);
	}

	props_frac.l2 /= R_dim;
	props_frac.w2 /= R_dim;
	props_frac.height /= R_dim;
	props_frac.p_init /= P_dim;
}
void AcidFrac::buildGrid()
{
	/*int counter = 0;
	const double hx = l2 / cellsNum_x;
	const double hy = w2 / cellsNum_y;
	const double hz = height / cellsNum_z;
	double x = 0.0, y = 0.0, z = 0.0;
	Type cur_type;
	for (int i = 0; i < cellsNum_x + 2; i++)
		{
		for (int j = 0; j < cellsNum_y + 1; j++)
			{
			for (int k = 0; k < cellsNum_z + 2; k++)
				{
				
					cells.push_back(Cell(counter++, x, y, z, hx, hy, hz, cur_type));
				if (k == 0 || k == cellsNum_z)
					z += hz / 2.0;
				else
					z += hz;
				}
			if (j == 0)
				y += hy / 2.0;
			else
				y += hy;
			}
		if (i == 0 || i == cellsNum_x)
			x += hx / 2.0;
		else
			x += hx;
		}*/
}