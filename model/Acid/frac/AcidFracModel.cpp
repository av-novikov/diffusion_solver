#include "model/Acid/frac/AcidFracModel.hpp"

using namespace acidfrac;

double acidfrac::Component::R = 8.3144598;
double acidfrac::Component::p_std = 101325.0;

AcidFrac::AcidFrac()
{
	grav = 9.8;
	Volume = 0.0;
}
AcidFrac::~AcidFrac()
{
	delete snapshotter;
}
void AcidFrac::setProps(Properties& props)
{
	props_frac = props.props_frac;
	props_sk = props.props_sk;
	skeletonsNum = props.props_sk.size();

	leftBoundIsRate = props.leftBoundIsRate;
	rightBoundIsPres = props.rightBoundIsPres;

	// Setting grid properties
	cellsNum_x = props.cellsNum_x;
	cellsNum_y = props.cellsNum_y;
	cellsNum_z = props.cellsNum_z;
	cellsNum = (cellsNum_x + 2) * (cellsNum_z + 2) * (cellsNum_y + 1);
	cellsNum_y_1d = props.cellsNum_y_1d;
	xe = props.xe;

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
	for (int i = 0; i < skeletonsNum; i++)
	{
		auto& sk = props_sk[i];
		sk.perm /= (R_dim * R_dim);
		sk.beta /= (1.0 / P_dim);
		sk.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
		sk.p_init /= P_dim;
		sk.p_out /= P_dim;
		sk.p_ref /= P_dim;
		sk.hx /= R_dim;
		sk.hz /= R_dim;

		xe[i] /= R_dim;
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
size_t AcidFrac::getInitId2OutCell(const FracCell& cell)
{
	assert(cell.type == FracType::FRAC_OUT);
	return 0;
}
void AcidFrac::build1dGrid(const FracCell& cell)
{
	assert(cell.type == FracType::FRAC_OUT);
	poro_grids.push_back(PoroGrid()); 
	PoroGrid& grid = poro_grids.back();
	frac2poro[&const_cast<FracCell&>(cell)] = &grid;
	grid.frac_nebr = &cell;
	grid.Volume = 0.0;
	const auto init_id = getInitId2OutCell(cell);
	grid.cellsNum = cellsNum_y_1d[init_id];
	grid.props_sk = &props_sk[init_id];
	auto& cells = grid.cells;

	cells.reserve(grid.cellsNum);
	int counter = 0;

	double y_prev = cell.y + cell.hy / 2.0;
	double logMax = log((xe[init_id] + y_prev) / y_prev);
	double logStep = logMax / (double)grid.cellsNum;

	grid.hx = cell.hx;	grid.hz = cell.hz;
	double hy = y_prev * (exp(logStep) - 1.0);
	double cm_y = y_prev;

	// Left border
	cells.push_back(PoroCell(counter++, cm_y, 0.0, grid.hx, grid.hz, PoroType::WELL_LAT));
	// Middle cells
	for (int j = 0; j < cellsNum_x; j++)
	{
		cm_y = y_prev * (exp(logStep) + 1.0) / 2.0;
		hy = y_prev * (exp(logStep) - 1.0);
		cells.push_back(PoroCell(counter++, cm_y, hy, grid.hx, grid.hz, PoroType::MIDDLE));
		grid.Volume += cells.back().V;
		y_prev *= exp(logStep);
	}
	// Right border
	cm_y = xe[init_id] + cell.y + cell.hy / 2.0;
	cells.push_back(PoroCell(counter++, cm_y, 0.0, grid.hx, grid.hz, PoroType::RIGHT));
}
void AcidFrac::buildGrid()
{
	int counter = 0;
	double hy = 0.0, hz = 0.0, init_dx = props_frac.w2;
	FracType cur_type;

	poro_grids.reserve(cellsNum_x * cellsNum_z);
	cells_frac.reserve(cellsNum);
	double x_prev = init_dx;
	double logMax = log((props_frac.l2 + init_dx) / init_dx);
	double logStep = logMax / (double)cellsNum_x;
	double x = x_prev, y = 0.0, z = 0.0;

	//double hz = props_sk.h2 - props_sk.h1;
	//double cm_z = props_sk[skel_idx].h1;
	double hx = x_prev * (exp(logStep) - 1.0);

	for (int i = 0; i < cellsNum_x + 2; i++)
	{
		if (i == 0)
		{
			cur_type = FracType::FRAC_IN;
			hx = 0.0;
		}
		else if (i == cellsNum_x + 1)
		{
			cur_type = FracType::FRAC_BORDER;
			hx = 0.0;
			x = props_frac.l2 + init_dx;
		}
		else
		{
			cur_type = FracType::FRAC_MID;
			x = x_prev * (exp(logStep) + 1.0) / 2.0;
			hx = x_prev * (exp(logStep) - 1.0);
		}

		z = 0.0;	hz = 0.0;
		for (int k = 0; k < cellsNum_z + 2; k++)
		{
			if (k == 0 || k == cellsNum_z + 1)
			{
				cur_type = FracType::FRAC_BORDER;
				hz = 0.0;
			}
			else
				hz = props_frac.height / cellsNum_z;

			y = 0.0;	hy = 0.0;
			for (int j = 0; j < cellsNum_y + 1; j++)
			{
				if (j == 0)
				{
					cur_type = FracType::FRAC_BORDER;
					hy = 0.0;
				}
				else if (j == cellsNum_y)
				{
					cur_type = FracType::FRAC_OUT;
					hy = init_dx / cellsNum_y;
				}
				else 
					hy = init_dx / cellsNum_y;

				if (cur_type == FracType::FRAC_OUT && hx * hz == 0.0)
					cur_type = FracType::FRAC_BORDER;

				cells_frac.push_back(FracCell(counter++, x - init_dx, y, z, hx, hy, hz, cur_type));
				Volume += cells_frac.back().V;

				if (j == 0)
					y += init_dx / cellsNum_y / 2.0;
				else
					y += hy;
			}

			if(i != 0 && k != 0 && i != cellsNum_x + 1 && k != cellsNum_z + 1)
				build1dGrid(cells_frac.back());

			if (k == 0 || k == cellsNum_z)
				z += props_frac.height / cellsNum_z / 2.0;
			else
				z += hz;
		}

		if(i != 0)
			x_prev *= exp(logStep);
	}
}
void AcidFrac::setPerforated()
{
	height_perf = 0.0;
	vector<pair<int, int> >::iterator it;
	for (const auto& cell : cells_frac)
	{
		if (cell.type == FracType::FRAC_IN)
		{
			Qcell[cell.num] = 0.0;
			height_perf += cell.hz;
		}
	}
};
void AcidFrac::setPeriod(int period)
{
	if (leftBoundIsRate)
	{
		Q_sum = rate[period];

		if (period == 0 || rate[period - 1] < EQUALITY_TOLERANCE) {
			std::map<int, double>::iterator it;
			for (it = Qcell.begin(); it != Qcell.end(); ++it)
				it->second = Q_sum * cells_frac[it->first].hz / height_perf;
		}
		else {
			std::map<int, double>::iterator it;
			for (it = Qcell.begin(); it != Qcell.end(); ++it)
				it->second = it->second * Q_sum / rate[period - 1];
		}
	}
	else
	{
		Pwf = pwf[period];
		Q_sum = 0.0;
	}
	c = cs[period];
}
void AcidFrac::setInitialState() 
{
	for (auto& cell : cells_frac)
	{
		cell.u_prev.p = cell.u_iter.p = cell.u_next.p = props_frac.p_init;
		cell.u_prev.c = cell.u_iter.c = cell.u_next.c = props_frac.c_init;

		if (cell.type == FracType::FRAC_OUT)
		{
			const auto& poro_grid = frac2poro[&cell];
			const auto& props = poro_grid->props_sk;
			for (auto& poro_cell : poro_grid->cells)
			{
				poro_cell.u_prev.m = poro_cell.u_iter.m = poro_cell.u_next.m = props->m_init;
				poro_cell.u_prev.p = poro_cell.u_iter.p = poro_cell.u_next.p = props->p_init;
				poro_cell.u_prev.sw = poro_cell.u_iter.sw = poro_cell.u_next.sw = props->sw_init;
				poro_cell.u_prev.xw = poro_cell.u_iter.xw = poro_cell.u_next.xw = props->xw_init;
				poro_cell.u_prev.xa = poro_cell.u_iter.xa = poro_cell.u_next.xa = props->xa_init;
				poro_cell.u_prev.xs = poro_cell.u_iter.xs = poro_cell.u_next.xs = 0.0;
			}
		}
	}
}