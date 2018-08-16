#include "model/Acid/recfrac/AcidRecFracModel.hpp"
#include <assert.h>

using namespace acidrecfrac;

double acidrecfrac::Component::R = 8.3144598;
double acidrecfrac::Component::p_std = 101325.0;

AcidRecFrac::AcidRecFrac()
{
	grav = 9.8;
	Volume_frac = Volume_poro = 0.0;
}
AcidRecFrac::~AcidRecFrac()
{
	delete snapshotter;
	//delete[] x_frac, x_poro, h;
}
void AcidRecFrac::setProps(Properties& props)
{
	props_frac = props.props_frac;
	props_sk = props.props_sk;
	skeletonsNum = props.props_sk.size();

	rightBoundIsPres = props.rightBoundIsPres;

	// Setting grid properties
	cellsNum_y_frac = props.cellsNum_y_frac;
	cellsNum_y_poro = props.cellsNum_y_poro;
	cellsNum_x = props.cellsNum_x;
	cellsNum_z = props.cellsNum_z;

	cellsNum_frac = (cellsNum_y_frac + 1) * (cellsNum_x + 2) * (cellsNum_z + 2);
	cellsNum_poro = (cellsNum_y_poro + 2) * (cellsNum_x + 2) * (cellsNum_z + 2);
	re = props.re;

	skeletonsNum = props.props_sk.size();
	for (int j = 0; j < skeletonsNum; j++)
	{
		props_sk[j].perm = MilliDarcyToM2(props_sk[j].perm);
	}

	periodsNum = props.timePeriods.size();
	rate.resize(periodsNum);
	pwf.resize(periodsNum);
	int rate_idx = 0, pres_idx = 0;
	for (int i = 0; i < periodsNum; i++)
	{
		LeftBoundIsRate.push_back(props.LeftBoundIsRate[i]);
		cs.push_back(props.cs[i]);
		period.push_back(props.timePeriods[i]);
		if (LeftBoundIsRate.back())
		{
			rate[i] = props.rates[rate_idx++] / 86400.0;
			pwf[i] = 0.0;
		}
		else
		{
			pwf[i] = props.pwf[pres_idx++];
			rate[i] = 0.0;
		}
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
	props_frac.w2_avg = props_frac.w2 * M_PI / 4.0;

	props_o.b = setDataset(props.B_oil, P_dim / BAR_TO_PA, 1.0);
	props_o.Rs = setDataset(props.Rs, P_dim / BAR_TO_PA, 1.0);
	props_g.rho = setDataset(props.rho_co2, P_dim / BAR_TO_PA, (P_dim * t_dim * t_dim / R_dim / R_dim));
}
void AcidRecFrac::makeDimLess()
{
	T_dim = props_sk[0].t_init;
	R_dim = props_frac.l2 * 10.0;
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
		sk.t_init /= T_dim;
	}

	Q_dim = R_dim * R_dim * R_dim / t_dim;
	for (int i = 0; i < periodsNum; i++)
	{
		period[i] /= t_dim;
		rate[i] /= Q_dim;
		pwf[i] /= P_dim;
	}

	grav /= (R_dim / t_dim / t_dim);
	Component::p_std /= P_dim;
	Component::R /= (P_dim * R_dim * R_dim * R_dim / T_dim);
	Component::T /= T_dim;

	props_w.visc /= (P_dim * t_dim);
	props_w.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_w.beta /= (1.0 / P_dim);
	props_w.D_e /= (R_dim * R_dim / t_dim);
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

	re /= R_dim;
	props_frac.l2 /= R_dim;
	props_frac.w2 /= R_dim;
	props_frac.height /= R_dim;
	props_frac.p_init /= P_dim;
}
void AcidRecFrac::buildFracGrid()
{
	int counter = 0;
	cells_frac.reserve(cellsNum_frac);
	Volume_frac = 0.0;
	
	double hx = 0.0, hy = 0.0, hz = 0.0;
	double cx = 0.0, cy = 0.0, cz = 0.0;
	FracType cur_type;

	// Left border
	cur_type = FracType::FRAC_IN;
	hx = 0.0;	 cx = 0.0;
	for (int k = 0; k < cellsNum_z + 2; k++)
	{
		if (k == 0)
			hz = cz = 0.0;
		else if (k == cellsNum_z + 1)
		{
			hz = 0.0;
			cz = props_frac.height;
		}
		else
		{
			hz = props_frac.height / cellsNum_z;
			cz = ((double)k - 0.5) * hz;
		}

		for (int i = 0; i < cellsNum_y_frac + 1; i++)
		{
			if (i == 0)
				hy = cy = 0.0;
			else
			{
				hy = props_frac.w2 / (double)cellsNum_y_frac;
				cy = ((double)i - 0.5) * hy;
			}

			cells_frac.push_back( FracCell(counter++, cx, cy, cz, hx, hy, hz, cur_type) );
			auto& cell = cells_frac.back();
			//cell.nebrs[0] = cell.num + (cellsNum_y_frac + 1) * (cellsNum_z + 2);
		}
	}
	// Middle
	hx = props_frac.l2 / cellsNum_x;	 cx = hx / 2;
	for (int j = 0; j < cellsNum_x; j++)
	{
		for (int k = 0; k < cellsNum_z + 2; k++)
		{
			if (k == 0)
				hz = cz = 0.0;
			else if (k == cellsNum_z + 1)
			{
				hz = 0.0;
				cz = props_frac.height;
			}
			else
			{
				hz = props_frac.height / cellsNum_z;
				cz = ((double)k - 0.5) * hz;
			}

			for (int i = 0; i < cellsNum_y_frac + 1; i++)
			{
				if (i == 0)
					hy = cy = 0.0;
				else
				{
					hy = props_frac.w2 / (double)cellsNum_y_frac;
					cy = ((double)i - 0.5) * hy;
				}

				if (i == cellsNum_y_frac)
					cur_type = FracType::FRAC_OUT;
				else if (i == 0)
					cur_type = FracType::FRAC_BORDER;
				else
					cur_type = FracType::FRAC_MID;

				if(k == 0 || k == cellsNum_z + 1)
					cur_type = FracType::FRAC_BORDER;

				cells_frac.push_back(FracCell(counter++, cx, cy, cz, hx, hy, hz, cur_type));
				auto& cell = cells_frac.back();
				Volume_frac += cell.V;
				/*if (cell.type != FracType::FRAC_BORDER)
				{
					cell.nebrs[0] = cell.num - 1;
					if(cell.type == FracType::FRAC_MID)
						cell.nebrs[1] = cell.num + 1;
					else if (cell.type == FracType::FRAC_OUT)
						cell.nebrs[1] = (j + 1) * (cellsNum_mu_poro + 2) * (cellsNum_z + 2) + k * (cellsNum_mu_poro + 2);
					cell.nebrs[2] = cell.num - (cellsNum_mu_frac + 1) * (cellsNum_z + 2);
					cell.nebrs[3] = cell.num + (cellsNum_mu_frac + 1) * (cellsNum_z + 2);
					cell.nebrs[4] = cell.num - (cellsNum_mu_frac + 1);
					cell.nebrs[5] = cell.num + (cellsNum_mu_frac + 1);
				}
				else
				{
					if (i == 0)
						cell.nebrs[0] = cell.num + 1;
					if (k == 0)
						cell.nebrs[0] = cell.num + (cellsNum_mu_frac + 1);
					else if (k == cellsNum_z + 1)
						cell.nebrs[0] = cell.num - (cellsNum_mu_frac + 1);
				}*/
			}
		}
		cx += hx;
	}
	// Right cells
	cur_type = FracType::FRAC_BORDER;
	cx = props_frac.l2;	 hx = 0.0;
	for (int k = 0; k < cellsNum_z + 2; k++)
	{
		if (k == 0)
			hz = cz = 0.0;
		else if (k == cellsNum_z + 1)
		{
			hz = 0.0;
			cz = props_frac.height;
		}
		else
		{
			hz = props_frac.height / cellsNum_z;
			cz = ((double)k - 0.5) * hz;
		}

		for (int i = 0; i < cellsNum_y_frac + 1; i++)
		{
			if (i == 0)
				hy = cy = 0.0;
			else
			{
				hy = props_frac.w2 / (double)cellsNum_y_frac;
				cy = ((double)i - 0.5) * hy;
			}

			cells_frac.push_back(FracCell(counter++, cx, cy, cz, hx, hy, hz, cur_type));
			auto& cell = cells_frac.back();
			//cell.nebrs[0] = cell.num - (cellsNum_y_frac + 1) * (cellsNum_z + 2);
		}
	}
}
void AcidRecFrac::buildPoroGrid()
{
	int counter = 0;
	cells_poro.reserve(cellsNum_poro);
	Volume_poro = 0.0;

	double r_prev = re;
	double logMax = log(re / props_frac.w2);
	double logStep = logMax / (double)cellsNum_y_poro;

	double hy = 0.0, hx = 0.0, hz = 0.0;
	double cy = 0.0, cx = 0.0, cz = 0.0;
	PoroType cur_type;

	// Left border
	cur_type = PoroType::SIDE_LEFT;
	hx = 0.0;	 cx = 0.0;
	for (int k = 0; k < cellsNum_z + 2; k++)
	{
		if (k == 0)
			hz = cz = 0.0;
		else if (k == cellsNum_z + 1)
		{
			hz = 0.0;
			cz = props_frac.height;
		}
		else
		{
			hz = props_frac.height / cellsNum_z;
			cz = ((double)k - 0.5) * hz;
		}

		for (int i = 0; i < cellsNum_y_poro + 2; i++)
		{
			if (i == 0)
			{
				hy = 0.0;
				r_prev = cy = props_frac.w2;
			}
			else if (i == cellsNum_y_poro + 1)
			{
				hy = 0.0;
				cy = re;
			}
			else
			{
				//hmu = (mu_e - mu_w) / (double)cellsNum_mu_poro;
				//cmu = mu_w + ((double)i - 0.5) * hmu;
				cy = r_prev * (exp(logStep) + 1.0) / 2.0;
				hy = r_prev * (exp(logStep) - 1.0);
				r_prev *= exp(logStep);
			}

			cells_poro.push_back(PoroCell(counter++, cx, cy, cz, hx, hy, hz, cur_type));
			auto& cell = cells_poro.back();
			//cell.nebrs[0] = cell.num + (cellsNum_y_poro + 2) * (cellsNum_z + 2);
		}
	}
	// Middle
	hx = props_frac.l2 / cellsNum_x;	 cx = hx / 2;
	for (int j = 0; j < cellsNum_x; j++)
	{
		for (int k = 0; k < cellsNum_z + 2; k++)
		{
			if (k == 0)
				hz = cz = 0.0;
			else if (k == cellsNum_z + 1)
			{
				hz = 0.0;
				cz = props_frac.height;
			}
			else
			{
				hz = props_frac.height / cellsNum_z;
				cz = ((double)k - 0.5) * hz;
			}

			for (int i = 0; i < cellsNum_y_poro + 2; i++)
			{
				if (i == 0)
				{
					hy = 0.0;
					r_prev = cy = props_frac.w2;
				}
				else if (i == cellsNum_y_poro + 1)
				{
					hy = 0.0;
					cy = re;
				}
				else
				{
					//hmu = (mu_e - mu_w) / (double)cellsNum_mu_poro;
					//cmu = mu_w + ((double)i - 0.5) * hmu;
					cy = r_prev * (exp(logStep) + 1.0) / 2.0;
					hy = r_prev * (exp(logStep) - 1.0);
					r_prev *= exp(logStep);
				}

				if (i == cellsNum_y_poro + 1)
					cur_type = PoroType::RIGHT;
				else if (i == 0)
					cur_type = PoroType::WELL_LAT;
				else
					cur_type = PoroType::MIDDLE;

				if (k == 0)
					cur_type = PoroType::BOTTOM;
				else if (k == cellsNum_z + 1)
					cur_type = PoroType::TOP;

				cells_poro.push_back(PoroCell(counter++, cx, cy, cz, hx, hy, hz, cur_type));
				auto& cell = cells_poro.back();
				Volume_poro += cell.V;
				/*if (cell.type == PoroType::MIDDLE)
				{
					cell.nebrs[0] = cell.num - 1;
					cell.nebrs[1] = cell.num + 1;
					cell.nebrs[2] = cell.num - (cellsNum_mu_poro + 2) * (cellsNum_z + 2);
					cell.nebrs[3] = cell.num + (cellsNum_mu_poro + 2) * (cellsNum_z + 2);
					cell.nebrs[4] = cell.num - (cellsNum_mu_poro + 2);
					cell.nebrs[5] = cell.num + (cellsNum_mu_poro + 2);
				}
				else if (cell.type == PoroType::WELL_LAT)
				{
					cell.nebrs[0] = (j + 1) * (cellsNum_mu_frac + 1) * (cellsNum_z + 2) + k * (cellsNum_mu_frac + 1) + cellsNum_mu_frac;
					cell.nebrs[1] = cell.num + 1;
					cell.nebrs[2] = cell.num - (cellsNum_mu_poro + 2) * (cellsNum_z + 2);
					cell.nebrs[3] = cell.num + (cellsNum_mu_poro + 2) * (cellsNum_z + 2);
					cell.nebrs[4] = cell.num - (cellsNum_mu_poro + 2);
					cell.nebrs[5] = cell.num + (cellsNum_mu_poro + 2);
				}
				else if (cell.type == PoroType::BOTTOM)
					cell.nebrs[0] = cell.num + (cellsNum_mu_poro + 2);
				else if (cell.type == PoroType::TOP)
					cell.nebrs[0] = cell.num - (cellsNum_mu_poro + 2);
				else if (cell.type == PoroType::RIGHT)
					cell.nebrs[0] = cell.num - 1;*/
			}
		}
		cx += hx;
	}
	// Right cells
	cur_type = PoroType::SIDE_RIGHT;
	cx = props_frac.l2;	 hx = 0.0;
	for (int k = 0; k < cellsNum_z + 2; k++)
	{
		if (k == 0)
			hz = cz = 0.0;
		else if (k == cellsNum_z + 1)
		{
			hz = 0.0;
			cz = props_frac.height;
		}
		else
		{
			hz = props_frac.height / cellsNum_z;
			cz = ((double)k - 0.5) * hz;
		}

		for (int i = 0; i < cellsNum_y_poro + 2; i++)
		{
			if (i == 0)
			{
				hy = 0.0;
				r_prev = cy = props_frac.w2;
			}
			else if (i == cellsNum_y_poro + 1)
			{
				hy = 0.0;
				cy = re;
			}
			else
			{
				//hmu = (mu_e - mu_w) / (double)cellsNum_mu_poro;
				//cmu = mu_w + ((double)i - 0.5) * hmu;
				cy = r_prev * (exp(logStep) + 1.0) / 2.0;
				hy = r_prev * (exp(logStep) - 1.0);
				r_prev *= exp(logStep);
			}

			cells_poro.push_back(PoroCell(counter++, cx, cy, cz, hx, hy, hz, cur_type));
			auto& cell = cells_poro.back();
			//cell.nebrs[0] = cell.num - (cellsNum_y_poro + 2) * (cellsNum_z + 2);
		}
	}

}
void AcidRecFrac::processGeometry()
{
/*	// Change FRAC_IN type among zero-volume cells
	for (auto& cell : cells_frac)
		if (cell.type == FracType::FRAC_IN)
			if (cell.h.mu * cell.h.z == 0.0)
				cell.type = FracType::FRAC_BORDER;
	// Neighbours check & set indices
	for (auto& cell : cells_frac)
	{
		for(size_t i = 0; i < NEBRS_NUM; i++)
			if (cell.nebrs[i] >= 0)
			{
				if (cell.type == FracType::FRAC_OUT && i == 1)
				{
					const auto& nebr_cell = cells_poro[cell.nebrs[i]];
					const auto& nebrs = nebr_cell.nebrs;
					const ptrdiff_t idx = std::find(nebrs.begin(), nebrs.end(), cell.num) - nebrs.begin();
					if (idx >= 0 && idx < nebrs.size())
					{
						cell.nebrs_idx[i] = idx;
						assert(nebrs[idx] == cell.num);
					}
				}
				else
				{
					const auto& nebr_cell = cells_frac[cell.nebrs[i]];
					const auto& nebrs = nebr_cell.nebrs;
					const ptrdiff_t idx = std::find(nebrs.begin(), nebrs.end(), cell.num) - nebrs.begin();
					if (idx >= 0 && idx < nebrs.size())
					{
						cell.nebrs_idx[i] = idx;
						assert(nebrs[idx] == cell.num);
					}
				}
			}
	}
	for (auto& cell : cells_poro)
	{
		for (size_t i = 0; i < NEBRS_NUM; i++)
			if (cell.nebrs[i] >= 0)
			{
				if (cell.type == PoroType::WELL_LAT && i == 0)
				{
					const auto& nebr_cell = cells_frac[cell.nebrs[i]];
					const auto& nebrs = nebr_cell.nebrs;
					const ptrdiff_t idx = std::find(nebrs.begin(), nebrs.end(), cell.num) - nebrs.begin();
					if (idx >= 0 && idx < nebrs.size())
					{
						cell.nebrs_idx[i] = idx;
						assert(nebrs[idx] == cell.num);
					}
				}
				else
				{
					const auto& nebr_cell = cells_poro[cell.nebrs[i]];
					const auto& nebrs = nebr_cell.nebrs;
					const ptrdiff_t idx = std::find(nebrs.begin(), nebrs.end(), cell.num) - nebrs.begin();
					if (idx >= 0 && idx < nebrs.size())
					{
						cell.nebrs_idx[i] = idx;
						assert(nebrs[idx] == cell.num);
					}
				}
			}
	}
	// Set lengths & faces
	for (auto& cell : cells_frac)
		for (size_t i = 0; i < NEBRS_NUM; i++)
			if (cell.nebrs[i] != -1)
			{
				if (cell.type == FracType::FRAC_OUT && i == 1)
				{
					const auto& nebr = cells_poro[cell.nebrs[i]];
					const Point p_mid = (cell.c + nebr.c) / 2.0;
					cell.faces_dist[i] = point::distance(cell.c, p_mid, cell.c);//(cell.h * (nebr.c - cell.c) / (nebr.c - cell.c).dist() / 2.0).norm();
					double S = (cell.V + nebr.V) / (2.0 * point::distance(cell.c, nebr.c, p_mid));
					fmap_inter[{cell.num, nebr.num}] = Face{ S };
				}
				else
				{
					const auto& nebr = cells_frac[cell.nebrs[i]];
					const Point p_mid = (cell.c + nebr.c) / 2.0;
					cell.faces_dist[i] = point::distance(cell.c, p_mid, cell.c);
					double tmp = point::distance(cell.c, nebr.c, p_mid);
					double S = (cell.V + nebr.V) / (2.0 * tmp);
					fmap_frac[{cell.num, nebr.num}] = Face{ S };
				}
			}
	for (auto& cell : cells_poro)
		for (size_t i = 0; i < NEBRS_NUM; i++)
			if (cell.nebrs[i] != -1)
			{
				if (cell.type == PoroType::WELL_LAT && i == 0)
				{
					const auto& nebr = cells_frac[cell.nebrs[i]];
					cell.faces_dist[i] = 0.0;
				}
				else
				{
					const auto& nebr = cells_poro[cell.nebrs[i]];
					const Point p_mid = (cell.c + nebr.c) / 2.0;
					cell.faces_dist[i] = point::distance(cell.c, p_mid, cell.c);
					double S = (cell.V + nebr.V) / (2.0 * point::distance(cell.c, nebr.c, p_mid));
					fmap_poro[{cell.num, nebr.num}] = Face{ S };
				}
			}*/
}
void AcidRecFrac::setPerforated()
{
	height_perf = 0.0;
	vector<pair<int, int> >::iterator it;
	for (const auto& cell : cells_frac)
	{
		if (cell.type == FracType::FRAC_IN)
		{
			Qcell[cell.num] = 0.0;
			height_perf += cell.y * cell.z;
		}
	}
};
void AcidRecFrac::setPeriod(int period)
{
	leftBoundIsRate = LeftBoundIsRate[period];
	if (leftBoundIsRate)
	{
		Q_sum = rate[period];

		if (period == 0 || rate[period - 1] < EQUALITY_TOLERANCE) {
			std::map<int, double>::iterator it;
			for (it = Qcell.begin(); it != Qcell.end(); ++it)
				it->second = Q_sum;// *cells_frac[it->first].hz * cells_frac[it->first].hy / height_perf;
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
void AcidRecFrac::setInitialState() 
{
	const auto& props = props_sk.back();
	for (auto& cell : cells_frac)
	{
		cell.u_prev.p = cell.u_iter.p = cell.u_next.p = props_frac.p_init - grav * props_w.dens_stc * cell.z;
		cell.u_prev.c = cell.u_iter.c = cell.u_next.c = props_frac.c_init;
	}
	for (auto& cell : cells_poro)
	{
		cell.u_prev.m = cell.u_iter.m = cell.u_next.m = props.m_init;
		cell.u_prev.p = cell.u_iter.p = cell.u_next.p = props.p_init - grav * props_w.dens_stc * cell.z;
		cell.u_prev.sw = cell.u_iter.sw = cell.u_next.sw = props.sw_init;
		cell.u_prev.xw = cell.u_iter.xw = cell.u_next.xw = props.xw_init;
		cell.u_prev.xa = cell.u_iter.xa = cell.u_next.xa = props.xa_init;
		cell.u_prev.xs = cell.u_iter.xs = cell.u_next.xs = 0.0;
	}

	x_frac = new FracTapeVariable[cellsNum_frac];
	x_poro = new PoroTapeVariable[cellsNum_poro];
	h = new adouble[ var_frac_size * cellsNum_frac + var_poro_size * cellsNum_poro ];
}
void AcidRecFrac::calculateTrans()
{
	/*double sum, k, k0, x0;
	double width = 0.0, cur_width;
	for (auto& grid : poro_grids)
	{
		k0 = grid.props_sk->perm;
		cur_width = 0.0;
		for (const auto& cell : grid.cells)
		{
			k = grid.props_sk->getPermCoseni(cell.u_next.m, cell.u_next.p).value();
			if (fabs(k - k0) / k0 < 10.0)
				break;
			else
				cur_width += cell.hx;
		}

		if (cur_width > width)
			width = cur_width;
	}
	std::cout << "width = " << width * R_dim << std::endl;
	for (auto& grid : poro_grids)
	{
		grid.width = width;
		k0 = grid.props_sk->perm;
		x0 = grid.cells[0].x;
		sum = 0.0;
		for (const auto& cell : grid.cells)
		{
			k = grid.props_sk->getPermCoseni(cell.u_next.m, cell.u_next.p).value();
			if (cell.x - x0 > width)
				break;
			else
				sum += k * cell.hx;
		}

		grid.trans = (sum > 0.0) ? sum / (k0 * width) : 1.0;
	}*/
}
double AcidRecFrac::getRate(int cur) const
{
	/*const FracCell& cell = cells_frac[cur];
	assert(cell.type == FracType::FRAC_IN);
	const FracCell& beta = cells_frac[cur + (cellsNum_z + 2) * (cellsNum_y + 1)];
	const FracVariable& next = cell.u_next;
	const FracVariable& nebr = beta.u_next;

	double alpha = -props_frac.w2 * props_frac.w2 / props_w.visc * (1.0 - (cell.y / props_frac.w2) * (cell.y / props_frac.w2));
	return alpha * cell.hy * cell.hz * (next.p - nebr.p) / (cell.hx + beta.hx);*/
	return 0;
};

PoroTapeVariable AcidRecFrac::solvePoro(const PoroCell& cell)
{
	if (cell.type == PoroType::MIDDLE)
		return solvePoroMid(cell);
	else if (cell.type == PoroType::WELL_LAT)
		return solvePoroLeft(cell);
	else if (cell.type == PoroType::RIGHT)
		return solvePoroRight(cell);
	else
		return solvePoroBorder(cell);
}
PoroTapeVariable AcidRecFrac::solvePoroMid(const PoroCell& cell)
{
	assert(cell.type == PoroType::MIDDLE);
	const auto& props = props_sk.back();
	const auto& next = x_poro[cell.num];
	const auto& prev = cell.u_prev;
	adouble rate = getReactionRate(next, props);

	PoroTapeVariable res;
	/*adouble m = props.getPoro(next.m, next.p);
	adouble m_prev = props.getPoro(prev.m, prev.p);
	res.m = (1.0 - m) * props.getDensity(next.p) - (1.0 - m_prev) * props.getDensity(prev.p) -
		ht * reac.indices[REACTS::CALCITE] * reac.comps[REACTS::CALCITE].mol_weight * rate;
	res.p = m * next.sw * props_w.getDensity(next.p, next.xa, next.xw, next.xs) -
		m_prev * prev.sw * props_w.getDensity(prev.p, prev.xa, prev.xw, prev.xs) -
		ht * (reac.indices[REACTS::ACID] * reac.comps[REACTS::ACID].mol_weight +
			reac.indices[REACTS::WATER] * reac.comps[REACTS::WATER].mol_weight +
			reac.indices[REACTS::SALT] * reac.comps[REACTS::SALT].mol_weight) * rate;
	res.sw = m * (1.0 - next.sw) * props_o.getDensity(next.p) -
		m_prev * (1.0 - prev.sw) * props_o.getDensity(prev.p);
	res.xw = m * next.sw * props_w.getDensity(next.p, next.xa, next.xw, next.xs) * next.xw -
		m_prev * prev.sw * props_w.getDensity(prev.p, prev.xa, prev.xw, prev.xw) * prev.xw -
		ht * reac.indices[REACTS::WATER] * reac.comps[REACTS::WATER].mol_weight * rate;
	res.xa = m * next.sw * props_w.getDensity(next.p, next.xa, next.xw, next.xs) * next.xa -
		m_prev * prev.sw * props_w.getDensity(prev.p, prev.xa, prev.xw, prev.xs) * prev.xa -
		ht * reac.indices[REACTS::ACID] * reac.comps[REACTS::ACID].mol_weight * rate;
	res.xs = m * next.sw * props_w.getDensity(next.p, next.xa, next.xw, next.xs) * next.xs -
		m_prev * prev.sw * props_w.getDensity(prev.p, prev.xa, prev.xw, prev.xs) * prev.xs -
		ht * reac.indices[REACTS::SALT] * reac.comps[REACTS::SALT].mol_weight * rate;

	size_t upwd_idx;

	int neighbor[NEBRS_NUM];
	getPoroNeighborIdx(cell.num, neighbor);
	for (size_t i = 0; i < NEBRS_NUM; i++)
	{
		const PoroCell& beta = cells_poro[neighbor[i]];
		const auto& nebr = x_poro[beta.num];
		upwd_idx = getUpwindIdx(cell, beta);
		const auto& upwd = x_poro[upwd_idx];

		adouble dens_w = getAverage(props_w.getDensity(next.p, next.xa, next.xw, next.xs), cell.faces_dist[i],
										props_w.getDensity(nebr.p, nebr.xa, nebr.xw, nebr.xs), beta.faces_dist[cell.nebrs_idx[i]]);
		adouble dens_o = getAverage(props_o.getDensity(next.p), cell.faces_dist[i],
										props_o.getDensity(nebr.p), beta.faces_dist[cell.nebrs_idx[i]]);
		adouble buf_w = ht / cell.V * getPoroTrans(cell, next, i, beta, nebr) * (next.p - nebr.p) *
			dens_w * props_w.getKr(upwd.sw, upwd.m, &props) / props_w.getViscosity(upwd.p, upwd.xa, upwd.xw, upwd.xs);
		adouble buf_o = ht / cell.V * getPoroTrans(cell, next, i, beta, nebr) * (next.p - nebr.p) *
			dens_o * props_o.getKr(upwd.sw, upwd.m, &props) / props_o.getViscosity(upwd.p);

		res.p += buf_w;
		res.sw += buf_o;
		res.xw += buf_w * upwd.xw;
		res.xa += buf_w * upwd.xa;
		res.xs += buf_w * upwd.xs;
	}

	//res.p *= sqrt(cell.V);
	//res.sw *= sqrt(cell.V);
	//res.xw *= sqrt(cell.V);
	//res.xa *= sqrt(cell.V);
	//res.xs *= sqrt(cell.V);

	if (cell.num % (cellsNum_y_poro + 2) > cellsNum_y_poro / 3)
	{
		res.p *= sqrt(cell.V);
		res.sw *= sqrt(cell.V);
		res.xw *= sqrt(cell.V);
		res.xa *= sqrt(cell.V);
		res.xs *= sqrt(cell.V);
	}*/

	return res;
}
PoroTapeVariable AcidRecFrac::solvePoroLeft(const PoroCell& cell)
{
	assert(cell.type == PoroType::WELL_LAT);
	const auto& props = props_sk.back();

	//const auto& beta1 = cells_frac[cell.nebrs[0]];
	//const auto& beta2 = cells_frac[beta1.num - 1];
	//const auto& beta3 = cells_frac[beta1.num - 2];
	//const auto& nebr1 = x_frac[beta1.num];
	//const auto& nebr2 = x_frac[beta2.num];
	//const auto& nebr3 = x_frac[beta3.num];

	const auto& next = x_poro[cell.num];
	const auto& prev = cell.u_prev;
	
	PoroTapeVariable res;
	/*res.m = (1.0 - props.getPoro(next.m, next.p)) * props.getDensity(next.p) - (1.0 - prev.m) * props.getDensity(prev.p) -
		ht * reac.indices[REACTS::CALCITE] * reac.comps[REACTS::CALCITE].mol_weight * getReactionRate(next, props);
	res.p = ((next.p - nebr1.p) - (nebr1.p - nebr2.p) * 
		(cell.faces_dist[0] + beta1.faces_dist[cell.nebrs_idx[0]]) / (beta1.faces_dist[0] + beta2.faces_dist[beta1.nebrs_idx[0]])) / P_dim;
	//res.p = (next.p - getQuadAppr({ nebr1.p, nebr2.p, nebr3.p }, { beta1.y, beta2.y, beta3.y }, cell.x)) / P_dim;
	res.sw = (next.sw - (1.0 - props.s_oc)) / P_dim;
	res.xw = ((next.xw - (1.0 - nebr1.c)) - (nebr2.c - nebr1.c) * 
		(cell.faces_dist[0] + beta1.faces_dist[cell.nebrs_idx[0]]) / (beta1.faces_dist[0] + beta2.faces_dist[beta1.nebrs_idx[0]])) / P_dim;
	res.xa = ((next.xa - nebr1.c) - (nebr1.c - nebr2.c) * 
		(cell.faces_dist[0] + beta1.faces_dist[cell.nebrs_idx[0]]) / (beta1.faces_dist[0] + beta2.faces_dist[beta1.nebrs_idx[0]])) / P_dim;
	res.xs = next.xs / P_dim;*/
	return res;
}
PoroTapeVariable AcidRecFrac::solvePoroRight(const PoroCell& cell)
{
	assert(cell.type == PoroType::RIGHT);
	//const auto& beta = cells_poro[cell.nebrs[0]];
	const auto& props = props_sk.back();

	const auto& next = x_poro[cell.num];
	//const auto& nebr = x_poro[beta.num];

	adouble rightIsPres = rightBoundIsPres;
	PoroTapeVariable res;
	/*res.m = (props.getPoro(next.m, next.p) - props.getPoro(nebr.m, nebr.p)) / P_dim;
	condassign(res.p, rightIsPres, (next.p - props.p_out + grav * props_w.dens_stc * cell.c.z) / P_dim, (next.p - nebr.p) / P_dim);
	res.sw = (next.sw - nebr.sw) / P_dim;
	res.xw = (next.xw - nebr.xw) / P_dim;
	res.xa = (next.xa - nebr.xa) / P_dim;
	res.xs = (next.xs - nebr.xs) / P_dim;*/
	return res;
}
PoroTapeVariable AcidRecFrac::solvePoroBorder(const PoroCell& cell)
{
	assert(cell.type == PoroType::SIDE_LEFT || cell.type == PoroType::SIDE_RIGHT || 
			cell.type == PoroType::TOP || cell.type == PoroType::BOTTOM);
	//const auto& beta = cells_poro[cell.nebrs[0]];
	const auto& props = props_sk.back();

	const auto& next = x_poro[cell.num];
	//const auto& nebr = x_poro[beta.num];

	adouble rightIsPres = rightBoundIsPres;
	PoroTapeVariable res;
	/*res.m = (props.getPoro(next.m, next.p) - props.getPoro(nebr.m, nebr.p)) / P_dim;
	res.p = (next.p - nebr.p + grav * props_w.dens_stc * (cell.c.z - beta.c.z)) / P_dim;
	res.sw = (next.sw - nebr.sw) / P_dim;
	res.xw = (next.xw - nebr.xw) / P_dim;
	res.xa = (next.xa - nebr.xa) / P_dim;
	res.xs = (next.xs - nebr.xs) / P_dim;*/
	return res;
}
FracTapeVariable AcidRecFrac::solveFrac(const FracCell& cell)
{
	if (cell.type == FracType::FRAC_MID || cell.type == FracType::FRAC_OUT)
		return solveFracMid(cell);
	else if (cell.type == FracType::FRAC_BORDER)
		return solveFracBorder(cell);
	else if (cell.type == FracType::FRAC_IN)
		return solveFracIn(cell);
}
FracTapeVariable AcidRecFrac::solveFracIn(const FracCell& cell)
{
	assert(cell.type == FracType::FRAC_IN);
	//const auto& beta = cells_frac[cell.nebrs[0]];
	const auto& next = x_frac[cell.num];
	//const auto& nebr = x_frac[beta.num];
	const adouble leftIsRate = leftBoundIsRate;
	FracTapeVariable res;
	/*const auto& out_cell = cells_poro[getFirstMuPoro(cell.num)];
	condassign(res.p, leftIsRate, ((next.p - beta.u_prev.p) + 
		2.0 * Qcell[cell.num] / (1.0 - (sinh(cell.c.mu) / sinh(out_cell.c.mu)) * (sinh(cell.c.mu) / sinh(out_cell.c.mu))) / props_frac.w2_avg / props_frac.w2_avg * props_w.visc
						/ fmap_frac.at({ cell.num, beta.num }).S * (cell.faces_dist[0] + beta.faces_dist[cell.nebrs_idx[0]])) / P_dim,
						(next.p - Pwf + grav * props_w.dens_stc * cell.c.z) / P_dim);
	condassign(res.c, leftIsRate, (next.c - nebr.c) / P_dim, (next.c - c) / P_dim);
	res.c = (next.c - c) / P_dim;*/
	return res;
}
FracTapeVariable AcidRecFrac::solveFracBorder(const FracCell& cell)
{
	assert(cell.type == FracType::FRAC_BORDER);
	const auto& next = x_frac[cell.num];
	//const auto& beta = cells_frac[cell.nebrs[0]];
	//const auto& nebr = x_frac[beta.num];

	FracTapeVariable res;

	//res.p = ((next.p - nebr.p) + (cell.c.z - beta.c.z) * grav * props_w.dens_stc) / P_dim;
	//res.c = (next.c - nebr.c) / P_dim;

	return res;
}
FracTapeVariable AcidRecFrac::solveFracMid(const FracCell& cell)
{
	assert(cell.type == FracType::FRAC_MID || FracType::FRAC_OUT);
	const auto& next = x_frac[cell.num];
	const auto& out_cell = cells_poro[getFirstMuPoro(cell.num)];

	FracTapeVariable res;
	/*res.p = 0.0;	
	res.c = (next.c - cell.u_prev.c) * cell.V;

	const double rat_cell = sinh(cell.c.mu) / sinh(out_cell.c.mu);
	const double alpha = -props_frac.w2_avg * props_frac.w2_avg / 2.0 / props_w.visc * (1.0 - rat_cell * rat_cell);

	adouble tmp, isFromHere;
	for (int i = 0; i < 4; i++)
	{
		if (cell.type == FracType::FRAC_OUT && i == 1)
		{ 
			const auto& beta = cells_poro[cell.nebrs[i]];
			assert(beta.type == PoroType::WELL_LAT);
			const auto& nebr = x_poro[beta.num];
			const double& s = fmap_inter.at({ cell.num, beta.num }).S;
			const auto pt = (cell.c + beta.c) / 2.0;
			const adouble vel = getVelocity(cell, beta);
			//const double vel_val = vel.value();
			//assert(vel_val * ((double)i-1.5) < 0.0);
			res.p += pow(-1, (double)i) * s * vel * cell.V;
			isFromHere = (cell.u_next.p < beta.u_next.p) ? false : true;
			condassign(tmp, isFromHere, x_frac[cell.num].c, x_poro[beta.num].xa);
			res.c -= pow(-1, (double)i) * ht * s * tmp * vel;
			adouble diff = props_w.D_e * (next.c - nebr.xa) / (cell.faces_dist[i] + beta.faces_dist[cell.nebrs_idx[i]]);
			res.c += ht * s * diff;
		}
		else
		{
			const auto& beta = cells_frac[cell.nebrs[i]];
			const auto& nebr = x_frac[beta.num];
			const double& s = fmap_frac.at({ cell.num, beta.num }).S;
			const auto pt = (cell.c + beta.c) / 2.0;
			const double rat_cur = sinh(pt.mu) / sinh(out_cell.c.mu);
			const adouble vel = getVelocity(cell, beta);
			const double vel_val = vel.value();
			//const double vel_val = vel.value();
			//assert(vel_val * ((double)i-1.5) < 0.0);
			res.p += pow(-1, (double)i) * s * vel * cell.V;
			isFromHere = (cell.u_next.p < beta.u_next.p) ? false : true;
			condassign(tmp, isFromHere, x_frac[cell.num].c, x_frac[beta.num].c);
			res.c -= pow(-1, (double)i) * ht * s * tmp * vel;
			if (i < 2)
			{
				adouble diff = props_w.D_e * (next.c - nebr.c) / (cell.faces_dist[i] + beta.faces_dist[cell.nebrs_idx[i]]);	
				res.c += ht * s * diff;
			}
		}
	}

	const auto& beta_z_minus = cells_frac[cell.nebrs[4]];
	const auto& beta_z_plus = cells_frac[cell.nebrs[5]];
	const auto& nebr_z_minus = x_frac[beta_z_minus.num];
	const auto& nebr_z_plus = x_frac[beta_z_plus.num];
	adouble vz_minus = alpha * ((next.p - nebr_z_minus.p) / (cell.h.z / 2.0 + beta_z_minus.h.z / 2.0) - grav * props_w.dens_stc);
	adouble vz_plus = alpha * ((next.p - nebr_z_plus.p) / (cell.h.z / 2.0 + beta_z_plus.h.z / 2.0) + grav * props_w.dens_stc);
	const double& sz_minus = fmap_frac.at({ cell.num, beta_z_minus.num }).S;
	const double& sz_plus = fmap_frac.at({ cell.num, beta_z_plus.num }).S;

	res.p += (sz_minus * vz_minus + sz_plus * vz_plus) * cell.V;
	
	isFromHere = (cell.u_next.p < beta_z_minus.u_next.p) ? false : true;
	condassign(tmp, isFromHere, x_frac[cell.num].c, x_frac[beta_z_minus.num].c);
	res.c -= ht * tmp * sz_minus * vz_minus;
	
	isFromHere = (cell.u_next.p < beta_z_plus.u_next.p) ? false : true;
	condassign(tmp, isFromHere, x_frac[cell.num].c, x_frac[beta_z_plus.num].c);
	res.c -= ht * tmp * sz_plus * vz_plus;*/
	
	return res;
}