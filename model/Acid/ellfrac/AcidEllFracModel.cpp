#include "model/Acid/ellfrac/AcidEllFracModel.hpp"
#include <assert.h>

using namespace acidellfrac;

double acidellfrac::Component::R = 8.3144598;
double acidellfrac::Component::p_std = 101325.0;

double Point::a;

AcidEllFrac::AcidEllFrac()
{
	grav = 9.8;
	Volume_frac = Volume_poro = 0.0;
}
AcidEllFrac::~AcidEllFrac()
{
	delete snapshotter;
	//delete[] x_frac, x_poro, h;
}
void AcidEllFrac::setProps(Properties& props)
{
	props_frac = props.props_frac;
	props_sk = props.props_sk;
	skeletonsNum = props.props_sk.size();

	rightBoundIsPres = props.rightBoundIsPres;

	// Setting grid properties
	cellsNum_mu_frac = props.cellsNum_mu_frac;
	cellsNum_mu_poro = props.cellsNum_mu_poro;
	cellsNum_nu = props.cellsNum_x;
	cellsNum_z = props.cellsNum_z;

	cellsNum_frac = (cellsNum_mu_frac + 1) * (cellsNum_nu + 2) * (cellsNum_z + 2);
	cellsNum_poro = (cellsNum_mu_poro + 2) * (cellsNum_nu + 2) * (cellsNum_z + 2);
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

	PoroCell::Point::a = FracCell::Point::a = props_frac.l2;

	props_o.b = setDataset(props.B_oil, P_dim / BAR_TO_PA, 1.0);
	props_o.Rs = setDataset(props.Rs, P_dim / BAR_TO_PA, 1.0);
	props_g.rho = setDataset(props.rho_co2, P_dim / BAR_TO_PA, (P_dim * t_dim * t_dim / R_dim / R_dim));
}
void AcidEllFrac::makeDimLess()
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
void AcidEllFrac::buildFracGrid()
{
	int counter = 0;
	cells_frac.reserve(cellsNum_frac);
	Volume_frac = 0.0;

	const double mu_w = asinh(props_frac.w2 / props_frac.l2);
	
	double hmu = 0.0, hnu = 0.0, hz = 0.0;
	double cmu = 0.0, cnu = 0.0, cz = 0.0;
	FracType cur_type;

	// Left border
	cur_type = FracType::FRAC_IN;
	hnu = 0.0;	 cnu = M_PI_2;
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

		for (int i = 0; i < cellsNum_mu_frac + 1; i++)
		{
			if (i == 0)
				hmu = cmu = 0.0;
			else
			{
				hmu = mu_w / (double)cellsNum_mu_frac;
				cmu = ((double)i - 0.5) * hmu;
			}

			cells_frac.push_back( FracCell(counter++, cmu, cnu, cz, hmu, hnu, hz, cur_type) );
			auto& cell = cells_frac.back();
			cell.nebrs[0] = cell.num + (cellsNum_mu_frac + 1) * (cellsNum_z + 2);
		}
	}
	// Middle border
	hnu = M_PI_2 / (cellsNum_nu + 1);	 cnu = M_PI_2 - hnu / 2;
	for (int j = 0; j < cellsNum_nu; j++)
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

			for (int i = 0; i < cellsNum_mu_frac + 1; i++)
			{
				if (i == 0)
					hmu = cmu = 0.0;
				else
				{
					hmu = mu_w / (double)cellsNum_mu_frac;
					cmu = ((double)i - 0.5) * hmu;
				}

				if (i == cellsNum_mu_frac)
					cur_type = FracType::FRAC_OUT;
				else if (i == 0)
					cur_type = FracType::FRAC_BORDER;
				else
					cur_type = FracType::FRAC_MID;

				if(k == 0 || k == cellsNum_z + 1)
					cur_type = FracType::FRAC_BORDER;

				cells_frac.push_back(FracCell(counter++, cmu, cnu, cz, hmu, hnu, hz, cur_type));
				Volume_frac += cells_frac.back().V;
				auto& cell = cells_frac.back();
				if (cell.type != FracType::FRAC_BORDER)
				{
					cell.nebrs[0] = cell.num - 1;
					if(cell.type == FracType::FRAC_MID)
						cell.nebrs[1] = cell.num + 1;
					else if (cell.type == FracType::FRAC_OUT)
						cell.nebrs[1] = j * (cellsNum_mu_poro + 2) * (cellsNum_z + 2) + k * (cellsNum_mu_poro + 2);
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
				}
			}
		}
		cnu -= hnu;
	}
	// Right cells
	cur_type = FracType::FRAC_BORDER;
	cnu = hnu;	 hnu = 0.0;
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

		for (int i = 0; i < cellsNum_mu_frac + 1; i++)
		{
			if (i == 0)
				hmu = cmu = 0.0;
			else
			{
				hmu = mu_w / (double)cellsNum_mu_frac;
				cmu = ((double)i - 0.5) * hmu;
			}

			auto& cell = cells_frac.back();
			cells_frac.push_back(FracCell(counter++, cmu, cnu, cz, hmu, hnu, hz, cur_type));
			cell.nebrs[0] = cell.num - (cellsNum_mu_frac + 1) * (cellsNum_z + 2);
		}
	}
}
void AcidEllFrac::buildPoroGrid()
{
	int counter = 0;
	cells_poro.reserve(cellsNum_poro);
	Volume_poro = 0.0;

	const double mu_w = asinh(props_frac.w2 / props_frac.l2);
	const double mu_e = asinh(re / props_frac.l2);

	double r_prev = mu_w;
	double logMax = log(mu_e / mu_w);
	double logStep = logMax / (double)cellsNum_mu_poro;

	double hmu = 0.0, hnu = 0.0, hz = 0.0;
	double cmu = 0.0, cnu = 0.0, cz = 0.0;
	PoroType cur_type;

	// Left border
	cur_type = PoroType::SIDE_LEFT;
	hnu = 0.0;	 cnu = M_PI_2;
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

		for (int i = 0; i < cellsNum_mu_poro + 2; i++)
		{
			if (i == 0)
			{
				hmu = 0.0;
				r_prev = cmu = mu_w;
			}
			else if (i == cellsNum_mu_poro + 1)
			{
				hmu = 0.0;
				cmu = mu_e;
			}
			else
			{
				//hmu = (mu_e - mu_w) / (double)cellsNum_mu_poro;
				//cmu = mu_w + ((double)i - 0.5) * hmu;
				cmu = r_prev * (exp(logStep) + 1.0) / 2.0;
				hmu = r_prev * (exp(logStep) - 1.0);
				r_prev *= exp(logStep);
			}

			cells_poro.push_back(PoroCell(counter++, cmu, cnu, cz, hmu, hnu, hz, cur_type));
			auto& cell = cells_poro.back();
			cell.nebrs[0] = cell.num + (cellsNum_mu_poro + 2) * (cellsNum_z + 2);
		}
	}
	// Middle border
	hnu = M_PI_2 / (cellsNum_nu + 1);	 cnu = M_PI_2 - hnu / 2;
	for (int j = 0; j < cellsNum_nu; j++)
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

			for (int i = 0; i < cellsNum_mu_poro + 2; i++)
			{
				if (i == 0)
				{
					hmu = 0.0;
					r_prev = cmu = mu_w;
				}
				else if (i == cellsNum_mu_poro + 1)
				{
					hmu = 0.0;
					cmu = mu_e;
				}
				else
				{
					//hmu = (mu_e - mu_w) / (double)cellsNum_mu_poro;
					//cmu = mu_w + ((double)i - 0.5) * hmu;
					cmu = r_prev * (exp(logStep) + 1.0) / 2.0;
					hmu = r_prev * (exp(logStep) - 1.0);
					r_prev *= exp(logStep);
				}

				if (i == cellsNum_mu_poro + 1)
					cur_type = PoroType::RIGHT;
				else if (i == 0)
					cur_type = PoroType::WELL_LAT;
				else
					cur_type = PoroType::MIDDLE;

				if (k == 0)
					cur_type = PoroType::BOTTOM;
				else if (k == cellsNum_z + 1)
					cur_type = PoroType::TOP;

				cells_poro.push_back(PoroCell(counter++, cmu, cnu, cz, hmu, hnu, hz, cur_type));
				Volume_poro += cells_poro.back().V;
				auto& cell = cells_poro.back();
				if (cell.type == PoroType::MIDDLE)
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
					cell.nebrs[0] = j * (cellsNum_mu_frac + 1) * (cellsNum_z + 2) + k * (cellsNum_mu_frac + 1) + cellsNum_mu_frac;
					cell.nebrs[1] = cell.num + 1;
				}
				else if (cell.type == PoroType::BOTTOM)
					cell.nebrs[0] = cell.num + (cellsNum_mu_poro + 2);
				else if (cell.type == PoroType::TOP)
					cell.nebrs[0] = cell.num - (cellsNum_mu_poro + 2);
				else if (cell.type == PoroType::RIGHT)
					cell.nebrs[0] = cell.num - 1;
			}
		}
		cnu -= hnu;
	}
	// Right cells
	cur_type = PoroType::SIDE_RIGHT;
	cnu = hnu;	 hnu = 0.0;
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

		for (int i = 0; i < cellsNum_mu_poro + 2; i++)
		{
			if (i == 0)
			{
				hmu = 0.0;
				r_prev = cmu = mu_w;
			}
			else if (i == cellsNum_mu_poro + 1)
			{
				hmu = 0.0;
				cmu = mu_e;
			}
			else
			{
				//hmu = (mu_e - mu_w) / (double)cellsNum_mu_poro;
				//cmu = mu_w + ((double)i - 0.5) * hmu;
				cmu = r_prev * (exp(logStep) + 1.0) / 2.0;
				hmu = r_prev * (exp(logStep) - 1.0);
				r_prev *= exp(logStep);
			}

			cells_poro.push_back(PoroCell(counter++, cmu, cnu, cz, hmu, hnu, hz, cur_type));
			auto& cell = cells_poro.back();
			cell.nebrs[0] = cell.num - (cellsNum_mu_poro + 2) * (cellsNum_z + 2);
		}
	}

}
void AcidEllFrac::processGeometry()
{
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
					cell.faces_dist[i] = (cell.h * (nebr.c - cell.c) / (nebr.c - cell.c).dist() / 2.0).norm();
					const Point p_mid = cell.c + cell.h * (nebr.c - cell.c) / (nebr.c - cell.c).dist() / 2.0;
					double S = (cell.V + nebr.V) / (2.0 * (nebr.c - cell.c).norm());
					fmap_inter[{cell.num, nebr.num}] = Face{ p_mid, S };
				}
				else
				{
					const auto& nebr = cells_frac[cell.nebrs[i]];
					cell.faces_dist[i] = (cell.h * (nebr.c - cell.c) / (nebr.c - cell.c).dist() / 2.0).norm();
					const Point p_mid = cell.c + cell.h * (nebr.c - cell.c) / (nebr.c - cell.c).dist() / 2.0;
					double S = (cell.V + nebr.V) / (2.0 * (nebr.c - cell.c).norm());
					fmap_frac[{cell.num, nebr.num}] = Face{ p_mid, S };
				}
			}
	for (auto& cell : cells_poro)
		for (size_t i = 0; i < NEBRS_NUM; i++)
			if (cell.nebrs[i] != -1)
			{
				if (cell.type == PoroType::WELL_LAT && i == 0)
				{
					const auto& nebr = cells_poro[cell.nebrs[i]];
					cell.faces_dist[i] = (cell.h * (nebr.c - cell.c) / (nebr.c - cell.c).dist() / 2.0).norm();
					const Point p_mid = cell.c + cell.h * (nebr.c - cell.c) / (nebr.c - cell.c).dist() / 2.0;
					double S = (cell.V + nebr.V) / (2.0 * (nebr.c - cell.c).norm());
					fmap_inter[{cell.num, nebr.num}] = Face{ p_mid, S };
				}
				else
				{
					const auto& nebr = cells_poro[cell.nebrs[i]];
					cell.faces_dist[i] = (cell.h * (nebr.c - cell.c) / (nebr.c - cell.c).dist() / 2.0).norm();
					const Point p_mid = cell.c + cell.h * (nebr.c - cell.c) / (nebr.c - cell.c).dist() / 2.0;
					double S = (cell.V + nebr.V) / (2.0 * (nebr.c - cell.c).norm());
					fmap_poro[{cell.num, nebr.num}] = Face{ p_mid, S };
				}
			}
}
void AcidEllFrac::setPerforated()
{
	height_perf = 0.0;
	vector<pair<int, int> >::iterator it;
	for (const auto& cell : cells_frac)
	{
		if (cell.type == FracType::FRAC_IN)
		{
			Qcell[cell.num] = 0.0;
			height_perf += cell.c.getH() * cell.h.mu * cell.h.z;
		}
	}
};
void AcidEllFrac::setPeriod(int period)
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
void AcidEllFrac::setInitialState() 
{
	const auto& props = props_sk.back();
	for (auto& cell : cells_frac)
	{
		cell.u_prev.p = cell.u_iter.p = cell.u_next.p = props_frac.p_init - grav * props_w.dens_stc * cell.c.z;
		cell.u_prev.c = cell.u_iter.c = cell.u_next.c = props_frac.c_init;
	}
	for (auto& cell : cells_poro)
	{
		cell.u_prev.m = cell.u_iter.m = cell.u_next.m = props.m_init;
		cell.u_prev.p = cell.u_iter.p = cell.u_next.p = props.p_init - grav * props_w.dens_stc * cell.c.z;
		cell.u_prev.sw = cell.u_iter.sw = cell.u_next.sw = props.sw_init;
		cell.u_prev.xw = cell.u_iter.xw = cell.u_next.xw = props.xw_init;
		cell.u_prev.xa = cell.u_iter.xa = cell.u_next.xa = props.xa_init;
		cell.u_prev.xs = cell.u_iter.xs = cell.u_next.xs = 0.0;
	}

	x_frac = new FracTapeVariable[cellsNum_frac];
	x_poro = new PoroTapeVariable[cellsNum_poro];
	h = new adouble[ var_frac_size * cellsNum_frac + var_poro_size * cellsNum_poro ];
}
void AcidEllFrac::calculateTrans()
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
double AcidEllFrac::getRate(int cur) const
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

PoroTapeVariable AcidEllFrac::solvePoro(const PoroCell& cell)
{
	if (cell.type == PoroType::MIDDLE)
		return solvePoroMid(cell);
	else if (cell.type == PoroType::WELL_LAT)
		return solvePoroLeft(cell);
	else if (cell.type == PoroType::RIGHT)
		return solvePoroRight(cell);
}
PoroTapeVariable AcidEllFrac::solvePoroMid(const PoroCell& cell)
{
	assert(cell.type == PoroType::MIDDLE);
	const auto& props = props_sk.back();
	const auto& next = x_poro[cell.num];
	const auto& prev = cell.u_prev;
	adouble rate = getReactionRate(next, props);

	PoroTapeVariable res;
	adouble m = props.getPoro(next.m, next.p);
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
	for (size_t i = 0; i < NEBRS_NUM; i++)
	{
		const PoroCell& beta = cells_poro[cell.nebrs[i]];
		const auto& nebr = x_poro[beta.num];
		upwd_idx = getUpwindIdx(cell, beta);
		const auto& upwd = x_poro[upwd_idx];

		adouble dens_w = getPoroAverage(props_w.getDensity(next.p, next.xa, next.xw, next.xs), cell.faces_dist[i],
										props_w.getDensity(nebr.p, nebr.xa, nebr.xw, nebr.xs), beta.faces_dist[cell.nebrs_idx[i]]);
		adouble dens_o = getPoroAverage(props_o.getDensity(next.p), cell.faces_dist[i],
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

	if (cell.num % (cellsNum_mu_poro + 2) > cellsNum_mu_poro / 3)
	{
		res.p *= sqrt(cell.V);
		res.sw *= sqrt(cell.V);
		res.xw *= sqrt(cell.V);
		res.xa *= sqrt(cell.V);
		res.xs *= sqrt(cell.V);
	}

	return res;
}
PoroTapeVariable AcidEllFrac::solvePoroLeft(const PoroCell& cell)
{
	assert(cell.type == PoroType::WELL_LAT);
	const auto& props = props_sk.back();

	const auto& beta1 = cells_frac[cell.nebrs[0]];
	const auto& beta2 = cells_frac[beta1.num - 1];
	//const auto& beta3 = cells_frac[grid.frac_nebr->num - 2];
	const auto& nebr1 = x_frac[beta1.num];
	const auto& nebr2 = x_frac[beta2.num];
	//const auto& nebr3 = x_frac[beta3.num];

	const auto& next = x_poro[cell.num];
	const auto& prev = cell.u_prev;
	
	PoroTapeVariable res;
	res.m = (1.0 - props.getPoro(next.m, next.p)) * props.getDensity(next.p) - (1.0 - prev.m) * props.getDensity(prev.p) -
		ht * reac.indices[REACTS::CALCITE] * reac.comps[REACTS::CALCITE].mol_weight * getReactionRate(next, props);
	res.p = ((next.p - nebr1.p) - (nebr1.p - nebr2.p) * (cell.c.mu - beta1.c.mu) / (beta1.c.mu - beta2.c.mu)) / P_dim;
	//res.p = (next.p - getQuadAppr({ nebr1.p, nebr2.p, nebr3.p }, { beta1.y, beta2.y, beta3.y }, cell.x)) / P_dim;
	res.sw = (next.sw - (1.0 - props.s_oc)) / P_dim;
	res.xw = ((next.xw - (1.0 - nebr1.c)) - (nebr2.c - nebr1.c) * (cell.c.mu - beta1.c.mu) / (beta1.c.mu - beta2.c.mu)) / P_dim;
	res.xa = ((next.xa - nebr1.c) - (nebr1.c - nebr2.c) * (cell.c.mu - beta1.c.mu) / (beta1.c.mu - beta2.c.mu)) / P_dim;
	res.xs = next.xs / P_dim;
	return res;
}
PoroTapeVariable AcidEllFrac::solvePoroRight(const PoroCell& cell)
{
	assert(cell.type == PoroType::RIGHT);
	const auto& beta = cells_poro[cell.nebrs[0]];
	const auto& props = props_sk.back();

	const auto& next = x_poro[cell.num];
	const auto& nebr = x_poro[beta.num];

	adouble rightIsPres = rightBoundIsPres;
	PoroTapeVariable res;
	res.m = (props.getPoro(next.m, next.p) - props.getPoro(nebr.m, nebr.p)) / P_dim;
	condassign(res.p, rightIsPres, (next.p - props.p_out + grav * props_w.dens_stc * cell.c.z) / P_dim, (next.p - nebr.p) / P_dim);
	res.sw = (next.sw - nebr.sw) / P_dim;
	res.xw = (next.xw - nebr.xw) / P_dim;
	res.xa = (next.xa - nebr.xa) / P_dim;
	res.xs = (next.xs - nebr.xs) / P_dim;
	return res;
}
PoroTapeVariable AcidEllFrac::solvePoroBorder(const PoroCell& cell)
{
	assert(cell.type == PoroType::SIDE_LEFT || cell.type == PoroType::SIDE_RIGHT || 
			cell.type == PoroType::TOP || cell.type == PoroType::BOTTOM);
	const auto& beta = cells_poro[cell.nebrs[0]];
	const auto& props = props_sk.back();

	const auto& next = x_poro[cell.num];
	const auto& nebr = x_poro[beta.num];

	adouble rightIsPres = rightBoundIsPres;
	PoroTapeVariable res;
	res.m = (props.getPoro(next.m, next.p) - props.getPoro(nebr.m, nebr.p)) / P_dim;
	res.p = (next.p - nebr.p + grav * props_w.dens_stc * cell.c.z) / P_dim;
	res.sw = (next.sw - nebr.sw) / P_dim;
	res.xw = (next.xw - nebr.xw) / P_dim;
	res.xa = (next.xa - nebr.xa) / P_dim;
	res.xs = (next.xs - nebr.xs) / P_dim;
	return res;
}
FracTapeVariable AcidEllFrac::solveFrac(const FracCell& cell)
{
	if (cell.type == FracType::FRAC_MID)
		return solveFracMid(cell);
	else if (cell.type == FracType::FRAC_BORDER)
		return solveFracBorder(cell);
	else if (cell.type == FracType::FRAC_IN)
		return solveFracIn(cell);
	else if (cell.type == FracType::FRAC_OUT)
		return solveFracOut(cell);
}
FracTapeVariable AcidEllFrac::solveFracIn(const FracCell& cell)
{
	assert(cell.type == FracType::FRAC_IN);
	const auto& beta = cells_frac[cell.nebrs[0]];
	const auto& next = x_frac[cell.num];
	const auto& nebr = x_frac[beta.num];
	const adouble leftIsRate = leftBoundIsRate;
	FracTapeVariable res;
	const auto car_point = cell.c.getCartesian();
	condassign(res.p, leftIsRate, ((next.p - beta.u_prev.p) + 
		2.0 * Qcell[cell.num] / (1.0 - (car_point.y / props_frac.w2) * (car_point.y / props_frac.w2)) / props_frac.w2 / props_frac.w2 * props_w.visc
						/ fmap_poro.at({ cell.num, beta.num }).S * (cell.faces_dist[0] + beta.faces_dist[cell.nebrs_idx[0]])) / P_dim,
						(next.p - Pwf + grav * props_w.dens_stc * cell.c.z) / P_dim);
	condassign(res.c, leftIsRate, (next.c - nebr.c) / P_dim, (next.c - c) / P_dim);
	res.c = (next.c - c) / P_dim;
	return res;
}
FracTapeVariable AcidEllFrac::solveFracBorder(const FracCell& cell)
{
	assert(cell.type == FracType::FRAC_BORDER);
	const auto& next = x_frac[cell.num];
	const auto& beta = cells_frac[cell.nebrs[0]];
	const auto& nebr = x_frac[beta.num];

	FracTapeVariable res;
	//if (cell.h.z == 0.0 && cell.h.x > 0.0 && cell.hy > 0.0)
	//{
		res.p = ((next.p - nebr.p) + (cell.c.z - beta.c.z) * grav * props_w.dens_stc) / P_dim;
		res.c = (next.c - nebr.c) / P_dim;
	/*}
	else
	{
		res.p = (next.p - nebr.p) / P_dim;
		res.c = (next.c - nebr.c) / P_dim;
	}*/

	return res;
}
FracTapeVariable AcidEllFrac::solveFracOut(const FracCell& cell)
{
	assert(cell.type == FracType::FRAC_OUT);
	const auto& next = x_frac[cell.num];
	
	FracTapeVariable res;
/*	int neighbor[6];
	getNeighborIdx(cell.num, neighbor);
	const auto& grid = poro_grids[frac2poro[cell.num]];
	const auto& beta_poro = grid.cells[0];
	const auto& nebr_x_minus = x_frac[neighbor[0]];
	const auto& nebr_x_plus = x_frac[neighbor[1]];
	const auto& nebr_y_minus = x_frac[neighbor[2]];
	const auto& nebr_y_plus = x_poro[grid.start_idx];
	const auto& nebr_z_minus = x_frac[neighbor[4]];
	const auto& nebr_z_plus = x_frac[neighbor[5]];

	const double y_minus = cell.y - cell.hy / 2.0;
	adouble vL = getFlowLeak(cell);
	double alpha = -props_frac.w2 * props_frac.w2 / props_w.visc * (1.0 - (cell.y / props_frac.w2) * (cell.y / props_frac.w2));
	adouble vx_minus = alpha * (next.p - nebr_x_minus.p) / (cell.hx + cells_frac[neighbor[0]].hx);
	adouble vx_plus = alpha * (next.p - nebr_x_plus.p) / (cell.hx + cells_frac[neighbor[1]].hx);
	adouble vz_minus = alpha * ((next.p - nebr_z_minus.p) / (cell.hz + cells_frac[neighbor[4]].hz) - grav * props_w.dens_stc / 2.0);
	adouble vz_plus = alpha * ((next.p - nebr_z_plus.p) / (cell.hz + cells_frac[neighbor[5]].hz) + grav * props_w.dens_stc / 2.0);
	adouble vy_minus = vL * (1.5 * (y_minus / props_frac.w2) - 0.5 * y_minus * y_minus * y_minus / props_frac.w2 / props_frac.w2 / props_frac.w2);
	adouble vy_plus = vL;

	adouble diff_minus = 2.0 * props_w.D_e * (next.c - nebr_y_minus.c) / (cell.hy + cells_frac[neighbor[2]].hy);
	adouble diff_plus = 2.0 * props_w.D_e * (next.c - nebr_y_plus.xa) / (cell.hy + beta_poro.hx);

	const double sx = cell.hy * cell.hz, sy = cell.hx * cell.hz, sz = cell.hx * cell.hy;
	res.p =	sx * (vx_minus + vx_plus) + sz * (vz_minus + vz_plus) - sy * (vy_plus - vy_minus);
	res.p *= cell.V;
	res.c = cell.V * (next.c - cell.u_prev.c) - ht *
		(sx * (x_frac[getUpwindIdx(cell, cells_frac[neighbor[0]])].c * vx_minus +
				x_frac[getUpwindIdx(cell, cells_frac[neighbor[1]])].c * vx_plus) +
		sy * (x_frac[getUpwindIdx(cell, cells_frac[neighbor[2]])].c * vy_minus -
				next.c * vy_plus - diff_plus - diff_minus) +
		sz * (x_frac[getUpwindIdx(cell, cells_frac[neighbor[4]])].c * vz_minus +
				x_frac[getUpwindIdx(cell, cells_frac[neighbor[5]])].c * vz_plus));*/
	return res;
}
FracTapeVariable AcidEllFrac::solveFracMid(const FracCell& cell)
{
	assert(cell.type == FracType::FRAC_MID);
	const auto& next = x_frac[cell.num];

	FracTapeVariable res;
/*	int neighbor[6];
	getNeighborIdx(cell.num, neighbor);
	const auto& nebr_x_minus = x_frac[neighbor[0]];
	const auto& nebr_x_plus = x_frac[neighbor[1]];
	const auto& nebr_y_minus = x_frac[neighbor[2]];
	const auto& nebr_y_plus = x_frac[neighbor[3]];
	const auto& nebr_z_minus = x_frac[neighbor[4]];
	const auto& nebr_z_plus = x_frac[neighbor[5]];

	double y_plus = cell.y + cell.hy / 2.0;
	double y_minus = cell.y - cell.hy / 2.0;
	adouble vL = getFlowLeak(cell);
	double alpha = -props_frac.w2 * props_frac.w2 / props_w.visc * (1.0 - (cell.y / props_frac.w2) * (cell.y / props_frac.w2));
	adouble vx_minus = alpha * (next.p - nebr_x_minus.p) / (cell.hx + cells_frac[neighbor[0]].hx);
	adouble vx_plus = alpha * (next.p - nebr_x_plus.p) / (cell.hx + cells_frac[neighbor[1]].hx);
	adouble vz_minus = alpha * ((next.p - nebr_z_minus.p) / (cell.hz + cells_frac[neighbor[4]].hz) - grav * props_w.dens_stc / 2.0);
	adouble vz_plus = alpha * ((next.p - nebr_z_plus.p) / (cell.hz + cells_frac[neighbor[5]].hz) + grav * props_w.dens_stc / 2.0);
	adouble vy_minus = vL * (1.5 * (y_minus / props_frac.w2) - 0.5 * y_minus * y_minus * y_minus / props_frac.w2 / props_frac.w2 / props_frac.w2);
	adouble vy_plus = vL * (1.5 * (y_plus / props_frac.w2) - 0.5 * y_plus * y_plus * y_plus / props_frac.w2 / props_frac.w2 / props_frac.w2);

	adouble diff_minus = 2.0 * props_w.D_e * (next.c - nebr_y_minus.c) / (cell.hy + cells_frac[neighbor[2]].hy);
	adouble diff_plus = 2.0 * props_w.D_e * (next.c - nebr_y_plus.c) / (cell.hy + cells_frac[neighbor[3]].hy);

	const double sx = cell.hy * cell.hz, sy = cell.hx * cell.hz, sz = cell.hx * cell.hy;
	res.p =	((sx * (vx_minus + vx_plus) + sz * (vz_minus + vz_plus) - sy * (vy_plus - vy_minus)) / alpha / sx * (cell.hx + cells_frac[neighbor[0]].hx));
	res.p *= cell.V;
	res.c = ((next.c - cell.u_prev.c) * cell.V - ht *
		(sx * (x_frac[getUpwindIdx(cell, cells_frac[neighbor[0]])].c * vx_minus +
				x_frac[getUpwindIdx(cell, cells_frac[neighbor[1]])].c * vx_plus) +
		sy * (x_frac[getUpwindIdx(cell, cells_frac[neighbor[2]])].c * vy_minus -
				x_frac[getUpwindIdx(cell, cells_frac[neighbor[3]])].c * vy_plus - diff_plus - diff_minus) +
		sz * (x_frac[getUpwindIdx(cell, cells_frac[neighbor[4]])].c * vz_minus +
				x_frac[getUpwindIdx(cell, cells_frac[neighbor[5]])].c * vz_plus)));*/
	return res;
}