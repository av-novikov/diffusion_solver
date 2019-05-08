#include "model/Acid/recfracmov/AcidRecFracMovModel.hpp"
#include <assert.h>
#include <numeric>

using namespace acidrecfracmov;

double acidrecfracmov::Component::R = 8.3144598;
double acidrecfracmov::Component::p_std = 101325.0;

AcidRecFracMov::AcidRecFracMov()
{
	grav = 9.8;
	Volume_frac = Volume_poro = 0.0;
	injected_sol_volume = injected_acid_volume = 0.0;
	max_vel_x = max_vel_y = max_vel_z = 0.0;
}
AcidRecFracMov::~AcidRecFracMov()
{
	delete snapshotter;
	delete[] x_frac, x_poro, h;
}
void AcidRecFracMov::setProps(Properties& props)
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
        props_sk[j].perm_max = MilliDarcyToM2(props_sk[j].perm_max);
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
	max_sol_volume = props.max_sol_volume;

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
void AcidRecFracMov::makeDimLess()
{
	T_dim = props_sk[0].t_init;
	//R_dim = props_frac.l2 / 8.0;
	R_dim = props_frac.l2 / 5.0;
	t_dim = 1.7 * 3600.0;
	//t_dim = 3600.0;
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
        sk.perm_max /= (R_dim * R_dim);
		sk.beta /= (1.0 / P_dim);
		sk.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
		sk.p_init /= P_dim;
		sk.p_out /= P_dim;
		sk.p_ref /= P_dim;
		sk.hx /= R_dim;
		sk.hz /= R_dim;
		sk.t_init /= T_dim;
		sk.height /= R_dim;
	}

	Q_dim = R_dim * R_dim * R_dim / t_dim;
	for (int i = 0; i < periodsNum; i++)
	{
		period[i] /= t_dim;
		rate[i] /= Q_dim;
		pwf[i] /= P_dim;
	}
	max_sol_volume /= R_dim * R_dim * R_dim;

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
	reac.reaction_const /= (pow(R_dim, 3 * reac.alpha - 2.0) / t_dim);
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
void AcidRecFracMov::buildFracGrid()
{
	int counter = 0;
	cells_frac.reserve(cellsNum_frac);
	Volume_frac = 0.0;
    const double w2 = props_frac.w2[0];
	
	double hx = 0.0, hy = 0.0, hz = 0.0;
	double cx = 0.0, cy = 0.0, cz = 0.0;
	FracType cur_type;

	//dist_x = props_frac.l2 / 100.0;
	//double x_prev = dist_x;
	//double x_logMax = log((props_frac.l2 + dist_x) / dist_x);
	//double x_logStep = x_logMax / (double)cellsNum_x;

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
			//hz = props_frac.height / cellsNum_z;
			cz += hz / 2.0;
			hz = props_sk[k - 1].height;
			cz += hz / 2.0;
		}

		for (int i = 0; i < cellsNum_y_frac + 1; i++)
		{
			if (i == 0)
				hy = cy = 0.0;
			else
			{
				hy = w2 / (double)cellsNum_y_frac;
				cy = ((double)i - 0.5) * hy;
			}

			cells_frac.push_back( FracCell(counter++, cx, cy, cz, hx, hy, hz, cur_type) );
			auto& cell = cells_frac.back();
		}
	}
	// Middle
	hx = props_frac.l2 / cellsNum_x;	 cx = hx / 2;
	for (int j = 0; j < cellsNum_x; j++)
	{
		//cx = x_prev * (exp(x_logStep) + 1.0) / 2.0 - dist_x;
		//hx = x_prev * (exp(x_logStep) - 1.0);
		//x_prev *= exp(x_logStep);

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
				//hz = props_frac.height / cellsNum_z;
				cz += hz / 2.0;
				hz = props_sk[k - 1].height;
				cz += hz / 2.0;
			}

			for (int i = 0; i < cellsNum_y_frac + 1; i++)
			{
				if (i == 0)
					hy = cy = 0.0;
				else
				{
					hy = w2 / (double)cellsNum_y_frac;
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
			//hz = props_frac.height / cellsNum_z;
			cz += hz / 2.0;
			hz = props_sk[k - 1].height;
			cz += hz / 2.0;
		}

		for (int i = 0; i < cellsNum_y_frac + 1; i++)
		{
			if (i == 0)
				hy = cy = 0.0;
			else
			{
				hy = w2 / (double)cellsNum_y_frac;
				cy = ((double)i - 0.5) * hy;
			}

			cells_frac.push_back(FracCell(counter++, cx, cy, cz, hx, hy, hz, cur_type));
			auto& cell = cells_frac.back();
			//cell.nebrs[0] = cell.num - (cellsNum_y_frac + 1) * (cellsNum_z + 2);
		}
	}
}
void AcidRecFracMov::buildPoroGrid()
{
	int counter = 0;
	cells_poro.reserve(cellsNum_poro);
	Volume_poro = 0.0;
    const double w2 = props_frac.w2[0];

	// Grid sparcity parameters
	int mult_num = 2 * cellsNum_y_poro / 3;
    double dist = 2.0 * w2;

    double delta2 = 0.95 * dist;
	double r_prev2 = re;
	double logMax2 = log((re - delta2) / (dist - delta2));
	double logStep2 = logMax2 / (double)(cellsNum_y_poro - mult_num);

    double delta1 = 0.95 * w2;
    double r_prev1 = re;
    double logMax1 = log((dist - delta1) / (w2 - delta1));
    double logStep1 = logMax1 / (double)mult_num;

	//double x_prev = dist_x;
	//double x_logMax = log((props_frac.l2 + dist_x) / dist_x);
	//double x_logStep = x_logMax / (double)cellsNum_x;

    double r_prev, delta = delta1, logStep = logStep1;

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
			//hz = props_frac.height / cellsNum_z;
			cz += hz / 2.0;
			hz = props_sk[k - 1].height;
			cz += hz / 2.0;
		}

		for (int i = 0; i < cellsNum_y_poro + 2; i++)
		{
            if (i == mult_num + 1)
            {
                delta = delta2;
                logStep = logStep2;
                r_prev = dist - delta;
            }

			if (i == 0)
			{
                delta = delta1;
                logStep = logStep1;
				hy = 0.0;
				cy = w2;
				r_prev = w2 - delta;
			}
			else if (i == cellsNum_y_poro + 1)
			{
				hy = 0.0;
				cy = re;
			}
			else
			{
				cy = delta + r_prev * (exp(logStep) + 1.0) / 2.0;
				hy = r_prev * (exp(logStep) - 1.0);
				r_prev *= exp(logStep);
			}

			cells_poro.push_back(PoroCell(counter++, cx, cy, cz, hx, hy, hz, cur_type));
			auto& cell = cells_poro.back();
		}
	}
	// Middle
	hx = props_frac.l2 / cellsNum_x;	 cx = hx / 2;
	for (int j = 0; j < cellsNum_x; j++)
	{
		//cx = x_prev * (exp(x_logStep) + 1.0) / 2.0 - dist_x;
		//hx = x_prev * (exp(x_logStep) - 1.0);
		//x_prev *= exp(x_logStep);

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
				//hz = props_frac.height / cellsNum_z;
				cz += hz / 2.0;
				hz = props_sk[k - 1].height;
				cz += hz / 2.0;
			}

			for (int i = 0; i < cellsNum_y_poro + 2; i++)
			{
                if (i == mult_num + 1)
                {
                    delta = delta2;
                    logStep = logStep2;
                    r_prev = dist - delta;
                }

				if (i == 0)
				{
                    delta = delta1;
                    logStep = logStep1;
					hy = 0.0;
					r_prev = w2 - delta;
					cy = w2;
				}
				else if (i == cellsNum_y_poro + 1)
				{
					hy = 0.0;
					cy = re;
				}
				else
				{
					cy = delta + r_prev * (exp(logStep) + 1.0) / 2.0;
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
			//hz = props_frac.height / cellsNum_z;
			cz += hz / 2.0;
			hz = props_sk[k - 1].height;
			cz += hz / 2.0;
		}

		for (int i = 0; i < cellsNum_y_poro + 2; i++)
		{
            if (i == mult_num + 1)
            {
                delta = delta2;
                logStep = logStep2;
                r_prev = dist - delta;
            }

			if (i == 0)
			{
                delta = delta1;
                logStep = logStep1;
				hy = 0.0;
				r_prev = w2 - delta;
				cy = w2;
			}
			else if (i == cellsNum_y_poro + 1)
			{
				hy = 0.0;
				cy = re;
			}
			else
			{
				cy = delta + r_prev * (exp(logStep) + 1.0) / 2.0;
				hy = r_prev * (exp(logStep) - 1.0);
				r_prev *= exp(logStep);
			}

			cells_poro.push_back(PoroCell(counter++, cx, cy, cz, hx, hy, hz, cur_type));
			auto& cell = cells_poro.back();
		}
	}

}
void AcidRecFracMov::processGeometry()
{
}
void AcidRecFracMov::setPerforated()
{
	height_perf = 0.0;
	vector<pair<int, int> >::iterator it;
	for (const auto& cell : cells_frac)
	{
		if (cell.type == FracType::FRAC_IN)
		{
			Qcell[cell.num] = 0.0;
			height_perf += cell.hy * cell.hz;
		}
	}
};
void AcidRecFracMov::setPeriod(int period)
{
	leftBoundIsRate = LeftBoundIsRate[period];
	if (leftBoundIsRate)
	{
		Q_sum = rate[period];
		if (Q_sum == 0.0)
		{
			for(auto& cell : cells_frac)
				cell.u_prev.c = cell.u_iter.c = cell.u_next.c = 0.0;
		}

		if (period == 0 || rate[period - 1] < EQUALITY_TOLERANCE) 
		{
			std::map<int, double>::iterator it;
            for (it = Qcell.begin(); it != Qcell.end(); ++it)
            {
                const double w2 = props_frac.w2[getWidthIdByFracId(it->first)];
                it->second = 3.0 * Q_sum * cells_frac[it->first].hy / 2.0 / w2 *
                    (1.0 - (cells_frac[it->first].y / w2) * (cells_frac[it->first].y / w2));
            }
		}
		else
		{
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
void AcidRecFracMov::setInitialState() 
{
	//const auto& props = props_sk.back();
	for (auto& cell : cells_frac)
	{
		cell.u_prev.p = cell.u_iter.p = cell.u_next.p = props_frac.p_init - grav * props_w.dens_stc * cell.z;
		cell.u_prev.c = cell.u_iter.c = cell.u_next.c = props_frac.c_init;
	}
	for (auto& cell : cells_poro)
	{
		cell.props = &props_sk[getSkeletonId(cell)];
		cell.u_prev.m = cell.u_iter.m = cell.u_next.m = cell.props->m_init;
		cell.u_prev.p = cell.u_iter.p = cell.u_next.p = cell.props->p_init - grav * props_w.dens_stc * cell.z;
		
		if (cell.type == PoroType::WELL_LAT)
		{
			cell.u_prev.sw = cell.u_iter.sw = cell.u_next.sw = 1.0 - cell.props->s_oc;
			cell.u_prev.xa = cell.u_iter.xa = cell.u_next.xa = 0.0;
			cell.u_prev.xw = cell.u_iter.xw = cell.u_next.xw = 1.0;
		}
		else
		{
			cell.u_prev.sw = cell.u_iter.sw = cell.u_next.sw = cell.props->sw_init;
			cell.u_prev.xa = cell.u_iter.xa = cell.u_next.xa = cell.props->xa_init;
			cell.u_prev.xw = cell.u_iter.xw = cell.u_next.xw = cell.props->xw_init;
		}
		cell.u_prev.xs = cell.u_iter.xs = cell.u_next.xs = 0.0;

	}

	x_frac = new FracTapeVariable[cellsNum_frac];
	x_poro = new PoroTapeVariable[cellsNum_poro];
	h = new adouble[ var_frac_size * cellsNum_frac + var_poro_size * cellsNum_poro ];
	trans.resize(cellsNum_x * cellsNum_z, 1.0);
	widths = props_frac.w2;
}
void AcidRecFracMov::calculateTrans()
{
	double sum, k, k0, y0, cur_width, m, m0;
	/*width = 0.0;
	for (const auto& cell : cells_poro)
	{
		if (cell.hx != 0.0 && cell.hz != 0.0)
		{
			k = cell.props->getPermCoseni(cell.u_next.m, cell.u_next.p).value();
			k0 = cell.props->perm;
			if (fabs(k - k0) / k0 < 10.0)
				continue;
			else
				width += cell.hy;
		}
	}
	width /= (cellsNum_z * cellsNum_x);
	std::cout << std::endl << "width = " << width * R_dim << std::endl << std::endl;*/
	int i_ind, k_ind;
	for (int tr_idx = 0; tr_idx < trans.size(); tr_idx++)
	{
		k_ind = int(tr_idx % cellsNum_z) + 1;
		i_ind = int(tr_idx / cellsNum_z) + 1;
		cur_width = sum = 0.0;
		for (int j = 0; j < cellsNum_y_poro + 2; j++)
		{
			const auto& pcell = cells_poro[j + (cellsNum_y_poro + 2) * (k_ind + (cellsNum_z + 2) * i_ind)];
			m = pcell.u_next.m;
			m0 = pcell.props->m_init;
			k = pcell.props->getPermCoseni(pcell.u_next.m, pcell.u_next.p).value();
			k0 = pcell.props->perm;
			//if (fabs(k - k0) / k0 >= 10.0)
			if(fabs(m - m0) > EQUALITY_TOLERANCE)
				cur_width += pcell.hy;
		}
		widths[tr_idx] = cur_width;
	}

    width = widths.sum() / widths.size();

	for (int tr_idx = 0; tr_idx < trans.size(); tr_idx++)
	{
		k_ind = int(tr_idx % cellsNum_z) + 1;
		i_ind = int(tr_idx / cellsNum_z) + 1;
		sum = 0.0;
		for (int j = 0; j < cellsNum_y_poro + 2; j++)
		{
			const auto& pcell = cells_poro[j + (cellsNum_y_poro + 2) * (k_ind + (cellsNum_z + 2) * i_ind)];
			assert(pcell.type == PoroType::WELL_LAT || pcell.type == PoroType::MIDDLE);
			k = pcell.props->getPermCoseni(pcell.u_next.m, pcell.u_next.p).value();
			if (pcell.y - props_frac.w2[tr_idx] > widths[tr_idx])
				break;
			else
				sum += (k / pcell.props->perm) * (pcell.hy / widths[tr_idx]);
		}
		trans[tr_idx] = (sum > 0.0) ? sum : 1.0;
	}
}
void AcidRecFracMov::checkPoroCells()
{
    for (auto& cell : cells_poro)
    {
        if (cell.type == PoroType::MIDDLE)
        {
            const auto& prop = *cell.props;
            if (cell.u_next.m > prop.m_max || prop.getPermCoseni(cell.u_next.m, cell.u_next.p).value() > prop.m_max)
            {
                //cell.type = PoroType::BE_FRAC;
                expandFracByCell(cell);
            }
        }
    }
}
void AcidRecFracMov::expandFracByCell(PoroCell &cell)
{
}
double AcidRecFracMov::getRate(const int idx) const
{
	const FracCell& cell = cells_frac[idx];
    const double w2 = props_frac.w2[getWidthIdByFracId(idx)];
	assert(cell.type == FracType::FRAC_IN);
	assert(cell.hy != 0.0 && cell.hz != 0.0);
	const FracCell& beta = cells_frac[cell.num + (cellsNum_z + 2) * (cellsNum_y_frac + 1)];
	const FracVariable& next = cell.u_next;
	const FracVariable& nebr = beta.u_next;
	return cell.hy * cell.hz * w2 * w2 / props_w.visc * (nebr.p - next.p) / (cell.hx + beta.hx) * (1.0 - (cell.y / w2) * (cell.y / w2));
};

PoroTapeVariable AcidRecFracMov::solvePoro(const PoroCell& cell, const Regime reg)
{
	if (cell.type == PoroType::MIDDLE)
		return solvePoroMid(cell);
	else if (cell.type == PoroType::WELL_LAT)
		return solvePoroLeft(cell, reg);
	else if (cell.type == PoroType::RIGHT)
		return solvePoroRight(cell);
	else
		return solvePoroBorder(cell);
}
PoroTapeVariable AcidRecFracMov::solvePoroMid(const PoroCell& cell)
{
	assert(cell.type == PoroType::MIDDLE);
	const auto& props = *cell.props;
	auto& next = x_poro[cell.num];
	const auto& prev = cell.u_prev;
	adouble rate = getReactionRate(next, props);// / reac.indices[REACTS::ACID] / reac.comps[REACTS::ACID].mol_weight;

	PoroTapeVariable res;
	adouble m = props.getPoro(next.m, next.p);
	adouble m_prev = props.getPoro(prev.m, prev.p);
	res.m = (1.0 - m) * props.getDensity(next.p) - (1.0 - m_prev) * props.getDensity(prev.p) -
		ht * reac.indices[REACTS::DOLOMITE] * reac.comps[REACTS::DOLOMITE].mol_weight * rate;
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
	double dist1, dist2;
	for (size_t i = 0; i < NEBRS_NUM; i++)
	{
		const PoroCell& beta = cells_poro[neighbor[i]];
		const auto& nebr = x_poro[beta.num];
		if (i < 2)
		{
			dist1 = cell.hy / 2.0;
			dist2 = beta.hy / 2.0;
		}
		else if (i < 4)
		{
			dist1 = cell.hz / 2.0;
			dist2 = beta.hz / 2.0;
		}
		else
		{
			dist1 = cell.hx / 2.0;
			dist2 = beta.hx / 2.0;
		}
		upwd_idx = getUpwindIdx(cell, beta);
		const auto& upwd = x_poro[upwd_idx];

		adouble dens_w = getAverage(props_w.getDensity(next.p, next.xa, next.xw, next.xs), dist1, props_w.getDensity(nebr.p, nebr.xa, nebr.xw, nebr.xs), dist2);
		adouble dens_o = getAverage(props_o.getDensity(next.p), dist1, props_o.getDensity(nebr.p), dist2);
		adouble buf_w = ht / cell.V * getPoroTrans(cell, next, beta, nebr) * (next.p - nebr.p) *
			dens_w * props_w.getKr(upwd.sw, upwd.m, cells_poro[upwd_idx].props) / props_w.getViscosity(upwd.p, upwd.xa, upwd.xw, upwd.xs);
		adouble buf_o = ht / cell.V * getPoroTrans(cell, next, beta, nebr) * (next.p - nebr.p) *
			dens_o * props_o.getKr(upwd.sw, upwd.m, cells_poro[upwd_idx].props) / props_o.getViscosity(upwd.p);

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

	/*if (cell.num % (cellsNum_y_poro + 2) > cellsNum_y_poro / 3)
	{
		res.p *= sqrt(cell.V);
		res.sw *= sqrt(cell.V);
		res.xw *= sqrt(cell.V);
		res.xa *= sqrt(cell.V);
		res.xs *= sqrt(cell.V);
	}*/
	res.m *= 2.0;
	res.p *= 2.0;
	res.sw *= 2.0;
	res.xa *= 2.0;
	res.xs *= 2.0;
	res.xw *= 2.0;

	return res;
}
PoroTapeVariable AcidRecFracMov::solvePoroLeft(const PoroCell& cell, const Regime reg)
{
	assert(cell.type == PoroType::WELL_LAT);
	const auto& props = *cell.props;

	const auto& beta_poro1 = cells_poro[cell.num + 1];
	const auto& beta_poro2 = cells_poro[cell.num + 2];
	const auto& nebr_poro1 = x_poro[cell.num + 1];
	const auto& nebr_poro2 = x_poro[cell.num + 2];

	const auto& beta1 = cells_frac[getFracNebr(cell.num)];
	const auto& beta2 = cells_frac[beta1.num - 1];
	//const auto& beta3 = cells_frac[beta1.num - 2];
	const auto& nebr1 = x_frac[beta1.num];
	const auto& nebr2 = x_frac[beta2.num];
	//const auto& nebr3 = x_frac[beta3.num];

	auto& next = x_poro[cell.num];
	const auto& prev = cell.u_prev;

	const double dist0 = cell.hy / 2.0;
	const double dist1 = beta1.hy / 2.0;
	const double dist2 = beta2.hy / 2.0;
	const double dist_poro1 = beta_poro1.hy / 2.0;
	const double dist_poro2 = beta_poro2.hy / 2.0;

	PoroTapeVariable res;
	res.m = ((next.m - nebr_poro1.m) /*- (nebr_poro1.m - nebr_poro2.m) *
		(dist0 + dist_poro1) / (dist_poro1 + dist_poro2)*/) / P_dim * 2.0;
	res.p = ((next.p - nebr1.p) - (nebr1.p - nebr2.p) *
		(dist0 + dist1) / (dist1 + dist2)) / P_dim * 2.0;
	res.sw = (next.sw - (1.0 - props.s_oc)) / P_dim * 2.0;
	res.xw = (next.xw - (1.0 - next.xa)) / P_dim * 2.0;
	res.xs = next.xs / P_dim * 2.0;

	if (reg == INJECTION || reg == STOP)
	{
		adouble vL = getFlowLeak(beta1);
		res.xa = (next.xa - nebr1.c) / P_dim / 5.0;
		//res.p = getPoroTrans(cell, next, beta_poro1, nebr_poro1) * props_w.getKr(next.sw, next.m, cell.props) / props_w.getViscosity(next.p, next.xa, next.xw, next.xs) * (next.p - nebr_poro1.p) 
		//	+ (cell.hx * cell.hz * (vL - 2.0 * props_w.D_e * (next.xa - nebr1.c) / beta1.hy));
		//res.xa = next.xa * getPoroTrans(cell, next, beta_poro1, nebr_poro1).value() * props_w.getKr(next.sw, next.m, cell.props).value() / props_w.getViscosity(next.p, next.xa, next.xw, next.xs).value() * (next.p - nebr_poro1.p).value();
		//if (reg == INJECTION)
		//	res.xa -= (next.m.value() * cell.hx * cell.hz * (nebr1.c * vL.value() /*- 2.0 * props_w.D_e * (next.xa - nebr1.c) / beta1.hy*/));
		//else
		//	res.xa -= (next.m * cell.hx * cell.hz * (next.xa * vL /*- 2.0 * props_w.D_e * (next.xa - nebr1.c) / beta1.hy*/));
	}

	return res;
}
PoroTapeVariable AcidRecFracMov::solvePoroRight(const PoroCell& cell)
{
	assert(cell.type == PoroType::RIGHT);
	const auto& beta = cells_poro[cell.num - 1];
	const auto& props = *cell.props;

	const auto& next = x_poro[cell.num];
	const auto& nebr = x_poro[beta.num];

	adouble rightIsPres = rightBoundIsPres;
	PoroTapeVariable res;
	res.m = (props.getPoro(next.m, next.p) - props.getPoro(nebr.m, nebr.p)) / P_dim;
	condassign(res.p, rightIsPres, (next.p - props.p_out + grav * props_w.dens_stc * cell.z) / P_dim, (next.p - nebr.p) / P_dim);
	res.sw = (next.sw - nebr.sw) / P_dim;
	res.xw = (next.xw - nebr.xw) / P_dim;
	res.xa = (next.xa - nebr.xa) / P_dim;
	res.xs = (next.xs - nebr.xs) / P_dim;
	return res;
}
PoroTapeVariable AcidRecFracMov::solvePoroBorder(const PoroCell& cell)
{
	assert(cell.type == PoroType::SIDE_LEFT || cell.type == PoroType::SIDE_RIGHT || 
			cell.type == PoroType::TOP || cell.type == PoroType::BOTTOM);
	int beta_idx = -1;
	if (cell.type == PoroType::SIDE_LEFT)
		beta_idx = cell.num + (cellsNum_z + 2) * (cellsNum_y_poro + 2);
	else if (cell.type == PoroType::SIDE_RIGHT)
		beta_idx = cell.num - (cellsNum_z + 2) * (cellsNum_y_poro + 2);
	else if (cell.type == PoroType::BOTTOM)
		beta_idx = cell.num + (cellsNum_y_poro + 2);
	else if (cell.type == PoroType::TOP)
		beta_idx = cell.num - (cellsNum_y_poro + 2);
	const auto& beta = cells_poro[beta_idx];
	const auto& props = *cell.props;

	const auto& next = x_poro[cell.num];
	const auto& nebr = x_poro[beta.num];

	adouble rightIsPres = rightBoundIsPres;
	PoroTapeVariable res;
	res.m = (props.getPoro(next.m, next.p) - props.getPoro(nebr.m, nebr.p)) / P_dim;
	res.p = (next.p - nebr.p + grav * props_w.dens_stc * (cell.z - beta.z)) / P_dim;
	res.sw = (next.sw - nebr.sw) / P_dim;
	res.xw = (next.xw - nebr.xw) / P_dim;
	res.xa = (next.xa - nebr.xa) / P_dim;
	res.xs = (next.xs - nebr.xs) / P_dim;
	return res;
}
FracTapeVariable AcidRecFracMov::solveFrac(const FracCell& cell, const Regime reg)
{
	if (cell.type == FracType::FRAC_MID || cell.type == FracType::FRAC_OUT)
		return solveFracMid(cell, reg);
	else if (cell.type == FracType::FRAC_BORDER)
		return solveFracBorder(cell);
	else if (cell.type == FracType::FRAC_IN)
		return solveFracIn(cell);
}
FracTapeVariable AcidRecFracMov::solveFracIn(const FracCell& cell)
{
	assert(cell.type == FracType::FRAC_IN);
	const auto& beta = cells_frac[cell.num + (cellsNum_z + 2) * (cellsNum_y_frac + 1)];
	const auto& next = x_frac[cell.num];
	const auto& nebr = x_frac[beta.num];
	const adouble leftIsRate = leftBoundIsRate;
	FracTapeVariable res;
    const double w2 = props_frac.w2[getWidthIdByFracId(cell.num)];
	if (leftBoundIsRate)
	{
		if (cell.hy * cell.hz != 0.0)
			res.p = ((next.p - beta.u_prev.p) +	Qcell[cell.num] / (1.0 - (cell.y / w2) * (cell.y / w2)) / w2 / w2 * props_w.visc
				/ (cell.hy * cell.hz) * beta.hx) / P_dim;
		else
			res.p = (next.p - nebr.p) / P_dim;
		res.c = (next.c - c) / P_dim;
	}
	else
	{
		res.p = (next.p - Pwf + grav * props_w.dens_stc * cell.z) / P_dim;
		res.c = (next.c - c) / P_dim;
	}
	return res;
}
FracTapeVariable AcidRecFracMov::solveFracBorder(const FracCell& cell)
{
	assert(cell.type == FracType::FRAC_BORDER);
	const auto& next = x_frac[cell.num];
	int beta_idx = -1;
	if (cell.hx == 0.0)
		beta_idx = cell.num - (cellsNum_z + 2) * (cellsNum_y_frac + 1);
	else if (cell.hy == 0.0)
		beta_idx = cell.num + 1;
	else
		beta_idx = (cell.z == 0.0) ? cell.num + (cellsNum_y_frac + 1) : cell.num - (cellsNum_y_frac + 1);
	const auto& beta = cells_frac[beta_idx];
	const auto& nebr = x_frac[beta.num];

	FracTapeVariable res;

	res.p = ((next.p - nebr.p) + (cell.z - beta.z) * grav * props_w.dens_stc) / P_dim;
	res.c = (next.c - nebr.c) / P_dim;

	return res;
}
FracTapeVariable AcidRecFracMov::solveFracMid(const FracCell& cell, const Regime reg)
{
	assert(cell.type == FracType::FRAC_MID || FracType::FRAC_OUT);
	const auto& next = x_frac[cell.num];

	FracTapeVariable res;
	int neighbor[6];
	getFracNeighborIdx(cell.num, neighbor);
	const auto& nebr_y_minus = x_frac[neighbor[0]];
	const auto& nebr_x_minus = x_frac[neighbor[2]];
	const auto& nebr_x_plus = x_frac[neighbor[3]];
	const auto& nebr_z_minus = x_frac[neighbor[4]];
	const auto& nebr_z_plus = x_frac[neighbor[5]];

    const auto& cell_x_minus = cells_frac[neighbor[2]];
    const auto& cell_x_plus = cells_frac[neighbor[3]];
    const auto& cell_y_minus = cells_frac[neighbor[0]];
    const auto& cell_z_minus = cells_frac[neighbor[4]];
    const auto& cell_z_plus = cells_frac[neighbor[5]];

	double y_plus = cell.y + cell.hy / 2.0;
	double y_minus = cell.y - cell.hy / 2.0;
	adouble vL = getFlowLeak(cell);
    const double w2 = props_frac.w2[getWidthIdByFracId(cell.num)];
	const double alpha = -w2 * w2 / props_w.visc * (1.0 - (cell.y / w2) * (cell.y / w2));	
	adouble vy_minus = vL * (1.5 * (y_minus / w2) - 0.5 * y_minus * y_minus * y_minus / w2 / w2 / w2);
	adouble vy_plus = vL * (1.5 * (y_plus / w2) - 0.5 * y_plus * y_plus * y_plus / w2 / w2 / w2);
	adouble vz_minus = alpha * ((next.p - nebr_z_minus.p) / (cell.hz + cell_z_minus.hz) - grav * props_w.dens_stc / 2.0);
	adouble vz_plus = alpha * ((next.p - nebr_z_plus.p) / (cell.hz + cell_z_plus.hz) + grav * props_w.dens_stc / 2.0);
	adouble vx_minus = alpha * (next.p - nebr_x_minus.p) / (cell.hx + cell_x_minus.hx);
	adouble vx_plus = alpha * (next.p - nebr_x_plus.p) / (cell.hx + cell_x_plus.hx);
	if (max_vel_x < fabs(vx_minus.value()))
		max_vel_x = fabs(vx_minus.value());
	if (max_vel_x < fabs(vx_plus.value()))
		max_vel_x = fabs(vx_plus.value());
	if (max_vel_y < fabs(vL.value()))
		max_vel_y = fabs(vL.value());
    if (max_vel_z < fabs(vz_plus.value()))
        max_vel_z = fabs(vz_plus.value());
    if (max_vel_z < fabs(vz_minus.value()))
        max_vel_z = fabs(vz_minus.value());

	/*adouble diff_minus = 2.0 * props_w.D_e * (next.c - nebr_y_minus.c) / (cell.hy + cells_frac[neighbor[0]].hy);
	adouble diff_plus, c_y_plus, c_y_minus;
	if (cell.type == FracType::FRAC_OUT)
	{
		diff_plus = 2.0 * props_w.D_e * (next.c - x_poro[neighbor[1]].xa) / cell.hy;
		c_y_plus = (reg == INJECTION) ? next.c : x_poro[neighbor[1]].xa;
	}
	else
	{
		diff_plus = 2.0 * props_w.D_e * (next.c - x_frac[neighbor[1]].c) / (cell.hy + cells_frac[neighbor[1]].hy);
		c_y_plus = (reg == INJECTION) ? next.c : x_frac[neighbor[1]].c;
	}
	c_y_minus = (reg == INJECTION) ? x_frac[neighbor[0]].c : next.c;*/
	const double sx = cell.hy * cell.hz, sy = cell.hx * cell.hz, sz = cell.hx * cell.hy;
	res.p = ((sx * (vx_minus + vx_plus) + sz * (vz_minus + vz_plus) - sy * (vy_plus - vy_minus)) / alpha / sx * (cell.hx + cells_frac[neighbor[4]].hx));
	res.p *= cell.V;
	/*res.c = ((next.c - cell.u_prev.c) * cell.V - ht *
		(sx * (x_frac[getUpwindIdx(cell, cells_frac[neighbor[4]])].c * vx_minus +
			x_frac[getUpwindIdx(cell, cells_frac[neighbor[5]])].c * vx_plus) +
			sy * (c_y_minus * vy_minus - c_y_plus * vy_plus - diff_plus - diff_minus) +
			sz * (x_frac[getUpwindIdx(cell, cells_frac[neighbor[2]])].c * vz_minus +
				x_frac[getUpwindIdx(cell, cells_frac[neighbor[3]])].c * vz_plus)));*/

    const double vx1 = vx_minus.value();
    const double vx2 = vx_plus.value();
    const double vy1 = vy_minus.value();
    const double vy2 = vy_plus.value();

    double rx_minus1, rx_minus2, rx_plus1, rx_plus2, ry_minus1, ry_minus2, ry_plus1, ry_plus2, c_y_plus;
    if (cell_x_minus.type == FracType::FRAC_IN)
    {
        rx_minus1 = cell.hx / 2.0;
        rx_minus2 = cell_x_minus.hx / 2.0;
    }
    else
    {
        rx_minus1 = (cell.hx + fabs(vx_minus.value()) * ht) / 2.0;
        rx_minus2 = (cell_x_minus.hx - fabs(vx_minus.value()) * ht) / 2.0;
    }
    rx_plus1 = (cell.hx - fabs(vx_plus.value()) * ht) / 2.0;
    rx_plus2 = (cell_x_plus.hx + fabs(vx_plus.value()) * ht) / 2.0;
    if (cell_y_minus.type == FracType::FRAC_BORDER)
    {
        ry_minus1 = cell.hy / 2.0;
        ry_minus2 = cell_y_minus.hy / 2.0;
    }
    else
    {
        ry_minus1 = (cell.hy + fabs(vy_minus.value()) * ht) / 2.0;
        ry_minus2 = (cell_y_minus.hy - fabs(vy_minus.value()) * ht) / 2.0;
    }
    ry_plus1 = (cell.hy - fabs(vy_plus.value()) * ht) / 2.0;
    if (cell.type == FracType::FRAC_OUT)
    {
        ry_plus2 = fabs(vy_plus.value()) * ht / 2.0;
        c_y_plus = cells_poro[neighbor[1]].u_prev.xa;
    }
    else
    {
        ry_plus2 = (cells_frac[neighbor[1]].hy + fabs(vy_plus.value()) * ht) / 2.0;
        c_y_plus = cells_frac[neighbor[1]].u_prev.c;
    }
	
	//adouble fx_plus = -vx_plus.value() * (rx_plus1 * cell_x_plus.u_prev.c + rx_plus2 * cell.u_prev.c) / (rx_plus1 + rx_plus2);
	//adouble fx_minus = vx_minus.value() * (rx_minus1 * cell_x_minus.u_prev.c + rx_minus2 * cell.u_prev.c) / (rx_minus1 + rx_minus2);
	//adouble fy_plus = vy_plus.value() * (ry_plus1 * c_y_plus + ry_plus2 * cell.u_prev.c) / (ry_plus1 + ry_plus2);
	//adouble fy_minus = vy_minus.value() * (ry_minus1 * cell_y_minus.u_prev.c + ry_minus2 * cell.u_prev.c) / (ry_minus1 + ry_minus2);
	adouble fx_plus = -next.c * vx_plus.value();
	adouble fx_minus = cell_x_minus.u_prev.c * vx_minus.value();
	adouble fy_plus = next.c * vy_plus.value();
	adouble fy_minus = cell_y_minus.u_prev.c * vy_minus.value();
    adouble fz_plus, fz_minus;
    if(cell.u_prev.p > cell_z_plus.u_prev.p)
        fz_plus = -next.c * vz_plus.value();
    else
        fz_plus = -cell_z_plus.u_prev.c * vz_plus.value();
    if (cell.u_prev.p > cell_z_minus.u_prev.p)
        fz_minus = next.c * vz_minus.value();
    else
        fz_minus = cell_z_minus.u_prev.c * vz_minus.value();

	res.c = (next.c - cell.u_prev.c) * cell.V + ht * (sx * (fx_plus - fx_minus) + sy * (fy_plus - fy_minus) + sz * (fz_plus - fz_minus));
    //res.c = (next.c - cell.u_prev.c) * cell.V + ht * (sx * vx_minus.value() * (cell.u_prev.c - cell_x_minus.u_prev.c) + 
    //                                                    sy * vy_minus.value() * (cell.u_prev.c - cell_y_minus.u_prev.c));


	res.p /= 2.5;
	res.c /= 2.5;
	return res;
}