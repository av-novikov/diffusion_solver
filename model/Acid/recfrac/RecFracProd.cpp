#include "model/Acid/recfrac/RecFracProd.hpp"
#include <assert.h>
#include <numeric>

using namespace acidrecfrac_prod;

RecFracProd::RecFracProd()
{
	grav = 9.8;
	Volume = 0.0;
}
RecFracProd::~RecFracProd()
{
	delete snapshotter;
	delete[] x, h;
}
void RecFracProd::setProps(Properties& props)
{
	props_sk = props.props_sk;
	rightBoundIsPres = props.rightBoundIsPres;

    prev_cellsNum_x = props.cellsNum_x;
    prev_cellsNum_y = props.cellsNum_y_poro;
    prev_cellsNum_z = props.cellsNum_z;
    cellsNum_x = props.prod_props.nx;
    cellsNum_y = props.prod_props.ny;
    cellsNum_z = props.cellsNum_z;
    prev_x_size = props.props_frac.l2;
    prev_y_size = props.re;
    prev_z_size = props.props_frac.height;
    x_size = props.prod_props.x_size;
    y_size = props.prod_props.y_size;
    z_size = props.prod_props.z_size;
    
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

	makeDimLess();

	props_o.b = setDataset(props.B_oil, P_dim / BAR_TO_PA, 1.0);
	props_o.Rs = setDataset(props.Rs, P_dim / BAR_TO_PA, 1.0);
}
void RecFracProd::makeDimLess()
{
	T_dim = props_sk[0].t_init;
	//R_dim = props_frac.l2 / 8.0;
    R_dim = 1000.0;
	t_dim = 1.7 * 3600.0;
	//t_dim = 3600.0;
	P_dim = props_sk[0].p_init;

    x_size /= R_dim;        y_size /= R_dim;        z_size /= R_dim;
    prev_x_size /= R_dim;   prev_y_size /= R_dim;   prev_z_size /= R_dim;

	// Temporal properties
	ht /= t_dim;
	ht_min /= t_dim;
	ht_max /= t_dim;
	// Skeleton properties
		auto& sk = props_sk[0];
		sk.perm /= (R_dim * R_dim);
		sk.beta /= (1.0 / P_dim);
		sk.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
		sk.p_init /= P_dim;
		sk.p_out /= P_dim;
		sk.p_ref /= P_dim;
		sk.hx /= R_dim;
		sk.hz /= R_dim;
		sk.t_init /= T_dim;
		sk.height /= R_dim;

	Q_dim = R_dim * R_dim * R_dim / t_dim;
	for (int i = 0; i < periodsNum; i++)
	{
		period[i] /= t_dim;
		rate[i] /= Q_dim;
		pwf[i] /= P_dim;
	}

	grav /= (R_dim / t_dim / t_dim);

	props_w.visc /= (P_dim * t_dim);
	props_w.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_w.beta /= (1.0 / P_dim);
	props_w.D_e /= (R_dim * R_dim / t_dim);
	props_o.visc /= (P_dim * t_dim);
	props_o.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_o.gas_dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_o.beta /= (1.0 / P_dim);
	props_o.p_ref /= P_dim;
}
void RecFracProd::buildGrid(std::vector<PoroCell>& cells_poro)
{
    int counter = 0;
    double cx, cy, cz, hx, hy, hz, cur_hx, cur_hy, cur_hz;
    hx = (x_size - prev_x_size) / (cellsNum_x - prev_cellsNum_x);
    hy = (y_size - prev_y_size) / (cellsNum_y - prev_cellsNum_y);
    hz = (z_size - prev_z_size) / (cellsNum_z - prev_cellsNum_z);

    for (auto& cell : cells_poro)
        cell.y -= cells_poro[0].y;

	// Left border
	auto cur_type = Type::SIDE_LEFT;
    cur_hx = 0.0;	 cx = 0.0;
	for (int k = 0; k < cellsNum_z + 2; k++)
	{
        cy = 0.0;
        for (int i = 0; i < cellsNum_y + 2; i++)
        {
            if (k < prev_cellsNum_z + 2 && i < prev_cellsNum_y)
            {
                cells.push_back(Cell(cells_poro[i + k * (cellsNum_y + 2)], counter, cur_type));
                cy = cells.back().y;    cz = cells.back().z;
            }
            else
            {
                cur_hy = (i == cellsNum_y + 1) ? 0.0 : hy;
                cells.push_back(Cell(counter++, cx, cy, cz, cur_hx, cur_hy, cur_hz, cur_type));
            }
        }
		/*if (k == 0)
			hz = cz = 0.0;
		else if (k == cellsNum_z + 1)
		{
			hz = 0.0;
			cz = z_size;
		}
		else
		{
			//hz = props_frac.height / cellsNum_z;
			cz += hz / 2.0;
			hz = props_sk[k - 1].height;
			cz += hz / 2.0;
		}

		for (int i = 0; i < cellsNum_y + 2; i++)
		{
			if (i == cellsNum_y + 1)
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

			cells.push_back(PoroCell(counter++, cx, cy, cz, hx, hy, hz, cur_type));
			auto& cell = cells_poro.back();
		}*/
	}
}
void RecFracProd::processGeometry()
{
}
void RecFracProd::setPerforated()
{
	/*height_perf = 0.0;
	vector<pair<int, int> >::iterator it;
	for (const auto& cell : cells_frac)
	{
		if (cell.type == FracType::FRAC_IN)
		{
			Qcell[cell.num] = 0.0;
			height_perf += cell.hy * cell.hz;
		}
	}*/
};
void RecFracProd::setPeriod(int period)
{
	/*leftBoundIsRate = LeftBoundIsRate[period];
	if (leftBoundIsRate)
	{
		Q_sum = rate[period];
		if (Q_sum == 0.0)
		{
			for(auto& cell : cells)
				cell.u_prev.p = cell.u_iter.c = cell.u_next.c = 0.0;
		}

		if (period == 0 || rate[period - 1] < EQUALITY_TOLERANCE) 
		{
			std::map<int, double>::iterator it;
			for (it = Qcell.begin(); it != Qcell.end(); ++it)
				it->second = 3.0 * Q_sum * cells[it->first].hy / 2.0 / props_frac.w2 * 
					(1.0 - (cells_frac[it->first].y / props_frac.w2) * (cells_frac[it->first].y / props_frac.w2));
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
	c = cs[period];*/
}
void RecFracProd::setInitialState()
{
	//const auto& props = props_sk.back();
/*	for (auto& cell : cells)
	{
		cell.u_prev.p = cell.u_iter.p = cell.u_next.p = props.p_init - grav * props_w.dens_stc * cell.z;
        cell.u_prev.p = cell.u_iter.p = cell.u_next.p = 1.0;
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
	h = new adouble[ var_frac_size * cellsNum_frac + var_poro_size * cellsNum_poro ];*/
}
double RecFracProd::getRate(const int idx) const
{
	/*const FracCell& cell = cells_frac[idx];
	assert(cell.type == FracType::FRAC_IN);
	assert(cell.hy != 0.0 && cell.hz != 0.0);
	const FracCell& beta = cells_frac[cell.num + (cellsNum_z + 2) * (cellsNum_y_frac + 1)];
	const FracVariable& next = cell.u_next;
	const FracVariable& nebr = beta.u_next;

	return cell.hy * cell.hz * props_frac.w2 * props_frac.w2 / props_w.visc * (nebr.p - next.p) / (cell.hx + beta.hx) * 
				(1.0 - (cell.y / props_frac.w2) * (cell.y / props_frac.w2));*/
    return 0.0;
};