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
void RecFracProd::setProps(Properties& props, const std::vector<double>& _widths)
{
	widths = _widths;
    prefix = props.prefix;
	props_sk = props.props_sk;
	rightBoundIsPres = props.rightBoundIsPres;
	int skeletonsNum = props.props_sk.size();
	for (int j = 0; j < skeletonsNum; j++)
	{
		props_sk[j].perm = MilliDarcyToM2(props_sk[j].perm);
	}

    prev_cellsNum_x = props.cellsNum_x;
    prev_cellsNum_y = props.cellsNum_y_poro;
    prev_cellsNum_z = props.cellsNum_z;
    cellsNum_x = props.prod_props.nx;
    cellsNum_y = props.prod_props.ny;
    cellsNum_z = props.cellsNum_z;
	cellsNum = (cellsNum_x + 2) * (cellsNum_y + 2) * cellsNum_z;
	num_input_cells = 0.7 * prev_cellsNum_y;
    prev_x_size = props.props_frac.l2;
    prev_y_size = props.re - props.props_frac.w2;
    prev_z_size = props.props_frac.height;
    x_size = props.prod_props.x_size;
    y_size = props.prod_props.y_size;
    z_size = props.prod_props.z_size;
    R_dim = props.prod_props.R_dim;

	period.clear();
    periodsNum = 1;// props.timePeriods.size();
	rate.resize(periodsNum);
	pwf.resize(periodsNum);
	int rate_idx = 0, pres_idx = 0;
	double max_t = 10 * 86400.0 * 365.0;
	for (int i = 0; i < periodsNum; i++)
	{
        LeftBoundIsRate.push_back(false);// props.LeftBoundIsRate[i]);
        period.push_back(max_t);// props.timePeriods[i]);
        pwf[i] = 0.5 * props_sk.back().p_out;// .pwf[pres_idx++];
        rate[i] = 0.0;
	}

	// Temporal properties
    //ht = 10000.0;// props.ht;
    ht = 1000.0;
    ht_min = 1000.0;// props.ht_min;
    ht_max = 1000 * 3600.0;// props.ht_max;

	props_w = props.props_w;
	props_w.visc = cPToPaSec(props_w.visc);
	props_o = props.props_o;
	props_o.visc = cPToPaSec(props_o.visc);

	makeDimLess();
}
void RecFracProd::makeDimLess()
{
	T_dim = props_sk[0].t_init;
	t_dim = 3600.0;
	//t_dim = 3600.0;
	P_dim = props_sk[0].p_init;

    x_size /= R_dim;        y_size /= R_dim;        z_size /= R_dim;
    prev_x_size /= R_dim;   prev_y_size /= R_dim;   prev_z_size /= R_dim;

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
		sk.p_init /= P_dim;
		sk.p_out /= P_dim;
		sk.p_ref /= P_dim;
		sk.hx /= R_dim;
		sk.hz /= R_dim;
		sk.t_init /= T_dim;
		sk.height /= R_dim;
		sk.s_compres /= P_dim;
		sk.s_failure /= P_dim;
	}
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
	props_w.p_ref /= P_dim;
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
    Volume = 0.0;

    const double tmp = cells_poro[0].y;
    for (auto& cell : cells_poro)
        cell.y -= tmp;

	const int k = 0;
    auto cur_type = Type::MIDDLE;
    for (int j = 0; j < cellsNum_x + 2; j++)
    {
        if (j < prev_cellsNum_x + 1)
        {
            cx = cells_poro[j * prev_cellsNum_z * (prev_cellsNum_y + 2)].x;
            cur_hx = cells_poro[j * prev_cellsNum_z * (prev_cellsNum_y + 2)].hx;
        }
        else
        {
            cx = cells_poro[(prev_cellsNum_x + 1) * prev_cellsNum_z * (prev_cellsNum_y + 2)].x +
                (double)(j - (prev_cellsNum_x + 1) + 0.5) * hx;
            cur_hx = (j == cellsNum_x + 1) ? 0.0 : hx;
        }                                            


        //for (int k = 0; k < cellsNum_z + 2; k++)
        //{
            cz = cells_poro[k * (prev_cellsNum_y + 2)].z;
            cur_hz = cells_poro[k * (prev_cellsNum_y + 2)].hz;
            cy = 0.0;

            for (int i = 0; i < cellsNum_y + 2; i++)
            {
                cur_type = (j == 0) ? Type::SIDE_LEFT : Type::MIDDLE;
                cur_type = (j == cellsNum_x + 1) ? Type::SIDE_RIGHT : cur_type;
               //cur_type = (k == 0) ? Type::BOTTOM : cur_type;
               // cur_type = (k == cellsNum_z + 1) ? Type::TOP : cur_type;
                cur_type = (i == 0 && cur_type == Type::MIDDLE) ? Type::WELL_LAT : cur_type;
                cur_type = (i == cellsNum_y + 1 && cur_type == Type::MIDDLE) ? Type::RIGHT : cur_type;
                cur_type = (i < num_input_cells + 1 && cur_type == Type::SIDE_LEFT) ? Type::FRAC_IN : cur_type;
                if (j < prev_cellsNum_x + 1 && i < prev_cellsNum_y + 1)
                {
                    cells.push_back(Cell(cells_poro[i + (prev_cellsNum_y + 2) * (k + j * prev_cellsNum_z)], counter++, cur_type));
                    const auto& cell = cells.back();
                    Volume += cell.V;
                    cy = cell.y + cell.hy / 2.0;
                }
                else
                {
                    cur_hy = (i < prev_cellsNum_y + 1) ? cells_poro[i].hy : hy;
                    cur_hy = (i == cellsNum_y + 1) ? 0.0 : cur_hy;
                    cells.push_back(Cell(counter++, cx, cy + cur_hy / 2.0, cz, cur_hx, cur_hy, cur_hz, cur_type));
                    Volume += cells.back().V;
                    cy += cur_hy;
                }
            }
        //}
    }
}
void RecFracProd::buildBetterGrid(std::vector<PoroCell>& cells_poro)
{
    int j_ind, counter = 0;
    double cx, cy, cz, hx, hz, cur_hx, cur_hy, cur_hz, width_cur;
    hx = (x_size - prev_x_size) / (cellsNum_x - add_cellsNum_x - prev_cellsNum_x + 1);
    hz = (z_size - prev_z_size) / (cellsNum_z - prev_cellsNum_z + 1);
    Volume = 0.0;

	/*std::fill(widths.begin(), widths.end(),
		cells_poro[num_input_cells].y + cells_poro[num_input_cells].hy / 2.0 - cells_poro[0].y);*/
	double max_width = *std::max_element(widths.begin(), widths.end());
	double max_prev_y = 1.5 * max_width;
	MIN_FRAC_WIDTH = 0.01 * fmax(max_width, cells_poro[num_input_cells].y + 
											cells_poro[num_input_cells].hy / 2.0 - cells_poro[0].y);

	double delta = 0.0;// -9.0 * max_prev_y;
	double r_prev = max_prev_y - delta;
	double logMax = log((y_size - delta) / (max_prev_y - delta));
	double logStep = logMax / (double)(cellsNum_y - 2);

    double delta_x = 0.0;
    double rw = 0.02 / R_dim;
    double x_prev = cells_poro[prev_cellsNum_z * (prev_cellsNum_y + 2)].hx - delta_x;
    double logMax_x = log((x_prev - delta_x) / (rw - delta_x));
    double logStep_x = logMax_x / (double)(add_cellsNum_x);
    //const double tmp = cells_poro[0].y;
    //for (auto& cell : cells_poro)
    //    cell.y -= tmp;

    auto cur_type = Type::MIDDLE;
    x_prev = rw - delta_x;
    for (int j = 0; j < cellsNum_x + 2; j++)
    {
        if (j < add_cellsNum_x + 1 && j > 0)
        {
            cx = delta_x + x_prev * (exp(logStep_x) + 1.0) / 2.0;
            cur_hx = x_prev * (exp(logStep_x) - 1.0);
            x_prev *= exp(logStep_x);
        }
        else if (j < add_cellsNum_x + prev_cellsNum_x)
        {
            int ind = (j == 0 ? j : j + 1 - add_cellsNum_x);
            cx = cells_poro[ind * prev_cellsNum_z * (prev_cellsNum_y + 2)].x;
            cur_hx = cells_poro[ind * prev_cellsNum_z * (prev_cellsNum_y + 2)].hx;
        }
        else
        {
            cx = cells_poro[(prev_cellsNum_x + 1) * prev_cellsNum_z * (prev_cellsNum_y + 2)].x +
                (double)(j - (add_cellsNum_x + prev_cellsNum_x) + 0.5) * hx;
            cur_hx = (j == cellsNum_x + 1) ? 0.0 : hx;
        }

		const int k = 0;
        //for (int k = 0; k < cellsNum_z + 2; k++)
        //{
            cz = cells_poro[k * (prev_cellsNum_y + 2)].z;
            cur_hz = cells_poro[k * (prev_cellsNum_y + 2)].hz;
            cy = cur_hy = 0.0;

			//width_cur = (0.95 - 0.9 * (double)(j) / (double)(cellsNum_x + 1)) * max_prev_y;
			j_ind = j > add_cellsNum_x ? j - add_cellsNum_x : 0;
			j_ind = j > add_cellsNum_x + prev_cellsNum_x - 1 ? prev_cellsNum_x - 1 : j_ind;
			width_cur = widths[j_ind];
			if (width_cur < MIN_FRAC_WIDTH)
				width_cur = MIN_FRAC_WIDTH;
            r_prev = max_prev_y - delta;
            for (int i = 0; i < cellsNum_y + 2; i++)
            {
                cur_type = (j == 0) ? Type::SIDE_LEFT : Type::MIDDLE;
                cur_type = (j == cellsNum_x + 1) ? Type::SIDE_RIGHT : cur_type;
                //cur_type = (k == 0) ? Type::BOTTOM : cur_type;
                //cur_type = (k == cellsNum_z + 1) ? Type::TOP : cur_type;
                cur_type = (i == 0 && cur_type == Type::MIDDLE) ? Type::WELL_LAT : cur_type;
                cur_type = (i == cellsNum_y + 1 && cur_type == Type::MIDDLE) ? Type::RIGHT : cur_type;
                cur_type = (i < 2 && cur_type == Type::SIDE_LEFT) ? Type::FRAC_IN : cur_type;

                if (j < add_cellsNum_x + prev_cellsNum_x && i < 3)
                {
                    //cells.push_back(Cell(cells_poro[i + (prev_cellsNum_y + 2) * (k + j * (prev_cellsNum_z + 2))], counter++, cur_type));
                    //const auto& old_cell = cells_poro[i + (prev_cellsNum_y + 2) * (k + j * (prev_cellsNum_z + 2))];
                    cells.push_back(Cell(counter++, cx, cy, cz, cur_hx, cur_hy, cur_hz, cur_type));
                    const auto& cell = cells.back();
                    Volume += cell.V;

					if (i == 0)
					{
						cy = width_cur / 2.0;
						cur_hy = width_cur;
					}
					else if (i == 1)
					{
						cy += (cur_hy + (max_prev_y - width_cur)) / 2.0;
						cur_hy = max_prev_y - width_cur;
					}
                    
                }
                else
                {
                    if (i == 0)
                    {
                        cy = 0.0;
                        cur_hy = 0.0;
                    }
                    else if (i == 1)
                    {
                        cy = width_cur / 2.0;
                        cur_hy = width_cur;
                    }
					else if (i == 2)
					{
						cy += (cur_hy + (max_prev_y - width_cur)) / 2.0;
						cur_hy = max_prev_y - width_cur;
					}
                    else if (i == cellsNum_y + 1)
                    {
                        cur_hy = 0.0;
                        cy = y_size;
                    }
                    else
                    {
                        cy = delta + r_prev * (exp(logStep) + 1.0) / 2.0;
                        cur_hy = r_prev * (exp(logStep) - 1.0);
                        r_prev *= exp(logStep);
                    }
                    cells.push_back(Cell(counter++, cx, cy, cz, cur_hx, cur_hy, cur_hz, cur_type));
                    Volume += cells.back().V;
                }
            }
        //}
    }
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
    leftBoundIsRate = false;
    Pwf = pwf[period];
    Q_sum = 0.0;
}
const std::array<double, 3> RecFracProd::calcAvgFracPerm(const std::vector<PoroCell>& cells_poro, const int j, const int k) const
{
	const int tr_idx = j > 0 ? j - 1 : 0;
	double m = 0.0, kx = 0.0, ky = 0.0, k0, vol = 0.0, s = 0, L = 0, y_start;
	y_start = cells_poro[0].y;
	for (int i = 0; i < prev_cellsNum_y; i++)
	{
		const auto& cell0 = cells_poro[i + (prev_cellsNum_y + 2) * (k + j * prev_cellsNum_z)];
		const auto& props = *cell0.props;
		
		if (props.isCompacted(cell0.u_next.m, cell0.u_next.p))
		{
			y_start += cell0.hy;
			continue;
		}

		if (cell0.y - y_start <= widths[tr_idx])
		{
			vol += cell0.V;
			s += cell0.hy * cell0.hz;
			L += cell0.hy;
			k0 = props.getPermCoseni(cell0.u_next.m, cell0.u_next.p).value();

			m += cell0.V * cell0.u_next.m;
			kx += cell0.hy * cell0.hz * k0;
			ky += cell0.hy / k0;
		}
	}
	if (vol == 0.0 || s == 0.0 || L == 0.0)
		return{ props_sk[0].m_init, props_sk[0].perm, props_sk[0].perm };
	else
		return{ m / vol, kx / s, L / ky };
	//return{ 0.3, props_sk[0].getPermCoseni(0.3, 0.0).value(), props_sk[0].getPermCoseni(0.3, 0.0).value() };
}
void RecFracProd::setInitialState(const std::vector<PoroCell>& cells_poro)
{
    int j_ind;
    for (int j = 0; j < cellsNum_x + 2; j++)
    {
        for (int k = 0; k < cellsNum_z; k++)
        {
            for (int i = 0; i < cellsNum_y + 2; i++)
            {
                auto& cell = cells[i + (cellsNum_y + 2) * (k + j * cellsNum_z)];
                cell.props = &props_sk[getSkeletonId(cell)];
                cell.u_prev.p = cell.u_iter.p = cell.u_next.p = cell.props->p_init - grav * props_w.dens_stc * cell.z;
                if (j < add_cellsNum_x + prev_cellsNum_x && i < 2)
                {
                    j_ind = j > add_cellsNum_x ? j - add_cellsNum_x + 1 : 1;
                    j_ind = j == 0 ? j : j_ind;
                    const auto pars = calcAvgFracPerm(cells_poro, j_ind, k);
                    cell.u_prev.m = cell.u_iter.m = cell.u_next.m = pars[0];
                    cell.u_prev.kx = cell.u_iter.kx = cell.u_next.kx = pars[1];
                    cell.u_prev.ky = cell.u_iter.ky = cell.u_next.ky = pars[2];
                }
                else
                {
                    cell.u_prev.m = cell.u_iter.m = cell.u_next.m = cell.props->m_init;
                    cell.u_prev.kx = cell.u_iter.kx = cell.u_next.kx = 
                    cell.u_prev.ky = cell.u_iter.ky = cell.u_next.ky = cell.props->perm;
                }
            }
        }
    }

	x = new TapeVariable[ cellsNum ];
	h = new adouble[ Variable::size * cellsNum ];
}
double RecFracProd::getRate(const Cell& cell) const
{
	assert(cell.hy != 0.0 && cell.hz != 0.0);
	const Cell& beta = cells[cell.num + cellsNum_z * (cellsNum_y + 2)];
	const Variable& next = cell.u_next;
	const Variable& nebr = beta.u_next;
    return getPoroTrans(cell, beta) * (nebr.p - next.p) / props_o.visc;
};

TapeVariable RecFracProd::solveFrac(const Cell& cell)
{
    if (cell.type == Type::MIDDLE)
        return solvePoroMid(cell);
    else if (cell.type == Type::FRAC_IN)
        return solvePoroLeft(cell);
    else if (cell.type == Type::RIGHT || cell.type == Type::SIDE_RIGHT)
        return solvePoroRight(cell);
    else
        return solvePoroBorder(cell);
}
TapeVariable RecFracProd::solvePoroMid(const Cell& cell)
{
    assert(cell.type == Type::MIDDLE);
    const auto& props = *cell.props;
    auto& next = x[cell.num];
    const auto& prev = cell.u_prev;

    TapeVariable res;
    res.p = props.getPoro(cell.u_next.m, next.p) * props_o.getRho(next.p) -
        props.getPoro(cell.u_prev.m, prev.p) * props_o.getRho(prev.p);// *(next.s - prev.s);
    //res.s = props.getPoro(cell.u_next.m, next.p) * props_o.getRho(next.p) * (prev.s - next.s);

    int upwd_idx;
    int neighbor[NEBRS_NUM];
    getPoroNeighborIdx(cell.num, neighbor);
    double ress = res.p.value();
    for (size_t i = 0; i < NEBRS_NUM; i++)
    {
        const auto& beta = cells[neighbor[i]];
        const auto& nebr = x[beta.num];
        upwd_idx = getUpwindIdx(cell, beta);
        const auto& upwd = x[upwd_idx];
        const double trans = getPoroTrans(cell, beta);

        res.p += ht / cell.V * trans * (next.p - nebr.p) * props_o.getRho(next.p) / props_o.visc;
        /*res.s += ht / cell.V * trans * (next.p - nebr.p) *
			props_o.getRho(next.p) * props_o.getKr(upwd.s, cell.u_next.m, cell.props) / props_o.visc;*/
    }

	res.p /= P_dim;
	//if(cell.y * R_dim < 0.01)
	//	res.p /= 1.0 / cell.hy;
    return res;
}
TapeVariable RecFracProd::solvePoroLeft(const Cell& cell)
{
    assert(cell.type == Type::FRAC_IN);
    const auto& beta = cells[cell.num + cellsNum_z  * (cellsNum_y + 2)];
    const auto& props = *cell.props;
    const auto& next = x[cell.num];
    const auto& nebr = x[beta.num];
    TapeVariable res;
    res.p = (next.p - Pwf + grav * props_w.dens_stc * cell.z) / P_dim;
    //res.s = (next.s - nebr.s) / P_dim;
    return res;
}
TapeVariable RecFracProd::solvePoroRight(const Cell& cell)
{
    assert(cell.type == Type::RIGHT || cell.type == Type::SIDE_RIGHT);
    int beta_idx = -1;
    if (cell.type == Type::SIDE_RIGHT)
        beta_idx = cell.num - cellsNum_z * (cellsNum_y + 2);
    else if (cell.type == Type::RIGHT)
        beta_idx = cell.num - 1;
    const auto& beta = cells[beta_idx];
    const auto& props = *cell.props;

    const auto& next = x[cell.num];
    const auto& nebr = x[beta.num];

    adouble rightIsPres = rightBoundIsPres;
    TapeVariable res;
    res.p = (next.p - props.p_out + grav * props_w.dens_stc * cell.z) / P_dim;
	//res.s = (next.s - props.sw_init) / P_dim;// (next.s - nebr.s) / P_dim;
    return res;
}
TapeVariable RecFracProd::solvePoroBorder(const Cell& cell)
{
    assert(cell.type == Type::SIDE_LEFT || cell.type == Type::TOP || cell.type == Type::BOTTOM || cell.type == Type::WELL_LAT);
    int beta_idx = -1;
    if (cell.type == Type::SIDE_LEFT)
        beta_idx = cell.num + cellsNum_z * (cellsNum_y + 2);
    else if (cell.type == Type::BOTTOM)
        beta_idx = cell.num + (cellsNum_y + 2);
    else if (cell.type == Type::TOP)
        beta_idx = cell.num - (cellsNum_y + 2);
    else if (cell.type == Type::WELL_LAT)
        beta_idx = cell.num + 1;

    const auto& beta = cells[beta_idx];
    const auto& props = *cell.props;
    const auto& next = x[cell.num];
    const auto& nebr = x[beta.num];

    TapeVariable res;
    res.p = (next.p - nebr.p + grav * props_w.dens_stc * (cell.z - beta.z)) / P_dim;
    //res.s = (next.s - nebr.s) / P_dim;
    return res;
}