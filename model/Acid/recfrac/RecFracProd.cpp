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
	cellsNum = (cellsNum_x + 2) * (cellsNum_y + 2) * (cellsNum_z + 2);
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
	double max_t = 86400.0 * 365.0;
	for (int i = 0; i < periodsNum; i++)
	{
        LeftBoundIsRate.push_back(false);// props.LeftBoundIsRate[i]);
        period.push_back(max_t);// props.timePeriods[i]);
        pwf[i] = 0.5 * props_sk.back().p_out;// .pwf[pres_idx++];
        rate[i] = 0.0;
	}

	// Temporal properties
    //ht = 10000.0;// props.ht;
    ht = 2000.0;
    ht_min = 10000.0;// props.ht_min;
    ht_max = 10000000.0;// props.ht_max;

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

    auto cur_type = Type::MIDDLE;
    for (int j = 0; j < cellsNum_x + 2; j++)
    {
        if (j < prev_cellsNum_x + 1)
        {
            cx = cells_poro[j * (prev_cellsNum_z + 2) * (prev_cellsNum_y + 2)].x;
            cur_hx = cells_poro[j * (prev_cellsNum_z + 2) * (prev_cellsNum_y + 2)].hx;
        }
        else
        {
            cx = cells_poro[(prev_cellsNum_x + 1) * (prev_cellsNum_z + 2) * (prev_cellsNum_y + 2)].x +
                (double)(j - (prev_cellsNum_x + 1) + 0.5) * hx;
            cur_hx = (j == cellsNum_x + 1) ? 0.0 : hx;
        }                                            


        for (int k = 0; k < cellsNum_z + 2; k++)
        {
            cz = cells_poro[k * (prev_cellsNum_y + 2)].z;
            cur_hz = cells_poro[k * (prev_cellsNum_y + 2)].hz;
            cy = 0.0;

            for (int i = 0; i < cellsNum_y + 2; i++)
            {
                cur_type = (j == 0) ? Type::SIDE_LEFT : Type::MIDDLE;
                cur_type = (j == cellsNum_x + 1) ? Type::SIDE_RIGHT : cur_type;
                cur_type = (k == 0) ? Type::BOTTOM : cur_type;
                cur_type = (k == cellsNum_z + 1) ? Type::TOP : cur_type;
                cur_type = (i == 0 && cur_type == Type::MIDDLE) ? Type::WELL_LAT : cur_type;
                cur_type = (i == cellsNum_y + 1 && cur_type == Type::MIDDLE) ? Type::RIGHT : cur_type;
                cur_type = (i < num_input_cells + 1 && cur_type == Type::SIDE_LEFT) ? Type::FRAC_IN : cur_type;
                if (j < prev_cellsNum_x + 1 && i < prev_cellsNum_y + 1)
                {
                    cells.push_back(Cell(cells_poro[i + (prev_cellsNum_y + 2) * (k + j * (prev_cellsNum_z + 2))], counter++, cur_type));
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
        }
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
void RecFracProd::setInitialState(const std::vector<PoroCell>& cells_poro)
{
    for (int j = 0; j < cellsNum_x + 2; j++)
    {
        for (int k = 0; k < cellsNum_z + 2; k++)
        {
            for (int i = 0; i < cellsNum_y + 2; i++)
            {
                auto& cell = cells[i + (cellsNum_y + 2) * (k + j * (cellsNum_z + 2))];
                cell.props = &props_sk[getSkeletonId(cell)];
                cell.u_prev.p = cell.u_iter.p = cell.u_next.p = cell.props->p_init - grav * props_w.dens_stc * cell.z;
				//cell.u_prev.s = cell.u_iter.s = cell.u_next.s = cell.props->sw_init;
                if (j < prev_cellsNum_x + 1 && i < prev_cellsNum_y + 1)
                {
                    auto& cell0 = cells_poro[i + (prev_cellsNum_y + 2) * (k + j * (prev_cellsNum_z + 2))];
                    cell.u_prev.m = cell.u_iter.m = cell.u_next.m = cell0.u_next.m;
                    if (i < num_input_cells + 1)
                        cell.u_prev.m = cell.u_iter.m = cell.u_next.m = 0.3;
                }
                else
                {
                    cell.u_prev.m = cell.u_iter.m = cell.u_next.m = cell.props->m_init;
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
	const Cell& beta = cells[cell.num + (cellsNum_z + 2) * (cellsNum_y + 2)];
	const Variable& next = cell.u_next;
	const Variable& nebr = beta.u_next;
	//double kr = props_o.getKr(next.s, next.m, cell.props).value();
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
    const auto& beta = cells[cell.num + (cellsNum_z + 2) * (cellsNum_y + 2)];
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
        beta_idx = cell.num - (cellsNum_z + 2) * (cellsNum_y + 2);
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
        beta_idx = cell.num + (cellsNum_z + 2) * (cellsNum_y + 2);
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