#include "model/Acid/2drec/Acid2dRecModel.hpp"
#include <assert.h>
#include <numeric>

using namespace acid2drec;

double acid2drec::Component::R = 8.3144598;
double acid2drec::Component::p_std = 101325.0;

Acid2dRecModel::Acid2dRecModel()
{
	grav = 9.8;
	Volume = 0.0;
	injected_sol_volume = injected_acid_volume = 0.0;
}
Acid2dRecModel::~Acid2dRecModel()
{
	delete snapshotter;
	delete[] x, h;
}
void Acid2dRecModel::setProps(Properties& props)
{
    prefix = props.prefix;
	props_sk = props.props_sk;
	skeletonsNum = props.props_sk.size();

	rightBoundIsPres = props.rightBoundIsPres;
	hx = props.hx;	hy = props.hy;	hz = props.hz;

	// Setting grid properties
	cellsNum_x = props.cellsNum_x;
	cellsNum_y = props.cellsNum_y;
    R_dim = props.R_dim;

	cellsNum = (cellsNum_y + 2) * (cellsNum_x + 2);

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
	max_sol_volume = props.max_sol_volume;
	max_matrix_acid_volume = props.max_matrix_acid_volume;

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
void Acid2dRecModel::makeDimLess()
{
	T_dim = props_sk[0].t_init;
	//R_dim = props_frac.l2 / 8.0;
	t_dim = 1.0 * 3600.0;
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
		sk.beta /= (1.0 / P_dim);
		sk.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
		sk.p_init /= P_dim;
		sk.p_out /= P_dim;
		sk.p_ref /= P_dim;
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
	max_matrix_acid_volume /= R_dim * R_dim * R_dim;

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

	hx /= R_dim;	hy /= R_dim;	hz /= R_dim;
}
void Acid2dRecModel::buildGrid()
{
	//int counter = 0;
	//cells.reserve(cellsNum);
	//Volume = 0.0;
	//
	//double cur_hx = 0.0, cur_hy = 0.0, cur_hz = 0.0;
	//double cx = 0.0, cy = 0.0, cz = 0.0;
	//Type cur_type;

	////dist_x = props_frac.l2 / 100.0;
	////double x_prev = dist_x;
	////double x_logMax = log((props_frac.l2 + dist_x) / dist_x);
	////double x_logStep = x_logMax / (double)cellsNum_x;

	//// Left border
	//cur_type = Type::FRAC_IN;
	//hx = 0.0;	 cx = 0.0;
	//const int k = 1;
	///*for (int k = 0; k < cellsNum_z + 2; k++)
	//{
	//	if (k == 0)
	//		hz = cz = 0.0;
	//	else if (k == cellsNum_z + 1)
	//	{
	//		hz = 0.0;
	//		cz = props_frac.height;
	//	}
	//	else
	//	{*/
	//		//hz = props_frac.height / cellsNum_z;
	//		//cz += hz / 2.0;
	//		hz = props_sk[k - 1].height;
	//		cz = hz / 2.0;
	//	//}

	//	for (int i = 0; i < cellsNum_y_frac + 1; i++)
	//	{
	//		if (i == 0)
	//			hy = cy = 0.0;
	//		else
	//		{
	//			hy = props_frac.w2 / (double)cellsNum_y_frac;
	//			cy = ((double)i - 0.5) * hy;
	//		}

	//		cells_frac.push_back( FracCell(counter++, cx, cy, cz, hx, hy, hz, cur_type) );
	//		auto& cell = cells_frac.back();
	//	}
	////}
	//// Middle
	//hx = props_frac.l2 / cellsNum_x;	 cx = hx / 2;
	//for (int j = 0; j < cellsNum_x; j++)
	//{
	//	//cx = x_prev * (exp(x_logStep) + 1.0) / 2.0 - dist_x;
	//	//hx = x_prev * (exp(x_logStep) - 1.0);
	//	//x_prev *= exp(x_logStep);

	//	/*for (int k = 0; k < cellsNum_z + 2; k++)
	//	{
	//		if (k == 0)
	//			hz = cz = 0.0;
	//		else if (k == cellsNum_z + 1)
	//		{
	//			hz = 0.0;
	//			cz = props_frac.height;
	//		}
	//		else
	//		{*/
	//			//hz = props_frac.height / cellsNum_z;
	//			//cz += hz / 2.0;
	//			hz = props_sk[k - 1].height;
	//			cz = hz / 2.0;
	//		//}

	//		for (int i = 0; i < cellsNum_y_frac + 1; i++)
	//		{
	//			if (i == 0)
	//				hy = cy = 0.0;
	//			else
	//			{
	//				hy = props_frac.w2 / (double)cellsNum_y_frac;
	//				cy = ((double)i - 0.5) * hy;
	//			}

	//			if (i == cellsNum_y_frac)
	//				cur_type = FracType::FRAC_OUT;
	//			else if (i == 0)
	//				cur_type = FracType::FRAC_BORDER;
	//			else
	//				cur_type = FracType::FRAC_MID;

	//			//if(k == 0 || k == cellsNum_z + 1)
	//			//	cur_type = FracType::FRAC_BORDER;

	//			cells_frac.push_back(FracCell(counter++, cx, cy, cz, hx, hy, hz, cur_type));
	//			auto& cell = cells_frac.back();
	//			Volume_frac += cell.V;
	//		}
	//	//}
	//	cx += hx;
	//}
	//// Right cells
	//cur_type = FracType::FRAC_BORDER;
	//cx = props_frac.l2;	 hx = 0.0;
	///*for (int k = 0; k < cellsNum_z + 2; k++)
	//{
	//	if (k == 0)
	//		hz = cz = 0.0;
	//	else if (k == cellsNum_z + 1)
	//	{
	//		hz = 0.0;
	//		cz = props_frac.height;
	//	}
	//	else
	//	{*/
	//		//hz = props_frac.height / cellsNum_z;
	//		//cz += hz / 2.0;
	//		hz = props_sk[k - 1].height;
	//		cz = hz / 2.0;
	//	//}

	//	for (int i = 0; i < cellsNum_y_frac + 1; i++)
	//	{
	//		if (i == 0)
	//			hy = cy = 0.0;
	//		else
	//		{
	//			hy = props_frac.w2 / (double)cellsNum_y_frac;
	//			cy = ((double)i - 0.5) * hy;
	//		}

	//		cells_frac.push_back(FracCell(counter++, cx, cy, cz, hx, hy, hz, cur_type));
	//		auto& cell = cells_frac.back();
	//		//cell.nebrs[0] = cell.num - (cellsNum_y_frac + 1) * (cellsNum_z + 2);
	//	}
	////}
}
void Acid2dRecModel::processGeometry()
{
}
void Acid2dRecModel::setPerforated()
{
	height_perf = 0.0;
	vector<pair<int, int> >::iterator it;
	for (const auto& cell : cells)
	{
		if (cell.type == Type::SIDE_LEFT)
		{
			Qcell[cell.num] = 0.0;
			height_perf += cell.hx * cell.hz;
		}
	}
};
void Acid2dRecModel::setPeriod(int period)
{
/*	leftBoundIsRate = LeftBoundIsRate[period];
	if (leftBoundIsRate)
	{
		Q_sum = rate[period];
		if (Q_sum == 0.0)
		{
			for(auto& cell : cells)
				cell.u_prev.c = cell.u_iter.c = cell.u_next.c = 0.0;
		}

		if (period == 0 || rate[period - 1] < EQUALITY_TOLERANCE) 
		{
			std::map<int, double>::iterator it;
			for (it = Qcell.begin(); it != Qcell.end(); ++it)
				it->second = 3.0 * Q_sum * cells_frac[it->first].hy / 2.0 / props_frac.w2 * 
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
void Acid2dRecModel::setInitialState() 
{
	for (auto& cell : cells)
	{
		cell.props = &props_sk[getSkeletonId(cell)];
		cell.u_prev.m = cell.u_iter.m = cell.u_next.m = cell.props->m_init;
		cell.u_prev.p = cell.u_iter.p = cell.u_next.p = cell.props->p_init - grav * props_w.dens_stc * cell.z;
		
		if (cell.type == Type::WELL_LAT)
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

	x = new TapeVariable[cellsNum];
	h = new adouble[ var_size * cellsNum ];
}
double Acid2dRecModel::getRate(const int idx) const
{
/*	const Cell& cell = cells[idx];
	assert(cell.type == FracType::FRAC_IN);
	assert(cell.hy != 0.0 && cell.hz != 0.0);
	const FracCell& beta = cells_frac[cell.num + cellsNum_z * (cellsNum_y_frac + 1)];
	const FracVariable& next = cell.u_next;
	const FracVariable& nebr = beta.u_next;

	return cell.hy * cell.hz * props_frac.w2 * props_frac.w2 / props_w.visc * (nebr.p - next.p) / (cell.hx + beta.hx) * 
				(1.0 - (cell.y / props_frac.w2) * (cell.y / props_frac.w2));*/
	return 0.0;
};

TapeVariable Acid2dRecModel::solve(const Cell& cell, const Regime reg)
{
	if (cell.type == Type::MIDDLE)
		return solveMid(cell);
	else if (cell.type == Type::WELL_LAT)
		return solveLeft(cell, reg);
	else if (cell.type == Type::RIGHT)
		return solveRight(cell);
	else
		return solveBorder(cell);
}
TapeVariable Acid2dRecModel::solveMid(const Cell& cell)
{
	assert(cell.type == PoroType::MIDDLE);
	const auto& props = *cell.props;
	auto& next = x[cell.num];
	const auto& prev = cell.u_prev;
	adouble rate = getReactionRate(next, props);// / reac.indices[REACTS::ACID] / reac.comps[REACTS::ACID].mol_weight;

	TapeVariable res;
/*	adouble m = props.getPoro(next.m, next.p);
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
	getNeighborIdx(cell.num, neighbor);
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
		adouble buf_w = ht / cell.V * getPoroTrans(cell, next, beta, nebr) * (next.p - nebr.p + (cell.z - beta.z) * grav * props_w.dens_stc) *
			dens_w * props_w.getKr(upwd.sw, upwd.m, cells_poro[upwd_idx].props) / props_w.getViscosity(upwd.p, upwd.xa, upwd.xw, upwd.xs);
		adouble buf_o = ht / cell.V * getPoroTrans(cell, next, beta, nebr) * (next.p - nebr.p + (cell.z - beta.z) * grav * props_w.dens_stc) *
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
TapeVariable Acid2dRecModel::solveLeft(const Cell& cell, const Regime reg)
{
	assert(cell.type == PoroType::WELL_LAT);
	const auto& props = *cell.props;

	const auto& beta_poro1 = cells[cell.num + 1];
	const auto& beta_poro2 = cells[cell.num + 2];
	const auto& nebr_poro1 = x[cell.num + 1];
	const auto& nebr_poro2 = x[cell.num + 2];

	//const auto& beta1 = cells[getFracNebr(cell.num)];
	//const auto& beta2 = cells[beta1.num - 1];
	//const auto& beta3 = cells_frac[beta1.num - 2];
	//const auto& nebr1 = x[beta1.num];
	//const auto& nebr2 = x[beta2.num];
	//const auto& nebr3 = x_frac[beta3.num];

	auto& next = x[cell.num];
	/*const auto& prev = cell.u_prev;

	const double dist0 = cell.hy / 2.0;
	const double dist1 = beta1.hy / 2.0;
	const double dist2 = beta2.hy / 2.0;
	const double dist_poro1 = beta_poro1.hy / 2.0;
	const double dist_poro2 = beta_poro2.hy / 2.0;*/

	TapeVariable res;
	res.m = ((next.m - nebr_poro1.m) /*- (nebr_poro1.m - nebr_poro2.m) *
		(dist0 + dist_poro1) / (dist_poro1 + dist_poro2)*/) / P_dim * 2.0;
	/*res.p = ((next.p - nebr1.p) - (nebr1.p - nebr2.p) *
		(dist0 + dist1) / (dist1 + dist2)) / P_dim * 2.0;
	res.sw = (next.sw - (1.0 - props.s_oc)) / P_dim * 2.0;
	res.xw = (next.xw - (1.0 - next.xa)) / P_dim * 2.0;
	res.xs = next.xs / P_dim * 2.0;*/

	return res;
}
TapeVariable Acid2dRecModel::solveRight(const Cell& cell)
{
	assert(cell.type == PoroType::RIGHT);
	const auto& beta = cells[cell.num - 1];
	const auto& props = *cell.props;

	const auto& next = x[cell.num];
	const auto& nebr = x[beta.num];

	adouble rightIsPres = rightBoundIsPres;
	TapeVariable res;
	res.m = (props.getPoro(next.m, next.p) - props.getPoro(nebr.m, nebr.p)) / P_dim;
	condassign(res.p, rightIsPres, (next.p - props.p_out + grav * props_w.dens_stc * cell.z) / P_dim, (next.p - nebr.p) / P_dim);
	res.sw = (next.sw - nebr.sw) / P_dim;
	res.xw = (next.xw - nebr.xw) / P_dim;
	res.xa = (next.xa - nebr.xa) / P_dim;
	res.xs = (next.xs - nebr.xs) / P_dim;
	return res;
}
TapeVariable Acid2dRecModel::solveBorder(const Cell& cell)
{
	assert(cell.type == PoroType::SIDE_LEFT || cell.type == PoroType::SIDE_RIGHT || 
			cell.type == PoroType::TOP || cell.type == PoroType::BOTTOM);
	int beta_idx = -1;
	if (cell.type == Type::SIDE_LEFT)
		beta_idx = cell.num + cellsNum_y + 2;
	else if (cell.type == Type::SIDE_RIGHT)
		beta_idx = cell.num - cellsNum_y + 2;
	/*else if (cell.type == PoroType::BOTTOM)
		beta_idx = cell.num + (cellsNum_y_poro + 2);
	else if (cell.type == PoroType::TOP)
		beta_idx = cell.num - (cellsNum_y_poro + 2);*/
	const auto& beta = cells[beta_idx];
	const auto& props = *cell.props;

	const auto& next = x[cell.num];
	const auto& nebr = x[beta.num];

	adouble rightIsPres = rightBoundIsPres;
	TapeVariable res;
	/*res.m = (props.getPoro(next.m, next.p) - props.getPoro(nebr.m, nebr.p)) / P_dim;
	res.p = (next.p - nebr.p + grav * props_w.dens_stc * (cell.z - beta.z)) / P_dim;
	res.sw = (next.sw - nebr.sw) / P_dim;
	res.xw = (next.xw - nebr.xw) / P_dim;
	res.xa = (next.xa - nebr.xa) / P_dim;
	res.xs = (next.xs - nebr.xs) / P_dim;*/
	return res;
}