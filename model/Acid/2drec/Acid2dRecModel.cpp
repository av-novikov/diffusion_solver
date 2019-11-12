#include "model/Acid/2drec/Acid2dRecModel.hpp"
#include <assert.h>
#include "method/mcmath.h"
#include <numeric>
#include <random>

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
	fieldData = props.fieldData;
    prefix = props.prefix;
	props_sk = props.props_sk;

	rightBoundIsPres = props.rightBoundIsPres;
	hx = props.hx;	hy = props.hy;	hz = props.hz;

	// Setting grid properties
	cellsNum_x = props.cellsNum_x;
	cellsNum_y = props.cellsNum_y;
    R_dim = props.R_dim;

	cellsNum = (cellsNum_y + 2) * (cellsNum_x + 2);

	if(!fieldData)
		skeletonsNum = props.props_sk.size();
	else
	{
		skeletonsNum = cellsNum_x * cellsNum_y;
		for (int i = 1; i < skeletonsNum; i++)
		{
			props_sk.push_back(props_sk[0]);
		}
	}
	for (int j = 0; j < skeletonsNum; j++)
	{
		props_sk[j].id = j;
		props_sk[j].perm = MilliDarcyToM2(props_sk[j].perm);
	}

	periodsNum = props.LeftBoundIsRate.size();
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
	max_acid_volume = props.max_acid_volume;

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
	P_dim = 100.0 * 1.E+5;// props_sk[0].p_init;
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

	hx /= R_dim;	hy /= R_dim;	hz /= R_dim;
	Q_dim = R_dim * R_dim * R_dim / t_dim;
	for (int i = 0; i < periodsNum; i++)
	{
		period[i] /= t_dim;
		rate[i] /= Q_dim;
		pwf[i] /= P_dim;
	}
	max_sol_volume /= R_dim * R_dim * R_dim;
	max_acid_volume /= R_dim * R_dim * R_dim;

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
}
void Acid2dRecModel::buildGrid()
{
	Volume = hx * hy * hz;
	double cx = 0.0, cy = 0.0, cz = hz / 2.0;
	const double cell_hx = hx / (double)cellsNum_x, cell_hy = hy / (double)cellsNum_y, cell_hz = hz;
	int counter = 0;

	// x = 0
	cells.push_back(Cell(counter++, cx, cy, cz, 0.0, 0.0, cell_hz, Type::LEFT));
	cy = cell_hy / 2.0;
	for (int j = 0; j < cellsNum_y; j++)
	{
		cells.push_back(Cell(counter++, cx, cy, cz, 0.0, cell_hy, cell_hz, Type::LEFT));
		cy += cell_hy;
	}
	cy -= cell_hy / 2.0;
	cells.push_back(Cell(counter++, cx, cy, cz, 0.0, 0.0, cell_hz, Type::LEFT));
	// middle x
	cx = -cell_hx / 2.0;
	for (int i = 0; i < cellsNum_x; i++)
	{
		cx += cell_hx;
		cy = 0.0;
		cells.push_back(Cell(counter++, cx, cy, cz, cell_hx, 0.0, cell_hz, Type::BOTTOM));
		cy += cell_hy / 2.0;
		for (int j = 0; j < cellsNum_y; j++)
		{
			cells.push_back(Cell(counter++, cx, cy, cz, cell_hx, cell_hy, cell_hz, Type::MIDDLE));
			cy += cell_hy;
		}
		cy -= cell_hy / 2.0;
		cells.push_back(Cell(counter++, cx, cy, cz, cell_hx, 0.0, cell_hz, Type::TOP));
	}
	// x = hx
	cx = hx;
	cy = 0.0;
	cells.push_back(Cell(counter++, cx, cy, cz, 0.0, 0.0, cell_hz, Type::RIGHT));
	cy += cell_hy / 2.0;
	for (int j = 0; j < cellsNum_y; j++)
	{
		cells.push_back(Cell(counter++, cx, cy, cz, 0.0, cell_hy, cell_hz, Type::RIGHT));
		cy += cell_hy;
	}
	cy -= cell_hy / 2.0;
	cells.push_back(Cell(counter++, cx, cy, cz, 0.0, 0.0, cell_hz, Type::RIGHT));
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
		if (cell.type == Type::BOTTOM)
		{
			Qcell[cell.num] = 0.0;
			height_perf += cell.hx * cell.hz;
		}
	}
};
void Acid2dRecModel::setPeriod(int period)
{
	leftBoundIsRate = LeftBoundIsRate[period];
	if (leftBoundIsRate)
	{
		Q_sum = rate[period];
		if (Q_sum == 0.0)
		{
			for(auto& cell : cells)
				cell.u_prev.xa = cell.u_iter.xa = cell.u_next.xa = 0.0;
		}

		if (period == 0 || rate[period - 1] < EQUALITY_TOLERANCE) 
		{
			std::map<int, double>::iterator it;
			for (it = Qcell.begin(); it != Qcell.end(); ++it)
				it->second = 0.0;
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
void Acid2dRecModel::calcCorrelatedPermeability(const double sigma, const double lam)
{
	std::random_device rd{};
	std::mt19937 gen{ rd() };
	std::normal_distribution<double> d{ 0, sigma };

	const int size = cellsNum_x * cellsNum_y;
	MCMatrix covar(size, size);
	MCVector ran(size), res(size);
	double dist, cur_cov;
	int i1, j1, i2, j2;
	for (int i = 0; i < size; i++)
	{
		i1 = i / cellsNum_y;
		j1 = i % cellsNum_y;
		const Cell& cell1 = cells[j1 + 1 + (i1 + 1) * (cellsNum_y + 2)];
		for (int j = i; j < size; j++)
		{
			i2 = j / cellsNum_y;
			j2 = j % cellsNum_y;
			const Cell& cell2 = cells[j2 + 1 + (i2 + 1) * (cellsNum_y + 2)];
			assert(cell1.V > 0.0 && cell2.V > 0.0);
 			dist = sqrt((cell2.x - cell1.x) * (cell2.x - cell1.x) +
						(cell2.y - cell1.y) * (cell2.y - cell1.y));
			cur_cov = exp(-dist / lam);
			covar[i][j] = covar[j][i] = cur_cov;
		}
		ran[i] = d(gen);
	}

	MC_Cholesky chol(covar);
	chol.Cholesky_Decomposition();
	assert(chol.isSPD());
	if (chol.isSPD())
		res = chol.Cholesky_Matrix.Transpose() * ran;

	for (int i = 0; i < size; i++)
	{
		auto& sk = props_sk[i];
		sk.perm *= exp(res[i]);
		if (sk.perm < 0.0)
			sk.perm = 0.0;
	}
}
void Acid2dRecModel::setInitialState() 
{
	if (fieldData)
	{
		const double perm_av = props_sk[0].perm;
		calcCorrelatedPermeability(0.5, 0.05 / R_dim);
	}

	for (auto& cell : cells)
	{
		const int prop_id = getSkeletonId(cell);
		auto& sk = props_sk[prop_id];
		cell.props = &props_sk[prop_id];
		cell.u_prev.m = cell.u_iter.m = cell.u_next.m = cell.props->m_init;
		cell.u_prev.p = cell.u_iter.p = cell.u_next.p = cell.props->p_init;
		cell.u_prev.sw = cell.u_iter.sw = cell.u_next.sw = cell.props->sw_init;
		cell.u_prev.xa = cell.u_iter.xa = cell.u_next.xa = cell.props->xa_init;
		cell.u_prev.xw = cell.u_iter.xw = cell.u_next.xw = cell.props->xw_init;
		cell.u_prev.xs = cell.u_iter.xs = cell.u_next.xs = 0.0;
	}

	x = new TapeVariable[ cellsNum ];
	h = new adouble[ var_size * cellsNum ];
}
double Acid2dRecModel::getRate(const int idx) const
{
	const Cell& cell = cells[idx];
	assert(cell.type == Type::BOTTOM);
	assert(cell.hx != 0.0 && cell.hz != 0.0);
	const Cell& beta = cells[cell.num + 1];
	const auto& next = x[cell.num];
	const auto& nebr = x[beta.num];
	return getTrans(cell, next, beta, nebr).value() / props_w.getViscosity(next.p, next.xa, next.xw, next.xs).value() * (next.p.value() - nebr.p.value());
};

TapeVariable Acid2dRecModel::solve(const Cell& cell, const Regime reg)
{
	if (cell.type == Type::MIDDLE)
		return solveMid(cell);
	else if (cell.type == Type::BOTTOM)
		return solveLeft(cell, reg);
	else if (cell.type == Type::TOP)
		return solveRight(cell);
	else
		return solveBorder(cell);
}
TapeVariable Acid2dRecModel::solveMid(const Cell& cell)
{
	assert(cell.type == Type::MIDDLE);
	const auto& props = *cell.props;
	auto& next = x[cell.num];
	const auto& prev = cell.u_prev;
	adouble rate = getReactionRate(next, props);

	TapeVariable res;
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
	getNeighborIdx(cell.num, neighbor);
	double dist1, dist2;
	for (size_t i = 0; i < NEBRS_NUM; i++)
	{
		const Cell& beta = cells[neighbor[i]];
		const auto& nebr = x[beta.num];
		if (i < 2)
		{
			dist1 = cell.hy / 2.0;
			dist2 = beta.hy / 2.0;
		}
		else
		{
			dist1 = cell.hx / 2.0;
			dist2 = beta.hx / 2.0;
		}
		upwd_idx = getUpwindIdx(cell, beta);
		const auto& upwd = x[upwd_idx];

		adouble dens_w = getAverage(props_w.getDensity(next.p, next.xa, next.xw, next.xs), dist1, 
									props_w.getDensity(nebr.p, nebr.xa, nebr.xw, nebr.xs), dist2);
		adouble dens_o = getAverage(props_o.getDensity(next.p), dist1,
									props_o.getDensity(nebr.p), dist2);
		adouble buf_w = ht / cell.V * getTrans(cell, next, beta, nebr) * (next.p - nebr.p) *
			dens_w * props_w.getKr(upwd.sw, upwd.m, cells[upwd_idx].props) / props_w.getViscosity(upwd.p, upwd.xa, upwd.xw, upwd.xs);
		adouble buf_o = ht / cell.V * getTrans(cell, next, beta, nebr) * (next.p - nebr.p) *
			dens_o * props_o.getKr(upwd.sw, upwd.m, cells[upwd_idx].props) / props_o.getViscosity(upwd.p);

		res.p += buf_w;
		res.sw += buf_o;
		res.xw += buf_w * upwd.xw;
		res.xa += buf_w * upwd.xa;
		res.xs += buf_w * upwd.xs;
	}
	return res;
}
TapeVariable Acid2dRecModel::solveLeft(const Cell& cell, const Regime reg)
{
	assert(cell.type == Type::BOTTOM);
	const auto& props = *cell.props;
	auto& next = x[cell.num];
	const auto& beta = cells[cell.num + 1];
	auto& nebr = x[beta.num];

	TapeVariable res;
	res.m = ((props.getPoro(next.m, next.p) - props.getPoro(nebr.m, nebr.p)) /*- (nebr_poro1.m - nebr_poro2.m) *
		(dist0 + dist_poro1) / (dist_poro1 + dist_poro2)*/) / P_dim;
	//res.p = ((next.p - nebr1.p) - (nebr1.p - nebr2.p) *
	//	(dist0 + dist1) / (dist1 + dist2)) / P_dim * 2.0;
	//res.sw = (next.sw - (1.0 - props.s_oc)) / P_dim * 2.0;
	//res.xw = (next.xw - (1.0 - next.xa)) / P_dim * 2.0;
	//res.xs = next.xs / P_dim * 2.0;
	if (leftBoundIsRate)
	{
		res.p = -Q_sum;
		for (const auto& it : Qcell)
		{
			const auto& p_cell = cells[it.first];
			auto& p_next = x[it.first];
			const auto& p_beta = cells[it.first + 1];
			auto& p_nebr = x[it.first + 1];
			res.p += getTrans(p_cell, p_next, p_beta, p_nebr) / props_w.getViscosity(next.p, next.xa, next.xw, next.xs) * (next.p - nebr.p);
		}
	}
	else
	{
		res.p = (next.p - Pwf) / P_dim;
	}
	res.sw = (next.sw - (1.0 - props.s_oc)) / P_dim;
	res.xw = (next.xw - (1.0 - next.xa)) / P_dim;
	res.xa = (next.xa - c) / P_dim;
	res.xs = next.xs / P_dim;

	return res;
}
TapeVariable Acid2dRecModel::solveRight(const Cell& cell)
{
	assert(cell.type == Type::TOP);
	const auto& beta = cells[cell.num - 1];
	const auto& props = *cell.props;

	const auto& next = x[cell.num];
	const auto& nebr = x[beta.num];

	adouble rightIsPres = rightBoundIsPres;
	TapeVariable res;
	res.m = (props.getPoro(next.m, next.p) - props.getPoro(nebr.m, nebr.p)) / P_dim;
	if (rightBoundIsPres)
		res.p = (next.p - props.p_out) / P_dim;
	else
		res.p = (next.p - nebr.p) / P_dim;
	res.sw = (next.sw - nebr.sw) / P_dim;
	res.xw = (next.xw - nebr.xw) / P_dim;
	res.xa = (next.xa - nebr.xa) / P_dim;
	res.xs = (next.xs - nebr.xs) / P_dim;
	return res;
}
TapeVariable Acid2dRecModel::solveBorder(const Cell& cell)
{
	assert(cell.type == Type::LEFT || cell.type == Type::RIGHT);
	int beta_idx = -1;
	if (cell.type == Type::LEFT)
		beta_idx = cell.num + cellsNum_y + 2;
	else if (cell.type == Type::RIGHT)
		beta_idx = cell.num - cellsNum_y - 2;
	const auto& beta = cells[beta_idx];
	const auto& props = *cell.props;

	const auto& next = x[cell.num];
	const auto& nebr = x[beta.num];

	TapeVariable res;
	res.m = (props.getPoro(next.m, next.p) - props.getPoro(nebr.m, nebr.p)) / P_dim;
	res.p = (next.p - nebr.p) / P_dim;
	res.sw = (next.sw - nebr.sw) / P_dim;
	res.xw = (next.xw - nebr.xw) / P_dim;
	res.xa = (next.xa - nebr.xa) / P_dim;
	res.xs = (next.xs - nebr.xs) / P_dim;
	return res;
}