#include "model/GasOil_RZ_NIT/GasOil_RZ_NIT.h"
#include "util/utils.h"

using namespace std;
using namespace gasOil_rz_NIT;

GasOil_RZ_NIT::GasOil_RZ_NIT()
{
}

GasOil_RZ_NIT::~GasOil_RZ_NIT()
{
}

void GasOil_RZ_NIT::setProps(Properties& props)
{
	// Setting initial condition
	varInit.t = props.T_init;
	varInit.p = props.p_init;
	varInit.s = props.s_init;
	varInit.p_bub = props.p_sat;

	if(props.p_sat > props.p_init)
		varInit.SATUR = true;
	else
		varInit.SATUR = false;

	// Setting grid properties
	r_w = props.r_w;
	r_e = props.r_e;
	cellsNum_r = props.cellsNum_r;
	cellsNum_z = props.cellsNum_z;
	cellsNum = (cellsNum_r + 2) * (cellsNum_z + 2);

	// Setting skeleton properties
	perfIntervals = props.perfIntervals;
	props_sk.m = props.m;
	for(int i = 0; i < cellsNum_z+2; i++)
		props_sk.perm_r.push_back ( MilliDarcyToM2( props.perm_r[i] ) );
	for(int i = 0; i < cellsNum_z+2; i++)
		props_sk.perm_z.push_back ( MilliDarcyToM2( props.perm_z[i] ) );
	props_sk.dens_stc = props.dens_sk_stc;
	props_sk.beta = props.beta_sk;
	props_sk.height = props.h2 - props.h1;
	props_sk.h1 = props.h1;
	props_sk.h2 = props.h2;

	props_sk.c = props.c_sk;
	props_sk.lambda_r = props.lambda_sk_r;
	props_sk.lambda_z = props.lambda_sk_z;

	periodsNum = props.timePeriods.size();
	props_sk.skin.resize(cellsNum_z+2);
	props_sk.r_eff.resize(cellsNum_z+2);
	props_sk.perm_eff.resize(cellsNum_z+2);
	for(int i = 0; i < periodsNum; i++)
	{
		period.push_back( props.timePeriods[i] );
		rate.push_back( props.rates[i] / 86400.0 );
		for(int j = 0; j < cellsNum_z+2; j++)
		{
			props_sk.skin[j].push_back( props.skins[i] );
			props_sk.r_eff[j].push_back( props.radius[i] );
			if(props.radius[i] > props.r_w)
				props_sk.perm_eff[j].push_back( MilliDarcyToM2(props.perm_r[j] * log(props.radius[i] / props.r_w) / (log(props.radius[i] / r_w) + props.skins[i]) ) );
			else
				props_sk.perm_eff[j].push_back( MilliDarcyToM2(props.perm_r[j]) );
		}
	}
	r_eff.resize(cellsNum_z+2);
	Perm_eff.resize(cellsNum_z+2);
	skin.resize(cellsNum_z+2);

	// Temporal properties
	ht = props.ht;
	ht_min = props.ht_min;
	ht_max = props.ht_max;

	// Oil properties
	props_oil.visc = cPToPaSec( props.visc_oil );
	props_oil.dens_stc = props.dens_oil_stc;
	props_oil.b_bore = props.b_oil_bore;
	props_oil.beta = props.beta_oil;
	
	props_oil.c = props.c_oil;
	props_oil.lambda = props.lambda_oil;
	props_oil.jt = props.jt_oil;
	props_oil.ad = props.ad_oil;

	// Gas properties
	props_gas.visc = cPToPaSec( props.visc_gas );
	props_gas.dens_stc = props.dens_gas_stc;
	
	props_gas.c = props.c_gas;
	props_gas.lambda = props.lambda_gas;
	props_gas.jt = props.jt_gas;
	props_gas.ad = props.ad_gas;

	alpha = props.alpha;
	depth_point = props.depth_point;
	L = props.L;

	makeDimLess();
	
	// Data sets
	props_oil.kr = setDataset(props.kr_oil, 1.0, 1.0);
	props_gas.kr = setDataset(props.kr_gas, 1.0, 1.0);

	props_oil.b = setDataset(props.B_oil, P_dim / BAR_TO_PA, 1.0);
	props_gas.b = setDataset(props.B_gas, P_dim / BAR_TO_PA, 1.0);
	Rs = setDataset(props.Rs, P_dim / BAR_TO_PA, 1.0);
	Prs = setInvDataset(props.Rs, 1.0, P_dim / BAR_TO_PA);
}

void GasOil_RZ_NIT::makeDimLess()
{
	// Main units
	R_dim = r_w;
	t_dim = 3600.0;
	P_dim = BAR_TO_PA;
	T_dim = varInit.t;

	// Temporal properties
	ht = ht / t_dim;
	ht_min = ht_min / t_dim;
	ht_max = ht_max / t_dim;

	// Initial condition
	varInit.t /= T_dim;
	varInit.p /= P_dim;
	varInit.p_bub /= P_dim;

	// Grid properties
	r_w = r_w / R_dim;
	r_e = r_e / R_dim;

	// Skeleton properties
	for(int i = 0; i < cellsNum_z+2; i++)
		props_sk.perm_r[i] = props_sk.perm_r[i] / R_dim / R_dim;
	for(int i = 0; i < cellsNum_z+2; i++)
		props_sk.perm_z[i] = props_sk.perm_z[i] / R_dim / R_dim;

	props_sk.dens_stc = props_sk.dens_stc / P_dim / t_dim / t_dim * R_dim * R_dim;
	props_sk.beta = props_sk.beta * P_dim;
	props_sk.height = props_sk.height / R_dim;
	props_sk.h1 = (props_sk.h1 - depth_point) / R_dim;
	props_sk.h2 = (props_sk.h2 - depth_point) / R_dim;

	props_sk.c = props_sk.c / R_dim / R_dim * T_dim * t_dim * t_dim;
	props_sk.lambda_r = props_sk.lambda_r * T_dim * t_dim / P_dim / R_dim / R_dim;
	props_sk.lambda_z = props_sk.lambda_z * T_dim * t_dim / P_dim / R_dim / R_dim;

	Q_dim = R_dim * R_dim * R_dim / t_dim;
	for(int i = 0; i < periodsNum; i++)
	{
		period[i] = period[i] / t_dim;
		rate[i] = rate[i] / Q_dim;
	}
	for(int j = 0; j < cellsNum_z+2; j++)
		for(int i = 0; i < periodsNum; i++)
		{
			props_sk.r_eff[j][i] = props_sk.r_eff[j][i] / R_dim;
			props_sk.perm_eff[j][i] = props_sk.perm_eff[j][i] / R_dim / R_dim;
		}

	// Oil properties
	props_oil.visc = props_oil.visc / (P_dim * t_dim);
	props_oil.dens_stc = props_oil.dens_stc / P_dim / t_dim / t_dim * R_dim * R_dim;
	props_oil.beta = props_oil.beta * P_dim;
	
	props_oil.c = props_oil.c / R_dim / R_dim * T_dim * t_dim * t_dim;
	props_oil.lambda = props_oil.lambda * T_dim * t_dim / P_dim / R_dim / R_dim;
	props_oil.jt = props_oil.jt * P_dim / T_dim;
	props_oil.ad = props_oil.ad * P_dim / T_dim;

	// Gas properties
	props_gas.visc = props_gas.visc / (P_dim * t_dim);
	props_gas.dens_stc = props_gas.dens_stc / P_dim / t_dim / t_dim * R_dim * R_dim;
	
	props_gas.c = props_gas.c / R_dim / R_dim * T_dim * t_dim * t_dim;
	props_gas.lambda = props_gas.lambda * T_dim * t_dim / P_dim / R_dim / R_dim;
	props_gas.jt = props_gas.jt * P_dim / T_dim;
	props_gas.ad = props_gas.ad * P_dim / T_dim;

	// Rest properties
	alpha = alpha / t_dim;
	//depth_point = 0.0;
	L = L / R_dim / R_dim * t_dim * t_dim;
}


void GasOil_RZ_NIT::buildGridLog()
{
	cells.reserve( cellsNum );

	Volume = 0.0;
	int counter = 0;

	double r_prev = r_w;
	double logMax = log(r_e / r_w);
	double logStep = logMax / (double)cellsNum_r;
	
	double hz = (props_sk.h2 - props_sk.h1) / (double)cellsNum_z;
	double cm_z = props_sk.h1;
	double hr = r_prev * (exp(logStep) - 1.0);
	double cm_r = r_w;

	// Left border
	cells.push_back( Cell(counter++, cm_r, cm_z, 0.0, 0.0) );
	cm_z += hz / 2.0;
	for(int i = 0; i < cellsNum_z; i++)
	{
		cells.push_back( Cell(counter++, cm_r, cm_z, 0.0, hz) );
		cm_z += hz;
	}
	cells.push_back( Cell(counter++, cm_r, cm_z-hz/2.0, 0.0, 0.0) );

	// Middle cells
	for(int j = 0; j < cellsNum_r; j++)
	{
		cm_z = props_sk.h1;
		cm_r = r_prev * (exp(logStep) + 1.0) / 2.0;
		hr = r_prev * (exp(logStep) - 1.0);

		cells.push_back( Cell(counter++, cm_r, cm_z, hr, 0.0) );
		cm_z += hz / 2.0;
		for(int i = 0; i < cellsNum_z; i++)
		{
			cells.push_back( Cell(counter++, cm_r, cm_z, hr, hz) );
			Volume += cells[cells.size()-1].V;
			cm_z += hz;
		}
		cells.push_back( Cell(counter++, cm_r, cm_z-hz/2.0, hr, 0.0) );

		r_prev = r_prev * exp(logStep);
	}

	// Right border
	cm_z = props_sk.h1;
	cm_r = r_e;
	
	cells.push_back( Cell(counter++, cm_r, cm_z, 0.0, 0.0) );
	cm_z += hz / 2.0;
	for(int i = 0; i < cellsNum_z; i++)
	{
		cells.push_back( Cell(counter++, cm_r, cm_z, 0.0, hz) );
		cm_z += hz;
	}
	cells.push_back( Cell(counter++, cm_r, cm_z-hz/2.0, 0.0, 0.0) );
}

void GasOil_RZ_NIT::setPerforated()
{
	height_perf = 0.0;
	vector<pair<int,int> >::iterator it;
	for(it = perfIntervals.begin(); it != perfIntervals.end(); ++it)
	{
		for(int i = it->first; i <= it->second; i++)
		{
			Qcell[i] = 0.0;
			height_perf += cells[i].hz;
		}
	}
}

void GasOil_RZ_NIT::setPeriod(int period)
{
	Q_sum = rate[period];

	if(period == 0 || rate[period-1] < EQUALITY_TOLERANCE ) {
		map<int,double>::iterator it;
		for(it = Qcell.begin(); it != Qcell.end(); ++it)
			it->second = Q_sum * cells[ it->first ].hz / height_perf;
	} else {
		map<int,double>::iterator it;
		for(it = Qcell.begin(); it != Qcell.end(); ++it)
			it->second = it->second * Q_sum / rate[period-1];
	}

	for(int i = 0; i < cellsNum_z+2; i++)
	{
		r_eff[i] = props_sk.r_eff[i][period];
		Perm_eff[i] = props_sk.perm_eff[i][period];
		skin[i] = props_sk.skin[i][period];
	}
}

void GasOil_RZ_NIT::setRateDeviation(int num, double ratio)
{
	Qcell[num] += Q_sum * ratio;
}

double GasOil_RZ_NIT::solve_eq1(int cur)
{
	int neighbor [4];
	getNeighborIdx(cur, neighbor);

	Cell& cell = cells[cur];
	Var2phaseNIT& next = cell.u_next;
	Var2phaseNIT& prev = cell.u_prev;
	double H = 0.0;
	
	H = ( getPoro(next.p) * next.s / getB_oil(next.p, next.p_bub, next.SATUR) - 
				getPoro(prev.p) * prev.s / getB_oil(prev.p, prev.p_bub, prev.SATUR) );

	for(int i = 0; i < 4; i++)
	{
		Var2phaseNIT& upwd = cells[ getUpwindIdx(cur, neighbor[i]) ].u_next;

		H += ht / cell.V * getTrans(cur, neighbor[i]) * (next.p - cells[ neighbor[i] ].u_next.p) *
			getKr_oil(upwd.s) / props_oil.visc / getB_oil(upwd.p, upwd.p_bub, upwd.SATUR);
	}

	return H;
}

double GasOil_RZ_NIT::solve_eq1_dp(int cur)
{
	double upwind;
	int neighbor [4];
	getNeighborIdx(cur, neighbor);

	Cell& cell = cells[cur];
	Var2phaseNIT& next = cell.u_next;
	double Boil_upwd;
	double Boil = getB_oil(next.p, next.p_bub, next.SATUR);
	double H = 0.0;

	H = (next.s * getPoro_dp(next.p) - 
		getPoro(next.p) * next.s * getB_oil_dp(next.p, next.p_bub, next.SATUR) / Boil ) / Boil;

	for(int i = 0; i < 4; i++)
	{
		upwind = upwindIsCur(cur, neighbor[i]);
		Var2phaseNIT& upwd = cells[ getUpwindIdx(cur, neighbor[i]) ].u_next;
		Boil_upwd = getB_oil(upwd.p, upwd.p_bub, upwd.SATUR);

		H += ht / cell.V * getTrans(cur, neighbor[i]) * 
			( getKr_oil(upwd.s) / props_oil.visc / Boil_upwd - 
			upwind * (next.p - cells[ neighbor[i] ].u_next.p) * getKr_oil(upwd.s) / props_oil.visc / Boil_upwd / Boil_upwd * getB_oil_dp(upwd.p, upwd.p_bub, upwd.SATUR) );
	}
	return H;
}

double GasOil_RZ_NIT::solve_eq1_ds(int cur)
{
	double upwind;	
	int neighbor [4];
	getNeighborIdx(cur, neighbor);

	Cell& cell = cells[cur];
	Var2phaseNIT& next = cell.u_next;
	double H = 0.0;

	H = getPoro(next.p) / getB_oil(next.p, next.p_bub, next.SATUR);

	for(int i = 0; i < 4; i++)
	{
		upwind = upwindIsCur(cur, neighbor[i]);
		Var2phaseNIT& upwd = cells[ getUpwindIdx(cur, neighbor[i]) ].u_next;

		H += ht / cell.V * getTrans(cur, neighbor[i]) * 
			upwind * (next.p - cells[ neighbor[i] ].u_next.p) * getKr_oil_ds(upwd.s) / props_oil.visc / getB_oil(upwd.p, upwd.p_bub, upwd.SATUR);
	}

	return H;
}

double GasOil_RZ_NIT::solve_eq1_dp_beta(int cur, int beta)
{
	Cell& cell = cells[cur];

	double upwind = upwindIsCur(cur, beta);
	Var2phaseNIT& upwd = cells[ getUpwindIdx(cur, beta) ].u_next;
	double Boil_upwd = getB_oil(upwd.p, upwd.p_bub, upwd.SATUR);

	return -ht / cell.V * getTrans(cur, beta) * 
			( getKr_oil(upwd.s) / props_oil.visc / Boil_upwd + 
			(1.0 - upwind) * (cell.u_next.p - cells[beta].u_next.p) * getKr_oil(upwd.s) / props_oil.visc / Boil_upwd / Boil_upwd * getB_oil_dp(upwd.p, upwd.p_bub, upwd.SATUR) );
}

double GasOil_RZ_NIT::solve_eq1_ds_beta(int cur, int beta)
{
	Cell& cell = cells[cur];

	double upwind = upwindIsCur(cur, beta);
	Var2phaseNIT& upwd = cells[ getUpwindIdx(cur, beta) ].u_next;

	return ht / cell.V * getTrans(cur, beta) * (1.0 - upwind) * (cell.u_next.p - cells[beta].u_next.p) *
			getKr_oil_ds(upwd.s) / props_oil.visc / getB_oil(upwd.p, upwd.p_bub, upwd.SATUR);
}

double GasOil_RZ_NIT::solve_eq2(int cur)
{
	int neighbor [4];
	getNeighborIdx(cur, neighbor);

	Cell& cell = cells[cur];
	Var2phaseNIT& next = cell.u_next;
	Var2phaseNIT& prev = cell.u_prev;
	double H = 0.0;

	H = getPoro(next.p) * ( (1.0 - next.s) / getB_gas(next.p) + next.s * getRs(next.p, next.p_bub, next.SATUR) / getB_oil(next.p, next.p_bub, next.SATUR) ) -
				getPoro(prev.p) * ( (1.0 - prev.s) / getB_gas(prev.p) + prev.s * getRs(prev.p, prev.p_bub, prev.SATUR) / getB_oil(prev.p, prev.p_bub, prev.SATUR) );

	for(int i = 0; i < 4; i++)
	{
		Var2phaseNIT& upwd = cells[ getUpwindIdx(cur, neighbor[i]) ].u_next;

		H += ht / cell.V * getTrans(cur, neighbor[i]) * (next.p - cells[ neighbor[i] ].u_next.p) * 
			( getKr_oil(upwd.s) * getRs(upwd.p, upwd.p_bub, upwd.SATUR) / props_oil.visc / getB_oil(upwd.p, upwd.p_bub, upwd.SATUR) +
			getKr_gas(upwd.s) / props_gas.visc / getB_gas(upwd.p) );
	}

	return H;
}

double GasOil_RZ_NIT::solve_eq2_dp(int cur)
{
	double upwind;	
	int neighbor [4];
	getNeighborIdx(cur, neighbor);

	Cell& cell = cells[cur];
	Var2phaseNIT& next = cell.u_next;
	double Boil_upwd, Bgas_upwd, rs_upwd;
	double Boil = getB_oil(next.p, next.p_bub, next.SATUR);
	double Bgas = getB_gas(next.p);
	double rs = getRs(next.p, next.p_bub, next.SATUR);
	double H = 0.0;

	H = ( (next.s * rs / Boil + (1.0 - next.s) / Bgas) * getPoro_dp(next.p) - 
		getPoro(next.p) * ( (1.0 - next.s) / Bgas / Bgas * getB_gas_dp(next.p) + 
		next.s * rs / Boil / Boil * getB_oil_dp(next.p, next.p_bub, next.SATUR) - 
		next.s / Boil * getRs_dp(next.p, next.p_bub, next.SATUR) ) );

	for(int i = 0; i < 4; i++)
	{
		upwind = upwindIsCur(cur, neighbor[i]);
		Var2phaseNIT& upwd = cells[ getUpwindIdx(cur, neighbor[i]) ].u_next;
		Boil_upwd = getB_oil(upwd.p, upwd.p_bub, upwd.SATUR);
		Bgas_upwd = getB_gas(upwd.p);
		rs_upwd = getRs(upwd.p, upwd.p_bub, upwd.SATUR);

		H += ht / cell.V * getTrans(cur, neighbor[i]) * 
			( getKr_oil(upwd.s) * rs_upwd / props_oil.visc / Boil_upwd + getKr_gas(upwd.s) / props_gas.visc / Bgas_upwd + 
			upwind * (next.p - cells[ neighbor[i] ].u_next.p) * 
			( getKr_oil(upwd.s) / props_oil.visc / Boil_upwd * (getRs_dp(upwd.p, upwd.p_bub, upwd.SATUR) - rs_upwd * getB_oil_dp(upwd.p, upwd.p_bub, upwd.SATUR) / Boil_upwd) - 
			getKr_gas(upwd.s) / props_gas.visc / Bgas_upwd / Bgas_upwd * getB_gas_dp(upwd.p) ));
	}

	return H;
}

double GasOil_RZ_NIT::solve_eq2_ds(int cur)
{
	double upwind;	
	int neighbor [4];
	getNeighborIdx(cur, neighbor);

	Cell& cell = cells[cur];
	Var2phaseNIT& next = cell.u_next;
	double H = 0.0;

	H = getPoro(next.p) * ( getRs(next.p, next.p_bub, next.SATUR) / getB_oil(next.p, next.p_bub, next.SATUR) - 1.0 / getB_gas(next.p) );
	for(int i = 0; i < 4; i++)
	{
		upwind = upwindIsCur(cur, neighbor[i]);
		Var2phaseNIT& upwd = cells[ getUpwindIdx(cur, neighbor[i]) ].u_next;

		H += ht / cell.V * getTrans(cur, neighbor[i]) * upwind * (next.p - cells[ neighbor[i] ].u_next.p) * 
			( getRs(upwd.p, upwd.p_bub, upwd.SATUR) * getKr_oil_ds(upwd.s) / props_oil.visc / getB_oil(upwd.p, upwd.p_bub, upwd.SATUR) + 
			getKr_gas_ds(upwd.s) / props_gas.visc / getB_gas(upwd.p) );
	}

	return H;
}

double GasOil_RZ_NIT::solve_eq2_dp_beta(int cur, int beta)
{
	Cell& cell = cells[cur];

	double upwind = upwindIsCur(cur, beta);
	Var2phaseNIT& upwd = cells[ getUpwindIdx(cur, beta) ].u_next;
	double Boil_upwd = getB_oil(upwd.p, upwd.p_bub, upwd.SATUR);
	double Bgas_upwd = getB_gas(upwd.p);
	double rs_upwd = getRs(upwd.p, upwd.p_bub, upwd.SATUR);

	return -ht / cell.V * getTrans(cur, beta) * 
			( getKr_oil(upwd.s) * rs_upwd / props_oil.visc / Boil_upwd + getKr_gas(upwd.s) / props_gas.visc / Bgas_upwd - 
			(1.0 - upwind) * (cells[cur].u_next.p - cells[beta].u_next.p) * 
			( getKr_oil(upwd.s) / props_oil.visc / Boil_upwd * (getRs_dp(upwd.p, upwd.p_bub, upwd.SATUR) - rs_upwd * getB_oil_dp(upwd.p, upwd.p_bub, upwd.SATUR) / Boil_upwd) - 
			getKr_gas(upwd.s) / props_gas.visc / Bgas_upwd / Bgas_upwd * getB_gas_dp(upwd.p) ));
}

double GasOil_RZ_NIT::solve_eq2_ds_beta(int cur, int beta)
{
	Cell& cell = cells[cur];

	double upwind = upwindIsCur(cur, beta);
	Var2phaseNIT& upwd = cells[ getUpwindIdx(cur, beta) ].u_next;

	return ht / cell.V * getTrans(cur, beta) * (1.0 - upwind) * (cell.u_next.p - cells[beta].u_next.p) *
		( getRs(upwd.p, upwd.p_bub, upwd.SATUR) * getKr_oil_ds(upwd.s) / props_oil.visc / getB_oil(upwd.p, upwd.p_bub, upwd.SATUR) +
		getKr_gas_ds(upwd.s) / props_gas.visc / getB_gas(upwd.p) );
}

double GasOil_RZ_NIT::solve_eqLeft(int cur)
{
	const int neighbor = cur + cellsNum_z + 2;
	Var2phaseNIT& next = cells[cur].u_next;
	Var2phaseNIT& upwd = cells[ getUpwindIdx(cur, neighbor) ].u_next;

	return getTrans(cur, neighbor) * getKr_oil(upwd.s) / props_oil.visc / getBoreB_oil(next.p, next.p_bub, next.SATUR) * (cells[neighbor].u_next.p - next.p) - Qcell[cur];
}

double GasOil_RZ_NIT::solve_eqLeft_dp(int cur)
{
	const int neighbor = cur + cellsNum_z + 2;
	Var2phaseNIT& next = cells[cur].u_next;
	Var2phaseNIT& upwd = cells[ getUpwindIdx(cur, neighbor) ].u_next;

	return -getTrans(cur, neighbor) * getKr_oil(upwd.s) / getBoreB_oil(next.p, next.p_bub, next.SATUR) / props_oil.visc;
}

double GasOil_RZ_NIT::solve_eqLeft_ds(int cur)
{
	const int neighbor = cur + cellsNum_z + 2;
	Var2phaseNIT& next = cells[cur].u_next;
	Var2phaseNIT& upwd = cells[ getUpwindIdx(cur, neighbor) ].u_next;

	return getTrans(cur, neighbor) * upwindIsCur(cur, neighbor) * getKr_oil_ds(upwd.s) / getBoreB_oil(next.p, next.p_bub, next.SATUR) / props_oil.visc * (cells[neighbor].u_next.p - next.p);
}

double GasOil_RZ_NIT::solve_eqLeft_dp_beta(int cur)
{
	const int neighbor = cur + cellsNum_z + 2;
	Var2phaseNIT& next = cells[cur].u_next;
	Var2phaseNIT& upwd = cells[ getUpwindIdx(cur, neighbor) ].u_next;

	return getTrans(cur, neighbor) * getKr_oil(upwd.s) / getBoreB_oil(next.p, next.p_bub, next.SATUR) / props_oil.visc;
}

double GasOil_RZ_NIT::solve_eqLeft_ds_beta(int cur)
{
	const int neighbor = cur + cellsNum_z + 2;
	Var2phaseNIT& next = cells[cur].u_next;
	Var2phaseNIT& upwd = cells[ getUpwindIdx(cur, neighbor) ].u_next;

	return getTrans(cur, neighbor) * (1.0-upwindIsCur(cur, neighbor)) * getKr_oil_ds(upwd.s) / getBoreB_oil(next.p, next.p_bub, next.SATUR) / props_oil.visc * (cells[neighbor].u_next.p - next.p);
}

double GasOil_RZ_NIT::solve_eq3(int cur)
{
	int neighbor [4];
	getNeighborIdx(cur, neighbor);

	Cell& cell = cells[cur];
	Var2phaseNIT& next = cell.u_next;
	Var2phaseNIT& prev = cell.u_prev;
	double H = 0.0;
	
	H = ( getPoro(next.p) * next.s * getRho_oil(next.p, next.p_bub, next.SATUR) - 
		getPoro(prev.p) * prev.s * getRho_oil(prev.p, prev.p_bub, prev.SATUR) ) / ht;

	for(int i = 0; i < 4; i++)
	{
		Var2phaseNIT& upwd = cells[ getUpwindIdx(cur, neighbor[i]) ].u_next;

		H += 1.0 / cell.V * getTrans(cur, neighbor[i]) * (next.p - cells[ neighbor[i] ].u_next.p) *
			getKr_oil(upwd.s) / props_oil.visc * getRho_oil(upwd.p, upwd.p_bub, upwd.SATUR);
	}

	return H;
}

double GasOil_RZ_NIT::solveH()
{
	double H = 0.0;
	double p1, p0;

	map<int,double>::iterator it = Qcell.begin();
	for(int i = 0; i < Qcell.size()-1; i++)	
	{
		p0 = cells[ it->first ].u_next.p;
		p1 = cells[ (++it)->first ].u_next.p;

		H += (p1 - p0) * (p1 - p0) / 2.0;
	}

	return H;
}