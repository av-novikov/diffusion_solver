#include "model/Oil1D_NIT/Oil1D_NIT.h"

using namespace std;
using namespace oil1D_NIT;

Oil1D_NIT::Oil1D_NIT()
{
	skeletonsNum = 1;
	cellsNum_z = 1;
}

Oil1D_NIT::~Oil1D_NIT()
{
}

void Oil1D_NIT::setProps(Properties& props)
{
	leftBoundIsRate = props.leftBoundIsRate;
	rightBoundIsPres = props.rightBoundIsPres;

	// Setting grid properties
	r_w = props.r_w;
	r_e = props.r_e;
	cellsNum_r = props.cellsNum_r;
	cellsNum = (cellsNum_r + 2) * cellsNum_z;

	// Setting skeleton properties
	perfIntervals = props.perfIntervals;
	
	props_sk = props.props_sk;
	props_sk[0].perm_r = MilliDarcyToM2( props_sk[0].perm_r );
	props_sk[0].perm_z = MilliDarcyToM2( props_sk[0].perm_z );

	periodsNum = props.timePeriods.size();
	for(int i = 0; i < periodsNum; i++)
	{
		period.push_back( props.timePeriods[i] );
		if(leftBoundIsRate)
			rate.push_back( props.rates[i] / 86400.0 );
		else
			pwf.push_back( props.pwf[i] );

		if(props_sk[0].radiuses_eff[i] > props.r_w)
			props_sk[0].perms_eff.push_back( MilliDarcyToM2(props.props_sk[0].perm_r * log(props.props_sk[0].radiuses_eff[i] / props.r_w) / (log(props.props_sk[0].radiuses_eff[i] / props.r_w) + props.props_sk[0].skins[i]) ) );
		else
			props_sk[0].perms_eff.push_back( MilliDarcyToM2(props.props_sk[0].perm_r) );
	}

	// Temporal properties
	ht = props.ht;
	ht_min = props.ht_min;
	ht_max = props.ht_max;

	// Oil properties
	props_oil = props.props_oil;
	props_oil.visc = cPToPaSec( props_oil.visc );

	alpha = props.alpha;
	depth_point = props.depth_point;

	makeDimLess();
}

void Oil1D_NIT::makeDimLess()
{
	R_dim = r_w;
	t_dim = 3600.0;
	P_dim = BAR_TO_PA;
	T_dim = props_sk[0].t_init;

	// Temporal properties
	ht /= t_dim;
	ht_min /= t_dim;
	ht_max /= t_dim;

	// Grid properties
	r_w = r_w / R_dim;
	r_e = r_e / R_dim;

	// Skeleton properties
	for(int i = 0; i < skeletonsNum; i++)
	{
		props_sk[i].perm_r /= (R_dim * R_dim);
		props_sk[i].perm_z /= (R_dim * R_dim);
		
		props_sk[i].beta /= (1.0 / P_dim);
		props_sk[i].dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
		props_sk[i].h1 = (props_sk[i].h1 - depth_point) / R_dim;
		props_sk[i].h2 = (props_sk[i].h2 - depth_point) / R_dim;
		props_sk[i].height /= R_dim;
		props_sk[i].p_init /= P_dim;
		props_sk[i].p_out /= P_dim;
		props_sk[i].t_init /= T_dim;

		props_sk[i].c = props_sk[i].c / R_dim / R_dim * T_dim * t_dim * t_dim;
		props_sk[i].lambda = props_sk[i].lambda * T_dim * t_dim / P_dim / R_dim / R_dim;

		for(int j = 0; j < periodsNum; j++)
		{
			props_sk[i].perms_eff[j] /= (R_dim * R_dim);
			props_sk[i].radiuses_eff[j] /= R_dim;
		}
	}

	Q_dim = R_dim * R_dim * R_dim / t_dim;
	for(int i = 0; i < periodsNum; i++)
	{
		period[i] /= t_dim;
		if(leftBoundIsRate)
			rate[i] /= Q_dim;
		else
			pwf[i] /= P_dim;
	}

	// Oil properties
	props_oil.visc /= (P_dim * t_dim);
	props_oil.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_oil.beta /= (1.0 / P_dim);
	props_oil.c = props_oil.c / R_dim / R_dim * T_dim * t_dim * t_dim;
	props_oil.lambda = props_oil.lambda * T_dim * t_dim / P_dim / R_dim / R_dim;
	props_oil.jt = props_oil.jt * P_dim / T_dim;
	props_oil.ad = props_oil.ad * P_dim / T_dim;

	// Rest properties
	alpha /= t_dim;
	//depth_point = 0.0;
}

void Oil1D_NIT::buildGridLog()
{
	cells.reserve( cellsNum );
	
	Volume = 0.0;
	int counter = 0;

	double r_prev = r_w;
	double logMax = log(r_e / r_w);
	double logStep = logMax / (double)cellsNum_r;
	
	double hr = r_prev * (exp(logStep) - 1.0);
	double cm = r_w;

	cells.push_back( Cell(counter++, cm, 0.0, props_sk[0].height) );
	cm += hr / 2.0;
	for(int i = 0; i < cellsNum_r; i++)
	{
		cm = r_prev * (1.0 + exp(logStep)) / 2.0;
		hr = r_prev * (exp(logStep) - 1.0);
		cells.push_back( Cell(counter++, cm, hr, props_sk[0].height) );

		Volume += cells[cells.size()-1].V;
		r_prev = r_prev * exp(logStep);
	}
	cells.push_back( Cell(counter++, r_e, 0.0, props_sk[0].height) );
}

void Oil1D_NIT::setInitialState()
{
	vector<Cell>::iterator it;
	for(it = cells.begin(); it != cells.end(); ++it)
	{
		it->u_prev.p = it->u_iter.p = it->u_next.p = props_sk[0].p_init;
		it->u_prev.t = it->u_iter.t = it->u_next.t = props_sk[0].t_init;
	}
}

void Oil1D_NIT::setPerforated()
{
	Qcell[0] = 0.0;
	height_perf = props_sk[0].height;
}

void Oil1D_NIT::setPeriod(int period)
{
	if(leftBoundIsRate)
		Q_sum = rate[period];
	else
	{
		Pwf = pwf[period];
		Q_sum = 0.0;
	}
	
	Qcell[0] = Q_sum;

	props_sk[0].radius_eff = props_sk[0].radiuses_eff[period];
	props_sk[0].perm_eff = props_sk[0].perms_eff[period];
	props_sk[0].skin = props_sk[0].skins[period];
}

double Oil1D_NIT::getRate()
{
	if(leftBoundIsRate)
		return Qcell[0];
	else {
		Cell& cell = cells[0];
		Cell& cell1 = cells[1];
		return getTrans(cell, cell1) / props_oil.b_bore / props_oil.visc * (cell1.u_next.p - cell.u_next.p);
	}
}

double Oil1D_NIT::solve_eq(int i)
{
	Cell& cell = cells[i];
	double H = getPoro(cell.u_next.p) * getRho(cell.u_next.p) - getPoro(cell.u_prev.p) * getRho(cell.u_prev.p);

	for(int k = i-1; k < i+2; k += 2)
	{
		Cell& cell1 = cells[k];
		H += ht / cell.V / props_oil.visc * getTrans(cell, cell1) * (cell.u_next.p - cell1.u_next.p) * getRho(cell, cell1);
	}
	
	return H;
}

double Oil1D_NIT::solve_eq_dp(int i)
{
	Cell& cell = cells[i];
	double H = getPoro(cell.u_next.p) * props_oil.dens_stc * props_oil.beta + getRho(cell.u_next.p) * props_sk[0].m * props_sk[0].beta;
	
	for(int k = i-1; k < i+2; k += 2)
	{
		Cell& cell1 = cells[k];
		H += ht / cell.V / props_oil.visc * getTrans(cell, cell1) * (getRho(cell, cell1) + (cell.u_next.p - cell1.u_next.p) * cell1.hr / (cell1.hr + cell.hr) * props_oil.dens_stc * props_oil.beta);
	}
	
	return H;
}

double Oil1D_NIT::solve_eq_dp_beta(int i, int beta)
{
	Cell& cell = cells[i];
	Cell& cell1 = cells[beta];
	
	return ht / cell.V / props_oil.visc * getTrans(cell, cell1) * ( (cell.u_next.p - cell1.u_next.p) * cell.hr / (cell.hr + cell1.hr) * props_oil.dens_stc * props_oil.beta - getRho(cell, cell1));
}

double Oil1D_NIT::solve_eqLeft()
{
	Cell& cell = cells[0];
	Cell& cell1 = cells[1];
	
	if( leftBoundIsRate )
		return getTrans(cell, cell1) / props_oil.b_bore / props_oil.visc * (cell1.u_next.p - cell.u_next.p) - Qcell[0];
	else
		return cell.u_next.p - Pwf;
}

double Oil1D_NIT::solve_eqLeft_dp()
{
	if( leftBoundIsRate )
		return -getTrans(cells[0], cells[1]) / props_oil.b_bore / props_oil.visc;
	else
		return 1.0;
}

double Oil1D_NIT::solve_eqLeft_dp_beta()
{
	if( leftBoundIsRate )
		return getTrans(cells[0], cells[1]) / props_oil.b_bore / props_oil.visc;
	else
		return 0.0;
}

double Oil1D_NIT::solve_eqRight()
{
	Cell& cell = cells[cellsNum-1];

	if( rightBoundIsPres )
		return cell.u_next.p - props_sk[0].p_out;
	else
		return cell.u_next.p - cells[cellsNum-2].u_next.p;
}

double Oil1D_NIT::solve_eqRight_dp()
{
	if( rightBoundIsPres )
		return 1.0;
	else
		return 1.0;
}

double Oil1D_NIT::solve_eqRight_dp_beta()
{
	if( rightBoundIsPres )
		return 0.0;
	else
		return -1.0;
}
