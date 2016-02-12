#include "model/Gas1D/Gas1D.h"
#include "util/utils.h"

#include <cassert>

using std::map;
using std::string;
using std::vector;
using namespace gas1D;

Gas1D::Gas1D()
{
	skeletonsNum = 1;
	cellsNum_z = 1;
}

Gas1D::~Gas1D()
{
}

void Gas1D::setProps(Properties& props)
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

	// Gas properties
	props_gas = props.props_gas;
	//props_gas.visc = cPToPaSec( props_gas.visc );

	alpha = props.alpha;
	depth_point = props.depth_point;

	makeDimLess();

	props_gas.visc_table = setDataset(props.visc_gas, P_dim / BAR_TO_PA, PaSec2cP(1.0) * P_dim * t_dim);
	props_gas.z = setDataset(props.z_factor, P_dim / BAR_TO_PA, 1.0);
}

void Gas1D::makeDimLess()
{
	// Main units
	R_dim = r_w;
	t_dim = 3600.0;
	P_dim = BAR_TO_PA;

	// Temporal properties
	ht /= t_dim;
	ht_min /= t_dim;
	ht_max /= t_dim;

	// Grid properties
	r_w /= R_dim;
	r_e /= R_dim;

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

	// Gas properties
	//props_gas.visc /= (P_dim * t_dim);
	props_gas.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);

	// Rest properties
	alpha /= t_dim;
	//depth_point = 0.0;
}

void Gas1D::buildGridLog()
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

void Gas1D::setInitialState()
{
	vector<Cell>::iterator it;
	for(it = cells.begin(); it != cells.end(); ++it)
		it->u_prev.p = it->u_iter.p = it->u_next.p = props_sk[0].p_init;
}

void Gas1D::setPerforated()
{
	Qcell[0] = 0.0;
	height_perf = props_sk[0].height;
}

void Gas1D::setPeriod(int period)
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

double Gas1D::getRate()
{
	if(leftBoundIsRate)
		return Qcell[0];
	else {
		Cell& cell = cells[0];
		Cell& cell1 = cells[1];
		return getTrans(cell, cell1) / P_ATM * getCoeff(cell, cell1) * (cell1.u_next.p - cell.u_next.p);
	}
}

double Gas1D::solve_eq(int i)
{
	Cell& cell = cells[i];
	double p_next = cell.u_next.p;
	double p_prev = cell.u_next.p;

	double H = getPoro(p_next) * p_next / getZ(p_next) - getPoro(p_prev) * p_prev / getZ(p_prev);

	for(int k = i-1; k < i+2; k += 2)
	{
		Cell& cell1 = cells[k];
		H += ht / cell.V * getTrans(cell, cell1) * (cell.u_next.p - cell1.u_next.p) * getCoeff(cell, cell1);
	}
	
	return H;
}

double Gas1D::solve_eq_dp(int i)
{
	Cell& cell = cells[i];
	double p = cell.u_next.p;
	double H = (p * getPoro_dp() + getPoro(p) * (1.0 - p / getZ(p) * getZ_dp(p))) / getZ(p);
	
	for(int k = i-1; k < i+2; k += 2)
	{
		Cell& cell1 = cells[k];
		H += ht / cell.V * getTrans(cell, cell1) * 
			( (cell.u_next.p - cell1.u_next.p) * getCoeff_dp(cell, cell1) + getCoeff(cell, cell1) );
	}
	
	return H;
}

double Gas1D::solve_eq_dp_beta(int i, int beta)
{
	Cell& cell = cells[i];
	Cell& cell1 = cells[beta];
	
	return ht / cell.V * getTrans(cell, cell1) * 
			( (cell.u_next.p - cell1.u_next.p) * getCoeff_dp_beta(cell, cell1) - getCoeff(cell, cell1) );
}

double Gas1D::solve_eqLeft()
{
	Cell& cell = cells[0];
	Cell& cell1 = cells[1];

	if( leftBoundIsRate )
		return getTrans(cell, cell1) / P_ATM * getCoeff(cell, cell1) * (cell1.u_next.p - cell.u_next.p) - Qcell[0];
	else
		return cell.u_next.p - Pwf;
}

double Gas1D::solve_eqLeft_dp()
{
	if( leftBoundIsRate )
	{
		Cell& cell = cells[0];
		Cell& cell1 = cells[1];
		return getTrans(cell, cell1) / P_ATM * 
				( (cell1.u_next.p - cell.u_next.p) * getCoeff_dp(cell, cell1) - getCoeff(cell, cell1) );
	} else
		return 1.0;
}

double Gas1D::solve_eqLeft_dp_beta()
{
	if( leftBoundIsRate )
	{
		Cell& cell = cells[0];
		Cell& cell1 = cells[1];
		return getTrans(cell, cell1) / P_ATM * 
				( (cell1.u_next.p - cell.u_next.p) * getCoeff_dp_beta(cell, cell1) + getCoeff(cell, cell1) );
	} else
		return 0.0;
}

double Gas1D::solve_eqRight()
{
	Cell& cell = cells[cellsNum-1];

	if( rightBoundIsPres )
		return cell.u_next.p - props_sk[0].p_out;
	else
		return cell.u_next.p - cells[cellsNum-2].u_next.p;
}

double Gas1D::solve_eqRight_dp()
{
	if( rightBoundIsPres )
		return 1.0;
	else
		return 1.0;
}

double Gas1D::solve_eqRight_dp_beta()
{
	if( rightBoundIsPres )
		return 0.0;
	else
		return -1.0;
}