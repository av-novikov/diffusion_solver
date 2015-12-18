#include "model/Oil1D/Oil1D.h"

using namespace std;
using namespace oil1D;

Oil1D::Oil1D()
{
}

Oil1D::~Oil1D()
{
}

void Oil1D::setProps(Properties& props)
{
	varInit.p = props.p_init;

	r_w = props.r_w;
	r_e = props.r_e;
	cellsNum_r = props.cellsNum_r;
	cellsNum = props.cellsNum_r + 2;

	perfIntervals = props.perfIntervals;
	props_sk.m = props.m;
	props_sk.perm = MilliDarcyToM2( props.perm );
	props_sk.dens_stc = props.dens_sk_stc;
	props_sk.beta = props.beta_sk;
	props_sk.height = props.height;

	periodsNum = props.timePeriods.size();	
	for(int i = 0; i < periodsNum; i++)
	{
		period.push_back( props.timePeriods[i] );
		rate.push_back( props.rates[i] / 86400.0 );
		props_sk.skin.push_back( props.skins[i] );
		props_sk.r_eff.push_back( props.radius[i] );

		if(props.radius[i] > props.r_w)
			props_sk.perm_eff.push_back( MilliDarcyToM2(props.perm * log(props.radius[i] / props.r_w) / (log(props.radius[i] / r_w) + props.skins[i]) ) );
		else
			props_sk.perm_eff.push_back( MilliDarcyToM2(props.perm) );
	}

	ht = props.ht;
	ht_min = props.ht_min;
	ht_max = props.ht_max;

	props_oil.visc = cPToPaSec( props.visc_oil );
	props_oil.dens_stc = props.dens_oil_stc;
	props_oil.beta = props.beta_oil;
	props_oil.b_bore = props.b_oil_bore;

	alpha = props.alpha;

	makeDimLess();
}

void Oil1D::makeDimLess()
{
	R_dim = r_w;
	t_dim = 3600.0;
	P_dim = BAR_TO_PA;

	// Temporal properties
	ht /= t_dim;
	ht_min /= t_dim;
	ht_max /= t_dim;

	// Initial condition
	varInit.p /= P_dim;

	// Grid properties
	r_w /= R_dim;
	r_e /= R_dim;

	props_sk.perm /= (R_dim * R_dim);
	props_sk.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_sk.beta /= ( 1.0 / P_dim);
	props_sk.height /= R_dim;

	Q_dim = R_dim * R_dim * R_dim / t_dim;
	for(int i = 0; i < periodsNum; i++)
	{
		period[i] /= t_dim;
		rate[i] /= Q_dim;
		props_sk.r_eff[i] /= R_dim;
		props_sk.perm_eff[i] /= (R_dim * R_dim);
	}

	props_oil.visc /= (P_dim * t_dim);
	props_oil.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_oil.beta /= ( 1.0 / P_dim);

	alpha /= t_dim;
}

void Oil1D::buildGridLog()
{
	cells.reserve( cellsNum );
	
	Volume = 0.0;
	int counter = 0;

	double r_prev = r_w;
	double logMax = log(r_e / r_w);
	double logStep = logMax / (double)cellsNum_r;
	
	double hr = r_prev * (exp(logStep) - 1.0);
	double cm = r_w;

	cells.push_back( Cell(counter++, cm, 0.0, props_sk.height) );
	cm += hr / 2.0;
	for(int i = 0; i < cellsNum_r; i++)
	{
		cm = r_prev * (1.0 + exp(logStep)) / 2.0;
		hr = r_prev * (exp(logStep) - 1.0);
		cells.push_back( Cell(counter++, cm, hr, props_sk.height) );

		Volume += cells[cells.size()-1].V;
		r_prev = r_prev * exp(logStep);
	}
	cells.push_back( Cell(counter++, r_e, 0.0, props_sk.height) );
}

void Oil1D::setPerforated()
{
	Qcell[0] = 0.0;
	height_perf = props_sk.height;
}

void Oil1D::setPeriod(int period)
{
	Q_sum = rate[period];
	Qcell[0] = Q_sum;

	r_eff = props_sk.r_eff[period];
	Perm_eff = props_sk.perm_eff[period];
	skin = props_sk.skin[period];
}

void Oil1D::setSnapshotter(string type)
{
	if(type == "VTK")
		snapshotter = new VTKSnapshotter<Oil1D>();
	else if(type == "GRDECL")
		snapshotter = new GRDECLSnapshotter<Oil1D>();
	else
		snapshotter = new GRDECLSnapshotter<Oil1D>();

	snapshotter->setModel(this);
}

void Oil1D::snapshot(int i)
{
	snapshotter->dump(i);
}

void Oil1D::snapshot_all(int i)
{
	snapshotter->dump_all(i);
}

double Oil1D::solve_eq(int i)
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

double Oil1D::solve_eq_dp(int i)
{
	Cell& cell = cells[i];
	double H = getPoro(cell.u_next.p) * props_oil.dens_stc * props_oil.beta + getRho(cell.u_next.p) * props_sk.m * props_sk.beta;
	
	for(int k = i-1; k < i+2; k += 2)
	{
		Cell& cell1 = cells[k];
		H += ht / cell.V / props_oil.visc * getTrans(cell, cell1) * (getRho(cell, cell1) + (cell.u_next.p - cell1.u_next.p) * cell1.hr / (cell1.hr + cell.hr) * props_oil.dens_stc * props_oil.beta);
	}
	
	return H;
}

double Oil1D::solve_eq_dp_beta(int i, int beta)
{
	Cell& cell = cells[i];
	Cell& cell1 = cells[beta];
	
	return ht / cell.V / props_oil.visc * getTrans(cell, cell1) * ( (cell.u_next.p - cell1.u_next.p) * cell.hr / (cell.hr + cell1.hr) * props_oil.dens_stc * props_oil.beta - getRho(cell, cell1));
}

double Oil1D::solve_left()
{
	return getTrans(cells[0], cells[1]) / props_oil.visc / props_oil.b_bore * (cells[1].u_next.p - cells[0].u_next.p) - Qcell[0];
}

double Oil1D::solve_left_dp()
{
	return -getTrans(cells[0], cells[1]) / props_oil.visc / props_oil.b_bore;
}

double Oil1D::solve_left_dp_beta()
{
	return getTrans(cells[0], cells[1]) / props_oil.visc / props_oil.b_bore;
}

double Oil1D::solve_right()
{
	return cells[cellsNum-1].u_next.p - varInit.p;
}

double Oil1D::solve_right_dp()
{
	return 1.0;
}

double Oil1D::solve_right_dp_beta()
{
	return 0.0;
}
