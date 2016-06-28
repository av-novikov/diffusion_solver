#include "model/3D/Perforation/GasOil_Perf_NIT.h"
#include "util/utils.h"

#include <cassert>

using namespace std;
using namespace std::placeholders;
using namespace gasOil_perf_nit;


GasOil_Perf_NIT::GasOil_Perf_NIT()
{
	isWriteSnaps = false;

	// Middle
	middleFoo.mat.resize(2);

	middleFoo.rhs.push_back(bind(&GasOil_Perf_NIT::solve_eq1, this, _1));
	middleFoo.mat[0].push_back(bind(&GasOil_Perf_NIT::solve_eq1_dp, this, _1, _2));
	middleFoo.mat[0].push_back(bind(&GasOil_Perf_NIT::solve_eq1_ds, this, _1, _2));
	middleFoo.mat[0].push_back(bind(&GasOil_Perf_NIT::solve_eq1_dp_beta, this, _1, _2));
	middleFoo.mat[0].push_back(bind(&GasOil_Perf_NIT::solve_eq1_ds_beta, this, _1, _2));
	middleFoo.mat[0].push_back(bind(&GasOil_Perf_NIT::solve_eq1_dp_beta, this, _1, _2));
	middleFoo.mat[0].push_back(bind(&GasOil_Perf_NIT::solve_eq1_ds_beta, this, _1, _2));
	middleFoo.mat[0].push_back(bind(&GasOil_Perf_NIT::solve_eq1_dp_beta, this, _1, _2));
	middleFoo.mat[0].push_back(bind(&GasOil_Perf_NIT::solve_eq1_ds_beta, this, _1, _2));
	middleFoo.mat[0].push_back(bind(&GasOil_Perf_NIT::solve_eq1_dp_beta, this, _1, _2));
	middleFoo.mat[0].push_back(bind(&GasOil_Perf_NIT::solve_eq1_ds_beta, this, _1, _2));
	middleFoo.mat[0].push_back(bind(&GasOil_Perf_NIT::solve_eq1_dp_beta, this, _1, _2));
	middleFoo.mat[0].push_back(bind(&GasOil_Perf_NIT::solve_eq1_ds_beta, this, _1, _2));
	middleFoo.mat[0].push_back(bind(&GasOil_Perf_NIT::solve_eq1_dp_beta, this, _1, _2));
	middleFoo.mat[0].push_back(bind(&GasOil_Perf_NIT::solve_eq1_ds_beta, this, _1, _2));

	middleFoo.rhs.push_back(bind(&GasOil_Perf_NIT::solve_eq2, this, _1));
	middleFoo.mat[1].push_back(bind(&GasOil_Perf_NIT::solve_eq2_dp, this, _1, _2));
	middleFoo.mat[1].push_back(bind(&GasOil_Perf_NIT::solve_eq2_ds, this, _1, _2));
	middleFoo.mat[1].push_back(bind(&GasOil_Perf_NIT::solve_eq2_dp_beta, this, _1, _2));
	middleFoo.mat[1].push_back(bind(&GasOil_Perf_NIT::solve_eq2_ds_beta, this, _1, _2));
	middleFoo.mat[1].push_back(bind(&GasOil_Perf_NIT::solve_eq2_dp_beta, this, _1, _2));
	middleFoo.mat[1].push_back(bind(&GasOil_Perf_NIT::solve_eq2_ds_beta, this, _1, _2));
	middleFoo.mat[1].push_back(bind(&GasOil_Perf_NIT::solve_eq2_dp_beta, this, _1, _2));
	middleFoo.mat[1].push_back(bind(&GasOil_Perf_NIT::solve_eq2_ds_beta, this, _1, _2));
	middleFoo.mat[1].push_back(bind(&GasOil_Perf_NIT::solve_eq2_dp_beta, this, _1, _2));
	middleFoo.mat[1].push_back(bind(&GasOil_Perf_NIT::solve_eq2_ds_beta, this, _1, _2));
	middleFoo.mat[1].push_back(bind(&GasOil_Perf_NIT::solve_eq2_dp_beta, this, _1, _2));
	middleFoo.mat[1].push_back(bind(&GasOil_Perf_NIT::solve_eq2_ds_beta, this, _1, _2));
	middleFoo.mat[1].push_back(bind(&GasOil_Perf_NIT::solve_eq2_dp_beta, this, _1, _2));
	middleFoo.mat[1].push_back(bind(&GasOil_Perf_NIT::solve_eq2_ds_beta, this, _1, _2));

	// Left
	leftFoo.mat.resize(2);

	leftFoo.rhs.push_back(bind(&GasOil_Perf_NIT::solve_eq1Left, this, _1));
	leftFoo.mat[0].push_back(bind(&GasOil_Perf_NIT::solve_eq1Left_dp, this, _1, _2));
	leftFoo.mat[0].push_back(bind(&GasOil_Perf_NIT::solve_eq1Left_ds, this, _1, _2));
	leftFoo.mat[0].push_back(bind(&GasOil_Perf_NIT::solve_eq1Left_dp_beta, this, _1, _2));
	leftFoo.mat[0].push_back(bind(&GasOil_Perf_NIT::solve_eq1Left_ds_beta, this, _1, _2));

	leftFoo.rhs.push_back(bind(&GasOil_Perf_NIT::solve_eq2Left, this, _1));
	leftFoo.mat[1].push_back(bind(&GasOil_Perf_NIT::solve_eq2Left_dp, this, _1, _2));
	leftFoo.mat[1].push_back(bind(&GasOil_Perf_NIT::solve_eq2Left_ds, this, _1, _2));
	leftFoo.mat[1].push_back(bind(&GasOil_Perf_NIT::solve_eq2Left_dp_beta, this, _1, _2));
	leftFoo.mat[1].push_back(bind(&GasOil_Perf_NIT::solve_eq2Left_ds_beta, this, _1, _2));
	leftFoo.mat[1].push_back(bind(&GasOil_Perf_NIT::solve_eq2Left_dp_beta, this, _1, _2));
	leftFoo.mat[1].push_back(bind(&GasOil_Perf_NIT::solve_eq2Left_ds_beta, this, _1, _2));

	// Right
	rightFoo.mat.resize(2);

	rightFoo.rhs.push_back(bind(&GasOil_Perf_NIT::solve_eq1Right, this, _1));
	rightFoo.mat[0].push_back(bind(&GasOil_Perf_NIT::solve_eq1Right_dp, this, _1, _2));
	rightFoo.mat[0].push_back(bind(&GasOil_Perf_NIT::solve_eq1Right_ds, this, _1, _2));
	rightFoo.mat[0].push_back(bind(&GasOil_Perf_NIT::solve_eq1Right_dp_beta, this, _1, _2));
	rightFoo.mat[0].push_back(bind(&GasOil_Perf_NIT::solve_eq1Right_ds_beta, this, _1, _2));

	rightFoo.rhs.push_back(bind(&GasOil_Perf_NIT::solve_eq2Right, this, _1));
	rightFoo.mat[1].push_back(bind(&GasOil_Perf_NIT::solve_eq2Right_dp, this, _1, _2));
	rightFoo.mat[1].push_back(bind(&GasOil_Perf_NIT::solve_eq2Right_ds, this, _1, _2));
	rightFoo.mat[1].push_back(bind(&GasOil_Perf_NIT::solve_eq2Right_dp_beta, this, _1, _2));
	rightFoo.mat[1].push_back(bind(&GasOil_Perf_NIT::solve_eq2Right_ds_beta, this, _1, _2));

	// Top
	topFoo.mat.resize(2);

	topFoo.rhs.push_back(bind(&GasOil_Perf_NIT::solve_eq1Top, this, _1));
	topFoo.mat[0].push_back(bind(&GasOil_Perf_NIT::solve_eq1Top_dp, this, _1, _2));
	topFoo.mat[0].push_back(bind(&GasOil_Perf_NIT::solve_eq1Top_ds, this, _1, _2));
	topFoo.mat[0].push_back(bind(&GasOil_Perf_NIT::solve_eq1Top_dp_beta, this, _1, _2));
	topFoo.mat[0].push_back(bind(&GasOil_Perf_NIT::solve_eq1Top_ds_beta, this, _1, _2));

	topFoo.rhs.push_back(bind(&GasOil_Perf_NIT::solve_eq2Top, this, _1));
	topFoo.mat[1].push_back(bind(&GasOil_Perf_NIT::solve_eq2Top_dp, this, _1, _2));
	topFoo.mat[1].push_back(bind(&GasOil_Perf_NIT::solve_eq2Top_ds, this, _1, _2));
	topFoo.mat[1].push_back(bind(&GasOil_Perf_NIT::solve_eq2Top_dp_beta, this, _1, _2));
	topFoo.mat[1].push_back(bind(&GasOil_Perf_NIT::solve_eq2Top_ds_beta, this, _1, _2));

	// Bot
	botFoo.mat.resize(2);

	botFoo.rhs.push_back(bind(&GasOil_Perf_NIT::solve_eq1Bot, this, _1));
	botFoo.mat[0].push_back(bind(&GasOil_Perf_NIT::solve_eq1Bot_dp, this, _1, _2));
	botFoo.mat[0].push_back(bind(&GasOil_Perf_NIT::solve_eq1Bot_ds, this, _1, _2));
	botFoo.mat[0].push_back(bind(&GasOil_Perf_NIT::solve_eq1Bot_dp_beta, this, _1, _2));
	botFoo.mat[0].push_back(bind(&GasOil_Perf_NIT::solve_eq1Bot_ds_beta, this, _1, _2));

	botFoo.rhs.push_back(bind(&GasOil_Perf_NIT::solve_eq2Bot, this, _1));
	botFoo.mat[1].push_back(bind(&GasOil_Perf_NIT::solve_eq2Bot_dp, this, _1, _2));
	botFoo.mat[1].push_back(bind(&GasOil_Perf_NIT::solve_eq2Bot_ds, this, _1, _2));
	botFoo.mat[1].push_back(bind(&GasOil_Perf_NIT::solve_eq2Bot_dp_beta, this, _1, _2));
	botFoo.mat[1].push_back(bind(&GasOil_Perf_NIT::solve_eq2Bot_ds_beta, this, _1, _2));
}

GasOil_Perf_NIT::~GasOil_Perf_NIT()
{
}

void GasOil_Perf_NIT::setProps(Properties& props)
{
	leftBoundIsRate = props.leftBoundIsRate;
	rightBoundIsPres = props.rightBoundIsPres;

	// Setting grid properties
	r_w = props.r_w;
	r_e = props.r_e;
	cellsNum_r = props.cellsNum_r;
	cellsNum_phi = props.cellsNum_phi;
	cellsNum_z = props.cellsNum_z;
	cellsNum = (cellsNum_r + 2) * (cellsNum_z + 2) * cellsNum_phi;

	// Setting skeleton properties
	perfTunnels = props.perfTunnels;
	
	skeletonsNum = props.props_sk.size();
	checkSkeletons(props.props_sk);
	props_sk = props.props_sk;
	for(int j = 0; j < skeletonsNum; j++)
	{
		props_sk[j].perm_r = MilliDarcyToM2( props_sk[j].perm_r );
		props_sk[j].perm_z = MilliDarcyToM2( props_sk[j].perm_z );
	}

	periodsNum = props.timePeriods.size();
	for(int i = 0; i < periodsNum; i++)
	{
		period.push_back( props.timePeriods[i] );
		if(leftBoundIsRate)
			rate.push_back( props.rates[i] / 86400.0 );
		else
			pwf.push_back( props.pwf[i] );
		for(int j = 0; j < skeletonsNum; j++)
		{
			if(props_sk[j].radiuses_eff[i] > props.r_w)
				props_sk[j].perms_eff.push_back( MilliDarcyToM2(props.props_sk[j].perm_r * log(props.props_sk[j].radiuses_eff[i] / props.r_w) / (log(props.props_sk[j].radiuses_eff[i] / props.r_w) + props.props_sk[j].skins[i]) ) );
			else
				props_sk[j].perms_eff.push_back( MilliDarcyToM2(props.props_sk[j].perm_r) );
		}
	}

	// Temporal properties
	ht = props.ht;
	ht_min = props.ht_min;
	ht_max = props.ht_max;

	// Oil properties
	props_oil = props.props_oil;
	props_oil.visc = cPToPaSec( props_oil.visc );

	// Gas properties
	props_gas = props.props_gas;
	props_gas.visc = cPToPaSec( props_gas.visc );

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

void GasOil_Perf_NIT::checkSkeletons(const vector<Skeleton_Props>& props)
{
	vector<Skeleton_Props>::const_iterator it = props.begin();
	double tmp;
	int indxs = 0;

//	assert( it->h2 - it->h1 == it->height );
	indxs += it->cellsNum_z;
	tmp = it->h2;
	++it;
	
	while(it != props.end())
	{
		//assert( it->h1 == tmp );
		//assert( it->h2 - it->h1 == it->height );
		indxs += it->cellsNum_z;
		tmp = it->h2;
		++it;
	}
	//assert( indxs == cellsNum_z );
}

void GasOil_Perf_NIT::makeDimLess()
{
	// Main units
	R_dim = r_w * 10.0;
	t_dim = 3600.0;
	P_dim = props_sk[0].p_init;
	if (props_sk[0].t_init != 0.0)
		T_dim = fabs(props_sk[0].t_init);
	else
		T_dim = 1.0;

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
		props_sk[i].p_bub /= P_dim;
		props_sk[i].t_init /= T_dim;

		props_sk[i].c = props_sk[i].c / R_dim / R_dim * T_dim * t_dim * t_dim;
		props_sk[i].lambda_r = props_sk[i].lambda_r * T_dim * t_dim / P_dim / R_dim / R_dim;
		props_sk[i].lambda_z = props_sk[i].lambda_z * T_dim * t_dim / P_dim / R_dim / R_dim;

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

	// Gas properties
	props_gas.visc /= (P_dim * t_dim);
	props_gas.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);

	props_gas.c = props_gas.c / R_dim / R_dim * T_dim * t_dim * t_dim;
	props_gas.lambda = props_gas.lambda * T_dim * t_dim / P_dim / R_dim / R_dim;
	props_gas.jt = props_gas.jt * P_dim / T_dim;
	props_gas.ad = props_gas.ad * P_dim / T_dim;

	// Rest properties
	alpha /= t_dim;
	//depth_point = 0.0;
	L = L / R_dim / R_dim * t_dim * t_dim;
}

void GasOil_Perf_NIT::buildGridLog()
{
	cells.reserve( cellsNum );

	Volume = 0.0;
	int counter = 0;
	int skel_idx = 0, cells_z = 0;

	double r_prev = r_w;
	double logMax = log(r_e / r_w);
	double logStep = logMax / (double)cellsNum_r;
	
	const double hphi = 2.0 * M_PI / (double)cellsNum_phi;
	double cm_phi = 0.0;
	double hz = 0.0;//(props_sk[skel_idx].h2 - props_sk[skel_idx].h1) / (double)(props_sk[skel_idx].cellsNum_z);
	double cm_z = props_sk[skel_idx].h1;
	double hr = r_prev * (exp(logStep) - 1.0);
	double cm_r = r_w;

	counter = 0;
	for(int k = 0; k < cellsNum_phi; k++)
	{
		skel_idx = 0;	cells_z = 0;

		r_prev = r_w;
		logMax = log(r_e / r_w);
		logStep = logMax / (double)cellsNum_r;

		hz = 0.0;
		cm_z = props_sk[0].h1;
		hr = r_prev * (exp(logStep) - 1.0);
		cm_r = r_w;

		// Left border
		cells.push_back( Cell(counter++, cm_r, cm_phi, cm_z, 0.0, hphi, 0.0) );
		for(int i = 0; i < cellsNum_z; i++)
		{
			hz = (props_sk[skel_idx].h2 - props_sk[skel_idx].h1) / (double)(props_sk[skel_idx].cellsNum_z);
			cm_z += (cells[cells.size()-1].hz + hz) / 2.0;

			cells.push_back( Cell(counter++, cm_r, cm_phi, cm_z, 0.0, hphi, hz) );
			cells_z++;

			if(cells_z >= props_sk[skel_idx].cellsNum_z)
			{
				cells_z = 0;
				skel_idx++;
			}
		}
		cells.push_back( Cell(counter++, cm_r, cm_phi, cm_z+hz/2.0, 0.0, hphi, 0.0) );

		// Middle cells
		for(int j = 0; j < cellsNum_r; j++)
		{
			skel_idx = 0;	cells_z = 0;
			cm_z = props_sk[0].h1;
			cm_r = r_prev * (exp(logStep) + 1.0) / 2.0;
			hr = r_prev * (exp(logStep) - 1.0);

			cells.push_back( Cell(counter++, cm_r, cm_phi, cm_z, hr, hphi, 0.0) );
			for(int i = 0; i < cellsNum_z; i++)
			{
				hz = (props_sk[skel_idx].h2 - props_sk[skel_idx].h1) / (double)(props_sk[skel_idx].cellsNum_z);
				cm_z += (cells[cells.size()-1].hz + hz) / 2.0;

				cells.push_back( Cell(counter++, cm_r, cm_phi, cm_z, hr, hphi, hz) );
				Volume += cells[cells.size()-1].V;
				cells_z++;

				if(cells_z >= props_sk[skel_idx].cellsNum_z)
				{
					cells_z = 0;
					skel_idx++;
				}
			}
			cells.push_back( Cell(counter++, cm_r, cm_phi, cm_z+hz/2.0, hr, hphi, 0.0) );

			r_prev = r_prev * exp(logStep);
		}

		// Right border
		cm_z = props_sk[0].h1;
		cm_r = r_e;
		
		cells.push_back( Cell(counter++, cm_r, cm_phi, cm_z, 0.0, hphi, 0.0) );
		skel_idx = 0;	cells_z = 0;
		for(int i = 0; i < cellsNum_z; i++)
		{
			hz = (props_sk[skel_idx].h2 - props_sk[skel_idx].h1) / (double)(props_sk[skel_idx].cellsNum_z);
			cm_z += (cells[cells.size()-1].hz + hz) / 2.0;

			cells.push_back( Cell(counter++, cm_r, cm_phi, cm_z, 0.0, hphi, hz) );
			cells_z++;

			if(cells_z >= props_sk[skel_idx].cellsNum_z)
			{
				cells_z = 0;
				skel_idx++;
			}
		}
		cells.push_back( Cell(counter++, cm_r, cm_phi, cm_z+hz/2.0, 0.0, hphi, 0.0) );
		cm_phi += hphi;
	}

	setUnused();
	buildTunnels();

	// Creating iterators
	midIter = new Iterator(&cells[cellsNum_z + 2], { 1, 0, 0 }, { cellsNum_r, cellsNum_phi - 1, cellsNum_z + 1 }, { cellsNum_r + 2, cellsNum_phi, cellsNum_z + 2 });
	midBegin = new Iterator( *midIter );
	midEnd = new Iterator(nullptr, { 1, 0, 0 }, { cellsNum_r, cellsNum_phi - 1, cellsNum_z + 1 }, { cellsNum_r + 2, cellsNum_phi, cellsNum_z + 2 });

	leftIter = new Iterator(&cells[0], { 0, 0, 0 }, { 0, cellsNum_phi - 1, cellsNum_z + 1 }, { cellsNum_r + 2, cellsNum_phi, cellsNum_z + 2 });
	leftBegin = new Iterator(*leftIter);
	leftEnd = new Iterator(nullptr, { 0, 0, 0 }, { 0, cellsNum_phi - 1, cellsNum_z + 1 }, { cellsNum_r + 2, cellsNum_phi, cellsNum_z + 2 });

	rightIter = new Iterator(&cells[(cellsNum_r+1)*(cellsNum_z+2)], { cellsNum_r+1, 0, 0 }, { cellsNum_r+1, cellsNum_phi - 1, cellsNum_z + 1 }, { cellsNum_r + 2, cellsNum_phi, cellsNum_z + 2 });
	rightBegin = new Iterator(*rightIter);
	rightEnd = new Iterator(nullptr, { cellsNum_r+1, 0, 0 }, { cellsNum_r+1, cellsNum_phi - 1, cellsNum_z + 1 }, { cellsNum_r + 2, cellsNum_phi, cellsNum_z + 2 });
}

void GasOil_Perf_NIT::setInitialState()
{
	vector<Cell>::iterator it;
	for(it = cells.begin(); it != cells.end(); ++it)
	{
		const Skeleton_Props& props = props_sk[ getSkeletonIdx(*it) ];
		it->u_prev.p = it->u_iter.p = it->u_next.p = props.p_init;
		it->u_prev.p_bub = it->u_iter.p_bub = it->u_next.p_bub = props.p_bub;
		it->u_prev.s = it->u_iter.s = it->u_next.s = props.s_init;
		it->u_prev.t = it->u_iter.t = it->u_next.t = props.t_init;
		if(props.p_bub > props.p_init)
			it->u_prev.SATUR = it->u_iter.SATUR = it->u_next.SATUR = true;
		else
			it->u_prev.SATUR = it->u_iter.SATUR = it->u_next.SATUR = false;
	}

	for (it = tunnelCells.begin(); it != tunnelCells.end(); ++it)
	{
		const Skeleton_Props& props = props_sk[getSkeletonIdx(*it)];
		it->u_prev.p = it->u_iter.p = it->u_next.p = props.p_init;
		it->u_prev.p_bub = it->u_iter.p_bub = it->u_next.p_bub = props.p_bub;
		it->u_prev.s = it->u_iter.s = it->u_next.s = props.s_init;
		it->u_prev.t = it->u_iter.t = it->u_next.t = props.t_init;
		if (props.p_bub > props.p_init)
			it->u_prev.SATUR = it->u_iter.SATUR = it->u_next.SATUR = true;
		else
			it->u_prev.SATUR = it->u_iter.SATUR = it->u_next.SATUR = false;
	}
}

void GasOil_Perf_NIT::setUnused()
{
	vector<pair<int, int> >::iterator it;
	int idx;

	for (it = perfTunnels.begin(); it != perfTunnels.end(); ++it)
	{
		idx = it->first;
		for (int i = 0; i <= it->second; i++)
		{
			cells[idx].isUsed = false;
			idx += (cellsNum_z + 2);
		}
	}
}

void GasOil_Perf_NIT::buildTunnels()
{
	int counter = 0;
	double r, phi, z, hr, hphi, hz;

	for (int k = 0; k < perfTunnels.size(); k++)
	{
		r = cells[perfTunnels[k].first].r;
		for (int i = 0; i < perfTunnels[k].second; i++)
		{
			Cell& cell = cells[perfTunnels[k].first + (i + 1) * (cellsNum_z + 2)];
			r = cell.r;
			hr = cell.hr;

			// Right
			phi = cell.phi - cell.hphi / 2.0;
			if (phi < 0.0)
				phi += 2.0 * M_PI;
			hphi = 0.0;
			z = cell.z;
			hz = cell.hz;
			tunnelCells.push_back(Cell(counter, r, phi, z, hr, hphi, hz, k));
			tunnelNebrMap[getIdx(cell.num - (cellsNum_z + 2)*(cellsNum_r + 2))] = counter;
			nebrMap[counter++] = make_pair<int,int>( getIdx(cell.num - (cellsNum_z + 2)*(cellsNum_r + 2)), getIdx(cell.num - 2*(cellsNum_z + 2)*(cellsNum_r + 2)));

			// Top
			phi = cell.phi;
			hphi = cell.hphi;
			z = cell.z + cell.hz / 2.0;
			hz = 0.0;
			tunnelCells.push_back(Cell(counter, r, phi, z, hr, hphi, hz, k));
			tunnelNebrMap[getIdx(cell.num - 1)] = counter;
			nebrMap[counter++] = make_pair<int,int>( getIdx(cell.num - 1), getIdx(cell.num - 2));

			// Left
			phi = cell.phi + cell.hphi / 2.0;
			if (phi > 2.0 * M_PI)
				phi -= 2.0 * M_PI;
			hphi = 0.0;
			z = cell.z;
			hz = cell.hz;
			tunnelCells.push_back(Cell(counter, r, phi, z, hr, hphi, hz, k));
			tunnelNebrMap[getIdx(cell.num + (cellsNum_z + 2)*(cellsNum_r + 2))] = counter;
			nebrMap[counter++] = make_pair<int, int>(getIdx(cell.num + (cellsNum_z + 2)*(cellsNum_r + 2)), getIdx(cell.num + 2*(cellsNum_z + 2)*(cellsNum_r + 2)));

			// Bottom
			phi = cell.phi;
			hphi = cell.hphi;
			z = cell.z - cell.hz / 2.0;
			hz = 0.0;
			tunnelCells.push_back(Cell(counter, r, phi, z, hr, hphi, hz, k));
			tunnelNebrMap[getIdx(cell.num + 1)] = counter;
			nebrMap[counter++] = make_pair<int, int>(getIdx(cell.num + 1), getIdx(cell.num + 2));

			r += cell.hr / 2.0;
		}

		Cell& cell = cells[perfTunnels[k].first];

		// Central
		hr = 0.0;
		phi = cell.phi;
		hphi = cell.hphi;
		z = cell.z;
		hz = cell.hz;
		tunnelCells.push_back(Cell(counter, r, phi, z, hr, hphi, hz, k));
		tunnelNebrMap[getIdx(cell.num + (perfTunnels[k].second + 1) * (cellsNum_z + 2))] = counter;
		nebrMap[counter++] = make_pair<int, int>(getIdx(cell.num + (perfTunnels[k].second + 1) * (cellsNum_z + 2)), getIdx(cell.num + (perfTunnels[k].second + 2) * (cellsNum_z + 2)));
	}
}

void GasOil_Perf_NIT::setPerforated()
{
	height_perf = 0.0;
	for (int i = 0; i < tunnelCells.size(); i++)
	{
			Qcell[i] = 0.0;
			Cell& cell = tunnelCells[i];
			height_perf += cell.hphi * cell.r * cell.hz;
	}
}

void GasOil_Perf_NIT::setPeriod(int period)
{
	if (leftBoundIsRate)
	{
		Q_sum = rate[period];

		if (period == 0 || rate[period - 1] < EQUALITY_TOLERANCE) {
			map<int, double>::iterator it;
			for (it = Qcell.begin(); it != Qcell.end(); ++it)
			{
				Cell& cell = tunnelCells[it->first];
				it->second = Q_sum * cell.hphi * cell.r * cell.hz / height_perf;
			}
		}
		else {
			map<int, double>::iterator it;
			for (it = Qcell.begin(); it != Qcell.end(); ++it)
				it->second = it->second * Q_sum / rate[period - 1];
		}
	}
	else
	{
		Pwf = pwf[period];
		Q_sum = 0.0;
	}

	for(int i = 0; i < skeletonsNum; i++)
	{
		props_sk[i].radius_eff = props_sk[i].radiuses_eff[period];
		props_sk[i].perm_eff = props_sk[i].perms_eff[period];
		props_sk[i].skin = props_sk[i].skins[period];
	}
}

void GasOil_Perf_NIT::setRateDeviation(int num, double ratio)
{
	Qcell[num] += Q_sum * ratio;
}

/*-------------- Middle cells ------------------*/

double GasOil_Perf_NIT::solve_eq1(int cur)
{
	Cell& cell = getCell(cur);
	Cell* neighbor [6];
	getNeighborIdx(cell, neighbor);

	Var2phaseNIT& next = cell.u_next;
	Var2phaseNIT& prev = cell.u_prev;
	
	double H = 0.0;
	H = ( getPoro(next.p, cell) * next.s / getB_oil(next.p, next.p_bub, next.SATUR) - 
				getPoro(prev.p, cell) * prev.s / getB_oil(prev.p, prev.p_bub, prev.SATUR) );

	for(int i = 0; i < 6; i++)
	{
		Cell& beta = *neighbor[i];
		const Var2phaseNIT& upwd = getUpwindIdx(&cell, neighbor[i])->u_next;

		H += ht / cell.V * getTrans(cell, beta) * (next.p - beta.u_next.p) *
			getKr_oil(upwd.s) / props_oil.visc / getB_oil(upwd.p, upwd.p_bub, upwd.SATUR);
	}

	return H;
}

double GasOil_Perf_NIT::solve_eq1_dp(int cur, int beta)
{
	double upwind;
	Cell& cell = getCell(cur);
	Cell* neighbor [6];
	getNeighborIdx(cell, neighbor);

	Var2phaseNIT& next = cell.u_next;
	double Boil_upwd;
	double Boil = getB_oil(next.p, next.p_bub, next.SATUR);
	
	double H = 0.0;
	H = (next.s * getPoro_dp(cell) - 
		getPoro(next.p, cell) * next.s * getB_oil_dp(next.p, next.p_bub, next.SATUR) / Boil ) / Boil;

	for(int i = 0; i < 6; i++)
	{
		upwind = upwindIsCur(&cell, neighbor[i]);
		const Var2phaseNIT& upwd = getUpwindIdx(&cell, neighbor[i])->u_next;
		Boil_upwd = getB_oil(upwd.p, upwd.p_bub, upwd.SATUR);
		Cell& beta = *neighbor[i];

		H += ht / cell.V * getTrans(cell, beta) * 
			( getKr_oil(upwd.s) / props_oil.visc / Boil_upwd - 
			upwind * (next.p - beta.u_next.p) * getKr_oil(upwd.s) / props_oil.visc / Boil_upwd / Boil_upwd * getB_oil_dp(upwd.p, upwd.p_bub, upwd.SATUR) );
	}
	return H;
}

double GasOil_Perf_NIT::solve_eq1_ds(int cur, int beta)
{
	double upwind;	
	Cell& cell = getCell(cur);
	Cell* neighbor [6];
	getNeighborIdx(cell, neighbor);

	Var2phaseNIT& next = cell.u_next;
	
	double H = 0.0;
	H = getPoro(next.p, cell) / getB_oil(next.p, next.p_bub, next.SATUR);

	for(int i = 0; i < 6; i++)
	{
		upwind = upwindIsCur(&cell, neighbor[i]);
		const Var2phaseNIT& upwd = getUpwindIdx(&cell, neighbor[i])->u_next;
		Cell& beta = *neighbor[i];

		H += ht / cell.V * getTrans(cell, beta) * 
			upwind * (next.p - beta.u_next.p) * getKr_oil_ds(upwd.s) / props_oil.visc / getB_oil(upwd.p, upwd.p_bub, upwd.SATUR);
	}

	return H;
}

double GasOil_Perf_NIT::solve_eq1_dp_beta(int cur, int beta)
{
	Cell& cell = getCell(cur);
	Cell& nebr = getCell(cur, beta);

	double upwind = upwindIsCur(&cell, &nebr);
	const Var2phaseNIT& upwd = getUpwindIdx(&cell, &nebr)->u_next;
	double Boil_upwd = getB_oil(upwd.p, upwd.p_bub, upwd.SATUR);

	return -ht / cell.V * getTrans(cell, nebr) * 
			( getKr_oil(upwd.s) / props_oil.visc / Boil_upwd + 
			(1.0 - upwind) * (cell.u_next.p - nebr.u_next.p) * getKr_oil(upwd.s) / props_oil.visc / Boil_upwd / Boil_upwd * getB_oil_dp(upwd.p, upwd.p_bub, upwd.SATUR) );
}

double GasOil_Perf_NIT::solve_eq1_ds_beta(int cur, int beta)
{
	Cell& cell = getCell(cur);
	Cell& nebr = getCell(cur, beta);

	double upwind = upwindIsCur(&cell, &nebr);
	const Var2phaseNIT& upwd = getUpwindIdx(&cell, &nebr)->u_next;

	return ht / cell.V * getTrans(cell, nebr) * (1.0 - upwind) * (cell.u_next.p - nebr.u_next.p) *
			getKr_oil_ds(upwd.s) / props_oil.visc / getB_oil(upwd.p, upwd.p_bub, upwd.SATUR);
}

double GasOil_Perf_NIT::solve_eq2(int cur)
{
	Cell& cell = getCell(cur);
	Cell* neighbor [6];
	getNeighborIdx(cell, neighbor);

	Var2phaseNIT& next = cell.u_next;
	Var2phaseNIT& prev = cell.u_prev;

	double H = 0.0;
	H = getPoro(next.p, cell) * ( (1.0 - next.s) / getB_gas(next.p) + next.s * getRs(next.p, next.p_bub, next.SATUR) / getB_oil(next.p, next.p_bub, next.SATUR) ) -
				getPoro(prev.p, cell) * ( (1.0 - prev.s) / getB_gas(prev.p) + prev.s * getRs(prev.p, prev.p_bub, prev.SATUR) / getB_oil(prev.p, prev.p_bub, prev.SATUR) );

	for(int i = 0; i < 6; i++)
	{
		const Var2phaseNIT& upwd = getUpwindIdx(&cell, neighbor[i])->u_next;
		Cell& beta = *neighbor[i];

		H += ht / cell.V * getTrans(cell, beta) * (next.p - beta.u_next.p) * 
			( getKr_oil(upwd.s) * getRs(upwd.p, upwd.p_bub, upwd.SATUR) / props_oil.visc / getB_oil(upwd.p, upwd.p_bub, upwd.SATUR) +
			getKr_gas(upwd.s) / props_gas.visc / getB_gas(upwd.p) );
	}

	return H;
}

double GasOil_Perf_NIT::solve_eq2_dp(int cur, int beta)
{
	double upwind;	
	Cell& cell = getCell(cur);
	Cell* neighbor [6];
	getNeighborIdx(cell, neighbor);

	Var2phaseNIT& next = cell.u_next;
	double Boil_upwd, Bgas_upwd, rs_upwd;
	double Boil = getB_oil(next.p, next.p_bub, next.SATUR);
	double Bgas = getB_gas(next.p);
	double rs = getRs(next.p, next.p_bub, next.SATUR);
	
	double H = 0.0;
	H = ( (next.s * rs / Boil + (1.0 - next.s) / Bgas) * getPoro_dp(cell) - 
		getPoro(next.p, cell) * ( (1.0 - next.s) / Bgas / Bgas * getB_gas_dp(next.p) + 
		next.s * rs / Boil / Boil * getB_oil_dp(next.p, next.p_bub, next.SATUR) - 
		next.s / Boil * getRs_dp(next.p, next.p_bub, next.SATUR) ) );

	for(int i = 0; i < 6; i++)
	{
		upwind = upwindIsCur(&cell, neighbor[i]);
		const Var2phaseNIT& upwd = getUpwindIdx(&cell, neighbor[i])->u_next;
		Boil_upwd = getB_oil(upwd.p, upwd.p_bub, upwd.SATUR);
		Bgas_upwd = getB_gas(upwd.p);
		rs_upwd = getRs(upwd.p, upwd.p_bub, upwd.SATUR);
		Cell& beta = *neighbor[i];

		H += ht / cell.V * getTrans(cell, beta) * 
			( getKr_oil(upwd.s) * rs_upwd / props_oil.visc / Boil_upwd + getKr_gas(upwd.s) / props_gas.visc / Bgas_upwd + 
			upwind * (next.p - beta.u_next.p) * 
			( getKr_oil(upwd.s) / props_oil.visc / Boil_upwd * (getRs_dp(upwd.p, upwd.p_bub, upwd.SATUR) - rs_upwd * getB_oil_dp(upwd.p, upwd.p_bub, upwd.SATUR) / Boil_upwd) - 
			getKr_gas(upwd.s) / props_gas.visc / Bgas_upwd / Bgas_upwd * getB_gas_dp(upwd.p) ));
	}

	return H;
}

double GasOil_Perf_NIT::solve_eq2_ds(int cur, int beta)
{
	double upwind;	
	Cell& cell = getCell(cur);
	Cell* neighbor [6];
	getNeighborIdx(cell, neighbor);

	Var2phaseNIT& next = cell.u_next;
	
	double H = 0.0;
	H = getPoro(next.p, cell) * ( getRs(next.p, next.p_bub, next.SATUR) / getB_oil(next.p, next.p_bub, next.SATUR) - 1.0 / getB_gas(next.p) );

	for(int i = 0; i < 6; i++)
	{
		upwind = upwindIsCur(&cell, neighbor[i]);
		const Var2phaseNIT& upwd = getUpwindIdx(&cell, neighbor[i])->u_next;
		Cell& beta = *neighbor[i];

		H += ht / cell.V * getTrans(cell, beta) * upwind * (next.p - beta.u_next.p) * 
			( getRs(upwd.p, upwd.p_bub, upwd.SATUR) * getKr_oil_ds(upwd.s) / props_oil.visc / getB_oil(upwd.p, upwd.p_bub, upwd.SATUR) + 
			getKr_gas_ds(upwd.s) / props_gas.visc / getB_gas(upwd.p) );
	}

	return H;
}

double GasOil_Perf_NIT::solve_eq2_dp_beta(int cur, int beta)
{
	Cell& cell = getCell(cur);
	Cell& nebr = getCell(cur, beta);

	double upwind = upwindIsCur(&cell, &nebr);
	const Var2phaseNIT& upwd = getUpwindIdx(&cell, &nebr)->u_next;
	double Boil_upwd = getB_oil(upwd.p, upwd.p_bub, upwd.SATUR);
	double Bgas_upwd = getB_gas(upwd.p);
	double rs_upwd = getRs(upwd.p, upwd.p_bub, upwd.SATUR);

	return -ht / cell.V * getTrans(cell, nebr) * 
			( getKr_oil(upwd.s) * rs_upwd / props_oil.visc / Boil_upwd + getKr_gas(upwd.s) / props_gas.visc / Bgas_upwd - 
			(1.0 - upwind) * (cell.u_next.p - nebr.u_next.p) * 
			( getKr_oil(upwd.s) / props_oil.visc / Boil_upwd * (getRs_dp(upwd.p, upwd.p_bub, upwd.SATUR) - rs_upwd * getB_oil_dp(upwd.p, upwd.p_bub, upwd.SATUR) / Boil_upwd) - 
			getKr_gas(upwd.s) / props_gas.visc / Bgas_upwd / Bgas_upwd * getB_gas_dp(upwd.p) ));
}

double GasOil_Perf_NIT::solve_eq2_ds_beta(int cur, int beta)
{
	Cell& cell = getCell(cur);
	Cell& nebr = getCell(cur, beta);

	double upwind = upwindIsCur(&cell, &nebr);
	const Var2phaseNIT& upwd = getUpwindIdx(&cell, &nebr)->u_next;

	return ht / cell.V * getTrans(cell, nebr) * (1.0 - upwind) * (cell.u_next.p - nebr.u_next.p) *
		( getRs(upwd.p, upwd.p_bub, upwd.SATUR) * getKr_oil_ds(upwd.s) / props_oil.visc / getB_oil(upwd.p, upwd.p_bub, upwd.SATUR) +
		getKr_gas_ds(upwd.s) / props_gas.visc / getB_gas(upwd.p) );
}

double GasOil_Perf_NIT::solveH()
{
	double H = 0.0;
	double p1, p0;

	map<int,double>::iterator it = Qcell.begin();
	for(int i = 0; i < Qcell.size()-1; i++)	
	{
		p0 = tunnelCells[ it->first ].u_next.p;
		p1 = tunnelCells[ (++it)->first ].u_next.p;

		H += (p1 - p0) * (p1 - p0) / 2.0;
	}

	return H;
}

double GasOil_Perf_NIT::getRate(int cur)
{
	Cell& cell = tunnelCells[cur];
	Cell& nebr = getCell(nebrMap[cur].first);
	const Var2phaseNIT& upwd = getUpwindIdx(&cell, &nebr)->u_next;
	return getTrans(cell, nebr) * getKr_oil(upwd.s) / props_oil.visc / getBoreB_oil(cell.u_next.p, cell.u_next.p_bub, cell.u_next.SATUR) * (nebr.u_next.p - cell.u_next.p);
}
