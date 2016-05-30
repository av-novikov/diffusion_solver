#include "model/3D/Perforation/Oil_Perf_NIT.h"
#include "util/utils.h"

#include <cassert>

using namespace std;
using namespace std::placeholders;
using namespace oil_perf_nit;


Oil_Perf_NIT::Oil_Perf_NIT()
{
	isWriteSnaps = false;

	// Middle
	middleFoo.mat.resize(1);

	middleFoo.rhs.push_back(bind(&Oil_Perf_NIT::solve_eq, this, _1));
	middleFoo.mat[0].push_back(bind(&Oil_Perf_NIT::solve_eq_dp, this, _1, _2));
	middleFoo.mat[0].push_back(bind(&Oil_Perf_NIT::solve_eq_dp_beta, this, _1, _2));
	middleFoo.mat[0].push_back(bind(&Oil_Perf_NIT::solve_eq_dp_beta, this, _1, _2));
	middleFoo.mat[0].push_back(bind(&Oil_Perf_NIT::solve_eq_dp_beta, this, _1, _2));
	middleFoo.mat[0].push_back(bind(&Oil_Perf_NIT::solve_eq_dp_beta, this, _1, _2));
	middleFoo.mat[0].push_back(bind(&Oil_Perf_NIT::solve_eq_dp_beta, this, _1, _2));
	middleFoo.mat[0].push_back(bind(&Oil_Perf_NIT::solve_eq_dp_beta, this, _1, _2));

	// Left
	leftFoo.mat.resize(1);

	leftFoo.rhs.push_back(bind(&Oil_Perf_NIT::solve_eqLeft, this, _1));
	leftFoo.mat[0].push_back(bind(&Oil_Perf_NIT::solve_eqLeft_dp, this, _1, _2));
	leftFoo.mat[0].push_back(bind(&Oil_Perf_NIT::solve_eqLeft_dp_beta, this, _1, _2));

	// Right
	rightFoo.mat.resize(1);

	rightFoo.rhs.push_back(bind(&Oil_Perf_NIT::solve_eqRight, this, _1));
	rightFoo.mat[0].push_back(bind(&Oil_Perf_NIT::solve_eqRight_dp, this, _1, _2));
	rightFoo.mat[0].push_back(bind(&Oil_Perf_NIT::solve_eqRight_dp_beta, this, _1, _2));

	// Top
	topFoo.mat.resize(1);

	topFoo.rhs.push_back(bind(&Oil_Perf_NIT::solve_eqTop, this, _1));
	topFoo.mat[0].push_back(bind(&Oil_Perf_NIT::solve_eqTop_dp, this, _1, _2));
	topFoo.mat[0].push_back(bind(&Oil_Perf_NIT::solve_eqTop_dp_beta, this, _1, _2));

	// Bot
	botFoo.mat.resize(1);

	botFoo.rhs.push_back(bind(&Oil_Perf_NIT::solve_eqBot, this, _1));
	botFoo.mat[0].push_back(bind(&Oil_Perf_NIT::solve_eqBot_dp, this, _1, _2));
	botFoo.mat[0].push_back(bind(&Oil_Perf_NIT::solve_eqBot_dp_beta, this, _1, _2));
}

Oil_Perf_NIT::~Oil_Perf_NIT()
{
}

void Oil_Perf_NIT::setProps(Properties& props)
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

	alpha = props.alpha;
	depth_point = props.depth_point;

	makeDimLess();
}

void Oil_Perf_NIT::checkSkeletons(const vector<Skeleton_Props>& props)
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

void Oil_Perf_NIT::makeDimLess()
{
	// Main units
	R_dim = r_w * 10.0;
	t_dim = 3600.0;
	P_dim = props_sk[0].p_init;
	if (props_sk[0].t_init != 0.0)
		T_dim = props_sk[0].t_init;
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

	// Rest properties
	alpha /= t_dim;
	//depth_point = 0.0;
}

void Oil_Perf_NIT::buildGridLog()
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

void Oil_Perf_NIT::setInitialState()
{
	vector<Cell>::iterator it;
	for(it = cells.begin(); it != cells.end(); ++it)
	{
		const Skeleton_Props& props = props_sk[ getSkeletonIdx(*it) ];
		it->u_prev.p = it->u_iter.p = it->u_next.p = props.p_init;
		it->u_prev.t = it->u_iter.t = it->u_next.t = props.t_init;
	}

	for (it = tunnelCells.begin(); it != tunnelCells.end(); ++it)
	{
		const Skeleton_Props& props = props_sk[getSkeletonIdx(*it)];
		it->u_prev.p = it->u_iter.p = it->u_next.p = props.p_init;
		it->u_prev.t = it->u_iter.t = it->u_next.t = props.t_init;
	}
}

void Oil_Perf_NIT::setUnused()
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

void Oil_Perf_NIT::buildTunnels()
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

void Oil_Perf_NIT::setPerforated()
{
	height_perf = 0.0;
	for (int i = 0; i < tunnelCells.size(); i++)
	{
			Qcell[i] = 0.0;
			Cell& cell = tunnelCells[i];
			height_perf += cell.hphi * cell.r * cell.hz;
	}
}

void Oil_Perf_NIT::setPeriod(int period)
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

void Oil_Perf_NIT::setRateDeviation(int num, double ratio)
{
	Qcell[num] += Q_sum * ratio;
}

/*-------------- Middle cells ------------------*/

double Oil_Perf_NIT::solve_eq(int cur)
{
	Cell& cell = getCell(cur);
	Cell* neighbor [6];
	getNeighborIdx(cell, neighbor);

	Var1phaseNIT& next = cell.u_next;
	Var1phaseNIT& prev = cell.u_prev;
	
	double H = 0.0;
	H = getPoro(next.p, cell) * getRho(next.p) - getPoro(prev.p, cell) * getRho(prev.p);

	for(int i = 0; i < 6; i++)
	{
		Cell& beta = *neighbor[i];
		const Var1phaseNIT& upwd = getUpwindIdx(&cell, neighbor[i])->u_next;

		H += ht / cell.V * getTrans(cell, beta) * (next.p - beta.u_next.p) / props_oil.visc * getRho(cell, beta);
	}

	return H;
}

double Oil_Perf_NIT::solve_eq_dp(int cur, int beta)
{
	double upwind;
	Cell& cell = getCell(cur);
	Cell* neighbor [6];
	getNeighborIdx(cell, neighbor);

	Var1phaseNIT& next = cell.u_next;
	
	double H = 0.0;
	H = getPoro(next.p, cell) * getRho_dp() + getRho(next.p) * getPoro_dp(cell);

	for(int i = 0; i < 6; i++)
	{
		upwind = upwindIsCur(&cell, neighbor[i]);
		const Var1phaseNIT& upwd = getUpwindIdx(&cell, neighbor[i])->u_next;
		Cell& beta = *neighbor[i];

		H += ht / cell.V / props_oil.visc * getTrans(cell, beta) *
			(getRho(cell, beta) + (next.p - beta.u_next.p) * getRho_dp(cell, beta));
	}
	return H;
}

double Oil_Perf_NIT::solve_eq_dp_beta(int cur, int beta)
{
	Cell& cell = getCell(cur);
	Cell& nebr = getCell(cur, beta);

	double upwind = upwindIsCur(&cell, &nebr);

	return ht / cell.V / props_oil.visc * getTrans(cell, nebr) *
		( (cell.u_next.p - nebr.u_next.p) * getRho_dp_beta(cell, nebr) - getRho(cell, nebr));
}

double Oil_Perf_NIT::solveH()
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

double Oil_Perf_NIT::getRate(int cur)
{
	Cell& cell = tunnelCells[cur];
	Cell& nebr = getCell(nebrMap[cur].first);
	const Var1phaseNIT& upwd = getUpwindIdx(&cell, &nebr)->u_next;
	return getTrans(cell, nebr) / props_oil.visc / getBoreB_oil(cell.u_next.p) * (nebr.u_next.p - cell.u_next.p);
}
