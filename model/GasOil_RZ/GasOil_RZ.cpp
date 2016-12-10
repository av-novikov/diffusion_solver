#include "model/GasOil_RZ/GasOil_RZ.h"
#include "util/utils.h"

#include <cassert>

using namespace std;
using namespace gasOil_rz;

GasOil_RZ::GasOil_RZ()
{
	x = new double[stencil * Variable::size];
	y = new double[Variable::size-1];
}
GasOil_RZ::~GasOil_RZ()
{
	delete x;
	delete y;
}
void GasOil_RZ::setProps(Properties& props)
{
	leftBoundIsRate = props.leftBoundIsRate;
	rightBoundIsPres = props.rightBoundIsPres;

	// Setting grid properties
	r_w = props.r_w;
	r_e = props.r_e;
	cellsNum_r = props.cellsNum_r;
	cellsNum_z = props.cellsNum_z;
	cellsNum = (cellsNum_r + 2) * (cellsNum_z + 2);

	// Setting skeleton properties
	perfIntervals = props.perfIntervals;
	
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

	makeDimLess();
	
	// Data sets
	props_oil.kr = setDataset(props.kr_oil, 1.0, 1.0);
	props_gas.kr = setDataset(props.kr_gas, 1.0, 1.0);

	props_oil.b = setDataset(props.B_oil, P_dim / BAR_TO_PA, 1.0);
	props_gas.b = setDataset(props.B_gas, P_dim / BAR_TO_PA, 1.0);
	props_oil.Rs = setDataset(props.Rs, P_dim / BAR_TO_PA, 1.0);
}
void GasOil_RZ::checkSkeletons(const vector<Skeleton_Props>& props)
{
	vector<Skeleton_Props>::const_iterator it = props.begin();
	double tmp;
	int indxs = 0;

	assert( it->h2 - it->h1 == it->height );
	indxs += it->cellsNum_z;
	tmp = it->h2;
	++it;
	
	while(it != props.end())
	{
		assert( it->h1 == tmp );
		assert( it->h2 - it->h1 == it->height );
		indxs += it->cellsNum_z;
		tmp = it->h2;
		++it;
	}
	assert( indxs == cellsNum_z );
}
void GasOil_RZ::makeDimLess()
{
	// Main units
	R_dim = r_w;
	t_dim = 3600.0;
	P_dim = props_sk[0].p_init;

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
	props_oil.p_sat /= P_dim;

	// Gas properties
	props_gas.visc /= (P_dim * t_dim);
	props_gas.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);

	// Rest properties
	alpha /= t_dim;
	//depth_point = 0.0;
}
void GasOil_RZ::buildGridLog()
{
	cells.reserve( cellsNum );

	Volume = 0.0;
	int counter = 0;
	int skel_idx = 0, cells_z = 0;

	double r_prev = r_w;
	double logMax = log(r_e / r_w);
	double logStep = logMax / (double)cellsNum_r;
	
	double hz = 0.0;//(props_sk[skel_idx].h2 - props_sk[skel_idx].h1) / (double)(props_sk[skel_idx].cellsNum_z);
	double cm_z = props_sk[skel_idx].h1;
	double hr = r_prev * (exp(logStep) - 1.0);
	double cm_r = r_w;

	// Left border
	cells.push_back( Cell(counter++, cm_r, cm_z, 0.0, 0.0) );
	for(int i = 0; i < cellsNum_z; i++)
	{
		hz = (props_sk[skel_idx].h2 - props_sk[skel_idx].h1) / (double)(props_sk[skel_idx].cellsNum_z);
		cm_z += (cells[cells.size()-1].hz + hz) / 2.0;

		cells.push_back( Cell(counter++, cm_r, cm_z, 0.0, hz) );
		cells_z++;

		if(cells_z >= props_sk[skel_idx].cellsNum_z)
		{
			cells_z = 0;
			skel_idx++;
		}
	}
	cells.push_back( Cell(counter++, cm_r, cm_z+hz/2.0, 0.0, 0.0) );

	// Middle cells
	for(int j = 0; j < cellsNum_r; j++)
	{
		skel_idx = 0;	cells_z = 0;
		cm_z = props_sk[0].h1;
		cm_r = r_prev * (exp(logStep) + 1.0) / 2.0;
		hr = r_prev * (exp(logStep) - 1.0);

		cells.push_back( Cell(counter++, cm_r, cm_z, hr, 0.0) );
		for(int i = 0; i < cellsNum_z; i++)
		{
			hz = (props_sk[skel_idx].h2 - props_sk[skel_idx].h1) / (double)(props_sk[skel_idx].cellsNum_z);
			cm_z += (cells[cells.size()-1].hz + hz) / 2.0;

			cells.push_back( Cell(counter++, cm_r, cm_z, hr, hz) );
			Volume += cells[cells.size()-1].V;
			cells_z++;

			if(cells_z >= props_sk[skel_idx].cellsNum_z)
			{
				cells_z = 0;
				skel_idx++;
			}
		}
		cells.push_back( Cell(counter++, cm_r, cm_z+hz/2.0, hr, 0.0) );

		r_prev = r_prev * exp(logStep);
	}

	// Right border
	cm_z = props_sk[0].h1;
	cm_r = r_e;
	
	cells.push_back( Cell(counter++, cm_r, cm_z, 0.0, 0.0) );
	skel_idx = 0;	cells_z = 0;
	for(int i = 0; i < cellsNum_z; i++)
	{
		hz = (props_sk[skel_idx].h2 - props_sk[skel_idx].h1) / (double)(props_sk[skel_idx].cellsNum_z);
		cm_z += (cells[cells.size()-1].hz + hz) / 2.0;

		cells.push_back( Cell(counter++, cm_r, cm_z, 0.0, hz) );
		cells_z++;

		if(cells_z >= props_sk[skel_idx].cellsNum_z)
		{
			cells_z = 0;
			skel_idx++;
		}
	}
	cells.push_back( Cell(counter++, cm_r, cm_z+hz/2.0, 0.0, 0.0) );
}
void GasOil_RZ::setInitialState()
{
	vector<Cell>::iterator it;
	for(it = cells.begin(); it != cells.end(); ++it)
	{
		Skeleton_Props* props = &props_sk[ getSkeletonIdx(*it) ];
		it->u_prev.p = it->u_iter.p = it->u_next.p = props->p_init;
		it->u_prev.p_bub = it->u_iter.p_bub = it->u_next.p_bub = props->p_bub;
		it->u_prev.s = it->u_iter.s = it->u_next.s = props->s_init;
		if (props->p_init > props_oil.p_sat)
		{
			it->u_prev.SATUR = it->u_iter.SATUR = it->u_next.SATUR = false;
			it->u_prev.s = it->u_iter.s = it->u_next.s = 1.0;
		}
		else
			it->u_prev.SATUR = it->u_iter.SATUR = it->u_next.SATUR = true;

		it->props = props;
	}
}
void GasOil_RZ::setPerforated()
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
void GasOil_RZ::setPeriod(int period)
{
	if (leftBoundIsRate)
	{
		Q_sum = rate[period];

		if (period == 0 || rate[period - 1] < EQUALITY_TOLERANCE) {
			map<int, double>::iterator it;
			for (it = Qcell.begin(); it != Qcell.end(); ++it)
				it->second = Q_sum * cells[it->first].hz / height_perf;
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
void GasOil_RZ::setRateDeviation(int num, double ratio)
{
	Qcell[num] += Q_sum * ratio;
}
double GasOil_RZ::solveH()
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
double GasOil_RZ::getRate(int cur)
{
	int neighbor = cur + cellsNum_z + 2;
	Variable& upwd = cells[ getUpwindIdx(cur, neighbor) ].u_next;
	Variable& next = cells[cur].u_next;
	return getTrans(cells[cur], cells[neighbor]) * props_oil.getKr(upwd.s).value() / 
			props_oil.getViscosity(next.p).value() / props_oil.getBoreB(next.p, next.p_bub, next.SATUR).value() * 
			(cells[neighbor].u_next.p - next.p);
}

void GasOil_RZ::solve_eqMiddle(int cur)
{
	const Cell& cell = cells[cur];
	const Skeleton_Props& props = *cell.props;

	trace_on(mid);
	adouble h[Variable::size-1];
	TapeVariable var[stencil];
	const Variable& prev = cell.u_prev;

	for (int i = 0; i < stencil; i++)
	{
		var[i].p <<= x[i * Variable::size];
		var[i].s <<= x[i * Variable::size + 1];
		var[i].p_bub <<= x[i * Variable::size + 2];
	}

	const TapeVariable& next = var[0];
	adouble satur = cell.u_next.SATUR;

	condassign(h[0], satur,
		props.getPoro(next.p) * next.s / props_oil.getB(next.p, next.p_bub, satur) -
			props.getPoro(prev.p) * prev.s / props_oil.getB(prev.p, prev.p_bub, prev.SATUR),
		props.getPoro(next.p) / props_oil.getB(next.p, next.p_bub, satur) -
			props.getPoro(prev.p) / props_oil.getB(prev.p, prev.p_bub, prev.SATUR));

	condassign(h[1], satur,
		props.getPoro(next.p) * ((1.0 - next.s) / props_gas.getB(next.p) + 
			next.s * props_oil.getRs(next.p, next.p_bub, satur) / props_oil.getB(next.p, next.p_bub, satur)) -
			props.getPoro(prev.p) * ((1.0 - prev.s) / props_gas.getB(prev.p) + 
			prev.s * props_oil.getRs(prev.p, prev.p_bub, prev.SATUR) / props_oil.getB(prev.p, prev.p_bub, prev.SATUR)),
		props.getPoro(next.p) * props_oil.getRs(next.p, next.p_bub, satur) / props_oil.getB(next.p, next.p_bub, satur) -
			props.getPoro(prev.p) * props_oil.getRs(prev.p, prev.p_bub, prev.SATUR) / props_oil.getB(prev.p, prev.p_bub, prev.SATUR));

	int neighbor[4];
	getNeighborIdx(cur, neighbor);
	adouble tmp[Variable::size-1];
	for (int i = 0; i < 4; i++)
	{
		const Cell& beta = cells[neighbor[i]];
		const int upwd_idx = (getUpwindIdx(cur, neighbor[i]) == cur) ? 0 : i + 1;
		const int nebr_idx = (i + 1) * Variable::size;
		const TapeVariable& nebr = var[i + 1];
		TapeVariable& upwd = var[upwd_idx];

		h[0] += ht / cell.V * getTrans(cell, beta) * (next.p - nebr.p) *
			props_oil.getKr(upwd.s) / props_oil.getViscosity(upwd.p) /
			props_oil.getB(upwd.p, upwd.p_bub, upwd.SATUR);

		h[1] += ht / cell.V * getTrans(cell, beta) * (next.p - nebr.p) *
			(props_oil.getKr(upwd.s) * props_oil.getRs(upwd.p, upwd.p_bub, upwd.SATUR)
				/ props_oil.getViscosity(upwd.p) / props_oil.getB(upwd.p, upwd.p_bub, upwd.SATUR) +
			props_gas.getKr(upwd.s) / props_gas.getViscosity(upwd.p) / props_gas.getB(upwd.p));
	}

	for (int i = 0; i < Variable::size-1; i++)
		h[i] >>= y[i];

	trace_off();
}
void GasOil_RZ::solve_eqLeft(int cur)
{
	const Cell& cell = cells[cur];
	const Cell& beta1 = cells[cur + cellsNum_z + 2];
	const Cell& beta2 = cells[cur + 2 * cellsNum_z + 4];

	trace_on(left);
	adouble h[Variable::size-1];
	TapeVariable var[3];
	adouble leftIsRate = leftBoundIsRate;
	for (int i = 0; i < 3; i++)
	{
		var[i].p <<= x[i * Variable::size];
		var[i].s <<= x[i * Variable::size + 1];
		var[i].p_bub <<= x[i * Variable::size + 2];
	}

	const TapeVariable& next = var[0];
	const TapeVariable& nebr1 = var[1];
	const TapeVariable& nebr2 = var[2];
	const int upwd_idx = (getUpwindIdx(cur, cur + cellsNum_z + 2) == cur) ? 0 : 1;
	TapeVariable& upwd = var[upwd_idx];

	condassign(h[0], leftIsRate,
		(adouble)(getTrans(cells[cur], cells[cur + cellsNum_z + 2])) * props_oil.getKr(upwd.s) /
			props_oil.getViscosity(next.p) / props_oil.getBoreB(next.p, next.p_bub, next.SATUR) *
			(nebr1.p - next.p) - Qcell[cur],
		next.p - Pwf);

	adouble satur = cell.u_next.SATUR;
	condassign(h[1], satur,
		(next.s - nebr1.s) / (adouble)(cell.r - beta1.r) - (nebr1.s - nebr2.s) / (adouble)(beta1.r - beta2.r),
		(next.p_bub - nebr1.p_bub) / (adouble)(cell.r - beta1.r) - (nebr1.p_bub - nebr2.p_bub) / (adouble)(beta1.r - beta2.r));

	for (int i = 0; i < Variable::size-1; i++)
		h[i] >>= y[i];

	trace_off();
}
void GasOil_RZ::solve_eqRight(int cur)
{
	const Cell& cell = cells[cur];

	trace_on(right);
	adouble h[Variable::size-1];
	TapeVariable var[2];
	adouble rightIsPres = rightBoundIsPres;
	for (int i = 0; i < 2; i++)
	{
		var[i].p <<= x[i * Variable::size];
		var[i].s <<= x[i * Variable::size + 1];
		var[i].p_bub <<= x[i * Variable::size + 2];
	}

	const TapeVariable& next = var[0];
	const TapeVariable& nebr = var[1];

	condassign(h[0], rightIsPres, next.p - (adouble)(cell.props->p_out), next.p - (adouble)(nebr.p));
	adouble satur = cell.u_next.SATUR;
	condassign(h[1], satur, next.s - nebr.s, next.p_bub - nebr.p_bub);

	for (int i = 0; i < Variable::size-1; i++)
		h[i] >>= y[i];

	trace_off();
}
void GasOil_RZ::solve_eqVertical(int cur)
{
	const Cell& cell = cells[cur];

	trace_on(vertical);
	adouble h[Variable::size-1];
	TapeVariable var[2];

	for (int i = 0; i < 2; i++)
	{
		var[i].p <<= x[i * Variable::size];
		var[i].s <<= x[i * Variable::size + 1];
		var[i].p_bub <<= x[i * Variable::size + 2];
	}

	const TapeVariable& next = var[0];
	const TapeVariable& nebr = var[1];

	h[0] = next.p - nebr.p;
	adouble satur = cell.u_next.SATUR;
	condassign(h[1], satur, next.s - nebr.s, next.p_bub - nebr.p_bub);

	for (int i = 0; i < Variable::size-1; i++)
		h[i] >>= y[i];

	trace_off();
}
void GasOil_RZ::setVariables(int cur)
{
	if (cur < cellsNum_z + 2)
	{
		// Left
		const Variable& next = cells[cur].u_next;
		const Variable& nebr1 = cells[cur + cellsNum_z + 2].u_next;
		const Variable& nebr2 = cells[cur + 2 * cellsNum_z + 4].u_next;

		for (int i = 0; i < Variable::size; i++)
		{
			x[i] = next.values[i];
			x[Variable::size + i] = nebr1.values[i];
			x[2 * Variable::size + i] = nebr2.values[i];
		}
	}
	else if (cur >= (cellsNum_z + 2) * (cellsNum_r + 1))
	{
		// Right
		const Variable& next = cells[cur].u_next;
		const Variable& nebr = cells[cur - cellsNum_z - 2].u_next;

		for (int i = 0; i < Variable::size; i++)
		{
			x[i] = next.values[i];
			x[Variable::size + i] = nebr.values[i];
		}
	}
	else if (cur % (cellsNum_z + 2) == 0)
	{
		// Top
		const Variable& next = cells[cur].u_next;
		const Variable& nebr = cells[cur + 1].u_next;

		for (int i = 0; i < Variable::size; i++)
		{
			x[i] = next.values[i];
			x[Variable::size + i] = nebr.values[i];
		}
	}
	else if ((cur + 1) % (cellsNum_z + 2) == 0)
	{
		// Bottom
		const Variable& next = cells[cur].u_next;
		const Variable& nebr = cells[cur - 1].u_next;

		for (int i = 0; i < Variable::size; i++)
		{
			x[i] = next.values[i];
			x[Variable::size + i] = nebr.values[i];
		}
	}
	else
	{
		// Middle
		const Variable& next = cells[cur].u_next;
		int neighbor[4];
		getNeighborIdx(cur, neighbor);

		for (int i = 0; i < Variable::size; i++)
		{
			x[i] = next.values[i];

			for (int j = 0; j < 4; j++)
			{
				const Variable& nebr = cells[neighbor[j]].u_next;
				x[(j + 1) * Variable::size + i] = nebr.values[i];
			}
		}
	}
}