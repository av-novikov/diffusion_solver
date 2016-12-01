#include "model/VPP2d/VPP2d.hpp"

#include <cassert>
#include <map>

#define P 0
#define S 1
#define C 2

using namespace vpp2d;
using std::vector;
using std::map;

VPP2d::VPP2d()
{
	x = new double[5 * Variable::size];
	y = new double[Variable::size];
}
VPP2d::~VPP2d()
{
	delete x;
	delete y;
}
void VPP2d::setProps(Properties& props)
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
	for (int j = 0; j < skeletonsNum; j++)
	{
		props_sk[j].perm_r = MilliDarcyToM2(props_sk[j].perm_r);
		props_sk[j].perm_z = MilliDarcyToM2(props_sk[j].perm_z);
	}

	periodsNum = props.timePeriods.size();
	for (int i = 0; i < periodsNum; i++)
	{
		period.push_back(props.timePeriods[i]);
		c.push_back(props.c[i]);
		if (leftBoundIsRate)
			rate.push_back(props.rates[i] / 86400.0);
		else
			pwf.push_back(props.pwf[i]);
		for (int j = 0; j < skeletonsNum; j++)
		{
			if (props_sk[j].radiuses_eff[i] > props.r_w)
				props_sk[j].perms_eff.push_back(MilliDarcyToM2(props.props_sk[j].perm_r * log(props.props_sk[j].radiuses_eff[i] / props.r_w) / (log(props.props_sk[j].radiuses_eff[i] / props.r_w) + props.props_sk[j].skins[i])));
			else
				props_sk[j].perms_eff.push_back(MilliDarcyToM2(props.props_sk[j].perm_r));
		}
	}

	// Temporal properties
	ht = props.ht;
	ht_min = props.ht_min;
	ht_max = props.ht_max;

	// Oil properties
	props_oil = props.props_o;
	props_oil.visc = cPToPaSec(props_oil.visc);

	// Water properties
	props_wat = props.props_w;
	props_wat.visc = cPToPaSec(props_wat.visc);

	depth_point = props.depth_point;

	makeDimLess();

	// Data sets
	props_oil.kr = setDataset(props.kr_o, 1.0, 1.0);
 	props_wat.kr = setDataset(props.kr_w, 1.0, 1.0);

//	props_oil.b = setDataset(props.B_o, P_dim / BAR_TO_PA, 1.0);
//	props_wat.b = setDataset(props.B_w, P_dim / BAR_TO_PA, 1.0);

//	props_wat.a = setDataset(props.a, P_dim / BAR_TO_PA, 1.0);
}
void VPP2d::makeDimLess()
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
	for (int i = 0; i < skeletonsNum; i++)
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

		for (int j = 0; j < periodsNum; j++)
		{
			props_sk[i].perms_eff[j] /= (R_dim * R_dim);
			props_sk[i].radiuses_eff[j] /= R_dim;
		}
	}

	Q_dim = R_dim * R_dim * R_dim / t_dim;
	for (int i = 0; i < periodsNum; i++)
	{
		period[i] /= t_dim;
		if (leftBoundIsRate)
			rate[i] /= Q_dim;
		else
			pwf[i] /= P_dim;
	}

	// Oil properties
	props_oil.visc /= (P_dim * t_dim);
	props_oil.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_oil.beta /= (1.0 / P_dim);
	props_oil.p_ref /= P_dim;

	// Water properties
	props_wat.visc /= (P_dim * t_dim);
	props_wat.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_wat.beta /= (1.0 / P_dim);
	props_wat.p_ref /= P_dim;
}
void VPP2d::checkSkeletons(const vector<Skeleton_Props>& props)
{
	vector<Skeleton_Props>::const_iterator it = props.begin();
	double tmp;
	int indxs = 0;

	assert(it->h2 - it->h1 == it->height);
	indxs += it->cellsNum_z;
	tmp = it->h2;
	++it;

	while (it != props.end())
	{
		assert(it->h1 == tmp);
		assert(it->h2 - it->h1 == it->height);
		indxs += it->cellsNum_z;
		tmp = it->h2;
		++it;
	}
	assert(indxs == cellsNum_z);
}
void VPP2d::buildGridLog()
{
	cells.reserve(cellsNum);

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
	cells.push_back(Cell(counter++, cm_r, cm_z, 0.0, 0.0));
	for (int i = 0; i < cellsNum_z; i++)
	{
		hz = (props_sk[skel_idx].h2 - props_sk[skel_idx].h1) / (double)(props_sk[skel_idx].cellsNum_z);
		cm_z += (cells[cells.size() - 1].hz + hz) / 2.0;

		cells.push_back(Cell(counter++, cm_r, cm_z, 0.0, hz));
		cells_z++;

		if (cells_z >= props_sk[skel_idx].cellsNum_z)
		{
			cells_z = 0;
			skel_idx++;
		}
	}
	cells.push_back(Cell(counter++, cm_r, cm_z + hz / 2.0, 0.0, 0.0));

	// Middle cells
	for (int j = 0; j < cellsNum_r; j++)
	{
		skel_idx = 0;	cells_z = 0;
		cm_z = props_sk[0].h1;
		cm_r = r_prev * (exp(logStep) + 1.0) / 2.0;
		hr = r_prev * (exp(logStep) - 1.0);

		cells.push_back(Cell(counter++, cm_r, cm_z, hr, 0.0));
		for (int i = 0; i < cellsNum_z; i++)
		{
			hz = (props_sk[skel_idx].h2 - props_sk[skel_idx].h1) / (double)(props_sk[skel_idx].cellsNum_z);
			cm_z += (cells[cells.size() - 1].hz + hz) / 2.0;

			cells.push_back(Cell(counter++, cm_r, cm_z, hr, hz));
			Volume += cells[cells.size() - 1].V;
			cells_z++;

			if (cells_z >= props_sk[skel_idx].cellsNum_z)
			{
				cells_z = 0;
				skel_idx++;
			}
		}
		cells.push_back(Cell(counter++, cm_r, cm_z + hz / 2.0, hr, 0.0));

		r_prev = r_prev * exp(logStep);
	}

	// Right border
	cm_z = props_sk[0].h1;
	cm_r = r_e;

	cells.push_back(Cell(counter++, cm_r, cm_z, 0.0, 0.0));
	skel_idx = 0;	cells_z = 0;
	for (int i = 0; i < cellsNum_z; i++)
	{
		hz = (props_sk[skel_idx].h2 - props_sk[skel_idx].h1) / (double)(props_sk[skel_idx].cellsNum_z);
		cm_z += (cells[cells.size() - 1].hz + hz) / 2.0;

		cells.push_back(Cell(counter++, cm_r, cm_z, 0.0, hz));
		cells_z++;

		if (cells_z >= props_sk[skel_idx].cellsNum_z)
		{
			cells_z = 0;
			skel_idx++;
		}
	}
	cells.push_back(Cell(counter++, cm_r, cm_z + hz / 2.0, 0.0, 0.0));
}
void VPP2d::setInitialState()
{
	vector<Cell>::iterator it;
	for (it = cells.begin(); it != cells.end(); ++it)
	{
		const Skeleton_Props& props = props_sk[getSkeletonIdx(*it)];
		it->props = const_cast<Skeleton_Props*>(&props);
		it->u_prev.p = it->u_iter.p = it->u_next.p = props.p_init;
		it->u_prev.s = it->u_iter.s = it->u_next.s = props.s_init;
		it->u_prev.c = it->u_iter.c = it->u_next.c = props.c_init;
	}
}
void VPP2d::setPerforated()
{
	height_perf = 0.0;
	vector<pair<int, int> >::iterator it;
	for (it = perfIntervals.begin(); it != perfIntervals.end(); ++it)
	{
		for (int i = it->first; i <= it->second; i++)
		{
			Qcell[i] = 0.0;
			height_perf += cells[i].hz;
		}
	}
}
void VPP2d::setPeriod(int period)
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

	Conc = c[period];

	for (int i = 0; i < skeletonsNum; i++)
	{
		props_sk[i].radius_eff = props_sk[i].radiuses_eff[period];
		props_sk[i].perm_eff = props_sk[i].perms_eff[period];
		props_sk[i].skin = props_sk[i].skins[period];
	}
}
void VPP2d::setRateDeviation(int num, double ratio)
{
	Qcell[num] += Q_sum * ratio;
}
double VPP2d::solveH()
{
	double H = 0.0;
	double p1, p0;

	map<int, double>::iterator it = Qcell.begin();
	for (int i = 0; i < Qcell.size() - 1; i++)
	{
		p0 = cells[it->first].u_next.p;
		p1 = cells[(++it)->first].u_next.p;

		H += (p1 - p0) * (p1 - p0) / 2.0;
	}

	return H;
}
double VPP2d::getRate(int cur)
{
	int neighbor = cur + cellsNum_z + 2;
	Variable& upwd = cells[getUpwindIdx(cur, neighbor)].u_next;
	Variable& next = cells[cur].u_next;
	return getTrans(cells[cur], cells[neighbor]) / props_wat.getViscosity(upwd.s).value() / props_wat.getBoreB(next.p).value() * (cells[neighbor].u_next.p - next.p);
}

void VPP2d::solve_eqMiddle(int cur)
{
	const Cell& cell = cells[cur];
	const Skeleton_Props& props = *cell.props;

	trace_on(mid);
	adouble h [Variable::size];
	//TapeVariable var[5];
	adouble var[5 * Variable::size];
	const Variable& prev = cell.u_prev;

	for (int i = 0; i < 5; i++)
	{
		var[i * Variable::size] <<= x[i * Variable::size];
		var[i * Variable::size + 1] <<= x[i * Variable::size + 1];
		var[i * Variable::size + 2] <<= x[i * Variable::size + 2];
	}

	//const TapeVariable& next = var[0];
	adouble* next = &var[0];

	/*h[0] = props.getPoro(next.p) * next.s / props_wat.getB(next.p) -
		props.getPoro(prev.p) * prev.s / props_wat.getB(prev.p);

	h[1] = props.getPoro(next.p) * (1.0 - next.s) / props_oil.getB(next.p) -
		props.getPoro(prev.p) * (1.0 - prev.s) / props_oil.getB(prev.p);

	h[2] = props.getPoro(next.p) * next.s * next.c / props_wat.getB(next.p) -
		props.getPoro(prev.p) * prev.s * prev.c / props_wat.getB(prev.p);*/

	h[0] = props.getPoro(next[P]) * next[S] / props_wat.getB(next[P]) -
		props.getPoro(prev.p) * prev.s / props_wat.getB(prev.p);

	h[1] = props.getPoro(next[P]) * (1.0 - next[S]) / props_oil.getB(next[P]) -
		props.getPoro(prev.p) * (1.0 - prev.s) / props_oil.getB(prev.p);

	h[2] = props.getPoro(next[P]) * next[S] * next[C] / props_wat.getB(next[P]) -
		props.getPoro(prev.p) * prev.s * prev.c / props_wat.getB(prev.p);

	int neighbor[4];
	getNeighborIdx(cur, neighbor);
	adouble tmp[Variable::size];
	for (int i = 0; i < 4; i++)
	{
		const Cell& beta = cells[neighbor[i]];
		//const int upwd_idx = (getUpwindIdx(cur, neighbor[i]) == cur) ? 0 : i + 1;
		const int nebr_idx = (i + 1) * Variable::size;
		adouble* nebr = &var[nebr_idx];
		//const TapeVariable& nebr = var[i + 1];
		//TapeVariable& upwd = var[upwd_idx];

		condassign(tmp[0], upwindIsCur(cur, neighbor[i]),
			(adouble)(ht / cell.V * getTrans(cell, beta)) * (next[P] - nebr[P]) *
			props_wat.getKr(next[S]) / props_wat.getViscosity(next[P]) / props_wat.getB(next[P]),
			(adouble)(ht / cell.V * getTrans(cell, beta)) * (next[P] - nebr[P]) *
			props_wat.getKr(nebr[S]) / props_wat.getViscosity(nebr[P]) / props_wat.getB(nebr[P]));
		/*tmp[0] = (adouble)(ht / cell.V * getTrans(cell, beta)) * (next[P] - nebr[P]) *
			props_wat.getKr(nebr[S]) / props_wat.getViscosity(nebr[P]) / props_wat.getB(nebr[P]);*/

		condassign(tmp[1], upwindIsCur(cur, neighbor[i]),
			(adouble)(ht / cell.V * getTrans(cell, beta)) * (next[P] - nebr[P]) *
			props_oil.getKr(next[S]) / props_oil.getViscosity(next[P]) / props_oil.getB(next[P]),
			(adouble)(ht / cell.V * getTrans(cell, beta)) * (next[P] - nebr[P]) *
			props_oil.getKr(nebr[S]) / props_oil.getViscosity(nebr[P]) / props_oil.getB(nebr[P]));

		condassign(tmp[2], upwindIsCur(cur, neighbor[i]),
			(adouble)(ht / cell.V * getTrans(cell, beta)) * (next[P] - nebr[P]) *
			props_wat.getKr(next[S]) / props_wat.getViscosity(next[P]) / props_wat.getB(next[P]) *
			next[C],
			(adouble)(ht / cell.V * getTrans(cell, beta)) * (next[P] - nebr[P]) *
			props_wat.getKr(nebr[S]) / props_wat.getViscosity(nebr[P]) / props_wat.getB(nebr[P]) *
			nebr[C]);

		h[0] += tmp[0];		h[1] += tmp[1];		h[2] += tmp[2];

		/*h[0] += (adouble)(ht / cell.V * getTrans(cell, beta)) * (next.p - nebr.p) *
			props_wat.getKr(upwd.s) / props_wat.getViscosity(upwd.p) / props_wat.getB(upwd.p);

		h[1] += (adouble)(ht / cell.V * getTrans(cell, beta)) * (next.p - nebr.p) *
			props_oil.getKr(upwd.s) / props_oil.getViscosity(upwd.p) / props_oil.getB(upwd.p);

		h[2] += (adouble)(ht / cell.V * getTrans(cell, beta)) * (next.p - nebr.p) *
			props_wat.getKr(upwd.s) / props_wat.getViscosity(upwd.p) / props_wat.getB(upwd.p) *
			upwd.c;*/
	}

	for (int i = 0; i < Variable::size; i++)
		h[i] >>= y[i];

	trace_off();
}
void VPP2d::solve_eqLeft(int cur)
{
	const Cell& cell = cells[cur];

	trace_on(left);
	adouble h[Variable::size];
	TapeVariable var[2];
	adouble leftIsRate = leftBoundIsRate;
	for (int i = 0; i < 2; i++)
	{
		var[i].p <<= x[i * Variable::size];
		var[i].s <<= x[i * Variable::size + 1];
		var[i].c <<= x[i * Variable::size + 2];
	}
	
	const TapeVariable& next = var[0];
	const TapeVariable& nebr = var[1];

	condassign(h[0], leftIsRate,
		(adouble)(getTrans(cells[cur], cells[cur+cellsNum_z+2])) / props_wat.getViscosity(next.p) / 
		props_wat.getBoreB(next.p) * (nebr.p - next.p) - (adouble)(Qcell[cur]),
		next.p - (adouble)(Pwf));

	h[1] = next.s - (adouble)(cell.props->s_oc);
	h[2] = next.c - (adouble)(Conc);

	for (int i = 0; i < Variable::size; i++)
		h[i] >>= y[i];

	trace_off();
}
void VPP2d::solve_eqRight(int cur)
{
	const Cell& cell = cells[cur];
	const Cell& beta1 = cells[cur - cellsNum_z - 2];
	const Cell& beta2 = cells[cur - 2 * cellsNum_z - 4];
	const Skeleton_Props& props = *cell.props;

	trace_on(right);
	adouble h[Variable::size];
	TapeVariable var[3];
	adouble rightIsPres = rightBoundIsPres;
	for (int i = 0; i < 3; i++)
	{
		var[i].p <<= x[i * Variable::size];
		var[i].s <<= x[i * Variable::size + 1];
		var[i].c <<= x[i * Variable::size + 2];
	}

	const TapeVariable& next = var[0];
	const TapeVariable& nebr1 = var[1];
	const TapeVariable& nebr2 = var[2];

	condassign(h[0], rightIsPres, next.p - (adouble)(props.p_out), next.p - (adouble)(nebr1.p));
	h[1] = (next.s - nebr1.s) / (adouble)(cell.r - beta1.r) - (nebr1.s - nebr2.s) / (adouble)(beta1.r - beta2.r);
	h[2] = (next.c - nebr1.c) / (adouble)(cell.r - beta1.r) - (nebr1.c - nebr2.c) / (adouble)(beta1.r - beta2.r);

	for (int i = 0; i < Variable::size; i++)
		h[i] >>= y[i];

	trace_off();
}
void VPP2d::solve_eqVertical(int cur)
{
	const Cell& cell = cells[cur];

	trace_on(vertical);
	adouble h[Variable::size];
	TapeVariable var[2];

	for (int i = 0; i < 2; i++)
	{
		var[i].p <<= x[i * Variable::size];
		var[i].s <<= x[i * Variable::size + 1];
		var[i].c <<= x[i * Variable::size + 2];
	}

	const TapeVariable& next = var[0];
	const TapeVariable& nebr = var[1];

	h[0] = next.p - nebr.p;
	h[1] = next.s - nebr.s;
	h[2] = next.c - nebr.c;

	for (int i = 0; i < Variable::size; i++)
		h[i] >>= y[i];

	trace_off();
}
void VPP2d::solveFixVar()
{
}
void VPP2d::setVariables(int cur)
{
	if (cur < cellsNum_z + 2)
	{
		// Left
		const Variable& next = cells[cur].u_next;
		const Variable& nebr = cells[cur + cellsNum_z + 2].u_next;
		
		for (int i = 0; i < Variable::size; i++)
		{
			x[i] = next.values[i];
			x[Variable::size + i] = nebr.values[i];
		}
	}
	else if (cur >= (cellsNum_z + 2) * (cellsNum_r + 1))
	{
		// Right
		const Variable& next = cells[cur].u_next;
		const Variable& nebr1 = cells[cur - cellsNum_z - 2].u_next;
		const Variable& nebr2 = cells[cur - 2 * cellsNum_z - 4].u_next;

		for (int i = 0; i < Variable::size; i++)
		{
			x[i] = next.values[i];
			x[Variable::size + i] = nebr1.values[i];
			x[2 * Variable::size + i] = nebr2.values[i];
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