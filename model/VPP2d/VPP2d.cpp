#include "model/VPP2d/VPP2d.hpp"

#include <cassert>
#include <map>

using namespace vpp2d;
using std::vector;
using std::map;

VPP2d::VPP2d()
{}

VPP2d::~VPP2d()
{}

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
	props_o = props.props_o;
	props_o.visc = cPToPaSec(props_o.visc);

	// Water properties
	props_w = props.props_w;
	props_w.visc = cPToPaSec(props_w.visc);

	depth_point = props.depth_point;

	makeDimLess();

	// Data sets
	props_o.kr = setDataset(props.kr_o, 1.0, 1.0);
	props_w.kr = setDataset(props.kr_w, 1.0, 1.0);

	props_o.b = setDataset(props.B_o, P_dim / BAR_TO_PA, 1.0);
	props_w.b = setDataset(props.B_w, P_dim / BAR_TO_PA, 1.0);

	a = setDataset(props.a, P_dim / BAR_TO_PA, 1.0);
}

void VPP2d::makeDimLess()
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
	props_o.visc /= (P_dim * t_dim);
	props_o.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_o.beta /= (1.0 / P_dim);

	// Water properties
	props_w.visc /= (P_dim * t_dim);
	props_w.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_w.beta /= (1.0 / P_dim);
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
	return getTrans(cells[cur], cells[neighbor]) * props_w.getKr(upwd.s) / props_o.visc / props_w.getBoreB(next.p) * (cells[neighbor].u_next.p - next.p);
}