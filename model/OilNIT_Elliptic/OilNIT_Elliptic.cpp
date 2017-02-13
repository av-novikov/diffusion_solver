#include "model/OilNIT_Elliptic/OilNIT_Elliptic.hpp"
#include "util/utils.h"

#include "adolc/drivers/drivers.h"
#include "adolc/adolc.h"

using namespace std;
using namespace oilnit_elliptic;

OilNIT_Elliptic::OilNIT_Elliptic()
{
	x = new double[stencil * (Variable::size - 1)];
	y = new double[Variable::size - 1];

	jac = new double*[Variable::size - 1];
	for (int i = 0; i < Variable::size - 1; i++)
		jac[i] = new double[stencil * Variable::size];
};
OilNIT_Elliptic::~OilNIT_Elliptic()
{
	delete x;
	delete y;

	for (int i = 0; i < Variable::size - 1; i++)
		delete[] jac[i];
	delete[] jac;
};
void OilNIT_Elliptic::setProps(Properties& props)
{
	leftBoundIsRate = props.leftBoundIsRate;
	rightBoundIsPres = props.rightBoundIsPres;

	// Setting grid properties
	r_w = props.r_w;
	r_e = props.r_e;
	l = props.l;
	cellsNum_mu = props.cellsNum_mu;
	cellsNum_nu = props.cellsNum_nu;
	cellsNum_z = props.cellsNum_z;
	cellsNum = (cellsNum_mu + 2) * (cellsNum_z + 2) * cellsNum_nu;

	// Setting skeleton properties
	depth_point = props.depth_point;

	skeletonsNum = props.props_sk.size();
	props_sk = props.props_sk;
	int idx_z = 1;
	for (auto sk = props_sk.begin(); sk != props_sk.end(); ++sk)
	{
		sk->start_z = idx_z;	 idx_z += sk->cellsNum_z;
		sk->perm_mu = MilliDarcyToM2(sk->perm_mu);
		sk->perm_z = MilliDarcyToM2(sk->perm_z);
		if (sk->isWellHere)
			sk_well = sk;
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
				props_sk[j].perms_eff.push_back(MilliDarcyToM2(props.props_sk[j].perm_mu * log(props.props_sk[j].radiuses_eff[i] / props.r_w) / (log(props.props_sk[j].radiuses_eff[i] / props.r_w) + props.props_sk[j].skins[i])));
			else
				props_sk[j].perms_eff.push_back(MilliDarcyToM2(props.props_sk[j].perm_mu));
		}
	}

	// Temporal properties
	ht = props.ht;
	ht_min = props.ht_min;
	ht_max = props.ht_max;

	// Oil properties
	props_oil = props.props_oil;
	props_oil.visc = cPToPaSec(props_oil.visc);

	makeDimLess();

	Cell::a = l / 2;
};
void OilNIT_Elliptic::makeDimLess()
{
	// Main units
	R_dim = r_e / 10.0;
	t_dim = 3600.0;
	P_dim = props_sk[0].p_init;

	// Temporal properties
	ht /= t_dim;
	ht_min /= t_dim;
	ht_max /= t_dim;

	// Grid properties
	r_w /= R_dim;
	r_e /= R_dim;
	l /= R_dim;

	// Skeleton properties
	for (int i = 0; i < skeletonsNum; i++)
	{
		props_sk[i].perm_mu /= (R_dim * R_dim);
		props_sk[i].perm_z /= (R_dim * R_dim);

		props_sk[i].beta /= (1.0 / P_dim);
		props_sk[i].dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
		props_sk[i].h1 = (props_sk[i].h1 - depth_point) / R_dim;
		props_sk[i].h2 = (props_sk[i].h2 - depth_point) / R_dim;
		props_sk[i].h_well = (props_sk[i].h_well - depth_point) / R_dim;
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
	props_oil.oil.rho_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_oil.gas.rho_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_oil.beta /= (1.0 / P_dim);
	props_oil.p_sat /= P_dim;
};
void OilNIT_Elliptic::setInitialState()
{
	vector<Cell>::iterator it;
	for (it = cells.begin(); it != cells.end(); ++it)
	{
		const Skeleton_Props& props = props_sk[getSkeletonIdx(*it)];
		it->props = const_cast<Skeleton_Props*>(&props);
		it->u_prev.p = it->u_iter.p = it->u_next.p = props.p_init;
		it->u_prev.t = it->u_iter.t = it->u_next.t = props.t_init;
	}

	for (it = wellCells.begin(); it != wellCells.end(); ++it)
	{
		const Skeleton_Props& props = props_sk[getSkeletonIdx(*it)];
		it->props = const_cast<Skeleton_Props*>(&props);
		it->u_prev.p = it->u_iter.p = it->u_next.p = props.p_init;
		it->u_prev.t = it->u_iter.t = it->u_next.t = props.t_init;
	}
};
void OilNIT_Elliptic::buildGridLog()
{
	cells.reserve(cellsNum);

	Volume = 0.0;
	int counter = 0;
	auto sk_it = props_sk.begin();	int cells_z = 0;

	const double mu_w = asinh(4.0 * r_w / l);
	const double mu_e = asinh(2.0 * r_e / l);
	double hmu = (mu_e - mu_w) / (double)cellsNum_mu;
	const double hnu = 2.0 * M_PI / (double)cellsNum_nu;
	double hz1, hz2, hz;

	double r_prev = mu_w;
	double logMax = log(mu_e / mu_w);
	double logStep = logMax / (double)cellsNum_mu;

	double upcoord;
	double z_prev1 = r_w;	double z_prev2 = r_w;
	double logMax_z1 = log((sk_well->h_well - sk_well->h1) / r_w);
	double logMax_z2 = log((sk_well->h2 - sk_well->h_well) / r_w);
	double logStep_z1 = 2.0 * logMax_z1 / (double)(sk_well->cellsNum_z - 1);
	double logStep_z2 = 2.0 * logMax_z2 / (double)(sk_well->cellsNum_z - 1);

	double cm_mu = mu_w / 2;
	double cm_nu = 0.0;
	double cm_z = sk_it->h1;

	counter = 0;
	for (int k = 0; k < cellsNum_nu; k++)
	{
		sk_it = props_sk.begin();	cells_z = 0;

		r_prev = mu_w;
		logMax = log(mu_e / mu_w);
		logStep = logMax / (double)cellsNum_mu;
		hmu = r_prev * (exp(logStep) - 1.0);
		cm_mu = mu_w / 2;

		hz = hz1 = hz2 = 0.0;
		cm_z = sk_it->h1;
		cm_nu = (double)k * 2.0 * M_PI / (double)cellsNum_nu;
		z_prev1 = (sk_well->h_well - sk_well->h1) * pow((sk_well->h_well - sk_well->h1) / r_w, -2.0 / (double)(sk_well->cellsNum_z - 1));
		z_prev2 = r_w;

		// Left border
		cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z, mu_w, hnu, 0.0, TOP));
		for (int i = 0; i < cellsNum_z; i++)
		{
			if (!sk_it->isWellHere)
			{
				hz = (sk_it->h2 - sk_it->h1) / (double)(sk_it->cellsNum_z);
				cm_z += (cells[cells.size() - 1].hz + hz) / 2.0;
			}
			else
			{
				upcoord = cm_z + hz / 2.0 + r_w;
				if (upcoord < sk_it->h_well - EQUALITY_TOLERANCE)
				{
					hz = z_prev1 * (exp(logStep_z1) - 1.0);
					z_prev1 *= exp(-logStep_z1);
				}
				else if (upcoord < sk_it->h_well + r_w - EQUALITY_TOLERANCE)
					hz = 2 * r_w;
				else
				{
					hz = z_prev2 * (exp(logStep_z2) - 1.0);
					z_prev2 *= exp(logStep_z2);
				}
				cm_z += (cells[cells.size() - 1].hz + hz) / 2.0;
			}

			if (k % (cellsNum_nu / 2) != 0)
				cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z, mu_w, hnu, hz, MIDDLE));
			else
				cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z, mu_w, hnu, hz, MIDDLE_SIDE));
			cells_z++;

			if (cells_z >= sk_it->cellsNum_z)
			{
				cells_z = 0;
				++sk_it;
			}
		}
		cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z + hz / 2.0, mu_w, hnu, 0.0, BOTTOM));

		// Middle cells
		for (int j = 0; j < cellsNum_mu; j++)
		{
			sk_it = props_sk.begin();	cells_z = 0;	hz = 0.0;
			cm_z = sk_it->h1;
			cm_mu = r_prev * (exp(logStep) + 1.0) / 2.0;
			hmu = r_prev * (exp(logStep) - 1.0);
			z_prev1 = (sk_well->h_well - sk_well->h1) * pow((sk_well->h_well - sk_well->h1) / r_w, -2.0 / (double)(sk_well->cellsNum_z - 1));
			z_prev2 = r_w;

			cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z, hmu, hnu, 0.0, TOP));
			for (int i = 0; i < cellsNum_z; i++)
			{
				if (!sk_it->isWellHere)
				{
					hz = (sk_it->h2 - sk_it->h1) / (double)(sk_it->cellsNum_z);
					cm_z += (cells[cells.size() - 1].hz + hz) / 2.0;
				}
				else
				{
					upcoord = cm_z + hz / 2.0 + r_w;
					if (upcoord < sk_it->h_well - EQUALITY_TOLERANCE)
					{
						hz = z_prev1 * (exp(logStep_z1) - 1.0);
						z_prev1 *= exp(-logStep_z1);
					}
					else if (upcoord < sk_it->h_well + r_w - EQUALITY_TOLERANCE)
						hz = 2 * r_w;
					else
					{
						hz = z_prev2 * (exp(logStep_z2) - 1.0);
						z_prev2 *= exp(logStep_z2);
					}
					cm_z += (cells[cells.size() - 1].hz + hz) / 2.0;
				}

				cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z, hmu, hnu, hz, MIDDLE));
				Volume += cells[cells.size() - 1].V;
				cells_z++;

				if (cells_z >= sk_it->cellsNum_z)
				{
					cells_z = 0;
					++sk_it;
				}
			}
			cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z + hz / 2.0, hmu, hnu, 0.0, BOTTOM));

			r_prev *= exp(logStep);
		}

		// Right border
		sk_it = props_sk.begin();	cells_z = 0;	hz = 0.0;
		cm_z = sk_it->h1;
		cm_mu = mu_e;
		z_prev1 = (sk_well->h_well - sk_well->h1) * pow((sk_well->h_well - sk_well->h1) / r_w, -2.0 / (double)(sk_well->cellsNum_z - 1));
		z_prev2 = r_w;

		cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z, 0.0, hnu, 0.0, RIGHT));
		for (int i = 0; i < cellsNum_z; i++)
		{
			if (!sk_it->isWellHere)
			{
				hz = (sk_it->h2 - sk_it->h1) / (double)(sk_it->cellsNum_z);
				cm_z += (cells[cells.size() - 1].hz + hz) / 2.0;
			}
			else
			{
				upcoord = cm_z + hz / 2.0 + r_w;
				if (upcoord < sk_it->h_well - EQUALITY_TOLERANCE)
				{
					hz = z_prev1 * (exp(logStep_z1) - 1.0);
					z_prev1 *= exp(-logStep_z1);
				}
				else if (upcoord < sk_it->h_well + r_w - EQUALITY_TOLERANCE)
					hz = 2 * r_w;
				else
				{
					hz = z_prev2 * (exp(logStep_z2) - 1.0);
					z_prev2 *= exp(logStep_z2);
				}
				cm_z += (cells[cells.size() - 1].hz + hz) / 2.0;

			}

			cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z, 0.0, hnu, hz, RIGHT));
			cells_z++;

			if (cells_z >= sk_it->cellsNum_z)
			{
				cells_z = 0;
				++sk_it;
			}
		}
		cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z + hz / 2.0, 0.0, hnu, 0.0, RIGHT));
	}

	setUnused();
	buildWellCells();
}
void OilNIT_Elliptic::setUnused()
{
	for (const auto& sk : props_sk)
	{
		if (sk.isWellHere)
		{
			for (int j = 0; j < cellsNum_nu; j++)
			{
				const int idx = sk.start_z + (sk.cellsNum_z - 1) / 2 +
					(cellsNum_mu + 2) * (cellsNum_z + 2) * j;
				cells[idx].isUsed = false;
			}
		}
	}
}
void OilNIT_Elliptic::buildWellCells()
{
	int counter = 0;
	for (const auto& sk : props_sk)
	{
		if (sk.isWellHere)
		{
			// lateral
			for (int j = 0; j < cellsNum_nu; j++)
			{
				const int idx = sk.start_z + (sk.cellsNum_z - 1) / 2 +
					(cellsNum_mu + 2) * (cellsNum_z + 2) * j;
				const Cell& cell = cells[idx];

				wellCells.push_back(Cell(counter, cell.mu + cell.hmu / 2.0, cell.nu, cell.z,
					0.0, cell.hnu, cell.hz, WELL_LAT));
				wellNebrMap[idx + cellsNum_z + 2] = counter;
				nebrMap[counter++] = make_pair<int, int>(idx + cellsNum_z + 2, idx + 2 * cellsNum_z + 4);
			}

			// top
			for (int j = 0; j < cellsNum_nu; j++)
			{
				const int idx = sk.start_z + (sk.cellsNum_z - 1) / 2 +
					(cellsNum_mu + 2) * (cellsNum_z + 2) * j;
				const Cell& cell = cells[idx];

				wellCells.push_back(Cell(counter, cell.mu, cell.nu, cell.z - cell.hz,
					cell.hmu, cell.hnu, 0.0, WELL_TOP));
				wellNebrMap[idx - 1] = counter;
				nebrMap[counter++] = make_pair<int, int>(idx - 1, idx - 2);
			}

			// bot 
			for (int j = 0; j < cellsNum_nu; j++)
			{
				const int idx = sk.start_z + (sk.cellsNum_z - 1) / 2 +
					(cellsNum_mu + 2) * (cellsNum_z + 2) * j;
				const Cell& cell = cells[idx];

				wellCells.push_back(Cell(counter, cell.mu, cell.nu, cell.z + cell.hz,
					cell.hmu, cell.hnu, 0.0, WELL_BOT));
				wellNebrMap[idx + 1] = counter;
				nebrMap[counter++] = make_pair<int, int>(idx + 1, idx + 2);
			}
		}
	}
}
void OilNIT_Elliptic::setPerforated()
{
	height_perf = 0.0;

	// lateral
	for (int i = 0; i < cellsNum_nu; i++)
	{
		Qcell[i] = 0.0;
		if(i < int(cellsNum_nu / 4) + 1)
			Qcell_ellipse[i] = 0.0;
		const Cell& cell = wellCells[i];
		height_perf += Cell::getH(cell.mu, cell.nu) * cell.hnu * cell.hz;
	}

	// top
	for (int i = cellsNum_nu; i < 2 * cellsNum_nu; i++)
	{
		Qcell[i] = 0.0;
		if (i - cellsNum_nu < int(cellsNum_nu / 4) + 1)
			Qcell_ellipse[i] = 0.0;
		const Cell& cell = wellCells[i];
		height_perf += Cell::getH(cell.mu, cell.nu) * Cell::getH(cell.mu, cell.nu) * cell.hnu * cell.hmu;
	}

	// bot
	for (int i = 2 * cellsNum_nu; i < 3 * cellsNum_nu; i++)
	{
		Qcell[i] = 0.0;
		const Cell& cell = wellCells[i];
		height_perf += Cell::getH(cell.mu, cell.nu) * Cell::getH(cell.mu, cell.nu) * cell.hnu * cell.hmu;
	}
}
void OilNIT_Elliptic::setPeriod(int period)	
{
	if (leftBoundIsRate)
	{
		Q_sum = rate[period];

		if (period == 0 || rate[period - 1] < EQUALITY_TOLERANCE)
		{
			Q_sum_quater = 0.0;
			map<int, double>::iterator it, it_ell;
			double S;
			for (it = Qcell.begin(); it != Qcell.end(); ++it)
			{
				const Cell& cell = wellCells[it->first];
				if (cell.type == WELL_LAT)
					S = Cell::getH(cell.mu, cell.nu) * cell.hnu * cell.hz;
				else
					S = Cell::getH(cell.mu, cell.nu) * Cell::getH(cell.mu, cell.nu) * cell.hnu * cell.hmu;

				it->second = Q_sum * S / height_perf;
				it_ell = Qcell_ellipse.find(it->first);
				if (it_ell != Qcell_ellipse.end())
				{
					it_ell->second = it->second;
					Q_sum_quater += it->second;
				}
			}
		}
		else {
			Q_sum_quater = 0.0;
			map<int, double>::iterator it, it_ell;
			for (it = Qcell.begin(); it != Qcell.end(); ++it)
			{
				it->second = it->second * Q_sum / rate[period - 1];

				it_ell = Qcell_ellipse.find(it->first);
				if (it_ell != Qcell_ellipse.end())
				{
					it_ell->second = it->second;
					Q_sum_quater += it->second;
				}
			}
		}
	}
	else
	{
		Pwf = pwf[period];
		Q_sum = Q_sum_quater = 0.0;
	}

	for (int i = 0; i < skeletonsNum; i++)
	{
		props_sk[i].radius_eff_mu = asinh(props_sk[i].radiuses_eff[period] / Cell::a);
		props_sk[i].perm_eff_mu = props_sk[i].perms_eff[period];
		props_sk[i].radius_eff_z = props_sk[i].radiuses_eff[period];
		props_sk[i].perm_eff_z = props_sk[i].perm_z / props_sk[i].perm_mu * props_sk[i].perms_eff[period];
		props_sk[i].skin = props_sk[i].skins[period];
	}
};
void OilNIT_Elliptic::setRateDeviation(int num, double ratio)
{
	const auto indices = getSymmetricalWellIndices(num);

	Qcell_ellipse[num] += Q_sum_quater * ratio;
	for (const auto& idx : indices)
		Qcell[idx] += Q_sum_quater * ratio;
}
double OilNIT_Elliptic::solveH()
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
double OilNIT_Elliptic::getRate(int cur) const
{
	const Cell& cell = wellCells[cur];
	const Cell& beta = cells[nebrMap.at(cur).first];
	const Variable& next = wellCells[cur].u_next;
	return getTrans(cell, beta) * props_oil.getDensity(next.p).value() / props_oil.oil.rho_stc / props_oil.getViscosity(next.p).value() * (beta.u_next.p - next.p);
}

void OilNIT_Elliptic::solve_eqMiddle(const Cell& cell, const int val)
{
	const Skeleton_Props& props = *cell.props;
	const int nebrNum = (cell.type == MIDDLE) ? 6 : 5;

	if (val == 0)
	{
		trace_on(tmid);

		adouble h[Variable::size - 1];
		TapeVariableNIT var[stencil];
		const Variable& prev = cell.u_prev;

		for (int i = 0; i < nebrNum + 1; i++)
			var[i].t <<= x[i];

		const TapeVariableNIT& next = var[0];

		h[0] = getCn(cell) * (next.t - prev.t) - getAd(cell) * (cell.u_next.p - prev.p);

		Cell* neighbor[6];
		getNeighbors(cell, neighbor);
		adouble tmp[Variable::size - 1];
		for (int i = 0; i < nebrNum; i++)
		{
			const Cell& beta = *neighbor[i];
			const TapeVariableNIT& nebr = var[i + 1];

			//h[0] += ht * ((max(0.0, getA(cell, neighbor, ) + getTherCond(cell, beta) / cell.V) * (next.t - nebr.t) + b * (cell.u_next.p - beta.u_next.p)) / getDistance(cell, beta);
		}

		h[0] >>= y[0];

		trace_off();
	}
	else if (val == 1)
	{
		trace_on(mid);

		adouble h[Variable::size - 1];
		TapeVariable var[stencil];
		const Variable& prev = cell.u_prev;

		for (int i = 0; i < nebrNum + 1; i++)
			var[i].p <<= x[i];

		const TapeVariable& next = var[0];

		h[0] = props.getPoro(next.p) * props_oil.getDensity(next.p) - props.getPoro(prev.p) * props_oil.getDensity(prev.p);

		Cell* neighbor[6];
		getNeighbors(cell, neighbor);
		adouble tmp[Variable::size - 1];
		for (int i = 0; i < nebrNum; i++)
		{
			const Cell& beta = *neighbor[i];
			const TapeVariable& nebr = var[i + 1];

			h[0] += ht / cell.V * getTrans(cell, beta) * (next.p - nebr.p) * getDensity(next.p, cell, nebr.p, beta) / props_oil.getViscosity(next.p);
		}

		h[0] >>= y[0];

		trace_off();
	}
}
void OilNIT_Elliptic::solve_eqWell(const Cell& cell, const int val)
{
	if (val == 0)
	{
		const Cell& beta1 = cells[nebrMap[cell.num].first];
		const Cell& beta2 = cells[nebrMap[cell.num].second];

		double dist1, dist2;
		if (cell.hz == 0)
		{
			dist1 = cell.z - beta1.z;
			dist2 = beta1.z - beta2.z;
		}
		else
		{
			dist1 = Cell::getH(cell.mu, cell.nu) * (cell.mu - beta1.mu);
			dist2 = (Cell::getH(beta1.mu, beta1.nu) * beta1.hmu + Cell::getH(beta2.mu, beta2.nu) * beta2.hnu) / 2.0;
		}

		trace_on(tleft);
		adouble h[Variable::size - 1];
		TapeVariableNIT var[TLstencil];
		for (int i = 0; i < TLstencil; i++)
			var[i].t <<= x[i];

		const TapeVariableNIT& next = var[0];
		const TapeVariableNIT& nebr1 = var[1];
		const TapeVariableNIT& nebr2 = var[2];

		h[0] = (next.t - nebr1.t) / (adouble)(dist1)-(nebr1.t - nebr2.t) / (adouble)(dist2);
		h[0] >>= y[0];

		trace_off();
	}
	else if (val == 1)
	{
		const Cell& beta = cells[nebrMap[cell.num].first];
		double dist;
		if (cell.hz == 0)
			dist = cell.z - beta.z;
		else
			dist = Cell::getH(cell.mu, cell.nu) * (cell.mu - beta.mu);

		trace_on(left);
		adouble h[Variable::size - 1];
		TapeVariable var[Lstencil];
		adouble leftIsRate = leftBoundIsRate;
		for (int i = 0; i < Lstencil; i++)
			var[i].p <<= x[i];

		const TapeVariable& next = var[0];
		const TapeVariable& nebr = var[1];

		condassign(h[0], leftIsRate,
			(adouble)(getTrans(cell, beta)) * props_oil.getDensity(next.p) / props_oil.oil.rho_stc /
			props_oil.getViscosity(next.p) * (nebr.p - next.p) - Qcell[cell.num],
			next.p - Pwf);

		h[0] >>= y[0];

		trace_off();
	}
}
void OilNIT_Elliptic::solve_eqRight(const Cell& cell, const int val)
{
	if (val == 0)
	{
		trace_on(tright);
		adouble h[Variable::size - 1];
		TapeVariableNIT var[2];
		for (int i = 0; i < Rstencil; i++)
			var[i].t <<= x[i];

		const TapeVariableNIT& next = var[0];
		const TapeVariableNIT& nebr = var[1];

		h[0] = next.t - (adouble)(cell.props->t_out);
		h[0] >>= y[0];

		trace_off();
	}
	else if (val == 1)
	{
		trace_on(right);
		adouble h[Variable::size - 1];
		TapeVariable var[2];
		adouble rightIsPres = rightBoundIsPres;
		for (int i = 0; i < Rstencil; i++)
			var[i].p <<= x[i];

		const TapeVariable& next = var[0];
		const TapeVariable& nebr = var[1];

		condassign(h[0], rightIsPres, next.p - (adouble)(cell.props->p_out), next.p - (adouble)(nebr.p));

		h[0] >>= y[0];

		trace_off();
	}
}
void OilNIT_Elliptic::solve_eqVertical(const Cell& cell, const int val)
{
	if(val == 0)
	{
		trace_on(tvertical);
		adouble h[Variable::size - 1];
		TapeVariableNIT var[Vstencil];

		for (int i = 0; i < Vstencil; i++)
			var[i].t <<= x[i];

		const TapeVariableNIT& next = var[0];
		const TapeVariableNIT& nebr = var[1];

		h[0] = next.t - nebr.t;
		h[0] >>= y[0];

		trace_off();
	}
	else if (val == 1)
	{
		trace_on(vertical);
		adouble h[Variable::size - 1];
		TapeVariable var[Vstencil];

		for (int i = 0; i < Vstencil; i++)
			var[i].p <<= x[i];

		const TapeVariable& next = var[0];
		const TapeVariable& nebr = var[1];

		h[0] = next.p - nebr.p;
		h[0] >>= y[0];

		trace_off();
	}
}
void OilNIT_Elliptic::setVariables(const Cell& cell, const int val)
{
	assert(cell.isUsed);

	if (cell.type == WELL_LAT || cell.type == WELL_TOP || cell.type == WELL_BOT) // Well
	{
		if (val == 0)
		{
			const Variable& next = cell.u_next;
			const Variable& nebr1 = cells[nebrMap[cell.num].first].u_next;
			const Variable& nebr2 = cells[nebrMap[cell.num].second].u_next;

			x[0] = next.values[val];
			x[1] = nebr1.values[val];
			x[2] = nebr2.values[val];

			solve_eqWell(cell, val);
			jacobian(left, Variable::size - 1, (Variable::size - 1) * TLstencil, x, jac);
		}
		else if (val == 1)
		{
			const Variable& next = cell.u_next;
			const Variable& nebr = cells[nebrMap[cell.num].first].u_next;

			x[0] = next.values[val];
			x[1] = nebr.values[val];

			solve_eqWell(cell, val);
			jacobian(left, Variable::size - 1, (Variable::size - 1) * Lstencil, x, jac);
		}
	}
	else if (cell.type == RIGHT) // Right
	{
		const Variable& next = cell.u_next;
		const Variable& nebr = cells[cell.num - cellsNum_z - 2].u_next;

		x[0] = next.values[ val ];
		x[1] = nebr.values[ val ];

		solve_eqRight(cell, val);
		jacobian(right, Variable::size - 1, (Variable::size - 1) * Rstencil, x, jac);
	}
	else if (cell.type == TOP) // Top
	{
		const Variable& next = cell.u_next;
		const Variable& nebr = cells[cell.num + 1].u_next;

		x[0] = next.values[ val ];
		x[1] = nebr.values[ val ];

		solve_eqVertical(cell, val);
		jacobian(vertical, Variable::size - 1, (Variable::size - 1) * Vstencil, x, jac);
	}
	else if (cell.type == BOTTOM) // Bottom
	{
		const Variable& next = cell.u_next;
		const Variable& nebr = cells[cell.num - 1].u_next;

		x[0] = next.values[ val ];
		x[1] = nebr.values[ val ];

		solve_eqVertical(cell, val);
		jacobian(vertical, Variable::size - 1, (Variable::size - 1) * Vstencil, x, jac);
	}
	else // Middle
	{
		const int nebrNum = (cell.type == MIDDLE) ? 6 : 5;
		const Variable& next = cell.u_next;
		Cell* neighbor[6];
		getNeighbors(cell, neighbor);

		x[0] = next.values[ val ];
		for (int j = 0; j < nebrNum; j++)
				x[(j + 1)] = neighbor[j]->u_next.values[ val ];

		solve_eqMiddle(cell, val);
		jacobian(mid, Variable::size - 1, (Variable::size - 1) * (nebrNum + 1), x, jac);
	}
}