#include "model/GasOil_Elliptic/GasOil_Elliptic.hpp"
#include "util/utils.h"

#include <cassert>

using namespace std;
using namespace gasOil_elliptic;

GasOil_Elliptic::GasOil_Elliptic()
{
};
GasOil_Elliptic::~GasOil_Elliptic()
{
};
void GasOil_Elliptic::setProps(Properties& props)
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
	perfIntervals = props.perfIntervals;
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

	makeDimLess();

	Cell::a = l / 2;
};
void GasOil_Elliptic::makeDimLess()
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
		props_sk[i].p_bub /= P_dim;

		for (int j = 0; j < periodsNum; j++)
		{
			props_sk[i].perms_eff[j] /= (R_dim * R_dim);
			props_sk[i].radiuses_eff[j] /= R_dim;
		}
	}
};
void GasOil_Elliptic::setInitialState()
{
	vector<Cell>::iterator it;
	for (it = cells.begin(); it != cells.end(); ++it)
	{
		const Skeleton_Props& props = props_sk[getSkeletonIdx(*it)];
		it->props = const_cast<Skeleton_Props*>(&props);
		it->u_prev.p = it->u_iter.p = it->u_next.p = props.p_init;
		it->u_prev.p_bub = it->u_iter.p_bub = it->u_next.p_bub = props.p_bub;
		it->u_prev.s = it->u_iter.s = it->u_next.s = props.s_init;
		if (props.p_init > props_oil.p_sat)
		{
			it->u_prev.SATUR = it->u_iter.SATUR = it->u_next.SATUR = false;
			it->u_prev.s = it->u_iter.s = it->u_next.s = 1.0;
		}
		else
			it->u_prev.SATUR = it->u_iter.SATUR = it->u_next.SATUR = true;
	}

	for (it = wellCells.begin(); it != wellCells.end(); ++it)
	{
		const Skeleton_Props& props = props_sk[getSkeletonIdx(*it)];
		it->props = const_cast<Skeleton_Props*>(&props);
		it->u_prev.p = it->u_iter.p = it->u_next.p = props.p_init;
		it->u_prev.p_bub = it->u_iter.p_bub = it->u_next.p_bub = props.p_bub;
		it->u_prev.s = it->u_iter.s = it->u_next.s = props.s_init;
		if (props.p_init > props_oil.p_sat)
		{
			it->u_prev.SATUR = it->u_iter.SATUR = it->u_next.SATUR = false;
			it->u_prev.s = it->u_iter.s = it->u_next.s = 1.0;
		}
		else
			it->u_prev.SATUR = it->u_iter.SATUR = it->u_next.SATUR = true;
	}
};
void GasOil_Elliptic::buildGridLog()
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

	double cm_mu = mu_w;
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
		cm_mu = mu_w;

		hz = hz1 = hz2 = 0.0;
		cm_z = sk_it->h1;
		cm_nu = (double)k * 2.0 * M_PI / (double)cellsNum_nu;
		z_prev1 = (sk_well->h_well - sk_well->h1) * pow((sk_well->h_well - sk_well->h1) / r_w, -2.0 / (double)(sk_well->cellsNum_z - 1));
		z_prev2 = r_w;

		// Left border
		cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z, mu_w, hnu, 0.0));
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

			cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z, mu_w, hnu, hz));
			cells_z++;

			if (cells_z >= sk_it->cellsNum_z)
			{
				cells_z = 0;
				++sk_it;
			}
		}
		cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z + hz / 2.0, mu_w, hnu, 0.0));

		// Middle cells
		for (int j = 0; j < cellsNum_mu; j++)
		{
			sk_it = props_sk.begin();	cells_z = 0;	hz = 0.0;
			cm_z = sk_it->h1;
			cm_mu = r_prev * (exp(logStep) + 1.0) / 2.0;
			hmu = r_prev * (exp(logStep) - 1.0);
			z_prev1 = (sk_well->h_well - sk_well->h1) * pow((sk_well->h_well - sk_well->h1) / r_w, -2.0 / (double)(sk_well->cellsNum_z - 1));
			z_prev2 = r_w;

			cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z, hmu, hnu, 0.0));
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

				cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z, hmu, hnu, hz));
				Volume += cells[cells.size() - 1].V;
				cells_z++;

				if (cells_z >= sk_it->cellsNum_z)
				{
					cells_z = 0;
					++sk_it;
				}
			}
			cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z + hz / 2.0, hmu, hnu, 0.0));

			r_prev *= exp(logStep);
		}

		// Right border
		sk_it = props_sk.begin();	cells_z = 0;	hz = 0.0;
		cm_z = sk_it->h1;
		cm_mu = mu_e;
		z_prev1 = (sk_well->h_well - sk_well->h1) * pow((sk_well->h_well - sk_well->h1) / r_w, -2.0 / (double)(sk_well->cellsNum_z - 1));
		z_prev2 = r_w;

		cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z, 0.0, hnu, 0.0));
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

			cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z, 0.0, hnu, hz));
			cells_z++;

			if (cells_z >= sk_it->cellsNum_z)
			{
				cells_z = 0;
				++sk_it;
			}
		}
		cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z + hz / 2.0, 0.0, hnu, 0.0));
	}

	setUnused();
	buildWellCells();
}
void GasOil_Elliptic::setUnused()
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
void GasOil_Elliptic::buildWellCells()
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
													0.0, cell.hnu, cell.hz));
				wellNebrMap[counter++] = idx + cellsNum_z + 2;
			}

			// top
			for (int j = 0; j < cellsNum_nu; j++)
			{
				const int idx = sk.start_z + (sk.cellsNum_z - 1) / 2 +
					(cellsNum_mu + 2) * (cellsNum_z + 2) * j;
				const Cell& cell = cells[idx];

				wellCells.push_back(Cell(counter, cell.mu, cell.nu, cell.z - cell.hz,
													cell.hmu, cell.hnu, 0.0));
				wellNebrMap[counter++] = idx - 1;
			}

			// bot 
			for (int j = 0; j < cellsNum_nu; j++)
			{
				const int idx = sk.start_z + (sk.cellsNum_z - 1) / 2 +
					(cellsNum_mu + 2) * (cellsNum_z + 2) * j;
				const Cell& cell = cells[idx];

				wellCells.push_back(Cell(counter, cell.mu, cell.nu, cell.z + cell.hz,
					cell.hmu, cell.hnu, 0.0));
				wellNebrMap[counter++] = idx + 1;
			}
		}
	}
}
void GasOil_Elliptic::setPerforated()
{
	height_perf = 0.0;

	// lateral
	for (int i = 0; i < cellsNum_nu; i++)
	{
		Qcell[i] = 0.0;
		const Cell& cell = wellCells[i];
		height_perf += Cell::a * sqrt(sinh(cell.mu) * sinh(cell.mu) + sin(cell.nu) * sin(cell.nu)) * cell.hnu * cell.hz;
	}

	// top
	for (int i = cellsNum_nu; i < 2 * cellsNum_nu; i++)
	{
		Qcell[i] = 0.0;
		const Cell& cell = wellCells[i];
		height_perf += Cell::a * Cell::a * (sinh(cell.mu) * sinh(cell.mu) + sin(cell.nu) * sin(cell.nu)) * cell.hnu * cell.hmu;
	}

	// bot
	for (int i = 2 * cellsNum_nu; i < 3 * cellsNum_nu; i++)
	{
		Qcell[i] = 0.0;
		const Cell& cell = wellCells[i];
		height_perf += Cell::a * Cell::a * (sinh(cell.mu) * sinh(cell.mu) + sin(cell.nu) * sin(cell.nu)) * cell.hnu * cell.hmu;
	}
}
void GasOil_Elliptic::setPeriod(int period)
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
		props_sk[i].radius_eff_mu = props_sk[i].radiuses_eff[period];
		props_sk[i].perm_eff_mu = props_sk[i].perms_eff[period];
		props_sk[i].radius_eff_z = props_sk[i].radiuses_eff[period];
		props_sk[i].perm_eff_z = props_sk[i].perms_eff[period];
		props_sk[i].skin = props_sk[i].skins[period];
	}
};

