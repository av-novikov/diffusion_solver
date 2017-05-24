#include "model/BlackOilNIT_Elliptic/BlackOilNIT_Elliptic.hpp"
#include "util/utils.h"

#include "adolc/drivers/drivers.h"
#include "adolc/adolc.h"

using namespace std;
using namespace blackoilnit_elliptic;

BlackOilNIT_Elliptic::BlackOilNIT_Elliptic()
{
	x = new double[stencil * var_size];
	y = new double[var_size - 1];

	jac = new double*[var_size - 1];
	for (int i = 0; i < var_size - 1; i++)
		jac[i] = new double[stencil * var_size];
};
BlackOilNIT_Elliptic::~BlackOilNIT_Elliptic()
{
	delete x;
	delete y;

	for (int i = 0; i < var_size - 1; i++)
		delete[] jac[i];
	delete[] jac;
};
void BlackOilNIT_Elliptic::setProps(Properties& props)
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

	perfIntervals = props.perfIntervals;

	// Setting skeleton properties
	depth_point = props.depth_point;

	skeletonsNum = props.props_sk.size();
	props_sk = props.props_sk;
	int idx_z = 1;
	for (auto sk = props_sk.begin(); sk != props_sk.end(); ++sk)
	{
		// TODO: some changes
		sk->start_z = 1;	 idx_z += sk->cellsNum_z;
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
				props_sk[j].perms_eff.push_back(getDamagedPerm(props_sk[j], i));
				//props_sk[j].perms_eff.push_back(MilliDarcyToM2(props.props_sk[0].perm_mu * log(props.props_sk[j].radiuses_eff[i] / props.r_w) / (log(props.props_sk[j].radiuses_eff[i] / props.r_w) + props.props_sk[j].skins[i])));
			else
				props_sk[j].perms_eff.push_back(MilliDarcyToM2(props.props_sk[j].perm_mu));
		}
	}

	// Temporal properties
	ht = props.ht;
	ht_min = props.ht_min;
	ht_max = props.ht_max;

	// Water properties
	props_wat = props.props_wat;
	props_wat.visc = cPToPaSec(props_wat.visc);
	// Oil properties
	props_oil = props.props_oil;
	props_oil.visc = cPToPaSec(props_oil.visc);
	// Gas properties
	props_gas = props.props_gas;
	props_gas.visc = cPToPaSec(props_gas.visc);

	L = props.L;

	makeDimLess();

	Cell::a = l / 2;

	props_oil.b = setDataset(props.B_oil, P_dim / BAR_TO_PA, 1.0);
	props_gas.b = setDataset(props.B_gas, P_dim / BAR_TO_PA, 1.0);
	props_oil.Rs = setDataset(props.Rs, P_dim / BAR_TO_PA, 1.0);
};
void BlackOilNIT_Elliptic::makeDimLess()
{
	// Main units
	R_dim = r_e / 20.0;
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
	l /= R_dim;

	// Skeleton properties
	for (int i = 0; i < skeletonsNum; i++)
	{
		props_sk[i].perm_mu /= (R_dim * R_dim);
		props_sk[i].perm_z /= (R_dim * R_dim);

		props_sk[i].beta /= (1.0 / P_dim);
		props_sk[i].p_ref /= P_dim;
		props_sk[i].dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
		props_sk[i].h1 = (props_sk[i].h1 - depth_point) / R_dim;
		props_sk[i].h2 = (props_sk[i].h2 - depth_point) / R_dim;
		props_sk[i].h_well = (props_sk[i].h_well - depth_point) / R_dim;
		props_sk[i].height /= R_dim;
		props_sk[i].p_init /= P_dim;
		props_sk[i].p_out /= P_dim;
		props_sk[i].p_sat /= P_dim;
		props_sk[i].t_init /= T_dim;

		props_sk[i].c = props_sk[i].c / R_dim / R_dim * T_dim * t_dim * t_dim;
		props_sk[i].lambda_r = props_sk[i].lambda_r * T_dim * t_dim / P_dim / R_dim / R_dim;
		props_sk[i].lambda_z = props_sk[i].lambda_z * T_dim * t_dim / P_dim / R_dim / R_dim;

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
	props_oil.c = props_oil.c / R_dim / R_dim * T_dim * t_dim * t_dim;
	props_oil.lambda = props_oil.lambda * T_dim * t_dim / P_dim / R_dim / R_dim;
	props_oil.jt = props_oil.jt * P_dim / T_dim;
	props_oil.ad = props_oil.ad * P_dim / T_dim;
	// Water
	props_wat.visc /= (P_dim * t_dim);
	props_wat.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_wat.beta /= (1.0 / P_dim);
	props_wat.p_ref /= P_dim;
	props_wat.c = props_wat.c / R_dim / R_dim * T_dim * t_dim * t_dim;
	props_wat.lambda = props_wat.lambda * T_dim * t_dim / P_dim / R_dim / R_dim;
	props_wat.jt = props_wat.jt * P_dim / T_dim;
	props_wat.ad = props_wat.ad * P_dim / T_dim;
	// Gas
	props_gas.visc /= (P_dim * t_dim);
	props_gas.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_gas.c = props_gas.c / R_dim / R_dim * T_dim * t_dim * t_dim;
	props_gas.lambda = props_gas.lambda * T_dim * t_dim / P_dim / R_dim / R_dim;
	props_gas.jt = props_gas.jt * P_dim / T_dim;
	props_gas.ad = props_gas.ad * P_dim / T_dim;

	L = L / R_dim / R_dim * t_dim * t_dim;
};
void BlackOilNIT_Elliptic::setInitialState()
{
	vector<Cell>::iterator it;
	for (it = cells.begin(); it != cells.end(); ++it)
	{
		const Skeleton_Props& props = props_sk[getSkeletonIdx(*it)];
		it->props = const_cast<Skeleton_Props*>(&props);
		it->u_prev.p = it->u_iter.p = it->u_next.p = props.p_init;
		it->u_prev.p_bub = it->u_iter.p_bub = it->u_next.p_bub = props.p_sat;
		it->u_prev.s_w = it->u_iter.s_w = it->u_next.s_w = props.sw_init;
		it->u_prev.s_o = it->u_iter.s_o = it->u_next.s_o = props.so_init;
		it->u_prev.t = it->u_iter.t = it->u_next.t = props.t_init;

		if (props.p_init > props.p_sat)
			it->u_prev.SATUR = it->u_iter.SATUR = it->u_next.SATUR = false;
		else
			it->u_prev.SATUR = it->u_iter.SATUR = it->u_next.SATUR = true;
	}

	for (it = wellCells.begin(); it != wellCells.end(); ++it)
	{
		const Skeleton_Props& props = props_sk[getSkeletonIdx(*it)];
		it->props = const_cast<Skeleton_Props*>(&props);
		it->u_prev.p = it->u_iter.p = it->u_next.p = props.p_init;
		it->u_prev.p_bub = it->u_iter.p_bub = it->u_next.p_bub = props.p_sat;
		it->u_prev.s_w = it->u_iter.s_w = it->u_next.s_w = props.sw_init;
		it->u_prev.s_o = it->u_iter.s_o = it->u_next.s_o = props.so_init;
		it->u_prev.t = it->u_iter.t = it->u_next.t = props.t_init;

		if (props.p_init > props.p_sat)
			it->u_prev.SATUR = it->u_iter.SATUR = it->u_next.SATUR = false;
		else
			it->u_prev.SATUR = it->u_iter.SATUR = it->u_next.SATUR = true;
	}
};
void BlackOilNIT_Elliptic::buildGridLog()
{
	cells.reserve(cellsNum);

	Volume = 0.0;
	int counter = 0;
	auto sk_it = props_sk.begin();	int cells_z = 0;

	const double mu_w = asinh(2.0 * r_w / l);
	mu_init = mu_w / 2.0;
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
		cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z, mu_w, hnu, 0.0, Type::TOP));
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
				cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z, mu_w, hnu, hz, Type::MIDDLE));
			else
				cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z, mu_w, hnu, hz, Type::MIDDLE_SIDE));
			cells_z++;

			if (cells_z >= sk_it->cellsNum_z)
			{
				cells_z = 0;
				++sk_it;
			}
		}
		cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z + hz / 2.0, mu_w, hnu, 0.0, Type::BOTTOM));

		// Middle cells
		for (int j = 0; j < cellsNum_mu; j++)
		{
			sk_it = props_sk.begin();	cells_z = 0;	hz = 0.0;
			cm_z = sk_it->h1;
			cm_mu = r_prev * (exp(logStep) + 1.0) / 2.0;
			hmu = r_prev * (exp(logStep) - 1.0);
			z_prev1 = (sk_well->h_well - sk_well->h1) * pow((sk_well->h_well - sk_well->h1) / r_w, -2.0 / (double)(sk_well->cellsNum_z - 1));
			z_prev2 = r_w;

			cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z, hmu, hnu, 0.0, Type::TOP));
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

				cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z, hmu, hnu, hz, Type::MIDDLE));
				Volume += cells[cells.size() - 1].V;
				cells_z++;

				if (cells_z >= sk_it->cellsNum_z)
				{
					cells_z = 0;
					++sk_it;
				}
			}
			cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z + hz / 2.0, hmu, hnu, 0.0, Type::BOTTOM));

			r_prev *= exp(logStep);
		}

		// Right border
		sk_it = props_sk.begin();	cells_z = 0;	hz = 0.0;
		cm_z = sk_it->h1;
		cm_mu = mu_e;
		z_prev1 = (sk_well->h_well - sk_well->h1) * pow((sk_well->h_well - sk_well->h1) / r_w, -2.0 / (double)(sk_well->cellsNum_z - 1));
		z_prev2 = r_w;

		cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z, 0.0, hnu, 0.0, Type::RIGHT));
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

			cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z, 0.0, hnu, hz, Type::RIGHT));
			cells_z++;

			if (cells_z >= sk_it->cellsNum_z)
			{
				cells_z = 0;
				++sk_it;
			}
		}
		cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z + hz / 2.0, 0.0, hnu, 0.0, Type::RIGHT));
	}

	setUnused();
	buildWellCells();
}
void BlackOilNIT_Elliptic::setUnused()
{
	/*for (const auto& sk : props_sk)
	{
		if (sk.isWellHere)
		{*/
			const auto& sk = *props_sk.begin();
			for (int j = 0; j < cellsNum_nu; j++)
			{
				const int idx = sk.start_z + (sk.cellsNum_z - 1) / 2 +
					(cellsNum_mu + 2) * (cellsNum_z + 2) * j;
				cells[idx].isUsed = false;
			}
	/*	}
	}*/
}
void BlackOilNIT_Elliptic::buildWellCells()
{
	int counter = 0;
	/*for (const auto& sk : props_sk)
	{
		if (sk.isWellHere)
		{*/
			const auto sk = *props_sk.begin();
			// lateral
			for (int j = 0; j < cellsNum_nu; j++)
			{
				const int idx = sk.start_z + (sk.cellsNum_z - 1) / 2 +
					(cellsNum_mu + 2) * (cellsNum_z + 2) * j;
				const Cell& cell = cells[idx];

				wellCells.push_back(Cell(counter, cell.mu + cell.hmu / 2.0, cell.nu, cell.z,
					0.0, cell.hnu, cell.hz, Type::WELL_LAT));
				wellNebrMap[idx + cellsNum_z + 2] = counter;
				nebrMap[counter++] = make_pair<int, int>(idx + cellsNum_z + 2, idx + 2 * cellsNum_z + 4);
			}

			// top
			for (int j = 0; j < cellsNum_nu; j++)
			{
				const int idx = sk.start_z + (sk.cellsNum_z - 1) / 2 +
					(cellsNum_mu + 2) * (cellsNum_z + 2) * j;
				const Cell& cell = cells[idx];

				wellCells.push_back(Cell(counter, cell.mu, cell.nu, cell.z - cell.hz / 2.0,
					cell.hmu, cell.hnu, 0.0, Type::WELL_TOP));
				wellNebrMap[idx - 1] = counter;
				nebrMap[counter++] = make_pair<int, int>(idx - 1, idx - 2);
			}

			// bot 
			for (int j = 0; j < cellsNum_nu; j++)
			{
				const int idx = sk.start_z + (sk.cellsNum_z - 1) / 2 +
					(cellsNum_mu + 2) * (cellsNum_z + 2) * j;
				const Cell& cell = cells[idx];

				wellCells.push_back(Cell(counter, cell.mu, cell.nu, cell.z + cell.hz / 2.0,
					cell.hmu, cell.hnu, 0.0, Type::WELL_BOT));
				wellNebrMap[idx + 1] = counter;
				nebrMap[counter++] = make_pair<int, int>(idx + 1, idx + 2);
			}
		/*}
	}*/
}
void BlackOilNIT_Elliptic::setPerforated()
{
	height_perf = 0.0;

	for (auto perfInt : perfIntervals)
	{
		for (int i = perfInt.first; i <= perfInt.second; i++)
		{
			const auto indices = getPerforationIndices(i);
			for (auto ind : indices)
			{
				Qcell[ind] = 0.0;
				const Cell& cell = wellCells[ind];
				if(cell.type == Type::WELL_LAT)
					height_perf += Cell::getH(cell.mu, cell.nu) * cell.hnu * cell.hz;
				else
					height_perf += Cell::getH(cell.mu, cell.nu) * Cell::getH(cell.mu, cell.nu) * cell.hnu * cell.hmu;
			}

			Qcell_ellipse[ indices[0] ] = 0.0;
			Qcell_ellipse[ indices[2] ] = 0.0;
		};
	};

	/*// lateral
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
	}*/
}
void BlackOilNIT_Elliptic::setPeriod(int period)	
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
				if (cell.type == Type::WELL_LAT)
					S = Cell::getH(cell.mu, cell.nu) * cell.hnu * cell.hz;
				else
					S = Cell::getH(cell.mu, cell.nu) * Cell::getH(cell.mu, cell.nu) * cell.hnu * cell.hmu;

				it->second = Q_sum * S / height_perf;
				it_ell = Qcell_ellipse.find(it->first);
				if (it_ell != Qcell_ellipse.end())
				{
					if(wellCells[it_ell->first].type == Type::WELL_LAT)
						it_ell->second = it->second;
					else
						it_ell->second = 2.0 * it->second;
					
					Q_sum_quater += it_ell->second;
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
					if (wellCells[it_ell->first].type == Type::WELL_LAT)
						it_ell->second = it->second;
					else
						it_ell->second = 2.0 * it->second;

					Q_sum_quater += it_ell->second;
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
		auto& props = props_sk[i];

		props.radius_eff_mu = asinh(props.radiuses_eff[period] / Cell::a);
		props.radius_eff_z = props.radiuses_eff[period];
		props.radius_eff = props.radiuses_eff[period];

		props.perm_eff_mu = props_sk[i].perms_eff[period];
		props.perm_eff_z = props.perm_z / props.perm_mu * props.perms_eff[period];

		props.skin = props.skins[period];
	}
};
void BlackOilNIT_Elliptic::setRateDeviation(int num, double ratio)
{
	const auto indices = getPerforationIndicesSameType(num);

	Qcell_ellipse[num] += Q_sum_quater * ratio;
	for (const auto& idx : indices)
		Qcell[idx] += Q_sum_quater * ratio / indices.size();
}
double BlackOilNIT_Elliptic::solveH()
{
	double H = 0.0;
	double p1, p0;

	map<int, double>::iterator it = Qcell.begin();
	for (int i = 0; i < Qcell.size() - 1; i++)
	{
		p0 = wellCells[it->first].u_next.p;
		p1 = wellCells[(++it)->first].u_next.p;

		H += (p1 - p0) * (p1 - p0) / 2.0;
	}

	return H;
}
double BlackOilNIT_Elliptic::getRate(int cur) const
{
	const Cell& cell = wellCells[cur];
	const Cell& beta = cells[nebrMap.at(cur).first];
	const Variable& next = cell.u_next;
	return getTrans(cell, beta) * props_oil.getKr(next.s_w, next.s_o, cell.props).value() / 
			props_oil.getB(next.p, next.p_bub, next.SATUR).value() / props_oil.getViscosity(next.p).value() * 
			(beta.u_next.p - next.p);
}
double BlackOilNIT_Elliptic::phaseTrans(const Cell& cell)
{
	const Skeleton_Props& props = *cell.props;
	const int nebrNum = (cell.type == Type::MIDDLE) ? 6 : 5;
	const Variable& next = cell.u_next;
	const Variable& prev = cell.u_prev;

	double H = (props.getPoro(next.p).value() * next.s_o * props_oil.getDensity(next.p, next.p_bub, next.SATUR).value() -
				props.getPoro(prev.p).value() * prev.s_o * props_oil.getDensity(prev.p, prev.p_bub, prev.SATUR).value()) / ht;

	Cell* neighbor[6];
	getNeighbors(cell, neighbor);
	for (int i = 0; i < nebrNum; i++)
	{
		const Cell& beta = *neighbor[i];
		const Cell& upwd_cell = getUpwindCell(cell, beta);
		const Variable& upwd = upwd_cell.u_next;
		const Variable& nebr = beta.u_next;

		H += 1.0 / cell.V * getTrans(cell, beta) * (next.p - nebr.p) *
			props_oil.getKr(upwd.s_w, upwd.s_o, upwd_cell.props).value() / props_oil.getViscosity(upwd.p).value() * 
			props_oil.getDensity(upwd.p, upwd.p_bub, upwd.SATUR).value();
	}

	return H;
}

void BlackOilNIT_Elliptic::solve_eqMiddle(const Cell& cell, const int val)
{
	const Skeleton_Props& props = *cell.props;
	const int nebrNum = (cell.type == Type::MIDDLE) ? 6 : 5;

	if (val == TEMP)
	{
		trace_on(tmid);

		adouble h[Variable::size - 1];
		TapeVariableNIT var[stencil];
		const Variable& prev = cell.u_prev;

		for (int i = 0; i < nebrNum + 1; i++)
			var[i].t <<= x[i];

		const TapeVariableNIT& next = var[0];

		h[0] = getCn(cell) * (next.t - prev.t) - getAd(cell) * (cell.u_next.p - prev.p) + ht * L * phaseTrans(cell);

		Cell* neighbor[6];
		getNeighbors(cell, neighbor);
		for (int i = 0; i < nebrNum; i++)
		{
			const Cell& beta = *neighbor[i];
			const TapeVariableNIT& nebr = var[i + 1];
			const auto mult = getDivCoeff(const_cast<Cell&>(cell), const_cast<Cell&>(beta), neighbor);

			h[0] += ht * (mult.ther * (next.t - nebr.t) + mult.pres * (cell.u_next.p - beta.u_next.p)) / fabs(getDistance(beta, cell));
		}

		adouble buf = h[0];
		adouble isMiddleSide = ((cell.nu == 0.0 || cell.nu == M_PI) && cell.mu < 10.0 * mu_init) ? true : false;
		condassign(h[0], isMiddleSide, buf / sqrt(cell.V));
		h[0] *= cell.V;
		h[0] >>= y[0];

		trace_off();
	}
	else if (val == PRES)
	{
		trace_on(mid);

		adouble h[Variable::size - 1], tmp;
		TapeVariable var[stencil];
		const Variable& prev = cell.u_prev;

		for (int i = 0; i < nebrNum + 1; i++)
		{
			var[i].p <<= x[i * var_size];
			var[i].s_w <<= x[i * var_size + 1];
			var[i].s_o <<= x[i * var_size + 2];
			var[i].p_bub <<= x[i * var_size + 3];
		}

		const TapeVariable& next = var[0];
		adouble satur = cell.u_next.SATUR;

		h[0] = props.getPoro(next.p) * next.s_w / props_wat.getB(next.p) -
			props.getPoro(prev.p) * prev.s_w / props_wat.getB(prev.p);
		h[1] = props.getPoro(next.p) * next.s_o / props_oil.getB(next.p, next.p_bub, next.SATUR) -
			props.getPoro(prev.p) * prev.s_o / props_oil.getB(prev.p, prev.p_bub, prev.SATUR);
		condassign(h[2], satur,
			props.getPoro(next.p) * ((1.0 - next.s_o - next.s_w) / props_gas.getB(next.p) +
				next.s_o * props_oil.getRs(next.p, next.p_bub, next.SATUR) / props_oil.getB(next.p, next.p_bub, next.SATUR)) -
			props.getPoro(prev.p) * ((1.0 - prev.s_o - prev.s_w) / props_gas.getB(prev.p) +
				prev.s_o * props_oil.getRs(prev.p, prev.p_bub, prev.SATUR) / props_oil.getB(prev.p, prev.p_bub, prev.SATUR)),
			props.getPoro(next.p) *	next.s_o * props_oil.getRs(next.p, next.p_bub, next.SATUR) / props_oil.getB(next.p, next.p_bub, next.SATUR) -
			props.getPoro(prev.p) *	prev.s_o * props_oil.getRs(prev.p, prev.p_bub, prev.SATUR) / props_oil.getB(prev.p, prev.p_bub, prev.SATUR));

		Cell* neighbor[6];
		getNeighbors(cell, neighbor);
		for (int i = 0; i < nebrNum; i++)
		{
			const Cell& beta = *neighbor[i];
			const Cell& upwd_cell = getUpwindCell(cell, beta);
			const int upwd_idx = (upwd_cell == cell) ? 0 : i + 1;
			TapeVariable& upwd = var[upwd_idx];
			const TapeVariable& nebr = var[i + 1];

			h[0] += ht / cell.V * getTrans(cell, beta) * (next.p - nebr.p) *
				props_wat.getKr(upwd.s_w, upwd.s_o, upwd_cell.props) / props_wat.getViscosity(upwd.p) /
				props_wat.getB(upwd.p);
			h[1] += ht / cell.V * getTrans(cell, beta) * (next.p - nebr.p) *
				props_oil.getKr(upwd.s_w, upwd.s_o, upwd_cell.props) / props_oil.getViscosity(upwd.p) /
				props_oil.getB(upwd.p, upwd.p_bub, upwd.SATUR);
			condassign(tmp, satur,
				ht / cell.V * getTrans(cell, beta) * (next.p - nebr.p) *
				(props_oil.getKr(upwd.s_w, upwd.s_o, upwd_cell.props) * props_oil.getRs(upwd.p, upwd.p_bub, upwd.SATUR) /
					props_oil.getViscosity(upwd.p) / props_oil.getB(upwd.p, upwd.p_bub, upwd.SATUR) +
					props_gas.getKr(upwd.s_w, upwd.s_o, upwd_cell.props) / props_gas.getViscosity(upwd.p) / props_gas.getB(upwd.p)),
				ht / cell.V * getTrans(cell, beta) * (next.p - nebr.p) *
				props_oil.getKr(upwd.s_w, upwd.s_o, upwd_cell.props) * props_oil.getRs(upwd.p, upwd.p_bub, upwd.SATUR) /
				props_oil.getViscosity(upwd.p) / props_oil.getB(upwd.p, upwd.p_bub, upwd.SATUR));
			h[2] += tmp;
		}

		h[0] *= sqrt(cell.V);		h[1] *= cell.V;		h[2] *= cell.V;
		adouble buf[3];
		buf[0] = h[0];		buf[1] = h[1];		buf[2] = h[2];
		adouble isMiddleSide = ((cell.nu == 0.0 || cell.nu == M_PI) && cell.mu < 2.0 * mu_init) ? true : false;
		condassign(h[0], isMiddleSide, buf[0] * 100.0);
		condassign(h[1], isMiddleSide, buf[1] * 100.0);
		condassign(h[2], isMiddleSide, buf[2] * 100.0);

		for (int i = 0; i < var_size-1; i++)
			h[i] >>= y[i];

		trace_off();
	}
}
void BlackOilNIT_Elliptic::solve_eqWell(const Cell& cell, const int val)
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
		dist1 = Cell::getH(cell.mu + sign(beta1.mu - cell.mu) * cell.hmu / 2.0, cell.nu) * (cell.mu - beta1.mu);
		dist2 = Cell::getH(beta1.mu + sign(beta2.mu - beta1.mu) * beta1.hmu / 2.0, beta1.nu) * (beta1.mu - beta2.mu);
	}

	if (val == TEMP)
	{
		trace_on(tleft);
		adouble h[Variable::size - 1];
		TapeVariableNIT var[TLstencil];
		for (int i = 0; i < TLstencil; i++)
			var[i].t <<= x[i];

		const TapeVariableNIT& next = var[0];
		const TapeVariableNIT& nebr1 = var[1];
		const TapeVariableNIT& nebr2 = var[2];

		adouble isPerforated;
		auto it = Qcell.find(cell.num);
		if (it != Qcell.end())
			isPerforated = true;
		else
			isPerforated = false;

		condassign(h[0], isPerforated, ((next.t - nebr1.t) / (adouble)(dist1) - (nebr1.t - nebr2.t) / (adouble)(dist2)) / P_dim, (next.t - nebr1.t) / P_dim);
		h[0] >>= y[0];

		trace_off();
	}
	else if (val == PRES)
	{
		trace_on(left);
		adouble h[Variable::size - 1];
		TapeVariable var[Lstencil];
		for (int i = 0; i < Lstencil; i++)
		{
			var[i].p <<= x[i * var_size];
			var[i].s_w <<= x[i * var_size + 1];
			var[i].s_o <<= x[i * var_size + 2];
			var[i].p_bub <<= x[i * var_size + 3];
		}

		const TapeVariable& next = var[0];
		const TapeVariable& nebr1 = var[1];
		const TapeVariable& nebr2 = var[1];

		double rate;
		adouble leftIsRate = leftBoundIsRate;
		adouble leftPresPerforated, leftPresNonPerforated;
		auto it = Qcell.find(cell.num);
		if (it != Qcell.end())
		{
			leftPresPerforated = !leftBoundIsRate;
			leftPresNonPerforated = false;
			rate = Qcell[cell.num];
		}
		else
		{
			rate = 0.0;
			leftPresPerforated = false;
			leftPresNonPerforated = !leftBoundIsRate;
		}

		condassign(h[0], leftIsRate,
			(adouble)(getTrans(cell, beta1)) * props_oil.getKr(next.s_w, next.s_o, cell.props) /
			props_oil.getViscosity(next.p) / props_oil.getB(next.p, next.p_bub, next.SATUR) * (nebr1.p - next.p) - rate,
			(next.p - Pwf) / P_dim);
		condassign(h[0], leftPresPerforated, (next.p - Pwf) / P_dim);
		condassign(h[0], leftPresNonPerforated, (next.p - nebr1.p) / P_dim);

		h[1] = ((next.s_w - nebr1.s_w) / (adouble)(dist1) - (nebr1.s_w - nebr2.s_w) / (adouble)(dist2)) / P_dim;

		adouble satur = cell.u_next.SATUR;
		condassign(h[2], satur,
			((next.s_o - nebr1.s_o) / (adouble)(dist1) - (nebr1.s_o - nebr2.s_o) / (adouble)(dist2)) / P_dim,
			((next.p_bub - nebr1.p_bub) / (adouble)(dist1) - (nebr1.p_bub - nebr2.p_bub) / (adouble)(dist2)) / P_dim);

		for (int i = 0; i < var_size-1; i++)
			h[i] >>= y[i];

		trace_off();
	}
}
void BlackOilNIT_Elliptic::solve_eqRight(const Cell& cell, const int val)
{
	if (val == TEMP)
	{
		trace_on(tright);
		adouble h[Variable::size - 1];
		TapeVariableNIT var[2];
		for (int i = 0; i < Rstencil; i++)
			var[i].t <<= x[i];

		const TapeVariableNIT& next = var[0];
		const TapeVariableNIT& nebr = var[1];

		h[0] = (next.t - (adouble)(cell.props->t_init)) / P_dim;
		h[0] >>= y[0];

		trace_off();
	}
	else if (val == PRES)
	{
		trace_on(right);
		adouble h[Variable::size - 1];
		TapeVariable var[2];
		adouble rightIsPres = rightBoundIsPres;
		for (int i = 0; i < Rstencil; i++)
		{
			var[i].p <<= x[i * var_size];
			var[i].s_w <<= x[i * var_size + 1];
			var[i].s_o <<= x[i * var_size + 2];
			var[i].p_bub <<= x[i * var_size + 3];
		}
		const TapeVariable& next = var[0];
		const TapeVariable& nebr = var[1];

		condassign(h[0], rightIsPres, (next.p - cell.props->p_out) / P_dim, (next.p - nebr.p) / P_dim);
		h[1] = (next.s_w - nebr.s_w) / P_dim;
		adouble satur = cell.u_next.SATUR;
		condassign(h[2], satur, (next.s_o - nebr.s_o) / P_dim, (next.p_bub - nebr.p_bub) / P_dim);

		for (int i = 0; i < var_size-1; i++)
			h[i] >>= y[i];

		trace_off();
	}
}
void BlackOilNIT_Elliptic::solve_eqVertical(const Cell& cell, const int val)
{
	if(val == TEMP)
	{
		trace_on(tvertical);
		adouble h[Variable::size - 1];
		TapeVariableNIT var[Vstencil];

		for (int i = 0; i < Vstencil; i++)
			var[i].t <<= x[i];

		const TapeVariableNIT& next = var[0];
		const TapeVariableNIT& nebr = var[1];

		h[0] = (next.t - nebr.t) / P_dim;
		h[0] >>= y[0];

		trace_off();
	}
	else if (val == PRES)
	{
		trace_on(vertical);
		adouble h[Variable::size - 1];
		TapeVariable var[Vstencil];

		for (int i = 0; i < Vstencil; i++)
		{
			var[i].p <<= x[i * var_size];
			var[i].s_w <<= x[i * var_size + 1];
			var[i].s_o <<= x[i * var_size + 2];
			var[i].p_bub <<= x[i * var_size + 3];
		}

		const TapeVariable& next = var[0];
		const TapeVariable& nebr = var[1];

		h[0] = (next.p - nebr.p) / P_dim;
		h[1] = (next.s_w - nebr.s_w) / P_dim;
		adouble satur = cell.u_next.SATUR;
		condassign(h[2], satur, (next.s_o - nebr.s_o) / P_dim, (next.p_bub - nebr.p_bub) / P_dim);

		for (int i = 0; i < var_size-1; i++)
			h[i] >>= y[i];

		trace_off();
	}
}
void BlackOilNIT_Elliptic::setVariables(const Cell& cell, const int val)
{
	assert(cell.isUsed);

	if (cell.type == Type::WELL_LAT || cell.type == Type::WELL_TOP || cell.type == Type::WELL_BOT) // Well
	{
		if (val == TEMP)
		{
			const Variable& next = cell.u_next;
			const Variable& nebr1 = cells[nebrMap[cell.num].first].u_next;
			const Variable& nebr2 = cells[nebrMap[cell.num].second].u_next;

			x[0] = next.values[val];
			x[1] = nebr1.values[val];
			x[2] = nebr2.values[val];

			solve_eqWell(cell, val);
			jacobian(tleft, 1, 1 * TLstencil, x, jac);
		}
		else if (val == PRES)
		{
			const Variable& next = cell.u_next;
			const Variable& nebr1 = cells[nebrMap[cell.num].first].u_next;
			const Variable& nebr2 = cells[nebrMap[cell.num].second].u_next;

			for (int i = 0; i < var_size; i++)
			{
				x[i] = next.values[i + 1];
				x[var_size + i] = nebr1.values[i + 1];
				x[2 * var_size + i] = nebr2.values[i + 1];
			}

			solve_eqWell(cell, val);
			jacobian(left, var_size - 1, var_size * Lstencil, x, jac);
		}
	}
	else if (cell.type == Type::RIGHT) // Right
	{
		const Variable& next = cell.u_next;
		const Variable& nebr1 = cells[cell.num - cellsNum_z - 2].u_next;
		const Variable& nebr2 = cells[cell.num - 2 * cellsNum_z - 4].u_next;

		if (val == TEMP) 
		{
			x[0] = next.values[val];
			x[1] = nebr1.values[val];
			solve_eqRight(cell, val);
			jacobian(tright, 1, 1 * Rstencil, x, jac);
		}
		else if (val == PRES)
		{
			for (int i = 0; i < var_size; i++)
			{
				x[i] = next.values[i+1];
				x[var_size + i] = nebr1.values[i+1];
			}
			solve_eqRight(cell, val);
			jacobian(right, var_size - 1, var_size * Rstencil, x, jac);
		}
	}
	else if (cell.type == Type::TOP) // Top
	{
		const Variable& next = cell.u_next;
		const Variable& nebr = cells[cell.num + 1].u_next;

		if (val == TEMP)
		{
			x[0] = next.values[val];
			x[1] = nebr.values[val];
			solve_eqVertical(cell, val);
			jacobian(tvertical, 1, 1 * Vstencil, x, jac);
		}
		else if (val == PRES)
		{
			for (int i = 0; i < var_size; i++)
			{
				x[i] = next.values[i + 1];
				x[var_size + i] = nebr.values[i + 1];
			}
			solve_eqVertical(cell, val);
			jacobian(vertical, var_size - 1, var_size * Vstencil, x, jac);
		}
	}
	else if (cell.type == Type::BOTTOM) // Bottom
	{
		const Variable& next = cell.u_next;
		const Variable& nebr = cells[cell.num - 1].u_next;

		if (val == TEMP)
		{
			x[0] = next.values[val];
			x[1] = nebr.values[val];
			solve_eqVertical(cell, val);
			jacobian(tvertical, 1, 1 * Vstencil, x, jac);
		}
		else if (val == PRES)
		{
			for (int i = 0; i < var_size; i++)
			{
				x[i] = next.values[i + 1];
				x[var_size + i] = nebr.values[i + 1];
			}
			solve_eqVertical(cell, val);
			jacobian(vertical, var_size - 1, var_size * Vstencil, x, jac);
		}
	}
	else // Middle
	{
		const int nebrNum = (cell.type == Type::MIDDLE) ? 6 : 5;
		const Variable& next = cell.u_next;
		Cell* neighbor[6];
		getNeighbors(cell, neighbor);

		if (val == TEMP)
		{
			x[0] = next.values[val];
			for (int j = 0; j < nebrNum; j++)
				x[(j + 1)] = neighbor[j]->u_next.values[val];

			solve_eqMiddle(cell, val);
			jacobian(tmid, 1, 1 * (nebrNum + 1), x, jac);
		}
		else if (val == PRES)
		{
			for (int i = 0; i < var_size; i++)
			{
				x[i] = next.values[i + 1];

				for (int j = 0; j < nebrNum; j++)
					x[(j + 1) * var_size + i] = neighbor[j]->u_next.values[i + 1];
			}

			solve_eqMiddle(cell, val);
			jacobian(mid, var_size - 1, var_size * (nebrNum + 1), x, jac);
		}
	}
}