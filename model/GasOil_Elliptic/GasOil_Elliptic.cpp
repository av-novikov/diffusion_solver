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
	for (int j = 0; j < skeletonsNum; j++)
	{
		props_sk[j].perm_mu = MilliDarcyToM2(props_sk[j].perm_mu);
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
};
/*void GasOil_Elliptic::buildGridLog()
{
	cells.reserve(cellsNum);

	Volume = 0.0;
	int counter = 0;
	int skel_idx = 0, cells_z = 0;

	const double mu_w = asinh(4.0 * r_w / l);
	const double mu_e = asinh(2.0 * r_e / l);
	const double hmu = (mu_e - mu_w) / (double)cellsNum_mu;
	const double hnu = 2.0 * M_PI / (double)cellsNum_nu;
	double hz = props_sk[skel_idx].height;

	double cm_mu = mu_w;
	double cm_nu = 0.0;
	double cm_z = props_sk[skel_idx].h1;

	counter = 0;
	
	// Top
	for (int j = 0; j < cellsNum_nu; j++)
	{
		cm_mu = mu_w;

		cells.push_back( Cell(counter++, cm_mu, cm_nu, cm_z, 0.0, hnu, 0.0) );
		for (int i = 0; i < cellsNum_mu; i++)
		{
			cells.push_back( Cell(counter++, cm_mu + hmu / 2.0, cm_nu, cm_z, hmu, hnu, 0.0) );
			cm_mu += hmu;
		}
		cells.push_back( Cell(counter++, cm_mu, cm_nu, cm_z, 0.0, hnu, 0.0) );

		cm_nu += hnu;
	}

	// Middle
	for (int k = 0; k < cellsNum_z; k++)
	{
		hz = (props_sk[skel_idx].h2 - props_sk[skel_idx].h1) / (double)(props_sk[skel_idx].cellsNum_z);
		cm_z += (cells[cells.size() - 1].hz + hz) / 2.0;

		cm_nu = 0.0;
		for (int j = 0; j < cellsNum_nu; j++)
		{
			cm_mu = mu_w;

			cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z, 0.0, hnu, hz));
			for (int i = 0; i < cellsNum_mu; i++)
			{
				cells.push_back( Cell(counter++, cm_mu + hmu / 2.0, cm_nu, cm_z, hmu, hnu, hz) );
				cm_mu += hmu;
			}
			cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z, 0.0, hnu, hz));

			cm_nu += hnu;
		}

		cells_z++;
		if (cells_z >= props_sk[skel_idx].cellsNum_z)
		{
			cells_z = 0;
			skel_idx++;
		}
	}

	// Bottom
	for (int j = 0; j < cellsNum_nu; j++)
	{
		cm_mu = mu_w;

		cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z, 0.0, hnu, 0.0));
		for (int i = 0; i < cellsNum_mu; i++)
		{
			cells.push_back(Cell(counter++, cm_mu + hmu / 2.0, cm_nu, cm_z, hmu, hnu, 0.0));
			cm_mu += hmu;
		}
		cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z, 0.0, hnu, 0.0));

		cm_nu += hnu;
	}
};*/

void GasOil_Elliptic::buildGridLog()
{
	cells.reserve(cellsNum);

	Volume = 0.0;
	int counter = 0;
	int skel_idx = 0, cells_z = 0;

	const double mu_w = asinh(4.0 * r_w / l);
	const double mu_e = asinh(2.0 * r_e / l);
	double hmu = (mu_e - mu_w) / (double)cellsNum_mu;
	const double hnu = 2.0 * M_PI / (double)cellsNum_nu;
	double hz = props_sk[skel_idx].height;

	double r_prev = mu_w;
	double logMax = log(mu_e / mu_w);
	double logStep = logMax / (double)cellsNum_mu;

	double cm_mu = mu_w;
	double cm_nu = 0.0;
	double cm_z = props_sk[skel_idx].h1;

	counter = 0;
	for (int k = 0; k < cellsNum_nu; k++)
	{
		skel_idx = 0;	cells_z = 0;
		
		r_prev = mu_w;
		logMax = log(mu_e / mu_w);
		logStep = logMax / (double)cellsNum_mu;
		hmu = r_prev * (exp(logStep) - 1.0);
		cm_mu = mu_w;

		hz = 0.0;
		cm_z = props_sk[0].h1;
		cm_nu = (double)k * 2.0 * M_PI / (double)cellsNum_nu;

		// Left border
		cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z, 0.0, hnu, 0.0));
		for (int i = 0; i < cellsNum_z; i++)
		{
			hz = (props_sk[skel_idx].h2 - props_sk[skel_idx].h1) / (double)(props_sk[skel_idx].cellsNum_z);
			cm_z += (cells[cells.size() - 1].hz + hz) / 2.0;

			cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z, 0.0, hnu, hz));
			cells_z++;

			if (cells_z >= props_sk[skel_idx].cellsNum_z)
			{
				cells_z = 0;
				skel_idx++;
			}
		}
		cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z + hz / 2.0, 0.0, hnu, 0.0));

		// Middle cells
		for (int j = 0; j < cellsNum_mu; j++)
		{
			skel_idx = 0;	cells_z = 0;
			cm_z = props_sk[0].h1;
			cm_mu = r_prev * (exp(logStep) + 1.0) / 2.0;
			hmu = r_prev * (exp(logStep) - 1.0);

			cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z, hmu, hnu, 0.0));
			for (int i = 0; i < cellsNum_z; i++)
			{
				hz = (props_sk[skel_idx].h2 - props_sk[skel_idx].h1) / (double)(props_sk[skel_idx].cellsNum_z);
				cm_z += (cells[cells.size() - 1].hz + hz) / 2.0;

				cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z, hmu, hnu, hz));
				Volume += cells[cells.size() - 1].V;
				cells_z++;

				if (cells_z >= props_sk[skel_idx].cellsNum_z)
				{
					cells_z = 0;
					skel_idx++;
				}
			}
			cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z + hz / 2.0, hmu, hnu, 0.0));

			r_prev *= exp(logStep);
		}

		// Right border
		cm_z = props_sk[0].h1;
		cm_mu = mu_e;

		cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z, 0.0, hnu, 0.0));
		skel_idx = 0;	cells_z = 0;
		for (int i = 0; i < cellsNum_z; i++)
		{
			hz = (props_sk[skel_idx].h2 - props_sk[skel_idx].h1) / (double)(props_sk[skel_idx].cellsNum_z);
			cm_z += (cells[cells.size() - 1].hz + hz) / 2.0;

			cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z, 0.0, hnu, hz));
			cells_z++;

			if (cells_z >= props_sk[skel_idx].cellsNum_z)
			{
				cells_z = 0;
				skel_idx++;
			}
		}
		cells.push_back(Cell(counter++, cm_mu, cm_nu, cm_z + hz / 2.0, 0.0, hnu, 0.0));
	}
}

void GasOil_Elliptic::setPerforated()
{}
void GasOil_Elliptic::setPeriod(int period) {}

