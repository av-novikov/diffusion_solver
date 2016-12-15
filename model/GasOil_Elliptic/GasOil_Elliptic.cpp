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

	skeletonsNum = props.props_sk.size();
	props_sk = props.props_sk;
	for (int j = 0; j < skeletonsNum; j++)
	{
		props_sk[j].perm_mu = MilliDarcyToM2(props_sk[j].perm_mu);
		props_sk[j].perm_z = MilliDarcyToM2(props_sk[j].perm_z);
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
void GasOil_Elliptic::buildGridLog()
{
	cells.reserve(cellsNum);

	Volume = 0.0;
	int counter = 0;
	int skel_idx = 0, cells_z = 0;

	const double mu_w = asinh(4.0 * r_w / l);
	const double mu_e = asinh(2.0 * r_e / l);
	const double hmu = (mu_e - mu_w) / (double)cellsNum_mu;
	const double hnu = 2.0 * M_PI / (double)cellsNum_nu;
	double hz = 0.0;

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
		//cm_z = 
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

		cm_z += hz;
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
};
void GasOil_Elliptic::setPerforated()
{}
void GasOil_Elliptic::setPeriod(int period) {}

