#include <new>
#include <vector>
#include "gtest/gtest.h"

#include "tests/oil1D-test.h"
#include "util/utils.h"

using namespace oil1D_NIT;
using std::make_pair;

Properties* Oil1D_Test::getProps()
{
	Properties* props = new oil1D_NIT::Properties();

	props->cellsNum_r = 100;
	props->cellsNum_z = 1;

	props->timePeriods.push_back(279216.0);
	props->timePeriods.push_back(1800000.0);
	
	props->leftBoundIsRate = true;
	props->rightBoundIsPres = true;

	props->rates.push_back(45.0);
	props->rates.push_back(0.0);

	props->ht = 100.0;
	props->ht_min = 100.0;
	props->ht_max  = 100000.0;
	
	props->alpha = 7200.0;

	props->perfIntervals.push_back( make_pair(0, 0) );

	props->r_w = 0.1524;
	props->r_e = 570.0;

	Skeleton_Props tmp;
	tmp.cellsNum_z = 1;
	tmp.m = 0.186;
	tmp.p_init = tmp.p_out = 250.807 * 100000.0;
	tmp.t_init = 302.058;
	tmp.h1 = 1500.0;
	tmp.h2 = 1514.0;
	tmp.height = 14.0;
	tmp.perm_r = 460.0;
	tmp.perm_z = 0.0;
	tmp.dens_stc = 2200.0;
	tmp.beta = 6.0 * 1.0e-10;
	
	tmp.skins.push_back(2.7);
	tmp.skins.push_back(2.7);

	tmp.radiuses_eff.push_back(0.5);
	tmp.radiuses_eff.push_back(0.5);

	tmp.c = 1800.0;
	tmp.kappa_eff = 0.0;
	tmp.lambda = 5.0;
	props->props_sk.push_back( tmp );

	props->depth_point = 1500.0;

	// Thermal properties
	props->props_oil.visc = 1.64;
	props->props_oil.b_bore = 1.245;
	props->props_oil.dens_stc = 736.0;
	props->props_oil.beta = 5.0 * 1.e-9;
	props->props_oil.jt = 4.0 * 1.e-7;
	props->props_oil.ad = 2.1 * 1.e-7;
	props->props_oil.c = 1880.0;
	props->props_oil.lambda = 0.16;

	return props;
}

void Oil1D_Test::test()
{
}

void Oil1D_Test::stat_test()
{
}

void Oil1D_Test::non_stat_test()
{
}

void Oil1D_Test::jt_temp_test()
{
}