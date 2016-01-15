#include <new>
#include <vector>
#include "gtest/gtest.h"

#include "tests/gas1Dsimple-test.h"
#include "util/utils.h"

using namespace gas1D;
using std::make_pair;

Properties* Gas1D_Simple_Test::getProps()
{
	Properties* props = new Properties();

	props->cellsNum_r = 100;

	props->timePeriods.push_back(10.0 * 365.0 * 86400.0);
	
	props->leftBoundIsRate = false;
	props->rightBoundIsPres = true;
	//props->rates.push_back(5000.0);
	props->pwf.push_back(100.0 * 100000.0);

	props->ht = 20000.0;
	props->ht_min = 10000.0;
	props->ht_max  = 10000000.0;
	
	props->alpha = 7200.0;

	props->perfIntervals.push_back( make_pair(0, 0) );

	props->r_w = 0.1;
	props->r_e = 1000.0;

	Skeleton_Props tmp;
	tmp.cellsNum_z = 1;
	tmp.m = 0.2;
	tmp.p_init = tmp.p_out = 150.0 * 1.0e+5;
	tmp.h1 = 1500.0;
	tmp.h2 = 1510.0;
	tmp.height = 10.0;
	tmp.perm_r = 10.0 / 0.986923;
	tmp.perm_z = 0.0;
	tmp.dens_stc = 2000.0;
	tmp.beta = 0.0;//3.0e-10;
	tmp.skins.push_back(0.0);
	tmp.radiuses_eff.push_back(props->r_w);
	props->props_sk.push_back( tmp );

	props->depth_point = 1500.0;

	props->props_gas.visc = 0.01;
	setDataFromFile(props->z_factor, "props/z.txt");

	return props;
}

double Gas1D_Simple_Test::getStatRate() const
{
	const Skeleton_Props& pr_sk = props->props_sk[0];
	double p_w = *(--props->pwf.end());

	return M_PI * pr_sk.perm_r * 0.986923 * 1.E-15 * pr_sk.height * (pr_sk.p_out * pr_sk.p_out - p_w * p_w) / 
			cPToPaSec(props->props_gas.visc) / log(props->r_e / props->r_w) / 1.E+5;
}

void Gas1D_Simple_Test::test()
{
	Gas1D_simple* model = scene.getModel();
	ASSERT_NEAR( model->getRate() * model->Q_dim * 86400.0, getStatRate() * 86400.0, getStatRate() * RATE_REL_TOL * 86400.0);
}