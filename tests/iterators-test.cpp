#include <new>
#include <vector>
#include "gtest/gtest.h"

#include "tests/iterators-test.h"
#include "util/utils.h"

using namespace gasOil_3d;
using std::make_pair;

Properties* Iterators_Test::getProps()
{
	gasOil_3d::Properties* props = new gasOil_3d::Properties();

	props->cellsNum_r = 30;
	props->cellsNum_phi = 1;
	props->cellsNum_z = 30;

	props->timePeriods.push_back(10.0 * 365.0 * 86400.0);

	props->leftBoundIsRate = false;
	props->rightBoundIsPres = true;
	//props->rates.push_back(1.0);
	props->pwf.push_back(200.0 * 1.E+5);

	props->ht = 1000.0;
	props->ht_min = 100.0;
	props->ht_max = 5000000.0;

	props->alpha = 7200.0;

	props->r_w = 0.1;
	props->r_e = 500.0;

	props->perfIntervals.push_back(make_pair(1, 1));

	gasOil_3d::Skeleton_Props tmp;
	tmp.cellsNum_z = 30;
	tmp.m = 0.1;
	tmp.p_init = tmp.p_out = tmp.p_bub = 250.0 * 1.0e+5;
	tmp.s_init = 0.999;
	tmp.h1 = 1500.0;
	tmp.h2 = 1500.25;
	tmp.height = 1.0;
	tmp.perm_r = 100.0;
	tmp.perm_z = 5.0;
	tmp.dens_stc = 2000.0;
	tmp.beta = 4.35113e-10;
	tmp.skins.push_back(0.0);
	tmp.radiuses_eff.push_back(props->r_w);
	props->props_sk.push_back(tmp);

	props->depth_point = 1500.0;

	props->props_oil.visc = 1.64;
	props->props_oil.b_bore = 1.245;
	props->props_oil.dens_stc = 736.0;
	props->props_oil.beta = 1.0 * 1.e-9;

	props->props_gas.visc = 0.02833;
	props->props_gas.dens_stc = 0.72275;

	props->props_gas.visc = 0.02833;
	props->props_gas.dens_stc = 0.8;

	// Defining relative permeabilities
	setDataFromFile(props->kr_oil, "props/koil.txt");
	setDataFromFile(props->kr_gas, "props/kgas.txt");

	// Defining volume factors
	//props->byDefault.B_oil = true;
	setDataFromFile(props->B_oil, "props/Boil.txt");
	//props->byDefault.B_gas = false;
	setDataFromFile(props->B_gas, "props/Bgas.txt");

	//props->byDefault.Rs = true;
	setDataFromFile(props->Rs, "props/Rs.txt");

	return props;
}

void Iterators_Test::run()
{
	props = getProps();
	scene.load(*props);
}

void Iterators_Test::test()
{
	GasOil_3D* model = scene.getModel();
	gasOil_3d::Iterator itr;

	for (itr = model->getLeftBegin(); itr != model->getLeftEnd(); ++itr)
		ASSERT_EQ(itr->num, itr.getIdx());

	for (itr = model->getMidBegin(); itr != model->getMidEnd(); ++itr)
		ASSERT_EQ(itr->num, itr.getIdx());

	for (itr = model->getRightBegin(); itr != model->getRightEnd(); ++itr)
		ASSERT_EQ(itr->num, itr.getIdx());
}
