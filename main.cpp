#include <new>
#include <vector>
#include <string>
#include <fstream>
#include <utility>
#include <iostream>
#include <valarray>

#include "util/utils.h"
#include "method/mcmath.h"
#include "Scene.h"

#include "model/GasOil_RZ/GasOil_RZ.h"
#include "model/VPP2d/VPP2d.hpp"

using namespace std;

gasOil_rz::Properties* getProps()
{
	gasOil_rz::Properties* props = new gasOil_rz::Properties();

	props->cellsNum_r = 50;
	props->cellsNum_z = 5;

	props->timePeriods.push_back(20.0 * 86400.0);
	props->timePeriods.push_back(40.0 * 86400.0);

	props->leftBoundIsRate = true;
	props->rightBoundIsPres = true;
	props->rates.push_back(10.0);
	props->rates.push_back(0.0);

	props->ht = 100.0;
	props->ht_min = 100.0;
	props->ht_max = 100000.0;

	props->alpha = 7200.0;

	props->r_w = 0.05;
	props->r_e = 1000.0;

	props->perfIntervals.push_back(make_pair(3, 3));

	gasOil_rz::Skeleton_Props tmp;
	tmp.cellsNum_z = 5;
	tmp.m = 0.1;
	tmp.p_init = tmp.p_out = tmp.p_bub = 200.0 * 1.0e+5;
	tmp.s_init = 1.0;
	tmp.h1 = 0.0;
	tmp.h2 = 10.0;
	tmp.height = 10.0;
	tmp.perm_r = 20.0;
	tmp.perm_z = 1.0;
	tmp.dens_stc = 2000.0;
	tmp.beta = 4.35113e-10;
	tmp.skins.push_back(0.0);
	tmp.skins.push_back(0.0);
	tmp.radiuses_eff.push_back(props->r_w);
	tmp.radiuses_eff.push_back(props->r_w);

	props->props_sk.push_back(tmp);

	props->depth_point = 0.0;

	props->props_oil.visc = 1.0;
	props->props_oil.b_bore = 1.0;
	props->props_oil.dens_stc = 736.0;
	props->props_oil.beta = 0.5 * 1.e-9;

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

/*vpp2d::Properties* getProps()
{
	vpp2d::Properties* props = new vpp2d::Properties();

	props->cellsNum_r = 100;
	props->cellsNum_z = 21;

	props->timePeriods.push_back(10.0 * 86400.0);
	props->timePeriods.push_back(20.0 * 86400.0);
	
	props->leftBoundIsRate = false;
	props->rightBoundIsPres = true;
	//props->rates.push_back(-100.0);
	//props->rates.push_back(0.0);
	props->c.push_back(0.1);
	props->c.push_back(0.0);
	props->pwf.push_back(250.0 * 1.E+5);
	props->pwf.push_back(200.0 * 1.E+5);

	props->ht = 100.0;
	props->ht_min = 100.0;
	props->ht_max = 100000.0;

	props->alpha = 7200.0;

	props->r_w = 0.05;
	props->r_e = 1000.0;

	props->perfIntervals.push_back(make_pair(6, 6));
	props->perfIntervals.push_back(make_pair(17, 17));

	vpp2d::Skeleton_Props tmp;
	tmp.cellsNum_z = 10;
	tmp.m = 0.1;
	tmp.p_init = tmp.p_out = 200.0 * 1.0e+5;
	tmp.s_init = 0.8;
	tmp.s_wc = 0.1;
	tmp.s_oc = 0.9;
	tmp.c_init = 0.0;
	tmp.h1 = 0.0;
	tmp.h2 = 10.0;
	tmp.height = 10.0;
	tmp.perm_r = 50.0;
	tmp.perm_z = 1.0;
	tmp.dens_stc = 2000.0;
	tmp.beta = 1.e-10;
	tmp.skins.push_back(0.0);
	tmp.skins.push_back(0.0);
	tmp.radiuses_eff.push_back(props->r_w);
	tmp.radiuses_eff.push_back(props->r_w);
	props->props_sk.push_back(tmp);

	tmp.cellsNum_z = 1;
	tmp.m = 0.1;
	tmp.p_init = tmp.p_out = 200.0 * 1.0e+5;
	tmp.s_init = 0.8;
	tmp.s_wc = 0.1;
	tmp.s_oc = 0.9;
	tmp.c_init = 0.0;
	tmp.h1 = 10.0;
	tmp.h2 = 20.0;
	tmp.height = 10.0;
	tmp.perm_r = 0.0;
	tmp.perm_z = 0.0;//	 50.0;
	tmp.dens_stc = 2000.0;
	tmp.beta = 1.e-10;
	tmp.skins.push_back(0.0);
	tmp.skins.push_back(0.0);
	tmp.radiuses_eff.push_back(props->r_w);
	tmp.radiuses_eff.push_back(props->r_w);
	props->props_sk.push_back(tmp);

	tmp.cellsNum_z = 10;
	tmp.m = 0.1;
	tmp.p_init = tmp.p_out = 200.0 * 1.0e+5;
	tmp.s_init = 0.8;
	tmp.s_wc = 0.1;
	tmp.s_oc = 0.9;
	tmp.c_init = 0.0;
	tmp.h1 = 20.0;
	tmp.h2 = 30.0;
	tmp.height = 10.0;
	tmp.perm_r = 500.0;
	tmp.perm_z = 10.0;//	 50.0;
	tmp.dens_stc = 2000.0;
	tmp.beta = 1.e-10;
	tmp.skins.push_back(0.0);
	tmp.skins.push_back(0.0);
	tmp.radiuses_eff.push_back(props->r_w);
	tmp.radiuses_eff.push_back(props->r_w);
	props->props_sk.push_back(tmp);

	props->depth_point = 0.0;

	props->props_o.visc = 1.0;
	props->props_o.dens_stc = 736.0;
	props->props_o.beta = 0.5 * 1.e-9;
	props->props_o.b_bore = 1.0;
	props->props_o.p_ref = 0.0;// tmp.p_out;

	props->props_w.visc = 1.0;
	props->props_w.dens_stc = 1000.0;
	props->props_w.beta = 0.5 * 1.e-9;
	props->props_w.b_bore = 1.0;
	props->props_w.p_ref = 0.0;// tmp.p_out;

	// Defining relative permeabilities
	setDataFromFile(props->kr_o, "props/vpp/kr_oil.txt");
	setDataFromFile(props->kr_w, "props/vpp/kr_wat.txt");

	// Defining volume factors
	//setDataFromFile(props->B_o, "props/vpp/Boil.txt");
	//setDataFromFile(props->B_w, "props/vpp/Bwat.txt");

	return props;
}*/

int main(int argc, char** argv)
{
	gasOil_rz::Properties* props = getProps();
	Scene<gasOil_rz::GasOil_RZ, gasOil_rz::GasOil2DSolver, gasOil_rz::Properties> scene;
	scene.load(*props);
	scene.setSnapshotterType("VTK");
	scene.start();

	/*vpp2d::Properties* props = getProps();
	Scene<vpp2d::VPP2d, vpp2d::VPPSolver, vpp2d::Properties> scene;
	scene.load(*props);
	scene.setSnapshotterType("VTK");
	scene.start();*/

	return 0;
}