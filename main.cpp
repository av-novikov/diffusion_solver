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
#include "model/GasOil_RZ/GasOil2DSolver.h"

#include "model/VPP2d/VPP2d.hpp"
#include "model/VPP2d/VPPSolver.hpp"

#include "model/Bingham1d/Bingham1d.hpp"
#include "model/Bingham1d/BingSolver.hpp"

#include "model/GasOil_Elliptic/GasOil_Elliptic.hpp"
#include "model/GasOil_Elliptic/GasOilEllipticSolver.hpp"

using namespace std;

/*gasOil_rz::Properties* getProps()
{
	gasOil_rz::Properties* props = new gasOil_rz::Properties();

	props->cellsNum_r = 20;
	props->cellsNum_z = 1;

	props->timePeriods.push_back(365.0 * 86400.0);
	props->timePeriods.push_back(400.0 * 86400.0);

	props->leftBoundIsRate = true;
	props->rightBoundIsPres = true;
	props->rates.push_back(170.0);
	props->rates.push_back(0.0);

	props->ht = 100.0;
	props->ht_min = 100.0;
	props->ht_max = 100000.0;

	props->alpha = 7200.0;

	props->r_w = 0.1;
	props->r_e = 1000.0;

	props->perfIntervals.push_back(make_pair(1, 1));

	gasOil_rz::Skeleton_Props tmp;
	tmp.cellsNum_z = 1;
	tmp.m = 0.1;
	tmp.p_init = tmp.p_out = tmp.p_bub = 200.0 * 1.0e+5;
	tmp.p_bub = 70.625 * 1.0e+5;
	tmp.s_init = 1.0;
	tmp.h1 = 0.0;
	tmp.h2 = 10.0;
	tmp.height = 10.0;
	tmp.perm_r = 20.0;
	tmp.perm_z = 20.0;
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
	props->props_oil.dens_stc = 887.261;
	props->props_oil.beta = 1.0 * 1.e-9;
	props->props_oil.p_sat = 70.625 * 1.0e+5;

	props->props_gas.visc = 0.03;
	props->props_gas.dens_stc = 0.8;

	// Defining relative permeabilities
	setDataFromFile(props->kr_oil, "props/koil_tempest.txt");
	setDataFromFile(props->kr_gas, "props/kgas_tempest.txt");

	// Defining volume factors
	//props->byDefault.B_oil = true;
	setDataFromFile(props->B_oil, "props/Boil_tempest.txt");
	//props->byDefault.B_gas = false;
	setDataFromFile(props->B_gas, "props/Bgas_tempest.txt");

	//props->byDefault.Rs = true;
	setDataFromFile(props->Rs, "props/Rs_tempest.txt");

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

gasOil_elliptic::Properties* getProps()
{
	gasOil_elliptic::Properties* props = new gasOil_elliptic::Properties();

	props->cellsNum_mu = 10;
	props->cellsNum_nu = 20;
	props->cellsNum_z = 9;

	props->timePeriods.push_back(10.0 * 86400.0);

	props->leftBoundIsRate = false;
	props->rightBoundIsPres = true;
	//props->rates.push_back(10.0);
	props->pwf.push_back(65.0 * 1.E+5);

	props->ht = 10.0;
	props->ht_min = 10.0;
	props->ht_max = 100000.0;

	props->alpha = 7200.0;

	props->r_w = 0.1;
	props->r_e = 1000.0;
	props->l = 100.0;

	props->depth_point = 0.0;

	gasOil_elliptic::Skeleton_Props tmp;
	tmp.isWellHere = true;
	tmp.cellsNum_z = 9;
	tmp.m = 0.1;
	tmp.p_init = tmp.p_out = tmp.p_bub = 70.625 * 1.0e+5;
	tmp.p_bub = 70.625 * 1.0e+5;
	tmp.s_init = 1.0;
	tmp.h1 = 0.0;
	tmp.h2 = 10.0;
	tmp.h_well = 5.0;
	tmp.height = 10.0;
	tmp.perm_mu = 20.0;
	tmp.perm_z = 20.0;
	tmp.dens_stc = 2000.0;
	tmp.beta = 4.35113e-10;
	tmp.skins.push_back(0.0);
	tmp.skins.push_back(0.0);
	tmp.radiuses_eff.push_back(props->r_w);

	props->props_sk.push_back(tmp);

	props->depth_point = 0.0;

	props->props_oil.visc = 1.0;
	props->props_oil.b_bore = 1.0;
	props->props_oil.dens_stc = 887.261;
	props->props_oil.beta = 1.0 * 1.e-9;
	props->props_oil.p_sat = 70.625 * 1.0e+5;

	props->props_gas.visc = 0.03;
	props->props_gas.dens_stc = 0.8;

	// Defining relative permeabilities
	setDataFromFile(props->kr_oil, "props/koil_tempest.txt");
	setDataFromFile(props->kr_gas, "props/kgas_tempest.txt");

	// Defining volume factors
	//props->byDefault.B_oil = true;
	setDataFromFile(props->B_oil, "props/Boil_tempest.txt");
	//props->byDefault.B_gas = false;
	setDataFromFile(props->B_gas, "props/Bgas_tempest.txt");

	//props->byDefault.Rs = true;
	setDataFromFile(props->Rs, "props/Rs_tempest.txt");

	return props;
}

int main(int argc, char** argv)
{
	/*gasOil_rz::Properties* props = getProps();
	Scene<gasOil_rz::GasOil_RZ, gasOil_rz::GasOil2DSolver, gasOil_rz::Properties> scene;
	scene.load(*props);
	scene.setSnapshotterType("VTK");
	scene.start();*/

	/*vpp2d::Properties* props = getProps();
	Scene<vpp2d::VPP2d, vpp2d::VPPSolver, vpp2d::Properties> scene;
	scene.load(*props);
	scene.setSnapshotterType("VTK");
	scene.start();*/

	gasOil_elliptic::Properties* props = getProps();
	Scene<gasOil_elliptic::GasOil_Elliptic, gasOil_elliptic::GasOilEllipticSolver, gasOil_elliptic::Properties> scene;
	scene.load(*props);
	scene.setSnapshotterType("VTK");
	//scene.getModel()->snapshot_all(0);
	scene.start();

	return 0;
}