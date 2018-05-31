#include <new>
#include <vector>
#include <string>
#include <fstream>
#include <utility>
#include <iostream>
#include <valarray>
#include <random>

#include "util/utils.h"
#include "method/mcmath.h"
#include "Scene.h"
#include "method/HypreInterface.hpp"

#include "model/GasOil_RZ/GasOil_RZ.h"
#include "model/GasOil_RZ/GasOil2DSolver.h"

#include "model/VPP2d/VPP2d.hpp"
#include "model/VPP2d/VPPSolver.hpp"

#include "model/Bingham1d/Bingham1d.hpp"
#include "model/Bingham1d/BingSolver.hpp"

#include "model/GasOil_Elliptic/GasOil_Elliptic.hpp"
#include "model/GasOil_Elliptic/GasOilEllipticSolver.hpp"

#include "model/OilNIT_Elliptic/OilNIT_Elliptic.hpp"
#include "model/OilNIT_Elliptic/OilNITEllipticSolver.hpp"

#include "model/BlackOilNIT_Elliptic/BlackOilNIT_Elliptic.hpp"
#include "model/BlackOilNIT_Elliptic/BlackOilNITEllipticSolver.hpp"

#include "model/BlackOil_RZ/BlackOil_RZ.hpp"
#include "model/BlackOil_RZ/BlackOil2DSolver.hpp"

#include "model/Acid/2d/Acid2d.hpp"
#include "model/Acid/2d/Acid2dSolver.hpp"

#include "model/Acid/1d/Acid1d.hpp"
#include "model/Acid/1d/Acid1dSolver.hpp"

#include "model/Acid/2dnit/Acid2dNIT.hpp"
#include "model/Acid/2dnit/Acid2dNITSolver.hpp"

#include "model/Acid/frac/AcidFracModel.hpp"

#include "model/WaxNIT/1d/WaxNIT1d.hpp"
#include "model/WaxNIT/1d/WaxNIT1dSolver.hpp"

#include "model/WaxNIT/2d/WaxNIT.hpp"
#include "model/WaxNIT/2d/WaxNITSolver.hpp"

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
/*gasOil_elliptic::Properties* getProps()
{
	gasOil_elliptic::Properties* props = new gasOil_elliptic::Properties();

	props->cellsNum_mu = 10;
	props->cellsNum_nu = 20;
	props->cellsNum_z = 9;

	props->timePeriods.push_back(10.0 * 86400.0);
	props->timePeriods.push_back(25.0 * 86400.0);

	props->leftBoundIsRate = false;
	props->rightBoundIsPres = true;
	//props->rates.push_back(10.0);
	props->pwf.push_back(65.0 * 1.E+5);
	props->pwf.push_back(70.625 * 1.E+5);

	props->ht = 10.0;
	props->ht_min = 10.0;
	props->ht_max = 1000000.0;

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
	tmp.perm_z = 2.0;
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
}*/
/*oilnit_elliptic::Properties* getProps()
{
	oilnit_elliptic::Properties* props = new oilnit_elliptic::Properties();

	props->cellsNum_mu = 15;
	props->cellsNum_nu = 40;
	props->cellsNum_z = 13;

	props->timePeriods.push_back(10.0 * 86400.0);
	props->timePeriods.push_back(20.0 * 86400.0);

	props->leftBoundIsRate = true;
	props->rightBoundIsPres = true;
	props->rates.push_back(50.0);
	props->rates.push_back(0.0);
	//props->pwf.push_back(150.0 * 1.E+5);
	//props->pwf.push_back(200.0 * 1.E+5);

	props->ht = 1000.0;
	props->ht_min = 100.0;
	props->ht_max = 100000.0;

	props->r_w = 0.1;
	props->r_e = 1000.0;
	props->l = 200.0;

	props->depth_point = 0.0;

	//props->perfIntervals.push_back(make_pair(3, 3));
	props->perfIntervals.push_back(make_pair(7, 7));
	props->perfIntervals.push_back(make_pair(10, 10));
	props->perfIntervals.push_back(make_pair(13, 13));

	oilnit_elliptic::Skeleton_Props tmp;
	tmp.isWellHere = true;
	tmp.cellsNum_z = 13;
	tmp.m = 0.15;
	tmp.p_init = tmp.p_out = 200.0 * 1.0e+5;
	tmp.t_init = 0.0;
	tmp.h1 = 0.0;
	tmp.h2 = 10.0;
	tmp.h_well = 5.0;
	tmp.height = 10.0;
	tmp.perm_mu = 20.0;
	tmp.perm_z = 2.0;
	tmp.dens_stc = 2000.0;
	tmp.beta = 4.35113e-10;
	tmp.lambda_r = 5.0;
	tmp.lambda_z = 5.0;
	tmp.c = 1800.0;

	//tmp.skins.push_back(2.68);
	//tmp.skins.push_back(2.68);
	tmp.skins.push_back(30.0);
	tmp.skins.push_back(30.0);
	tmp.radiuses_eff.push_back(1.0);
	tmp.radiuses_eff.push_back(1.0);

	props->props_sk.push_back(tmp);
	
	tmp.m = 0.10;
	tmp.perm_mu = 10.0;
	tmp.perm_z = 1.0;
	tmp.skins[0] = 50.0;
	tmp.skins[1] = 50.0;
	tmp.radiuses_eff[0] = 1.0;
	tmp.radiuses_eff[1] = 1.0;
	props->props_sk.push_back(tmp);

	tmp.m = 0.15;
	tmp.perm_mu = 20.0;
	tmp.perm_z = 2.0;
	tmp.skins[0] = 100.0;
	tmp.skins[1] = 100.0;
	tmp.radiuses_eff[0] = 1.0;
	tmp.radiuses_eff[1] = 1.0;
	props->props_sk.push_back(tmp);

	tmp.m = 0.20;
	tmp.perm_mu = 30.0;
	tmp.perm_z = 3.0;
	tmp.skins[0] = 0.0;
	tmp.skins[1] = 0.0;
	tmp.radiuses_eff[0] = 1.0;
	tmp.radiuses_eff[1] = 1.0;
	props->props_sk.push_back(tmp);

	tmp.m = 0.15;
	tmp.perm_mu = 20.0;
	tmp.perm_z = 2.0;
	tmp.skins[0] = 30.0;
	tmp.skins[1] = 30.0;
	tmp.radiuses_eff[0] = 1.0;
	tmp.radiuses_eff[1] = 1.0;
	props->props_sk.push_back(tmp);

	props->depth_point = 0.0;

	props->props_oil.visc = 1.0;
	props->props_oil.oil.rho_stc = 887.261;
	props->props_oil.beta = 1.0 * 1.e-9;
	props->props_oil.p_sat = 70.625 * 1.0e+5;
	props->props_oil.jt = 4.0 * 1.e-7;
	props->props_oil.ad = 2.0 * 1.e-7;
	props->props_oil.c = 1880.0;
	props->props_oil.lambda = 0.16;
	
	props->alpha = 3600.0;
	props->wellboreDuration = 86400.0;

	return props;
}*/
/*blackoil_rz::Properties* getProps()
{
	blackoil_rz::Properties* props = new blackoil_rz::Properties();

	props->cellsNum_r = 20;
	props->cellsNum_z = 1;

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

	props->r_w = 0.1;
	props->r_e = 1000.0;

	props->perfIntervals.push_back(make_pair(1, 1));

	blackoil_rz::Skeleton_Props tmp;
	tmp.cellsNum_z = 1;
	tmp.m0 = 0.1;
	tmp.p_init = tmp.p_out = tmp.p_sat = tmp.p_ref = 70.625 * 1.0e+5;
	tmp.sw_init = 0.2;	tmp.so_init = 0.8;
	tmp.s_wc = 0.1;		tmp.s_oc = 0.1;		tmp.s_gc = 0.05;
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

	props->props_oil.visc = 5.0;
	props->props_oil.dens_stc = 887.261;
	props->props_oil.beta = 1.0 * 1.e-9;
	props->props_oil.p_ref = tmp.p_ref;

	props->props_wat.visc = 1.0;
	props->props_wat.dens_stc = 1000.0;
	props->props_wat.beta = 1.0 * 1.e-9;
	props->props_wat.p_ref = tmp.p_ref;

	props->props_gas.visc = 0.03;
	props->props_gas.dens_stc = 0.8;

	// Defining relative permeabilities
	//setDataFromFile(props->kr_oil, "props/koil_tempest.txt");
	//setDataFromFile(props->kr_gas, "props/kgas_tempest.txt");

	// Defining volume factors
	//props->byDefault.B_oil = true;
	setDataFromFile(props->B_oil, "props/Boil_tempest.txt");
	//props->byDefault.B_gas = false;
	setDataFromFile(props->B_gas, "props/Bgas_tempest.txt");

	//props->byDefault.Rs = true;
	setDataFromFile(props->Rs, "props/Rs_tempest.txt");

	return props;
}*/
/*acid2d::Properties* getProps()
{
	acid2d::Properties* props = new acid2d::Properties;

	props->cellsNum_r = 100;
	props->cellsNum_z = 1;

	props->timePeriods.push_back(5.0 * 3600.0);
	props->timePeriods.push_back(30.0 * 3600.0);
	props->leftBoundIsRate = false;
	props->rightBoundIsPres = true;
	//props->rates.push_back(0.0);
	props->pwf.push_back(250.0 * 1.0e+5);
	props->pwf.push_back(200.0 * 1.0e+5);
	props->xa.push_back(0.0);
	props->xa.push_back(0.0);

	props->ht = 0.1;
	props->ht_min = 0.1;
	props->ht_max = 10000.0;

	props->alpha = 7200.0;

	props->r_w = 0.1;
	props->r_e = 150.0;

	props->perfIntervals.push_back(make_pair(1, 1));
	//props->perfIntervals.push_back(make_pair(15, 15));

	acid2d::Skeleton_Props tmp;
	tmp.cellsNum_z = 1;
	tmp.m_init = 0.1;
	tmp.p_init = tmp.p_out = tmp.p_ref = tmp.p_sat = 200.0 * 1.0e+5;
	tmp.sw_init = 0.2;	tmp.so_init = 0.8;
	tmp.xa_init = 0.0;	tmp.xw_init = 1.0;
	tmp.s_wc = 0.0;		tmp.s_oc = 0.0;		tmp.s_gc = 0.0;
	tmp.xa_eqbm = 0.0;
	tmp.h1 = 0.0;
	tmp.h2 = 10.0;
	tmp.height = 10.0;
	tmp.perm_r = 100.0;
	tmp.perm_z = 10.0;
	tmp.dens_stc = 2000.0;
	tmp.beta = 4.35113e-10;
	tmp.skins.push_back(0.0);
	tmp.skins.push_back(0.0);
	tmp.radiuses_eff.push_back(props->r_w);
	tmp.radiuses_eff.push_back(props->r_w);
	props->props_sk.push_back(tmp);

	/*tmp.cellsNum_z = 10;
	tmp.m_init = 0.1;
	tmp.p_init = tmp.p_out = tmp.p_ref = 200.0 * 1.0e+5;
	tmp.sw_init = 0.2;	tmp.so_init = 0.8;
	tmp.xa_init = 0.0;	tmp.xw_init = 1.0;
	tmp.s_wc = 0.0;		tmp.s_oc = 0.0;		tmp.s_gc = 0.00;
	tmp.xa_eqbm = 0.0;
	tmp.h1 = 5.0;
	tmp.h2 = 10.0;
	tmp.height = 5.0;
	tmp.perm_r = 20.0;
	tmp.perm_z = 6.0;
	tmp.dens_stc = 2000.0;
	tmp.beta = 4.35113e-10;
	tmp.skins.push_back(0.0);
	tmp.skins.push_back(0.0);
	tmp.radiuses_eff.push_back(props->r_w);
	tmp.radiuses_eff.push_back(props->r_w);
	props->props_sk.push_back(tmp);

	props->depth_point = 0.0;

	props->props_o.visc = 1.0;
	props->props_o.dens_stc = 887.261;
	props->props_o.beta = 1.0 * 1.e-9;
	props->props_o.p_ref = tmp.p_ref;

	props->props_w.visc = 1.0;
	props->props_w.dens_stc = 1000.0;
	props->props_w.beta = 1.0 * 1.e-9;
	props->props_w.p_ref = tmp.p_ref;

	props->props_g.visc = 0.06;
	props->props_g.dens_stc = 0.8;
	props->props_g.co2 = acid2d::getCO2();

	setDataFromFile(props->B_oil, "props/Boil_tempest.txt");
	setDataFromFile(props->Rs, "props/Rs_tempest.txt");
	setDataFromFile(props->rho_co2, "props/acid/co2_dens.txt");
	
	return props;
}*/
/*acid2dnit::Properties* getProps()
{
	acid2dnit::Properties* props = new acid2dnit::Properties;

	props->cellsNum_r = 100;
	props->cellsNum_z = 5;

	props->timePeriods.push_back(5.0 * 3600.0);
	props->timePeriods.push_back(30.0 * 3600.0);
	props->leftBoundIsRate = false;
	props->rightBoundIsPres = true;
	//props->rates.push_back(0.0);
	props->pwf.push_back(250.0 * 1.0e+5);
	props->pwf.push_back(200.0 * 1.0e+5);
	props->xa.push_back(0.15);
	props->xa.push_back(0.0);
	props->temps.push_back(300.0);
	props->temps.push_back(300.0);

	props->ht = 0.1;
	props->ht_min = 0.1;
	props->ht_max = 10000.0;

	props->alpha = 7200.0;

	props->r_w = 0.1;
	props->r_e = 150.0;

	props->perfIntervals.push_back(make_pair(3, 3));
	//props->perfIntervals.push_back(make_pair(15, 15));

	acid2dnit::Skeleton_Props tmp;
	tmp.cellsNum_z = 5;
	tmp.m_init = 0.1;
	tmp.p_init = tmp.p_out = tmp.p_ref = tmp.p_sat = 200.0 * 1.0e+5;
	tmp.t_init = 300.0;
	tmp.sw_init = 0.2;	tmp.so_init = 0.8;
	tmp.xa_init = 0.0;	tmp.xw_init = 1.0;
	tmp.s_wc = 0.0;		tmp.s_oc = 0.0;		tmp.s_gc = 0.0;
	tmp.xa_eqbm = 0.0;
	tmp.h1 = 0.0;
	tmp.h2 = 10.0;
	tmp.height = 10.0;
	tmp.perm_r = 100.0;
	tmp.perm_z = 10.0;
	tmp.dens_stc = 2000.0;
	tmp.beta = 4.35113e-10;
	tmp.lambda_r = 5.0;
	tmp.lambda_z = 5.0;
	tmp.c = 1800.0;
	tmp.skins.push_back(0.0);
	tmp.skins.push_back(0.0);
	tmp.radiuses_eff.push_back(props->r_w);
	tmp.radiuses_eff.push_back(props->r_w);
	props->props_sk.push_back(tmp);

	/*tmp.cellsNum_z = 10;
	tmp.m_init = 0.1;
	tmp.p_init = tmp.p_out = tmp.p_ref = 200.0 * 1.0e+5;
	tmp.sw_init = 0.2;	tmp.so_init = 0.8;
	tmp.xa_init = 0.0;	tmp.xw_init = 1.0;
	tmp.s_wc = 0.0;		tmp.s_oc = 0.0;		tmp.s_gc = 0.00;
	tmp.xa_eqbm = 0.0;
	tmp.h1 = 5.0;
	tmp.h2 = 10.0;
	tmp.height = 5.0;
	tmp.perm_r = 20.0;
	tmp.perm_z = 6.0;
	tmp.dens_stc = 2000.0;
	tmp.beta = 4.35113e-10;
	tmp.skins.push_back(0.0);
	tmp.skins.push_back(0.0);
	tmp.radiuses_eff.push_back(props->r_w);
	tmp.radiuses_eff.push_back(props->r_w);
	props->props_sk.push_back(tmp);

	props->depth_point = 0.0;

	props->props_o.visc = 1.0;
	props->props_o.dens_stc = 887.261;
	props->props_o.beta = 1.0 * 1.e-9;
	props->props_o.p_ref = tmp.p_ref;
	props->props_o.jt = 4.0 * 1.e-7;
	props->props_o.ad = 2.0 * 1.e-7;
	props->props_o.c = 1880.0;
	props->props_o.lambda = 0.16;

	props->props_w.visc = 1.0;
	props->props_w.dens_stc = 1000.0;
	props->props_w.beta = 1.0 * 1.e-9;
	props->props_w.p_ref = tmp.p_ref;
	props->props_w.jt = 2.0 * 1.e-7;
	props->props_w.ad = 2.0 * 1.e-7;
	props->props_w.c = 1880.0;
	props->props_w.lambda = 0.16;

	props->props_g.visc = 0.06;
	props->props_g.dens_stc = 0.8;
	props->props_g.co2 = acid2dnit::getCO2();
	props->props_g.jt = -1.7 * 1.e-6;
	props->props_g.ad = 3.6 * 1.e-6;
	props->props_g.c = 3400.0;
	props->props_g.lambda = 0.06;

	setDataFromFile(props->B_oil, "props/Boil_tempest.txt");
	setDataFromFile(props->Rs, "props/Rs_tempest.txt");
	setDataFromFile(props->rho_co2, "props/acid/co2_dens.txt");

	return props;
}*/
/*acid1d::Properties* getProps()
{
	acid1d::Properties* props = new acid1d::Properties;

	props->cellsNum_x = 100;

	props->timePeriods.push_back(5.0 * 3600.0);
	props->timePeriods.push_back(30.0 * 3600.0);
	props->leftBoundIsRate = false;
	props->rightBoundIsPres = true;
	//props->rates.push_back(0.0);
	props->pwf.push_back(210.0 * 1.0e+5);
	props->pwf.push_back(200.0 * 1.0e+5);
	props->xa.push_back(0.15);
	props->xa.push_back(0.0);

	props->ht = 0.1;
	props->ht_min = 0.1;
	props->ht_max = 10000.0;

	props->alpha = 7200.0;

	props->r_w = 0.1;
	props->r_e = 150.0;

	auto& tmp = props->props_sk;
	tmp.m_init = 0.1;
	tmp.p_init = tmp.p_out = tmp.p_ref = tmp.p_sat = 200.0 * 1.0e+5;
	tmp.sw_init = 0.2;	tmp.so_init = 0.8;
	tmp.xa_init = 0.0;	tmp.xw_init = 1.0;
	tmp.s_wc = 0.0;		tmp.s_oc = 0.0;		tmp.s_gc = 0.0;
	tmp.xa_eqbm = 0.0;
	tmp.h1 = 0.0;
	tmp.h2 = 10.0;
	tmp.height = 10.0;
	tmp.perm = 100.0;
	tmp.dens_stc = 2000.0;
	tmp.beta = 4.35113e-10;
	tmp.skins.push_back(0.0);
	tmp.skins.push_back(0.0);
	tmp.radiuses_eff.push_back(props->r_w);
	tmp.radiuses_eff.push_back(props->r_w);

	props->depth_point = 0.0;

	props->props_o.visc = 1.0;
	props->props_o.dens_stc = 887.261;
	props->props_o.beta = 1.0 * 1.e-9;
	props->props_o.p_ref = tmp.p_ref;

	props->props_w.visc = 1.0;
	props->props_w.dens_stc = 1000.0;
	props->props_w.beta = 1.0 * 1.e-9;
	props->props_w.p_ref = tmp.p_ref;

	props->props_g.visc = 0.06;
	props->props_g.dens_stc = 0.8;
	props->props_g.co2 = acid1d::getCO2();

	setDataFromFile(props->B_oil, "props/Boil_tempest.txt");
	//setDataFromFile(props->Rs, "props/Rs_tempest.txt");
	setDataFromFile(props->rho_co2, "props/acid/co2_dens.txt");

	return props;
}*/
/*blackoilnit_elliptic::Properties* getProps()
{
	blackoilnit_elliptic::Properties* props = new blackoilnit_elliptic::Properties();

	props->cellsNum_mu = 15;
	props->cellsNum_nu = 28;
	props->cellsNum_z = 7;

	props->timePeriods.push_back(100.0 * 86400.0);
	//props->timePeriods.push_back(20.0 * 86400.0);

	props->leftBoundIsRate = false;
	props->rightBoundIsPres = true;
	//props->rates.push_back(10.0);
	//props->rates.push_back(0.0);
	props->pwf.push_back(50.0 * 1.E+5);
	//props->pwf.push_back(70.625 * 1.E+5);

	props->ht = 1000.0;
	props->ht_min = 100.0;
	props->ht_max = 100000.0;

	props->alpha = 7200.0;

	props->r_w = 0.1;
	props->r_e = 1000.0;
	props->l = 100.0;

	props->depth_point = 0.0;

	//props->perfIntervals.push_back(make_pair(3, 3));
	props->perfIntervals.push_back(make_pair(4, 4));
	props->perfIntervals.push_back(make_pair(7, 7));
	props->perfIntervals.push_back(make_pair(10, 10));

	blackoilnit_elliptic::Skeleton_Props tmp;
	tmp.isWellHere = true;
	tmp.cellsNum_z = 7;
	tmp.m = 0.15;
	tmp.p_init = tmp.p_out = tmp.p_sat = tmp.p_ref = 70.625 * 1.0e+5;
	tmp.sw_init = 0.2;	tmp.so_init = 0.8;
	tmp.s_wc = 0.1;		tmp.s_oc = 0.1;		tmp.s_gc = 0.05;
	tmp.t_init = 0.0;
	tmp.h1 = 0.0;
	tmp.h2 = 10.0;
	tmp.h_well = 5.0;
	tmp.height = 10.0;
	tmp.perm_mu = 20.0;
	tmp.perm_z = 2.0;
	tmp.dens_stc = 2000.0;
	tmp.beta = 4.35113e-10;
	tmp.lambda_r = 0.0;// 5.0;
	tmp.lambda_z = 0.0;// 5.0;
	tmp.c = 1800.0;

	//tmp.skins.push_back(2.68);
	//tmp.skins.push_back(2.68);
	tmp.skins.push_back(30.0);
	tmp.skins.push_back(30.0);
	tmp.radiuses_eff.push_back(1.0);
	tmp.radiuses_eff.push_back(1.0);

	props->props_sk.push_back(tmp);

	tmp.m = 0.10;
	tmp.perm_mu = 10.0;
	tmp.perm_z = 1.0;
	tmp.skins[0] = 30.0;
	tmp.skins[1] = 30.0;
	tmp.radiuses_eff[0] = 1.0;
	tmp.radiuses_eff[1] = 1.0;
	props->props_sk.push_back(tmp);

	tmp.m = 0.15;
	tmp.perm_mu = 20.0;
	tmp.perm_z = 2.0;
	tmp.skins[0] = 30.0;
	tmp.skins[1] = 30.0;
	tmp.radiuses_eff[0] = 1.0;
	tmp.radiuses_eff[1] = 1.0;
	props->props_sk.push_back(tmp);

	tmp.m = 0.2;
	tmp.perm_mu = 30.0;
	tmp.perm_z = 3.0;
	tmp.skins[0] = 30.0;
	tmp.skins[1] = 30.0;
	tmp.radiuses_eff[0] = 1.0;
	tmp.radiuses_eff[1] = 1.0;
	props->props_sk.push_back(tmp);

	tmp.m = 0.15;
	tmp.perm_mu = 20.0;
	tmp.perm_z = 2.0;
	tmp.skins[0] = 30.0;
	tmp.skins[1] = 30.0;
	tmp.radiuses_eff[0] = 1.0;
	tmp.radiuses_eff[1] = 1.0;
	props->props_sk.push_back(tmp);

	props->depth_point = 0.0;

	props->props_oil.visc = 1.0;
	props->props_oil.dens_stc = 887.261;
	props->props_oil.beta = 1.0 * 1.e-9;
	props->props_oil.jt = 0.0;// 4.0 * 1.e-7;
	props->props_oil.ad = 2.0 * 1.e-7;
	props->props_oil.c = 1880.0;
	props->props_oil.lambda = 0.0;// 0.16;

	props->props_wat.visc = 0.6;
	props->props_wat.dens_stc = 1000.0;
	props->props_wat.beta = 1.0 * 1.e-9;
	props->props_wat.p_ref = tmp.p_ref;
	props->props_wat.jt = 0.0;// 2.0 * 1.e-7;
	props->props_wat.ad = 1.5 * 1.e-7;
	props->props_wat.c = 4181.0;
	props->props_wat.lambda = 0.0;// 0.65;

	props->props_gas.visc = 0.03;
	props->props_gas.dens_stc = 0.8;
	props->props_gas.jt = 0.0;// -1.7 * 1.e-6;
	props->props_gas.ad = 3.6 * 1.e-6;
	props->props_gas.c = 3400.0;
	props->props_gas.lambda = 0.0;// 0.06;

	props->L = 0.0;//  -50.0 * 1.e+3;
	// Defining relative permeabilities
	//setDataFromFile(props->kr_oil, "props/koil_tempest.txt");
	//setDataFromFile(props->kr_gas, "props/kgas_tempest.txt");

	// Defining volume factors
	//props->byDefault.B_oil = true;
	setDataFromFile(props->B_oil, "props/Boil_tempest.txt");
	//props->byDefault.B_gas = false;
	setDataFromFile(props->B_gas, "props/Bgas_tempest.txt");

	//props->byDefault.Rs = true;
	setDataFromFile(props->Rs, "props/Rs_tempest.txt");

	return props;
}*/
/*bing1d::Properties* getProps()
{
	bing1d::Properties* props = new bing1d::Properties;

	props->cellsNum_r = 100;

	props->timePeriods.push_back(5000.0 * 24 * 3600.0);
	//props->timePeriods.push_back(30.0 * 3600.0);
	props->leftBoundIsRate = false;
	props->rightBoundIsPres = true;
	//props->rates.push_back(1.0);
	props->pwf.push_back(150.0 * 1.0e+5);
	//props->pwf.push_back(200.0 * 1.0e+5);

	props->ht = 100.0;
	props->ht_min = 100.0;
	props->ht_max = 1000000.0;

	props->r_w = 0.1;
	props->r_e = 1000.0;

	props->perfIntervals.push_back(make_pair(0, 0));

	bing1d::Skeleton_Props tmp;
	tmp.m = 0.1;
	tmp.cellsNum_z = 1;
	tmp.p_init = tmp.p_out = 200.0 * 1.0e+5;
	tmp.height = 10.0;
	tmp.perm_r = 100.0;
	tmp.perm_z = 100.0;
	tmp.dens_stc = 2000.0;
	tmp.beta = 4.35113e-10;
	tmp.skins.push_back(0.0);
	//tmp.skins.push_back(0.0);
	tmp.radiuses_eff.push_back(props->r_w);
	//tmp.radiuses_eff.push_back(props->r_w);
	props->props_sk = tmp;

	props->depth_point = 0.0;

	props->props_oil.visc = 10.0;
	props->props_oil.dens_stc = 887.261;
	props->props_oil.beta = 1.0 * 1.e-9;
	props->props_oil.p_ref = tmp.p_init;
	props->props_oil.tau0 = 500;
	props->props_oil.d = 0.01;
	props->props_oil.m = 0.002;
	setDataFromFile(props->u_dp_dimless, "props/u_dp_500_002.txt");

	return props;
}*/
/*wax_nit::Properties* getProps()
{
	wax_nit::Properties* props = new wax_nit::Properties();

	props->cellsNum_r = 200;
	props->cellsNum_z = 1;

	props->timePeriods.push_back(500.0 * 3600.0);
	props->timePeriods.push_back(1.0 * 370.0 * 86400.0);

	props->leftBoundIsRate = false;
	props->rightBoundIsPres = true;
	props->pwf.push_back(160.625 * 1.0e+5);
	props->pwf.push_back(100.625 * 1.0e+5);
	//props->rates.push_back(0.1);
	//props->rates.push_back(0.0);

	props->ht = 1.0;
	props->ht_min = 1.0;
	props->ht_max = 500000.0;

	props->alpha = 7200.0;

	props->r_w = 0.1;
	props->r_e = 1000.0;

	props->perfIntervals.push_back(make_pair(1, 1));

	wax_nit::Skeleton_Props tmp;
	tmp.cellsNum_z = 1;
	tmp.m_init = 0.1;
	tmp.p_init = tmp.p_out = tmp.p_ref = 180.625 * 1.0e+5;
	tmp.t_init = 291.0;
	tmp.p_sat = 150.625 * 1.0e+5;
	tmp.t_sat = tmp.t_init;
	tmp.sw_init = 0.02;	tmp.so_init = 0.9799;	tmp.sg_init = 0.0;
	tmp.s_wc = 0.0;		tmp.s_oc = 0.0;		tmp.s_gc = 0.0;
	tmp.h1 = 0.0;
	tmp.h2 = 6.5;
	tmp.height = 6.5;
	tmp.perm_r = 450.0;
	tmp.perm_z = 45.0;
	tmp.dens_stc = 2000.0;
	tmp.beta = 4.35113e-10;
	tmp.lambda_r = 5.0;
	tmp.lambda_z = 5.0;
	tmp.c = 1800.0;
	tmp.skins.push_back(0.0);
	tmp.skins.push_back(0.0);

	tmp.radiuses_eff.push_back(props->r_w);
	tmp.radiuses_eff.push_back(props->r_w);

	props->props_sk.push_back(tmp);

	props->depth_point = 0.0;

	props->props_oil.visc = 8.36;
	props->props_oil.dens_stc = 855.6;
	props->props_oil.beta = 1.22 * 1.e-9;
	props->props_oil.p_ref = tmp.p_ref;
	props->props_oil.gamma = 0.5;

	props->props_oil.jt = 1.0 * 1.e-7;
	props->props_oil.ad = 1.0 * 1.e-7;
	props->props_oil.c = 1233.0;
	props->props_oil.lambda = 0.16;

	props->props_wat.visc = 4.0;
	props->props_wat.dens_stc = 1000.0;
	props->props_wat.beta = 1.0 * 1.e-9;
	props->props_wat.p_ref = tmp.p_ref;

	props->props_wat.jt = 2.0 * 1.e-7;
	props->props_wat.ad = 2.0 * 1.e-7;
	props->props_wat.c = 1880.0;
	props->props_wat.lambda = 0.16;

	props->props_gas.visc = 0.03;
	props->props_oil.dens_gas_stc = props->props_gas.dens_stc = 1.45;
	props->props_oil.dens_wax_stc = props->props_wax.dens_stc = 900.0;

	props->props_gas.jt = -1.7 * 1.e-6;
	props->props_gas.ad = 3.6 * 1.e-6;
	props->props_gas.c = 3400.0;
	props->props_gas.lambda = 0.06;

	props->L = -150.0 * 1.e+3;

	// Defining relative permeabilities
	//setDataFromFile(props->kr_oil, "props/koil_tempest.txt");
	//setDataFromFile(props->kr_gas, "props/kgas_tempest.txt");

	// Defining volume factors
	//props->byDefault.B_oil = true;
	setDataFromFile(props->B_oil, "props/new/Boil50.txt");
	//props->byDefault.B_gas = false;
	setDataFromFile(props->B_gas, "props/Bgas.txt");

	//props->byDefault.Rs = true;
	setDataFromFile(props->Rs, "props/new/Rs50.txt");
	setDataFromFile(props->lp, "props/lpx05.txt");

	return props;
}*/
/*wax_nit1d::Properties* getProps()
{
	wax_nit1d::Properties* props = new wax_nit1d::Properties();

	props->cellsNum_x = 100;
	props->timePeriods.push_back(1.0 * 370.0 * 86400.0);

	props->leftBoundIsRate = false;
	props->rightBoundIsPres = true;
	props->pwf.push_back(160.625 * 1.0e+5);

	props->ht = 1.0;
	props->ht_min = 1.0;
	props->ht_max = 500000.0;

	props->alpha = 7200.0;

	props->r_w = 0.1;
	props->r_e = 1.0;

	auto& tmp = props->props_sk;
	tmp.m_init = 0.1;
	tmp.p_init = tmp.p_out = tmp.p_ref = 150.625 * 1.0e+5;
	tmp.p_sat = 150.625 * 1.0e+5;
	tmp.so_init = 0.99999;
	tmp.s_wc = 0.0;		tmp.s_oc = 0.0;		tmp.s_gc = 0.0;
	tmp.h1 = 0.0;
	tmp.h2 = 0.1;
	tmp.height = 0.1;
	tmp.perm = 450.0;
	tmp.dens_stc = 2000.0;
	tmp.beta = 4.35113e-10;
	tmp.skins.push_back(0.0);
	tmp.radiuses_eff.push_back(props->r_w);

	props->depth_point = 0.0;

	props->props_oil.visc = 8.36;
	props->props_oil.dens_stc = 855.6;
	props->props_oil.beta = 1.22 * 1.e-9;
	props->props_oil.p_ref = tmp.p_ref;
	props->props_oil.gamma = 0.5;

	props->props_oil.dens_gas_stc = 1.45;
	props->props_oil.dens_wax_stc = props->props_wax.dens_stc = 900.0;

	// Defining relative permeabilities
	//setDataFromFile(props->kr_oil, "props/koil_tempest.txt");
	//setDataFromFile(props->kr_gas, "props/kgas_tempest.txt");

	// Defining volume factors
	//props->byDefault.B_oil = true;
	setDataFromFile(props->B_oil, "props/new/Boil50.txt");
	//props->byDefault.B_gas = false;
	setDataFromFile(props->B_gas, "props/Bgas.txt");

	//props->byDefault.Rs = true;
	setDataFromFile(props->Rs, "props/new/Rs50.txt");
	setDataFromFile(props->lp, "props/lpx05.txt");

	return props;
}*/
acidfrac::Properties* getProps()
{
	typedef acidfrac::Properties Properties;
	Properties* props = new Properties;

	props->ht = 0.01;
	props->ht_min = 0.01;
	props->ht_max = 1000.0;

	props->timePeriods.push_back(3600.0 / 3);
	props->timePeriods.push_back(5.0 * 3600.0);
	//props->timePeriods.push_back(10.0 * 3600.0);
	//props->leftBoundIsRate = false;
	props->LeftBoundIsRate.push_back(false);
	props->LeftBoundIsRate.push_back(true);
	//props->LeftBoundIsRate.push_back(true);
	props->rightBoundIsPres = true;
	props->pwf.push_back(300.0 * 1.0e+5);
	props->rates.push_back(0.0);
	props->cs.push_back(0.15);
	props->cs.push_back(0.0);

	props->props_frac.l2 = 20.0;
	props->props_frac.w2 = 0.01;

	props->props_frac.p_init = 200.0 * BAR_TO_PA;
	props->props_frac.c_init = 0.0;
	props->props_frac.height = 10.0;

	props->cellsNum_x = 20;
	props->cellsNum_y = 20;
	props->cellsNum_z = 1;

	acidfrac::Skeleton_Props props_sk;
	props_sk.m_init = 0.1;
	props_sk.t_init = 300.0;
	props_sk.p_init = props_sk.p_out = props_sk.p_ref = props->props_frac.p_init;
	props_sk.sw_init = 0.2;		props_sk.so_init = 0.8;
	props_sk.xa_eqbm = 0.0;
	props_sk.xa_init = 0.0;	props_sk.xw_init = 1.0;
	props_sk.s_wc = 0.0;		props_sk.s_oc = 0.0;		props_sk.s_gc = 0.0;
	props_sk.perm = 100.0;
	props_sk.dens_stc = 2000.0;
	props_sk.beta = 4.35113e-10;
	props_sk.height = props->props_frac.height;

	//default_random_engine generator;
	//normal_distribution<double> distribution(100.0, 30.0);
	/*vector<double> perms;
	ifstream file;
	file.open("props/perm.txt", ifstream::in);
	double temp1;
	while (!file.eof())
	{
		file >> temp1;
		perms.push_back(temp1);
		if (file.eof())
			break;
	}
	file.close();*/
	for (int i = 0; i < props->cellsNum_x; i++)
	{
		props->xe.push_back(200.0);
		props->cellsNum_y_1d.push_back(100);
		props->props_sk.push_back(props_sk);
		props->props_sk.back().perm = 100.0;// perms[i];// distribution(generator);
	}

	props->props_o.visc = 1.0;
	props->props_o.dens_stc = 887.261;
	props->props_o.beta = 1.0 * 1.e-9;
	props->props_o.p_ref = props_sk.p_ref;

	props->props_w.visc = 1.0;
	props->props_w.dens_stc = 1000.0;
	props->props_w.beta = 1.0 * 1.e-9;
	props->props_w.p_ref = props_sk.p_ref;
	props->props_w.D_e = 1.E-8;

	props->props_g.visc = 0.06;
	props->props_g.dens_stc = 0.8;
	props->props_g.co2 = acidfrac::getCO2();

	return props;
}

double acid1d::Component::T = 300.0;
double acid2d::Component::T = 300.0;
double acid2dnit::Component::T = 300.0;
double acidfrac::Component::T = 300.0;
double blackoilnit_elliptic::Water_Props::dens_stc = 1000.0;
double blackoilnit_elliptic::Oil_Props::dens_stc = 887.261;
double blackoilnit_elliptic::Gas_Props::dens_stc = 0.8;

int main(int argc, char* argv[])
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
	/*oilnit_elliptic::Properties* props = getProps();
	Scene<oilnit_elliptic::OilNIT_Elliptic, oilnit_elliptic::OilNITEllipticSolver<ParSolver>, oilnit_elliptic::Properties> scene;
	scene.load(*props/*, argc, argv);
	scene.setSnapshotterType("VTK");
	scene.start();*/
	/*blackoilnit_elliptic::Properties* props = getProps();
	Scene<blackoilnit_elliptic::BlackOilNIT_Elliptic, blackoilnit_elliptic::BlackOilNITEllipticSolver<ParSolver>, blackoilnit_elliptic::Properties> scene;
	scene.load(*props);// , argc, argv);
	scene.setSnapshotterType("VTK");
	scene.start();*/
	/*gasOil_elliptic::Properties* props = getProps();
	Scene<gasOil_elliptic::GasOil_Elliptic, gasOil_elliptic::GasOilEllipticSolver, gasOil_elliptic::Properties> scene;
	scene.load(*props);
	scene.setSnapshotterType("VTK");
	scene.start();*/
	/*blackoil_rz::Properties* props = getProps();
	Scene<blackoil_rz::BlackOil_RZ, blackoil_rz::BlackOil2dSolver, blackoil_rz::Properties> scene;
	scene.load(*props);
	scene.setSnapshotterType("VTK");
	scene.start();*/
	/*acid2d::Properties* props = getProps();
	Scene<acid2d::Acid2d, acid2d::Acid2dSolver, acid2d::Properties> scene;
	scene.load(*props);
	scene.setSnapshotterType("VTK");
	scene.start();*/
	/*acid2dnit::Properties* props = getProps();
	Scene<acid2dnit::Acid2dNIT, acid2dnit::Acid2dNITSolver, acid2dnit::Properties> scene;
	scene.load(*props);
	scene.setSnapshotterType("VTK");
	scene.start();*/
	/*bing1d::Properties* props = getProps();
	Scene<bing1d::Bingham1d, bing1d::Bing1dSolver, bing1d::Properties> scene;
	scene.load(*props);
	scene.setSnapshotterType("VTK");
	scene.start();*/
	/*wax_nit1d::Properties* props = getProps();
	Scene<wax_nit1d::WaxNIT1d, wax_nit1d::WaxNIT1dSolver, wax_nit1d::Properties> scene;
	scene.load(*props);
	scene.setSnapshotterType("VTK");
	scene.start();*/
	acidfrac::Properties* props = getProps();
	Scene<acidfrac::AcidFrac, acidfrac::AcidFracSolver, acidfrac::Properties> scene;
	scene.load(*props);
	scene.start();
	/*acid1d::Properties* props = getProps();
	Scene<acid1d::Acid1d, acid1d::Acid1dSolver, acid1d::Properties> scene;
	scene.load(*props);
	scene.setSnapshotterType("VTK");
	scene.start();*/
	return 0;
}