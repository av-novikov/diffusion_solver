#ifndef CONFIG_HPP_
#define CONFIG_HPP_

#include "model/Bingham1d/Bingham1d.hpp"
#include "model/Bingham1d/BingSolver.hpp"
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
#include "model/Acid/ellfrac/AcidEllFracModel.hpp"
#include "model/Acid/recfrac/AcidRecFracModel.hpp"
#include "model/Acid/recfrac/RecFracProdSolver.hpp"
#include "model/Acid/recfrac/RecFracProd.hpp"
#include "model/Acid/recfracmov/AcidRecFracMovModel.hpp"
#include "model/WaxNIT/1d/WaxNIT1d.hpp"
#include "model/WaxNIT/1d/WaxNIT1dSolver.hpp"
#include "model/WaxNIT/2d/WaxNIT.hpp"
#include "model/WaxNIT/2d/WaxNITSolver.hpp"
#include "Scene.h"

#include <exception>
#include <boost/program_options.hpp>

using std::exception;
namespace po = boost::program_options;

namespace issues
{
	template<typename TProp>
	TProp* getProps() { return NULL; };
	template<> blackoil_rz::Properties* getProps<blackoil_rz::Properties>()
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
		tmp.m_init = 0.1;
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
	}
	template<> acid2d::Properties* getProps<acid2d::Properties>()
	{
		acid2d::Properties* props = new acid2d::Properties;

		props->cellsNum_r = 100;
		props->cellsNum_z = 1;

		props->timePeriods.push_back(0.2 * 3600.0);
		props->timePeriods.push_back(5.0 * 3600.0);
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
		props->props_sk.push_back(tmp);*/

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
	}
	template<> acid1d::Properties* getProps<acid1d::Properties>()
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
	}
	template<> blackoilnit_elliptic::Properties* getProps<blackoilnit_elliptic::Properties>()
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
	}
	template<> bing1d::Properties* getProps<bing1d::Properties>()
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
	}
	template<> wax_nit::Properties* getProps<wax_nit::Properties>()
	{
		wax_nit::Properties* props = new wax_nit::Properties();

		props->cellsNum_r = 200;
		props->cellsNum_z = 1;

		props->timePeriods.push_back(120.0 * 3600.0);
		props->timePeriods.push_back(1.0 * 370.0 * 86400.0);

		props->leftBoundIsRate = false;
		props->rightBoundIsPres = true;
		props->pwf.push_back(150.725 * 1.0e+5);
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
		props->props_oil.gamma = 0.25;

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
		setDataFromFile(props->B_oil, "props/new/Boil150.txt");
		//props->byDefault.B_gas = false;
		setDataFromFile(props->B_gas, "props/Bgas.txt");

		//props->byDefault.Rs = true;
		setDataFromFile(props->Rs, "props/new/Rs150.txt");
		setDataFromFile(props->lp, "props/lpx05.txt");

		return props;
	}
	template<> wax_nit1d::Properties* getProps<wax_nit1d::Properties>()
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
	}
	template<> acidfrac::Properties* getProps<acidfrac::Properties>()
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
		props->props_sk.push_back(props_sk);
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
		file.close();
		for (int i = 0; i < props->cellsNum_x; i++)
		{
		props->xe.push_back(200.0);
		props->cellsNum_y_1d.push_back(100);
		props->props_sk.push_back(props_sk);
		props->props_sk.back().perm = 100.0;// perms[i];// distribution(generator);
		}*/

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
	template<> acidellfrac::Properties* getProps<acidellfrac::Properties>()
	{
		typedef acidellfrac::Properties Properties;
		Properties* props = new Properties;

		props->ht = 0.0001;
		props->ht_min = 0.1;
		props->ht_max = 500.0;

		props->timePeriods.push_back(3600.0 / 3.0);
		props->timePeriods.push_back(3.0 * 3600.0);
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
		props->re = props->props_frac.l2;

		props->props_frac.p_init = 200.0 * BAR_TO_PA;
		props->props_frac.c_init = 0.0;
		props->props_frac.height = 10.0;

		props->cellsNum_x = 5;
		props->cellsNum_mu_frac = 10;
		props->cellsNum_mu_poro = 20;
		props->cellsNum_z = 1;

		acidellfrac::Skeleton_Props props_sk;
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
		props->props_sk.push_back(props_sk);
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
		file.close();
		for (int i = 0; i < props->cellsNum_x; i++)
		{
		props->xe.push_back(200.0);
		props->cellsNum_y_1d.push_back(100);
		props->props_sk.push_back(props_sk);
		props->props_sk.back().perm = 100.0;// perms[i];// distribution(generator);
		}*/

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
		props->props_g.co2 = acidellfrac::getCO2();

		return props;
	}
	template<> acidrecfrac::Properties* getProps<acidrecfrac::Properties>()
	{
		typedef acidrecfrac::Properties Properties;
		Properties* props = new Properties;
        props->prefix = "snaps/";

		props->ht = 0.0001;
		props->ht_min = props->ht;
		props->ht_max = 10.0;
 
		props->timePeriods.push_back(0.1 * 3600.0);
		props->timePeriods.push_back(1.0 * 3600.0);
		//props->timePeriods.push_back(10.0 * 3600.0);
		//props->leftBoundIsRate = false;
		props->LeftBoundIsRate.push_back(false);
		//props->LeftBoundIsRate.push_back(false);
		props->LeftBoundIsRate.push_back(true);
		props->rightBoundIsPres = true;
		props->pwf.push_back(420.0 * 1.0e+5);
		//props->pwf.push_back(300.0 * 1.0e+5);
		//props->rates.push_back(-10.0);
		props->rates.push_back(0.0);
		props->cs.push_back(0.15);
		props->cs.push_back(0.0);
		props->max_sol_volume = 50.0;
        
		props->props_frac.l2 = 240.0;
		props->props_frac.w2 = 0.001;
		props->props_frac.height = 18.87;
		props->re = 200.0;

		props->props_frac.p_init = 200.0 * BAR_TO_PA;
		props->props_frac.c_init = 0.0;

		props->cellsNum_x = 20;
		props->cellsNum_y_frac = 20;
		props->cellsNum_y_poro = 75;
		props->cellsNum_z = 1;

        props->prod_props.x_size = props->prod_props.y_size = 1000.0;
        props->prod_props.z_size = props->props_frac.height;
        props->prod_props.nx = 110;
        props->prod_props.ny = 100;
        props->R_dim = props->prod_props.R_dim = props->props_frac.l2 / 5.0;

		acidrecfrac::Skeleton_Props props_sk;
		props_sk.m_init = 0.09;
		props_sk.m_max = 0.4;
		props_sk.A = 60.0;
		props_sk.t_init = 300.0;
		props_sk.p_init = props_sk.p_out = props_sk.p_ref = props->props_frac.p_init;
		props_sk.sw_init = 0.01;					props_sk.so_init = 0.99;
		props_sk.xa_eqbm = 0.0;
		props_sk.xa_init = 0.0;					props_sk.xw_init = 1.0;
		props_sk.xa_init = props_sk.xa_eqbm;	props_sk.xw_init = 1.0 - props_sk.xa_eqbm;
		props_sk.s_wc = 0.0;					props_sk.s_oc = 0.0;		props_sk.s_gc = 0.0;
		props_sk.perm = 0.5;
		props_sk.dens_stc = 2000.0;
		props_sk.beta = 4.35113e-10;
		props_sk.height = props->props_frac.height;
		props->props_sk.push_back(props_sk);
		/*default_random_engine generator;
		normal_distribution<double> distribution(0.15, 0.02);
		for (int i = 0; i < props->cellsNum_x * props->cellsNum_z; i++)
		{
		acidrecfrac::Skeleton_Props prop = props_sk;
		double tmp = distribution(generator);
		prop.m_init = (tmp > 0.01) ? tmp : 0.01;
		prop.perm = props_sk.getPermCoseni(prop.m_init, 0.0).value();
		props->props_sk.push_back(prop);
		}*/
		/*props_sk.height = 3.05;
		props_sk.m_init = 0.08;
		props_sk.perm = 0.3;
		props->props_sk.push_back(props_sk);

		props_sk.height = 0.91;
		props_sk.m_init = 0.1;
		props_sk.perm = 1;
		props->props_sk.push_back(props_sk);

		props_sk.height = 0.61;
		props_sk.m_init = 0.12;
		props_sk.perm = 1;
		props->props_sk.push_back(props_sk);*/

		props->props_o.visc = 4.75;
		props->props_o.dens_stc = 887.261;
		props->props_o.beta = 1.0 * 1.e-9;
		props->props_o.p_ref = props_sk.p_ref;

		props->props_w.visc = 12.0;
		props->props_w.dens_stc = 1000.0;
		props->props_w.beta = 1.0 * 1.e-9;
		props->props_w.p_ref = props_sk.p_ref;
		props->props_w.D_e = 0.0;// 1.E-8;

		props->props_g.visc = 0.06;
		props->props_g.dens_stc = 0.8;
		props->props_g.co2 = acidrecfrac::getCO2();

		return props;
	}
	template<> acidrecfracmov::Properties* getProps<acidrecfracmov::Properties>()
	{
		typedef acidrecfracmov::Properties Properties;
		Properties* props = new Properties;

		props->ht = 0.0005;
		props->ht_min = props->ht;
		props->ht_max = 10.0;

		props->timePeriods.push_back(0.5 * 3600.0);
		props->timePeriods.push_back(2.0 * 3600.0);
		//props->timePeriods.push_back(10.0 * 3600.0);
		//props->leftBoundIsRate = false;
		props->LeftBoundIsRate.push_back(false);
		//props->LeftBoundIsRate.push_back(false);
		props->LeftBoundIsRate.push_back(true);
		props->rightBoundIsPres = true;
		props->pwf.push_back(300.0 * 1.0e+5);
		//props->pwf.push_back(300.0 * 1.0e+5);
		//props->rates.push_back(-10.0);
		props->rates.push_back(0.0);
		props->cs.push_back(0.15);
		props->cs.push_back(0.0);
		props->max_sol_volume = 50.0;

        props->cellsNum_x = 20;
        props->cellsNum_y_frac = 20;
        props->cellsNum_y_poro = 75;
        props->cellsNum_z = 1;

		props->props_frac.l2 = 240.0;
        props->props_frac.w2.resize(props->cellsNum_x * props->cellsNum_z);
        props->props_frac.w2 = props->props_frac.w2_init = 0.001;
		props->props_frac.height = 18.87;
		props->re = 200.0;

		props->props_frac.p_init = 200.0 * BAR_TO_PA;
		props->props_frac.c_init = 0.0;

		acidrecfracmov::Skeleton_Props props_sk;
		props_sk.m_init = 0.09;
		props_sk.m_max = 0.4;
        props_sk.perm_max = 1.E+9;
		props_sk.A = 60.0;
		props_sk.t_init = 300.0;
		props_sk.p_init = props_sk.p_out = props_sk.p_ref = props->props_frac.p_init;
		props_sk.sw_init = 0.1;					props_sk.so_init = 0.9;
		props_sk.xa_eqbm = 0.0;
		props_sk.xa_init = 0.0;					props_sk.xw_init = 1.0;
		props_sk.xa_init = props_sk.xa_eqbm;	props_sk.xw_init = 1.0 - props_sk.xa_eqbm;
		props_sk.s_wc = 0.0;					props_sk.s_oc = 0.0;		props_sk.s_gc = 0.0;
		props_sk.perm = 5.0;
		props_sk.dens_stc = 2000.0;
		props_sk.beta = 4.35113e-10;
		props_sk.height = props->props_frac.height;
		props->props_sk.push_back(props_sk);
		/*default_random_engine generator;
		normal_distribution<double> distribution(0.15, 0.02);
		for (int i = 0; i < props->cellsNum_x * props->cellsNum_z; i++)
		{
		acidrecfrac::Skeleton_Props prop = props_sk;
		double tmp = distribution(generator);
		prop.m_init = (tmp > 0.01) ? tmp : 0.01;
		prop.perm = props_sk.getPermCoseni(prop.m_init, 0.0).value();
		props->props_sk.push_back(prop);
		}*/
		/*props_sk.height = 3.05;
		props_sk.m_init = 0.08;
		props_sk.perm = 0.3;
		props->props_sk.push_back(props_sk);

		props_sk.height = 0.91;
		props_sk.m_init = 0.1;
		props_sk.perm = 1;
		props->props_sk.push_back(props_sk);

		props_sk.height = 0.61;
		props_sk.m_init = 0.12;
		props_sk.perm = 1;
		props->props_sk.push_back(props_sk);*/

		props->props_o.visc = 4.75;
		props->props_o.dens_stc = 887.261;
		props->props_o.beta = 1.0 * 1.e-9;
		props->props_o.p_ref = props_sk.p_ref;

		props->props_w.visc = 12.0;
		props->props_w.dens_stc = 1000.0;
		props->props_w.beta = 1.0 * 1.e-9;
		props->props_w.p_ref = props_sk.p_ref;
		props->props_w.D_e = 0.0;// 1.E-8;

		props->props_g.visc = 0.06;
		props->props_g.dens_stc = 0.8;
		props->props_g.co2 = acidrecfracmov::getCO2();

		return props;
	}
	template<> acid2dnit::Properties* getProps<acid2dnit::Properties>()
	{
		acid2dnit::Properties* props = new acid2dnit::Properties;

		props->cellsNum_r = 100;
		props->cellsNum_z = 1;

		props->timePeriods.push_back(0.5 * 3600.0);
		props->timePeriods.push_back(2.0 * 3600.0);
		//props->leftBoundIsRate = false;
		props->LeftBoundIsRate.push_back(false);
		props->LeftBoundIsRate.push_back(true);
		props->rightBoundIsPres = true;
		props->pwf.push_back(400.0 * 1.0e+5);
		//props->pwf.push_back(200.0 * 1.0e+5);
		props->rates.push_back(0.0);
		props->xa.push_back(0.15);
		props->xa.push_back(1.E-5);
		props->temps.push_back(300.0);
		props->temps.push_back(300.0);

		props->ht = 0.01;
		props->ht_min = 0.01;
		props->ht_max = 10.0;

		props->alpha = 7200.0;

		props->r_w = 0.1;
		props->r_e = 200.0;

		props->perfIntervals.push_back(make_pair(1, 1));
		//props->perfIntervals.push_back(make_pair(15, 15));
		props->max_sol_volume = 10.0;
		acid2dnit::Skeleton_Props tmp;
		tmp.cellsNum_z = 1;
		tmp.m_init = 0.1;
		tmp.m_max = 0.4;
		tmp.A = 60.0;
		tmp.p_init = tmp.p_out = tmp.p_ref = tmp.p_sat = 200.0 * 1.0e+5;
		tmp.t_init = 300.0;
		tmp.xa_eqbm = 1.E-5;
		tmp.sw_init = 0.2;	tmp.so_init = 0.8;
		tmp.xa_init = tmp.xa_eqbm;	tmp.xw_init = 1.0 - tmp.xa_eqbm;
		tmp.s_wc = 0.0;		tmp.s_oc = 0.0;		tmp.s_gc = 0.0;
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
		props->props_sk.push_back(tmp);*/

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
	}

    template <typename Properties>
    int setOptions(Properties& props, int ac, char* av[])
    {
        try
        {
            po::options_description desc("Allowed options");
            desc.add_options()
                ("help", "produce help message")
                ("dir", po::value<std::string>(), "set output directory")
                ("pwf", po::value<double>(), "set injection pressure")
                ("perm", po::value<double>(), "set initial permeability")
                ;

            po::variables_map vm;
            po::store(po::parse_command_line(ac, av, desc), vm);
            po::notify(vm);

            if (vm.count("help")) { cout << desc << "\n"; return 0; }
            if (vm.count("dir"))        props->prefix = vm["dir"].as<std::string>();
            if (vm.count("pwf"))        props->pwf[0] = vm["pwf"].as<double>() * BAR_TO_PA;
            if (vm.count("perm"))       props->props_sk[0].perm = vm["perm"].as<double>();
        }
        catch (exception& e)
        {
            cerr << "error: " << e.what() << "\n";
            return 1;
        }
        catch (...)
        {
            cerr << "Exception of unknown type!\n";
        }
    }

	template <typename TProps, typename TModel, typename TSolver>
	int run(int ac, char* av[])
	{
		TProps* props = issues::getProps<TProps>();
		Scene<TModel, TSolver, TProps> scene;
		scene.load(*props);
		scene.start();
	};
    template <>
    int run<acidrecfrac::Properties, acidrecfrac::AcidRecFrac, acidrecfrac::AcidRecFracSolver>(int ac, char* av[])
    { 
        using namespace acidrecfrac;
        auto props = issues::getProps<Properties>();
        int res = setOptions(props, ac, av);
        Scene<AcidRecFrac, AcidRecFracSolver, Properties> scene;
        scene.load(*props);
		scene.start();

        auto model0 = scene.getModel();
        acidrecfrac_prod::RecFracProd model;
        model.load(*props, model0->getPoroMesh());
        acidrecfrac_prod::RecFracProdSolver method(&model);
        method.start();

		return res;
    };
};

double acid1d::Component::T = 300.0;
double acid2d::Component::T = 300.0;
double acid2dnit::Component::T = 300.0;
double acidfrac::Component::T = 300.0;
double acidellfrac::Component::T = 300.0;
double acidrecfrac::Component::T = 300.0;
double acidrecfracmov::Component::T = 300.0;
double blackoilnit_elliptic::Water_Props::dens_stc = 1000.0;
double blackoilnit_elliptic::Oil_Props::dens_stc = 887.261;
double blackoilnit_elliptic::Gas_Props::dens_stc = 0.8;

#endif /* CONFIG_HPP_ */
