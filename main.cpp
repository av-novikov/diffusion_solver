#include <new>
#include <vector>
#include <string>
#include <fstream>
#include <utility>
#include <iostream>
#include <valarray>

#include "gtest/gtest.h"

#include "util/utils.h"
#include "method/mcmath.h"
#include "Scene.h"

#include "model/GasOil_RZ/GasOil_RZ.h"

using namespace std;

/*gasOil_rz::Properties* getProps()
{
	gasOil_rz::Properties* props = new gasOil_rz::Properties();

	props->cellsNum_r = 100;
	props->cellsNum_z = 5;

	props->timePeriods.push_back(1296000.0);
	props->timePeriods.push_back(2000000.0);
	
	props->leftBoundIsRate = true;
	props->rightBoundIsPres = true;
	props->rates.push_back(90.0);
	props->rates.push_back(0.0);
 
	props->skins.push_back(0.0); 
	props->skins.push_back(0.0); 

	props->radius.push_back(0.1524);
	props->radius.push_back(0.1524);

	props->ht = 10.0;
	props->ht_min = 1.0;
	props->ht_max  = 100000.0;
	
	props->alpha = 7200.0;

	props->perfIntervals.push_back( make_pair(1, 3) );

	props->r_w = 0.1524;
	props->r_e = 3000.0;

	props->p_init = 160.0 * 1.0e+5;
	props->p_out = 160.0 * 1.0e+5;
	props->p_sat = 160.0 * 1.0e+5;
	props->s_init = 0.999;
	props->h1 = 1500.0;
	props->h2 = 1505.0;
	props->depth_point = 1500.0;

	props->perm_r.reserve(props->cellsNum_z+2);
	props->perm_z.reserve(props->cellsNum_z+2);
	double perm = 50.0;
	for(int i = 0; i < props->cellsNum_z+2; i++)
	{
		/*if(i < 10)
			perm = 50.0;
		else
			perm = 1000.0;
		props->perm_r.push_back( perm );
		props->perm_z.push_back( 0.0 );//perm / 50.0 );
	}

	props->m = 0.01;
	props->visc_gas = 0.02833;
	props->visc_oil = 0.25137;
	props->b_oil_bore = 1.56;
	props->dens_oil_stc = 855.0;
	props->dens_gas_stc = 0.72275;
	props->dens_sk_stc = 2000.0;
	props->beta_oil = 1.282*1.E-9;
	props->beta_sk = 6.41*1.E-10;
	
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
}*/

gasOil_rz::Properties* getProps()
{
	gasOil_rz::Properties* props = new gasOil_rz::Properties();

	props->cellsNum_r = 50;
	props->cellsNum_z = 5;

	props->timePeriods.push_back(20.0 * 86400.0);
	props->timePeriods.push_back(40.0 * 86400.0);
	//props->timePeriods.push_back(15.0 * 86400.0);
	//props->timePeriods.push_back(30.0 * 86400.0);

	props->leftBoundIsRate = true;
	props->rightBoundIsPres = true;
	props->rates.push_back(50.0);
	props->rates.push_back(0.0);
	//props->pwf.push_back(100.0 * 1.E+5);
	//props->pwf.push_back(200.0 * 1.E+5);
	//props->pwf.push_back(100.0 * 1.E+5);
	//props->pwf.push_back(200.0 * 1.E+5);

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
	tmp.perm_r = 50.0;
	tmp.perm_z = 50.0;
	tmp.dens_stc = 2000.0;
	tmp.beta = 4.35113e-10;

	tmp.skins.push_back(0.0);
	tmp.skins.push_back(0.0);
	//tmp.skins.push_back(0.0);
	//tmp.skins.push_back(0.0);

	tmp.radiuses_eff.push_back(props->r_w);
	tmp.radiuses_eff.push_back(props->r_w);
	//tmp.radiuses_eff.push_back(props->r_w);
	//tmp.radiuses_eff.push_back(props->r_w);

	props->props_sk.push_back(tmp);

	props->depth_point = 0.0;

	props->props_oil.visc = 1.0;
	props->props_oil.b_bore = 1.0;
	props->props_oil.dens_stc = 736.0;
	props->props_oil.beta = 1.0 * 1.e-9;

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

/*gasOil_rz::Properties* getProps()
{
	gasOil_rz::Properties* props = new gasOil_rz::Properties();

	props->cellsNum_r = 100;
	props->cellsNum_z = 5;

	props->timePeriods.push_back(121993.0);
	props->timePeriods.push_back(321813.0);
	props->timePeriods.push_back(494420.0);
	props->timePeriods.push_back(674042.0);
	props->timePeriods.push_back(1965038.0);
	
	props->rates.push_back(0.0);
	props->rates.push_back(33.0);
	props->rates.push_back(10.7);
	props->rates.push_back(46.8599);
	props->rates.push_back(0.0);
 
	props->skins.push_back(0.0); 
	props->skins.push_back(3.5); 
	props->skins.push_back(11.0);
	props->skins.push_back(5.6);
	props->skins.push_back(5.6);

	props->radius.push_back(3.0);
	props->radius.push_back(3.0);
	props->radius.push_back(3.0);
	props->radius.push_back(1.2);
	props->radius.push_back(1.2);

	props->ht = 10.0;
	props->ht_min = 1.0;
	props->ht_max  = 100000.0;
	
	props->alpha = 7200.0;

	props->perfIntervals.push_back( make_pair(2, 4) );

	props->r_w = 0.1524;
	props->r_e = 765.0;

	props->p_init = 169.0 * 101325.0;
	props->p_sat = 169.0 * 101325.0;
	props->s_init = 0.999;
	props->h1 = 1500.0;
	props->h2 = 1513.33;
	props->depth_point = 1500.0;

	props->perm_r.reserve(props->cellsNum_z+2);
	props->perm_z.reserve(props->cellsNum_z+2);
	double perm = 380.0;
	for(int i = 0; i < props->cellsNum_z+2; i++)
	{
		/*if(i < 10)
			perm = 50.0;
		else
			perm = 1000.0;
		props->perm_r.push_back( perm );
		props->perm_z.push_back( perm / 5.0 );
	}

	props->m = 0.28;
	props->visc_gas = 0.02833;
	props->visc_oil = 9.614;
	props->b_oil_bore = 1.1;
	props->dens_oil_stc = 855.0;
	props->dens_gas_stc = 0.72275;
	props->dens_sk_stc = 2000.0;
	props->beta_oil = 1.0*1.E-9;
	props->beta_sk = 1.0*1.E-9;

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
}*/

int main(int argc, char** argv)
{
	gasOil_rz::Properties* props = getProps();
	Scene<gasOil_rz::GasOil_RZ, gasOil_rz::GasOil2DSolver, gasOil_rz::Properties> scene;
	scene.load(*props);
	scene.setSnapshotterType("VTK");
	scene.start();

	return 0;
}