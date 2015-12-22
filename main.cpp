#include <new>
#include <vector>
#include <string>
#include <fstream>
#include <utility>
#include <iostream>

#include "gtest/gtest.h"

#include "method/mcmath.h"
#include "Scene.h"
#include "model/Oil1D/Oil1D.h"
#include "model/Gas1D/Gas1D.h"
#include "model/Oil1D_NIT/Oil1D_NIT.h"
#include "model/Oil_RZ/Oil_RZ.h"
#include "model/GasOil_RZ/GasOil_RZ.h"
#include "model/GasOil_RZ_NIT/GasOil_RZ_NIT.h"

using namespace std;

void setDataFromFile(vector< pair<double,double> >& vec, string fileName)
{
	ifstream file;
	file.open(fileName.c_str(), ifstream::in);
	
	double temp1, temp2;
	while( !file.eof() )
	{
		file >> temp1;
		if( file.eof() )
			break;
		file >> temp2;
		vec.push_back(make_pair(temp1, temp2));
	}

	file.close();
}

/*gasOil_rz_NIT::Properties* getProps()
{
	gasOil_rz_NIT::Properties* props = new gasOil_rz_NIT::Properties();

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

	props->perfIntervals.push_back( make_pair(3, 3) );

	props->r_w = 0.1524;
	props->r_e = 765.0;

	props->p_init = 169.0 * 101325.0;
	props->p_sat = 169.0 * 101325.0;
	props->T_init = 302.058;
	props->s_init = 0.999;
	props->h1 = 1500.0;
	props->h2 = 1540.0;
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
		props->perm_z.push_back( 0.0 );
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

	// Thermal properties
	props->jt_oil = 4.3*1.e-6;
	props->jt_gas = -1.4*1.e-5;
	props->ad_oil = 1.0e-6;
	props->ad_gas = 1.0e-6;
	props->c_oil = 1600.0;
	props->c_gas = 2000.0;
	props->c_sk = 800.0;
	props->kappa_eff = 0.0;//5.5*1e-6;
	props->L = -50.0*1.e3;

	props->lambda_sk_r = 41.0;
	props->lambda_sk_z = 41.0;
	props->lambda_oil = 5.4;
	props->lambda_gas = 0.05;
	
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

/*gasOil_rz_NIT::Properties* getProps()
{
	gasOil_rz_NIT::Properties* props = new gasOil_rz_NIT::Properties();

	props->cellsNum_r = 100;
	props->cellsNum_z = 25;

	props->timePeriods.push_back(1296000.0);
	props->timePeriods.push_back(2000000.0);
	
	props->rates.push_back(30.0);
	props->rates.push_back(0.0);
 
	props->skins.push_back(0.0); 
	props->skins.push_back(0.0); 

	props->radius.push_back(0.1524);
	props->radius.push_back(0.1524);

	props->ht = 10.0;
	props->ht_min = 1.0;
	props->ht_max  = 100000.0;
	
	props->alpha = 7200.0;

	props->perfIntervals.push_back( make_pair(9, 9) );
	props->perfIntervals.push_back( make_pair(11, 11) );
	props->perfIntervals.push_back( make_pair(13, 13) );
	props->perfIntervals.push_back( make_pair(15, 15) );

	props->r_w = 0.1524;
	props->r_e = 3000.0;

	props->p_init = 160.0 * 1.0e+5;
	props->p_sat = 160.0 * 1.0e+5;
	props->T_init = 302.058;
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
		props->perm_z.push_back( perm / 50.0 );
	}

	props->m = 0.01;
	props->visc_gas = 0.02833;
	props->visc_oil = 0.25137;
	props->b_oil_bore = 1.1;
	props->dens_oil_stc = 800.026;
	props->dens_gas_stc = 0.72275;
	props->dens_sk_stc = 2000.0;
	props->beta_oil = 4.0*1.E-9;
	props->beta_sk = 0.0;

	// Thermal properties
	props->jt_oil = 4.5*1.e-6;
	props->jt_gas = -1.4*1.e-5;
	props->ad_oil = 2.e-6;
	props->ad_gas = 2.e-6;
	props->c_oil = 1600.0;
	props->c_gas = 2000.0;
	props->c_sk = 800.0;
	props->kappa_eff = 0.0;//5.5*1e-6;
	props->L = -50.0*1.e3;

	props->lambda_sk_r = 7.0;
	props->lambda_sk_z = 5.0;
	props->lambda_oil = 0.4;
	props->lambda_gas = 0.05;
	
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

/*gasOil_rz::Properties* getProps()
{
	gasOil_rz::Properties* props = new gasOil_rz::Properties();

	props->cellsNum_r = 100;
	props->cellsNum_z = 3;

	props->timePeriods.push_back(10.0 * 365.0 * 86400.0);
	
	props->leftBoundIsRate = true;
	props->rightBoundIsPres = false;
	props->rates.push_back(5.0);

	props->ht = 1000.0;
	props->ht_min = 100.0;
	props->ht_max  = 5000000.0;
	
	props->alpha = 7200.0;

	props->perfIntervals.push_back( make_pair(1, 3) );

	props->r_w = 0.1;
	props->r_e = 2500.0;

	gasOil_rz::Skeleton_Props tmp;
	tmp.cellsNum_z = 3;
	tmp.m = 0.1;
	tmp.p_init = tmp.p_out = tmp.p_bub = 160.0 * 1.0e+5;
	tmp.s_init = 0.999;
	tmp.h1 = 1500.0;
	tmp.h2 = 1503.0;
	tmp.height = 3.0;
	tmp.perm_r = 100.0;
	tmp.perm_z = 0.0;
	tmp.dens_stc = 2000.0;
	tmp.beta = 4.35113e-10;
	tmp.skins.push_back(0.0);
	tmp.radiuses_eff.push_back(props->r_w);
	props->props_sk.push_back( tmp );

	props->depth_point = 1500.0;

	props->props_oil.visc = 0.62807;
	props->props_oil.b_bore = 1.09754;
	props->props_oil.dens_stc = 800.026;
	props->props_oil.beta = 1.24703e-09;

	props->props_gas.visc = 0.02833;
	props->props_gas.dens_stc = 0.72275;

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

/*oil_rz::Properties* getProps()
{
	oil_rz::Properties* props = new oil_rz::Properties();

	props->cellsNum_r = 100;
	props->cellsNum_z = 5;

	props->timePeriods.push_back(10.0 * 365.0 * 86400.0);
	
	props->leftBoundIsRate = true;
	props->rightBoundIsPres = false;
	props->rates.push_back(100.0);

	props->ht = 1000.0;
	props->ht_min = 100.0;
	props->ht_max  = 5000000.0;
	
	props->alpha = 7200.0;

	props->perfIntervals.push_back( make_pair(1, 1) );
	props->perfIntervals.push_back( make_pair(3, 3) );
	props->perfIntervals.push_back( make_pair(5, 5) );

	props->r_w = 0.1;
	props->r_e = 2500.0;

	oil_rz::Skeleton_Props tmp;
	tmp.cellsNum_z = 1;
	tmp.m = 0.1;
	tmp.p_init = tmp.p_out = 160.0 * 1.0e+5;
	tmp.h1 = 1500.0;
	tmp.h2 = 1530.0;
	tmp.height = 30.0;
	tmp.perm_r = 100.0;
	tmp.perm_z = 10.0;
	tmp.dens_stc = 2000.0;
	tmp.beta = 4.35113e-10;
	tmp.skins.push_back(0.0);
	tmp.radiuses_eff.push_back(props->r_w);
	props->props_sk.push_back( tmp );

	tmp.cellsNum_z = 1;
	tmp.m = 0.1;
	tmp.p_init = tmp.p_out = 0.0 * 1.0e+5;
	tmp.h1 = 1530.0;
	tmp.h2 = 1540.0;
	tmp.height = 10.0;
	tmp.perm_r = 0.0;
	tmp.perm_z = 0.0;
	tmp.dens_stc = 2000.0;
	tmp.beta = 4.35113e-10;
	tmp.skins.push_back(0.0);
	tmp.radiuses_eff.push_back(props->r_w);
	props->props_sk.push_back( tmp );

	tmp.cellsNum_z = 1;
	tmp.m = 0.1;
	tmp.p_init = tmp.p_out = 180.0 * 1.0e+5;
	tmp.h1 = 1540.0;
	tmp.h2 = 1560.0;
	tmp.height = 20.0;
	tmp.perm_r = 100.0;
	tmp.perm_z = 10.0;
	tmp.dens_stc = 2000.0;
	tmp.beta = 4.35113e-10;
	tmp.skins.push_back(0.0);
	tmp.radiuses_eff.push_back(props->r_w);
	props->props_sk.push_back( tmp );

	tmp.cellsNum_z = 1;
	tmp.m = 0.1;
	tmp.p_init = tmp.p_out = 0.0 * 1.0e+5;
	tmp.h1 = 1560.0;
	tmp.h2 = 1570.0;
	tmp.height = 10.0;
	tmp.perm_r = 0.0;
	tmp.perm_z = 0.0;
	tmp.dens_stc = 2000.0;
	tmp.beta = 4.35113e-10;
	tmp.skins.push_back(0.0);
	tmp.radiuses_eff.push_back(props->r_w);
	props->props_sk.push_back( tmp );

	tmp.cellsNum_z = 1;
	tmp.m = 0.1;
	tmp.p_init = tmp.p_out = 200.0 * 1.0e+5;
	tmp.h1 = 1570.0;
	tmp.h2 = 1600.0;
	tmp.height = 30.0;
	tmp.perm_r = 100.0;
	tmp.perm_z = 10.0;
	tmp.dens_stc = 2000.0;
	tmp.beta = 4.35113e-10;
	tmp.skins.push_back(0.0);
	tmp.radiuses_eff.push_back(props->r_w);
	props->props_sk.push_back( tmp );

	props->depth_point = 1500.0;

	props->props_oil.visc = 0.62807;
	props->props_oil.b_bore = 1.09754;
	props->props_oil.dens_stc = 800.026;
	props->props_oil.beta = 1.24703e-09;

	return props;
}*/

/*oil_rz::Properties* getProps()
{
	oil_rz::Properties* props = new oil_rz::Properties();

	props->cellsNum_r = 100;
	props->cellsNum_z = 3;

	props->timePeriods.push_back(10.0 * 365.0 * 86400.0);
	
	props->leftBoundIsRate = true;
	props->rightBoundIsPres = false;
	props->rates.push_back(5.0);

	props->ht = 1000.0;
	props->ht_min = 100.0;
	props->ht_max  = 5000000.0;
	
	props->alpha = 7200.0;

	props->perfIntervals.push_back( make_pair(1, 3) );

	props->r_w = 0.1;
	props->r_e = 2500.0;

	oil_rz::Skeleton_Props tmp;
	tmp.cellsNum_z = 3;
	tmp.m = 0.1;
	tmp.p_init = tmp.p_out = 160.0 * 1.0e+5;
	tmp.h1 = 1500.0;
	tmp.h2 = 1503.0;
	tmp.height = 3.0;
	tmp.perm_r = 100.0;
	tmp.perm_z = 0.0;
	tmp.dens_stc = 2000.0;
	tmp.beta = 4.35113e-10;
	tmp.skins.push_back(0.0);
	tmp.radiuses_eff.push_back(props->r_w);
	props->props_sk.push_back( tmp );

	props->depth_point = 1500.0;

	props->props_oil.visc = 0.62807;
	props->props_oil.b_bore = 1.09754;
	props->props_oil.dens_stc = 800.026;
	props->props_oil.beta = 1.24703e-09;

	return props;
}*/

/*oil1D::Properties* getProps()
{
	oil1D::Properties* props = new oil1D::Properties();

	props->cellsNum_r = 100;

	props->timePeriods.push_back(1296000.0);
	props->timePeriods.push_back(2000000.0);
	
	props->rates.push_back(30.0);
	props->rates.push_back(0.0);
	 
	props->skins.push_back(0.0); 
	props->skins.push_back(0.0); 

	props->radius.push_back(0.1524);
	props->radius.push_back(0.1524);

	props->ht = 10.0;
	props->ht_min = 1.0;
	props->ht_max  = 100000.0;
	
	props->alpha = 7200.0;

	props->perfIntervals.push_back( make_pair(0, 0) );

	props->r_w = 0.1524;
	props->r_e = 3000.0;

	props->p_init = 160.0 * 100000.0;
	props->height = 1.0;

	props->perm = 50.0;

	props->m = 0.01;
	props->visc_oil = 0.25137;
	props->dens_oil_stc = 855.0;
	props->dens_sk_stc = 2000.0;
	props->beta_oil = 1.282*1.E-9;
	props->beta_sk = 6.41*1.E-10;
	props->b_oil_bore = 1.56;
	return props;
}

*/

gas1D::Properties* getProps()
{
	gas1D::Properties* props = new gas1D::Properties();

	props->cellsNum_r = 100;

	props->timePeriods.push_back(10.0 * 365.0 * 86400.0);
	
	props->leftBoundIsRate = false;
	props->rightBoundIsPres = true;
	//props->rates.push_back(5000.0);
	props->pwf.push_back(140.0 * 100000.0);

	props->ht = 20000.0;
	props->ht_min = 10000.0;
	props->ht_max  = 10000000.0;
	
	props->alpha = 7200.0;

	props->perfIntervals.push_back( make_pair(0, 0) );

	props->r_w = 0.1;
	props->r_e = 2500.0;

	gas1D::Skeleton_Props tmp;
	tmp.cellsNum_z = 1;
	tmp.m = 0.1;
	tmp.p_init = tmp.p_out = 160.0 * 1.0e+5;
	tmp.h1 = 1500.0;
	tmp.h2 = 1503.0;
	tmp.height = 3.0;
	tmp.perm_r = 1.0;
	tmp.perm_z = 0.0;
	tmp.dens_stc = 2000.0;
	tmp.beta = 4.35113e-10;
	tmp.skins.push_back(0.0);
	tmp.radiuses_eff.push_back(props->r_w);
	props->props_sk.push_back( tmp );

	props->depth_point = 1500.0;

	props->props_gas.visc = 0.02833;
	props->props_gas.b_bore = 1.0 / 160.0;

	// Defining relative permeabilities
	setDataFromFile(props->z_factor, "props/z.txt");

	return props;
}

int main(int argc, char** argv)
{
	/*gasOil_rz_NIT::Properties* props = getProps();
	Scene<gasOil_rz_NIT::GasOil_RZ_NIT, gasOil_rz_NIT::GasOil2DNITSolver, gasOil_rz_NIT::Properties> scene;
	scene.load(*props);
	scene.setSnapshotterType("VTK");
	scene.start();*/

	/*oil1D::Properties* props = getProps();
	Scene<oil1D::Oil1D, oil1D::Oil1DSolver, oil1D::Properties> scene;
	scene.load(*props);
	scene.setSnapshotterType("VTK");
	scene.start();*/

	/*gasOil_rz::Properties* props = getProps();
	Scene<gasOil_rz::GasOil_RZ, gasOil_rz::GasOil2DSolver, gasOil_rz::Properties> scene;
	scene.load(*props);
	scene.setSnapshotterType("VTK");
	scene.start();*/

	gas1D::Properties* props = getProps();
	const double p_bhp [7] = {80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0};

	for(int i = 0; i < 7; i++)
	{
		props->pwf.clear();
		props->pwf.push_back( p_bhp[i] * 100000.0 );
		Scene<gas1D::Gas1D, gas1D::Gas1DSolver, gas1D::Properties> scene;	
		scene.load(*props, i);
		scene.setSnapshotterType("VTK");
		scene.start();
	}


	/*oil_rz::Properties* props = getProps();
	Scene<oil_rz::Oil_RZ, oil_rz::OilRZSolver, oil_rz::Properties> scene;
	scene.load(*props);
	scene.setSnapshotterType("VTK");
	scene.start();*/

	//testing::InitGoogleTest(&argc, argv);
	//int res = RUN_ALL_TESTS();

	return 0;
}