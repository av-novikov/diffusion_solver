#ifndef GASOILRZ_PROPERTIES_HPP_
#define GASOILRZ_PROPERTIES_HPP_

#include "util/utils.h"

#include "adolc/adouble.h"
#include "adolc/taping.h"

namespace gasOil_rz
{
	const int mid = 1;
	const int left = 2;
	const int right = 3;
	const int vertical = 4;

	struct Skeleton_Props
	{
		// Porosity in STC
		double m;
		inline adouble getPoro(adouble p) const
		{
			return (adouble)(m)* ((adouble)(1.0) + (adouble)(beta)* (p /*- cell.props->p_init*/));
		};
		// Density of skeleton matter in STC [kg/m3]
		double dens_stc;
		// Compessibility [1/Pa]
		double beta;
		// Permeability along radial direction [mD]
		double perm_r;
		// Permeability along vertical direction [mD]
		double perm_z;
		inline double getPerm_r(const double r) const
		{
			return (r > radius_eff ? perm_r : perm_eff);
		};

		// Permeability of colmatage zone [mD]
		std::vector<double> perms_eff;
		// Radius of colmatage zone [m]
		std::vector<double> radiuses_eff;
		// Vector of skins
		std::vector<double> skins;
		double perm_eff;
		double radius_eff;
		double skin;

		// Top and bottom depth of perforation
		double h1, h2;
		// Height of formation [m]
		double height;

		int cellsNum_z;

		// Connate saturations
		double s_wc;
		double s_oc;

		double p_out;
		double p_init;
		double p_bub;
		double s_init;
	};

	struct Oil_Props
	{
		double p_sat;

		// Viscosity [cP]
		double visc;
		Interpolate* visc_table;
		inline adouble getViscosity(const adouble p) const
		{
			return (adouble)(visc);
		};
		// Density of fluid in STC [kg/m3]
		double dens_stc;
		// Volume factor for well bore
		double b_bore;
		// Compessibility [1/Pa]
		double beta;
		// Relative fluid permeability
		Interpolate* kr;
		inline adouble getKr(adouble sat_oil) const
		{
			return kr->Solve(sat_oil);
		};

		// Fluid volume factor
		Interpolate* b;
		inline adouble getB(adouble p, adouble p_bub, adouble SATUR) const
		{
			adouble tmp;
			//adouble isAboveSat = (p > (adouble)p_sat) ? 1.0 : 0.0;
			condassign(tmp, SATUR, b->Solve(p), b->Solve(p_bub) * (1.0 + beta * (p_bub - p)));
			//condassign(tmp, isAboveSat, (adouble)(b->Solve(p_sat) * (1.0 + beta * (p_sat - p))));
			return tmp;
		};
		inline adouble getBoreB(adouble p, adouble p_bub, adouble SATUR) const
		{
			//return props_oil.b_bore;
			return getB(p, p_bub, SATUR);
		};

		// Gas-oil ratio
		Interpolate* Rs;
		inline adouble getRs(adouble p, adouble p_bub, adouble SATUR) const
		{
			adouble tmp;
			//adouble isAboveSat = (p > (adouble)p_sat) ? 1.0 : 0;
			condassign(tmp, SATUR, Rs->Solve(p), Rs->Solve(p_bub));
			//condassign(tmp, isAboveSat, (adouble)(Rs->Solve(p_sat)));
			return tmp;
		};
	};
	
	struct Gas_Props
	{
		// Viscosity [cP]
		double visc;
		Interpolate* visc_table;
		inline adouble getViscosity(const adouble p) const
		{
			return (adouble)(visc);
		};
		// Density of fluid in STC [kg/m3]
		double dens_stc;
		// Volume factor for well bore
		double b_bore;
		// Compessibility [1/Pa]
		double beta;
		// Relative fluid permeability
		Interpolate* kr;
		inline adouble getKr(adouble sat_oil) const
		{
			return kr->Solve(sat_oil);
		};
		// Fluid volume factor
		Interpolate* b;
		inline adouble getB(adouble p) const
		{
			return b->Solve(p);
		};
	};

	struct Properties
	{
		// Vector of start times of periods [sec]
		std::vector<double> timePeriods;
		// Vector of rates [m3/day]
		std::vector<double> rates;
		// Vector of BHPs [Pa]
		std::vector<double> pwf;

		// If left boundary condition would be 2nd type
		bool leftBoundIsRate;
		// If right boundary condition would be 1st type
		bool rightBoundIsPres;

		// Perforated intervals
		std::vector<std::pair<int, int> > perfIntervals;
		// Time step limits
		// Initial time step [sec]
		double ht;
		// Minimal time step [sec]
		double ht_min;
		// Maximum time step [sec]
		double ht_max;
		// During the time flow rate decreases 'e' times in well test [sec] 
		double alpha;

		// Inner radius of well [m]
		double r_w;
		// Radius of formation [m]
		double r_e;

		// Number of cells in radial direction
		int cellsNum_r;
		// Number of cells in vertical direction
		int cellsNum_z;

		std::vector<Skeleton_Props> props_sk;
		Oil_Props props_oil;
		Gas_Props props_gas;

		double depth_point;

		// Data set (saturation, relative oil permeability)
		std::vector< std::pair<double, double> > kr_oil;
		// Data set (saturation, relative gas permeability)
		std::vector< std::pair<double, double> > kr_gas;

		// Data set (pressure, oil volume factor) ([Pa], [m3/m3])
		std::vector< std::pair<double, double> > B_oil;
		// Data set (pressure, gas volume factor) ([Pa], [m3/m3])
		std::vector< std::pair<double, double> > B_gas;

		// Data set (pressure, oil viscosity) ([Pa], [cP])
		std::vector< std::pair<double, double> > visc_oil;
		// Data set (pressure, gas viscosity) ([Pa], [cP])
		std::vector< std::pair<double, double> > visc_gas;

		// Data set (pressure, gas content in oil) ([Pa], [m3/m3])
		std::vector< std::pair<double, double> > Rs;
	};
}

#endif /* GASOILRZ_PROPERTIES_HPP_ */