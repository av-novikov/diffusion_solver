#ifndef GASOILRZ_PROPERTIES_HPP_
#define GASOILRZ_PROPERTIES_HPP_

#include "model/cells/Variables.hpp"
#include "model/cells/CylCell2D.h"
#include "util/utils.h"

#include "adolc/adouble.h"
#include "adolc/taping.h"

namespace gasOil_rz
{
	const int mid1 = 1;
	const int mid2 = 2;
	const int left = 3;
	const int right = 4;

	const int neighbors = 4;
	const int mid_stenc = 3;
	const int b_stenc = 2;
	const int schemeVarNum = mid_stenc * (neighbors + 1);
	const int boundVarNum = b_stenc * (neighbors + 1);

	struct Skeleton_Props
	{
		// Porosity in STC
		double m;
		inline double getPoro(double p) const
		{
			return m * (1.0 + beta * (p /*- cell.props->p_init*/));
		};
		inline adouble getPoro(adouble p) const
		{
			return m * (1.0 + beta * (p /*- cell.props->p_init*/));
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

		double p_out;
		double p_init;
		double p_bub;
		double s_init;
	};

	struct Oil_Props
	{
		// Viscosity [cP]
		double visc;
		// Density of fluid in STC [kg/m3]
		double dens_stc;
		// Volume factor for well bore
		double b_bore;
		// Compessibility [1/Pa]
		double beta;
		// Relative fluid permeability
		Interpolate* kr;
		inline double getKr(double sat_oil) const
		{
			if (sat_oil > 1.0)
				return 1.0;
			else
				return kr->Solve(sat_oil);
		};
		inline adouble getKr(adouble sat_oil) const
		{
			if (sat_oil > 1.0)
				return 1.0;
			else
				return kr->Solve(sat_oil);
		};

		// Fluid volume factor
		Interpolate* b;
		inline double getB(double p, double p_bub, bool SATUR) const
		{
			if (SATUR)
				return b->Solve(p);
			else
				return b->Solve(p_bub) * (1.0 + beta * (p_bub - p));
		};
		inline adouble getB(adouble p, adouble p_bub, bool SATUR) const
		{
			if (SATUR)
				return b->Solve(p);
			else
				return b->Solve(p_bub) * (1.0 + beta * (p_bub - p));
		};
		inline double getBoreB(double p, double p_bub, bool SATUR) const
		{
			//return props_oil.b_bore;
			return getB(p, p_bub, SATUR);
		};
		inline adouble getBoreB(adouble p, adouble p_bub, bool SATUR) const
		{
			//return props_oil.b_bore;
			return getB(p, p_bub, SATUR);
		};

		// Gas-oil ratio
		Interpolate* Rs;
		inline double getRs(double p, double p_bub, bool SATUR) const
		{
			if (SATUR)
				return Rs->Solve(p);
			else
				return Rs->Solve(p_bub);
		};
		inline adouble getRs(adouble p, adouble p_bub, bool SATUR) const
		{
			if (SATUR)
				return Rs->Solve(p);
			else
				return Rs->Solve(p_bub);
		};
	};
	
	struct Gas_Props
	{
		// Viscosity [cP]
		double visc;
		// Density of fluid in STC [kg/m3]
		double dens_stc;
		// Volume factor for well bore
		double b_bore;
		// Compessibility [1/Pa]
		double beta;
		// Relative fluid permeability
		Interpolate* kr;
		inline double getKr(double sat_oil) const
		{
			if (sat_oil > 1.0)
				return 0.0;
			else
				return kr->Solve(sat_oil);
		};
		inline adouble getKr(adouble sat_oil) const
		{
			if (sat_oil > 1.0)
				return 0.0;
			else
				return kr->Solve(sat_oil);
		};
		// Fluid volume factor
		Interpolate* b;
		inline double getB(double p) const
		{
			return b->Solve(p);
		};
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

		// Data set (pressure, gas content in oil) ([Pa], [m3/m3])
		std::vector< std::pair<double, double> > Rs;
	};
}

#endif /* GASOILRZ_PROPERTIES_HPP_ */