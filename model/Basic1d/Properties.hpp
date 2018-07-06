#ifndef BASIC1D_PROPERTIES_HPP_
#define BASIC1D_PROPERTIES_HPP_

#include "adolc/adouble.h"
#include "adolc/taping.h"
#include "util/utils.h"

namespace basic1d
{
	// ADOLC stencil ids
	const int mid = 1;
	const int left = 2;
	const int right = 3;
		
	struct Skeleton_Props
	{
		double beta;
		double dens_stc;
		double perm;

		double m_init;
		double p_ref;
		inline adouble getPoro(adouble p) const
		{
			return (adouble)(m_init)* ((adouble)(1.0) + (adouble)(beta)* (p - p_ref));
		};
		inline double getPerm(const double r) const
		{
			return (r > radius_eff ? perm : perm_eff);
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

		// Connate saturations
		double s_wc, s_oc, s_gc;
		
		double p_out;

		double hx, hz;
		int cellsNum;
	};
	struct Liquid_Props
	{
		double p_ref;

		// Viscosity [cP]
		double visc;
		Interpolate* visc_table;
		// Density of fluid in STC [kg/m3]
		double dens_stc;
		// Volume factor for well bore
		double b_bore;
		// Compessibility [1/Pa]
		double beta;

		// Gas-oil ratio
		Interpolate* Rs;
		inline adouble getRs(adouble p, adouble p_bub, adouble SATUR) const
		{
			adouble tmp;
			condassign(tmp, SATUR, Rs->Solve(p), Rs->Solve(p_bub));
			return tmp;
		};
	};
	struct Gas_Props
	{
		// Viscosity [cP]
		double visc;
		Interpolate* visc_table;
		// Density of fluid in STC [kg/m3]
		double dens_stc;
		// Volume factor for well bore
		double b_bore;
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
		size_t cellsNum_x;
		double depth_point;

		std::vector< std::pair<double, double> > kr_wat;
		std::vector< std::pair<double, double> > kr_oil;
		std::vector< std::pair<double, double> > kr_gas;
		
		std::vector< std::pair<double, double> > B_wat;
		std::vector< std::pair<double, double> > B_oil;
		std::vector< std::pair<double, double> > B_gas;

		std::vector< std::pair<double, double> > visc_wat;
		std::vector< std::pair<double, double> > visc_oil;
		std::vector< std::pair<double, double> > visc_gas;

		std::vector< std::pair<double, double> > Rs;
	};
}

#endif /* BASIC1D_PROPERTIES_HPP_ */
