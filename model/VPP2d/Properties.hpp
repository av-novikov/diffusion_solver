#ifndef VPP2D_PROPERTIES_HPP_
#define VPP2D_PROPERTIES_HPP_

#include "adolc/adouble.h"
#include "adolc/taping.h"

namespace vpp2d
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

		// Connate saturations
		double s_wc;
		double s_oc;

		// Initial values
		double p_init;
		double s_init;
		double c_init;
	};

	struct Water_Props
	{
		Interpolate* a;

		// Viscosity [cP]
		double visc;
		Interpolate* visc_table;
		inline adouble getViscosity(const adouble p) const
		{
			return visc;
		};

		// Fluid volume factor
		Interpolate* b;
		double b_bore;
		// Density of fluid in STC [kg/m3]
		double dens_stc;
		// Compessibility [1/Pa]
		double beta;
		inline adouble getB(const adouble p) const
		{
			//return b->Solve(p);
			return b_bore * (1 - beta * (p - p_ref));
		};
		inline adouble getBoreB(const adouble p) const
		{
			return getB(p);
		};
		// Reference pressure [Pa]
		double p_ref;
		inline adouble getDensity(const adouble p) const
		{
			return dens_stc * (1.0 + beta * (p - p_ref));
		};

		// Relative fluid permeability
		Interpolate* kr;
		inline adouble getKr(const adouble s_l) const
		{
			if (s_l > 1.0)
				return kr->Solve(1.0);
			else
				return kr->Solve(s_l);
		};
		
	};
	struct Oil_Props
	{
		// Viscosity [cP]
		double visc;
		Interpolate* visc_table;
		inline adouble getViscosity(const adouble p) const
		{
			return visc;
		};

		// Fluid volume factor
		Interpolate* b;
		double b_bore;
		inline adouble getB(const adouble p) const
		{
			//return b->Solve(p);
			return b_bore * (1 - beta * (p - p_ref));
		};
		inline adouble getBoreB(const adouble p) const
		{
			return getB(p);
		};
		// Density of fluid in STC [kg/m3]
		double dens_stc;
		// Compessibility [1/Pa]
		double beta;
		// Reference pressure [Pa]
		double p_ref;
		inline adouble getDensity(const adouble p) const
		{
			return dens_stc * (1.0 + beta * (p - p_ref));
		};

		// Relative fluid permeability
		Interpolate* kr;
		inline adouble getKr(const adouble s_l) const
		{
			if (s_l > 1.0)
				return 0.0;
			else
				return kr->Solve(s_l);
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
		// Vector of concentration of addition in water
		std::vector<double> c;

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
		Water_Props props_w;
		Oil_Props props_o;

		double depth_point;

		// Data set (saturation, relative oil permeability)
		std::vector< std::pair<double, double> > kr_o;
		// Data set (saturation, relative gas permeability)
		std::vector< std::pair<double, double> > kr_w;

		// Data set (pressure, oil volume factor) ([Pa], [m3/m3])
		std::vector< std::pair<double, double> > B_o;
		// Data set (pressure, gas volume factor) ([Pa], [m3/m3])
		std::vector< std::pair<double, double> > B_w;

		// Data set (pressure, oil volume factor) ([Pa], [cP])
		std::vector< std::pair<double, double> > visc_o;
		// Data set (pressure, gas volume factor) ([Pa], [cP])
		std::vector< std::pair<double, double> > visc_w;

		// Data set (pressure, gas volume factor) ([Pa], [1])
		std::vector< std::pair<double, double> > a;
	};

};



#endif /* VPP2D_PROPERTIES_HPP_ */