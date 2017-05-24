#ifndef BLACKOILNITELLIPTIC_PROPERTIES_HPP_
#define BLACKOILNITELLIPTIC_PROPERTIES_HPP_

#include "util/utils.h"

#include "adolc/adouble.h"
#include "adolc/taping.h"

namespace blackoilnit_elliptic
	{
		const int mid = 1;
		const int left = 2;
		const int right = 3;
		const int vertical = 4;

		const int tmid = 5;
		const int tleft = 6;
		const int tright = 7;
		const int tvertical = 8;

		struct Skeleton_Props
		{
			bool isWellHere = false;

			// Porosity in STC
			double m;
			double p_ref;
			inline adouble getPoro(adouble p) const
			{
				return (adouble)(m)* ((adouble)(1.0) + (adouble)(beta)* (p - p_ref));
			};
			// Density of skeleton matter in STC [kg/m3]
			double dens_stc;
			inline adouble getDensity(adouble p) const
			{
				return (adouble)(dens_stc);
			};
			// Compessibility [1/Pa]
			double beta;
			// Permeability along radial direction [mD]
			double perm_mu;
			// Permeability along vertical direction [mD]
			double perm_z;
			inline double getPerm_mu(const double mu) const
			{
				return (mu > radius_eff_mu ? perm_mu : perm_eff_mu);
			};
			inline double getPerm_z(const double z) const
			{
				return (z > radius_eff_z ? perm_z : perm_eff_z);
			};

			// Permeability of colmatage zone [mD]
			std::vector<double> perms_eff;
			// Radius of colmatage zone [m]
			std::vector<double> radiuses_eff;
			// Vector of skins
			std::vector<double> skins;
			double perm_eff_mu, perm_eff_z;
			double radius_eff_mu, radius_eff_z, radius_eff;
			double skin;

			// Top and bottom depth of perforation
			double h1, h2;
			double h_well;
			// Height of formation [m]
			double height;

			int start_z;
			int cellsNum_z;

			double p_out;
			double p_init;
			double p_sat;
			double so_init;
			double sw_init;
			double t_init;
			// Connate saturations
			double s_wc, s_oc, s_gc;

			// Mass heat capacity [J/kg/K]
			double c;
			// Thermal conductivity coefficient [W/m/K]
			double lambda_r;
			// Thermal conductivity coefficient [W/m/K]
			double lambda_z;
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
			static double dens_stc;
			// Volume factor for well bore
			double b_bore;
			// Fluid volume factor
			Interpolate* b;
			inline adouble getB(adouble p) const
			{
				return b->Solve(p);
			};
			inline adouble getDensity(adouble p) const
			{
				return dens_stc / getB(p);
			};
			// Relative fluid permeability
			Interpolate* kr;
			inline adouble getKr(adouble s_w, adouble s_o, const Skeleton_Props* props) const
			{
				adouble isAboveZero = (1.0 - s_o - s_w - props->s_gc > 0.0) ? true : false;
				adouble isAboveCritical = (1.0 - s_o - s_w > 1.0 - props->s_wc - props->s_oc) ? true : false;
				adouble tmp;
				condassign(tmp, isAboveZero, 0.5 * pow(((adouble)(1.0 - props->s_gc) - s_w - s_o) / (adouble)(1.0 - props->s_wc - props->s_oc - props->s_gc), 3.0), (adouble)0.0);
				condassign(tmp, isAboveCritical, (adouble)0.5);
				return tmp;
				//return kr->Solve(s_w);
			};

			// Mass heat capacity [J/kg/K]
			double c;
			// Thermal conductivity coefficient [W/m/K]
			double lambda;
			// Joule-thompson coefficient [K/Pa]
			double jt;
			// Adiabatic coefficient [K/Pa]
			double ad;
		};
		struct Water_Props
		{
			double p_ref;
			// Density of fluid in STC [kg/m3]
			static double dens_stc;
			// Viscosity [cP]
			double visc;
			Interpolate* visc_table;
			inline adouble getViscosity(const adouble p) const
			{
				return (adouble)(visc);
			};
			// Density of fluid in STC [kg/m3]
			inline adouble getDensity(adouble p) const
			{
				return dens_stc / getB(p);
			};
			// Fluid volume factor
			Interpolate* b;
			inline adouble getB(adouble p) const
			{
				//adouble tmp;
				//condassign(tmp, SATUR, b->Solve(p), b->Solve(p_bub));
				return exp((adouble)beta * (p - p_ref));
			};
			// Relative fluid permeability
			Interpolate* kr;
			inline adouble getKr(adouble s_w, adouble s_o, const Skeleton_Props* props) const
			{
				adouble isAboveZero = (s_w - props->s_wc > 0.0) ? true : false;
				adouble isAboveCritical = (s_w > 1.0 - props->s_oc - props->s_gc) ? true : false;
				adouble tmp;
				condassign(tmp, isAboveZero, pow((s_w - (adouble)props->s_wc) / (adouble)(1.0 - props->s_wc - props->s_oc - props->s_gc), 3.0), (adouble)0.0);
				condassign(tmp, isAboveCritical, (adouble)1.0);
				return tmp;
				//return kr->Solve(s_w);
			};
			// Compessibility [1/Pa]
			double beta;

			// Mass heat capacity [J/kg/K]
			double c;
			// Thermal conductivity coefficient [W/m/K]
			double lambda;
			// Joule-thompson coefficient [K/Pa]
			double jt;
			// Adiabatic coefficient [K/Pa]
			double ad;
		};
		struct Oil_Props
		{
			double beta;
			double p_ref;
			// Density of fluid in STC [kg/m3]
			static double dens_stc;
			// Viscosity [cP]
			double visc;
			// Fluid volume factor
			Interpolate* b;
			inline adouble getB(adouble p, adouble p_bub, adouble SATUR) const
			{
				adouble tmp;
				condassign(tmp, SATUR, b->Solve(p), b->Solve(p_bub) * exp((adouble)beta * (p_bub - p)));
				return tmp;
			};
			// Gas-oil ratio
			Interpolate* Rs;
			inline adouble getRs(adouble p, adouble p_bub, adouble SATUR) const
			{
				adouble tmp;
				condassign(tmp, SATUR, Rs->Solve(p), Rs->Solve(p_bub));
				return tmp;
			};
			// Density of fluid in STC [kg/m3]
			inline adouble getDensity(adouble p, adouble p_bub, adouble SATUR) const
			{
				return (dens_stc + getRs(p, p_bub, SATUR) * Gas_Props::dens_stc) / getB(p, p_bub, SATUR);
			};
			inline adouble getViscosity(const adouble p) const
			{
				return (adouble)(visc);
			};
			// Relative fluid permeability
			Interpolate* kr;
			inline adouble getKr(adouble s_w, adouble s_o, const Skeleton_Props* props) const
			{
				adouble isAboveZero = (s_o - props->s_oc > 0.0) ? true : false;
				adouble isAboveCritical = (s_o > 1.0 - props->s_wc - props->s_gc) ? true : false;
				adouble tmp;
				condassign(tmp, isAboveZero, pow((s_o - (adouble)props->s_oc) / (adouble)(1.0 - props->s_wc - props->s_oc - props->s_gc), 3.0), (adouble)0.0);
				condassign(tmp, isAboveCritical, (adouble)1.0);
				return tmp;
				//return kr->Solve(s_w);
			};

			// Mass heat capacity [J/kg/K]
			double c;
			// Thermal conductivity coefficient [W/m/K]
			double lambda;
			// Joule-thompson coefficient [K/Pa]
			double jt;
			// Adiabatic coefficient [K/Pa]
			double ad;
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
			// Length of well [m]
			double l;
			// Radius of formation [m]
			double r_e;

			int cellsNum_mu;
			int cellsNum_nu;
			int cellsNum_z;

			std::vector<Skeleton_Props> props_sk;
			Oil_Props props_oil;
			Water_Props props_wat;
			Gas_Props props_gas;
			double L;

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
	};

#endif /* BLACKOILNITELLIPTIC_PROPERTIES_HPP_ */