#ifndef BLACKOILNITELLIPTIC_PROPERTIES_HPP_
#define BLACKOILNITELLIPTIC_PROPERTIES_HPP_

#include "util/utils.h"
#include "model/Acid/2d/Reactions.hpp"

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
			inline adouble getPoro(adouble p) const
			{
				return (adouble)(m)* ((adouble)(1.0) + (adouble)(beta)* (p /*- cell.props->p_init*/));
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
			double t_init;

			// Mass heat capacity [J/kg/K]
			double c;
			// Thermal conductivity coefficient [W/m/K]
			double lambda_r;
			// Thermal conductivity coefficient [W/m/K]
			double lambda_z;
		};
		struct Oil_Props
		{
			acid2d::Component oil;
			acid2d::Component gas;

			double p_sat;

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
				return (adouble)(oil.rho_stc)* ((adouble)(1.0) + (adouble)(beta)* p);
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

			double depth_point;

			// Data set (pressure, oil viscosity) ([Pa], [cP])
			std::vector< std::pair<double, double> > visc_oil;
			// Data set (pressure, gas viscosity) ([Pa], [cP])
			std::vector< std::pair<double, double> > visc_gas;
		};
	};

#endif /* BLACKOILNITELLIPTIC_PROPERTIES_HPP_ */