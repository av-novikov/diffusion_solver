#ifndef ACID2D_PROPERTIES_HPP_
#define ACID2D_PROPERTIES_HPP_

#include <array>

#include "model/Basic2d/Properties.hpp"

#include "adolc/adouble.h"
#include "adolc/taping.h"

namespace acid2d
{
	// ADOLC stencil ids
	const int mid = basic2d::mid;
	const int left = basic2d::left;
	const int right = basic2d::right;
	const int vertical = basic2d::vertical;

	struct Component
	{
		static double R;
		static double T;

		double mol_weight;
		double rho_stc;

		double stehiometr_mult;

	};
	struct Mineral
	{
		double mol_weight;
		double rho_stc;
		double beta;
		inline adouble getDensity(adouble p) const
		{
			return (adouble)(rho_stc)* ((adouble)(1.0) + (adouble)(beta)* p);
		};
	};
	struct Skeleton_Props : public basic2d::Skeleton_Props
	{
		int min_id;
		static const Mineral calcite;
		static const Mineral dolomite;
		Mineral& cur_mineral;

		// Initial values
		double m_init;
		double p_init;
		double s_init;
		double Ya_init;
		double Ys_init;

		/*Skeleton_Props(int _min_idx)
		{ 
			cur_mineral = 
		};
		Skeleton_Props(const Skeleton_Props& props) : basic2d::Skeleton_Props(props), mineral(props.mineral) {};
		Liquid_Props& operator=(const Liquid_Props& props)
		{
			basic2d::Liquid_Props::operator=(props);
			components = props.components;
			return *this;
		}*/
	};
	struct LiquidComponent
	{
		double rho_stc;
		double beta;
		inline adouble getDensity(adouble p) const
		{
			return (adouble)(rho_stc)* ((adouble)(1.0) + (adouble)(beta)* p);
		};
	};
	struct Liquid_Props : public basic2d::Liquid_Props
	{
		std::array<LiquidComponent, 3> components;
		const LiquidComponent& acid = components[0];
		const LiquidComponent& salt = components[1];
		const LiquidComponent& water = components[2];

		inline adouble getDensity(adouble p, adouble Ya, adouble Ys) const
		{
			return	acid.getDensity(p) * Ya + 
					salt.getDensity(p) * Ys + 
					water.getDensity(p) * ((adouble)(1.0) - Ya - Ys);
		};

		Liquid_Props() {};
		Liquid_Props(const Liquid_Props& props) : basic2d::Liquid_Props(props), components(props.components) {};
		Liquid_Props& operator=(const Liquid_Props& props) 
		{
			basic2d::Liquid_Props::operator=(props);
			components = props.components;
			return *this;
		}
	};
	struct GasComponent
	{
		static double R;
		static double T;
		double rho_stc;
		double mol_weight;
		double z;
		Interpolate* z_table;
		inline adouble getDensity(adouble p) const
		{
			return p * (adouble)(mol_weight / (z * R * T));
		};
	};
	struct Gas_Props : public basic2d::Gas_Props
	{
		std::array<GasComponent, 1> components;
		const GasComponent& co2 = components[0];

		inline adouble getDensity(adouble p) const
		{
			return co2.getDensity(p);
		};

		Gas_Props() {};
		Gas_Props(const Gas_Props& props) : basic2d::Gas_Props(props), components(props.components) {};
		Gas_Props& operator=(const Gas_Props& props)
		{
			basic2d::Gas_Props::operator=(props);
			components = props.components;
			return *this;
		}
	};
	struct Properties : public basic2d::Properties
	{
		std::vector<Skeleton_Props> props_sk;
		Liquid_Props props_l;
		Gas_Props props_g;

		// Data set (saturation, relative oil permeability)
		std::vector< std::pair<double, double> > kr_l;
		// Data set (saturation, relative gas permeability)
		std::vector< std::pair<double, double> > kr_g;
	};
};

#endif /* ACID2D_PROPERTIES_HPP_ */