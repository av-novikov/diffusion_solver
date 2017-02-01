#ifndef ACID2D_PROPERTIES_HPP_
#define ACID2D_PROPERTIES_HPP_

#include <array>

#include "model/Basic2d/Properties.hpp"
#include "model/Acid/2d/Reactions.hpp"

namespace acid2d
{
	// ADOLC stencil ids
	const int mid = basic2d::mid;
	const int left = basic2d::left;
	const int right = basic2d::right;
	const int vertical = basic2d::vertical;

	struct SolidComponent : Component
	{
		double beta;
		inline adouble getDensity(adouble p) const
		{
			return (adouble)(rho_stc)* ((adouble)(1.0) + (adouble)(beta)* p);
		};
	};
	struct LiquidComponent : Component
	{
		double beta;
		inline adouble getDensity(adouble p) const
		{
			return (adouble)(rho_stc)* ((adouble)(1.0) + (adouble)(beta)* p);
		};
	};
	struct GasComponent : Component
	{
		double z;
		Interpolate* z_table;
		inline adouble getDensity(adouble p) const
		{
			return p * (adouble)(mol_weight / (z * R * T));
		};
	};

	struct Skeleton_Props : public basic2d::Skeleton_Props
	{
		SolidComponent& cur_mineral;

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
	struct Liquid_Props : public basic2d::Liquid_Props
	{
		LiquidComponent& acid;
		LiquidComponent& salt;
		LiquidComponent& water;

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
	struct Gas_Props : public basic2d::Gas_Props
	{
		GasComponent& co2;

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