#ifndef ACID2D_REACTIONS_HPP_
#define ACID2D_REACTIONS_HPP_

#define KKAL_2_J 4186.8

#include <array>

#include "adolc/adouble.h"
#include "adolc/taping.h"

#include "model/Basic2d/Properties.hpp"

namespace acid2d
{
	struct Component
	{
		static double p_std;
		static double R;
		static double T;

		double mol_weight;
		double rho_stc;
	};
	struct SolidComponent : Component
	{
		inline adouble getMolarDensity(adouble p) const
		{
			return rho_stc / mol_weight;
		};
		double beta;
		inline adouble getDensity(adouble p) const
		{
			return rho_stc;
		};
	};
	struct LiquidComponent : Component
	{
		inline adouble getMolarDensity(adouble p) const
		{
			return getDensity(p) / mol_weight;
		};
		double beta;
		inline adouble getDensity(adouble p) const
		{
			return (adouble)(rho_stc)* ((adouble)(1.0) + (adouble)(beta)* (p - p_std));
		};
	};
	struct GasComponent : Component
	{
		inline adouble getMolarDensity(adouble p) const
		{
			return getDensity(p) / mol_weight;
		};
		double z;
		Interpolate* z_table;
		inline adouble getDensity(adouble p) const
		{
			//return p * (adouble)(mol_weight / (z * R * T));
			return rho_stc * p / p_std;
		};
	};

	// Laboratory
	inline SolidComponent getCaCO3()
	{
		SolidComponent comp;
		comp.mol_weight = 100.0;
		comp.rho_stc = 2710.0;
		return comp;
	};
	inline LiquidComponent getHCl()
	{
		LiquidComponent comp;
		comp.mol_weight = 36.0;
		return comp;
	};
	inline SolidComponent getCaCl2()
	{
		SolidComponent comp;
		comp.mol_weight = 110.0;
		return comp;
	};
	inline LiquidComponent getH2O()
	{
		LiquidComponent comp;
		comp.mol_weight = 18.0;
		return comp;
	};
	inline GasComponent getCO2()
	{
		GasComponent comp;
		comp.mol_weight = 44.0;
		comp.z = 1.0;
		comp.rho_stc = 1.223;
		return comp;
	};

	template <int N>
	struct Reaction
	{
		static const int comp_num = N;
		std::array<Component, N> comps;
		std::array<double, N> indices;

		double alpha;
		double surf_init;
		double activation_energy;
		double reaction_const;
		inline adouble getReactionRate(double m0, adouble m) const
		{
			return 1.e-5 * /*reaction_const * surf_init **/ (1.0 - m); /*/ (1 - m0) *
				exp(-activation_energy / Component::R / Component::T)*/;
		}
	};

	static const int calcite_components_num = 5;
	struct CalciteReaction : Reaction<calcite_components_num>
	{
		enum REACTS {CALCITE, ACID, SALT, WATER, CO2};

		CalciteReaction()
		{
			comps[REACTS::CALCITE	] = getCaCO3();		indices[REACTS::CALCITE	] = -1.0;
			comps[REACTS::ACID		] = getHCl();		indices[REACTS::ACID	] = -2.0;
			comps[REACTS::SALT		] = getCaCl2();		indices[REACTS::SALT	] = 1.0;
			comps[REACTS::WATER		] = getH2O();		indices[REACTS::WATER	] = 1.0;
			comps[REACTS::CO2		] = getCO2();		indices[REACTS::CO2		] = 1.0;

			activation_energy = 13.0 * KKAL_2_J;
			reaction_const = 1.51 * 1.e+5;
			surf_init = 0.175;
			alpha = 1.0;
		};
	};

	static const int dolomite_components_num = 5;
	struct DolomiteReaction : Reaction<dolomite_components_num>
	{
	};

};

#endif /* ACID2D_REACTIONS_HPP_ */