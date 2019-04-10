#ifndef ACIDRECFRACMOV_REACTIONS_HPP_
#define ACIDRECFRACMOV_REACTIONS_HPP_

#define KKAL_2_J 4186.8

#include <array>

#include "adolc/adouble.h"
#include "adolc/taping.h"

#include "model/Basic2d/Properties.hpp"

namespace acidrecfracmov
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
		comp.mol_weight = 110;
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
		comp.rho_stc = 1.98;
		return comp;
	};
	inline SolidComponent getDolomite()
	{
		SolidComponent comp;
		comp.mol_weight = 184.4;
		comp.rho_stc = 2860.0;
		return comp;
	};
	inline SolidComponent getMgCl2CaCl2()
	{
		SolidComponent comp;
		comp.mol_weight = 111 + 95.2;
		comp.rho_stc = (2150 * 111 + 2316 * 95.2) / (111 + 95.2);
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
		inline adouble getReactionRate(const double m0, const double m_max, const adouble m) const
		{
			return reaction_const * surf_init /* pow((m_max - m) / (m_max - m0), 2)*/ *
				exp(-activation_energy / Component::R / Component::T);
		}		
		inline adouble getSpecificReactionRate() const
		{
			return reaction_const * exp(-activation_energy / Component::R / Component::T);
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

			activation_energy = 15.0 * KKAL_2_J;
			reaction_const = 7.29 * 1.e+7;
			surf_init = 100000.0;
			alpha = 0.63;
		};
	};
	static const int dolomite_components_num = 5;
	struct DolomiteReaction : Reaction<dolomite_components_num>
	{
		enum REACTS { DOLOMITE, ACID, SALT, WATER, CO2 };

		DolomiteReaction()
		{
			comps[REACTS::DOLOMITE] = getDolomite();	indices[REACTS::DOLOMITE] = -1.0;
			comps[REACTS::ACID] = getHCl();		indices[REACTS::ACID] = -4.0;
			comps[REACTS::SALT] = getMgCl2CaCl2(); indices[REACTS::SALT] = 1.0;
			comps[REACTS::WATER] = getH2O();		indices[REACTS::WATER] = 2.0;
			comps[REACTS::CO2] = getCO2();		indices[REACTS::CO2] = 2.0;

			activation_energy = 8.31 * 11320.0;
			alpha = 0.618 / 1.5;
			reaction_const = 9.4 * pow(10, 11 - 3.0 * alpha) * 2.E-6;
			surf_init = 100000.0;
		};
	};

};

#endif /* ACIDRECFRACMOV_REACTIONS_HPP_ */