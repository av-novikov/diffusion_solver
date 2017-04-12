#ifndef ACID2D_REACTIONS_HPP_
#define ACID2D_REACTIONS_HPP_

#include <array>

#include "adolc/adouble.h"
#include "adolc/taping.h"

namespace acid2d
{
	struct Component
	{
		static const double R;
		static double T;

		double mol_weight;
		double rho_stc;

		inline adouble getMolarDensity() const
		{
			return rho_stc / mol_weight;
		};
	};

	// Laboratory
	inline Component getCaCO3()
	{
		Component comp;
		comp.mol_weight = 100.0;
		return comp;
	};
	inline Component getHCl()
	{
		Component comp;
		comp.mol_weight = 36.0;
		return comp;
	};
	inline Component getCaCl2()
	{
		Component comp;
		comp.mol_weight = 110;
		return comp;
	};
	inline Component getH2O()
	{
		Component comp;
		comp.mol_weight = 18.0;
		return comp;
	};
	inline Component getCO2()
	{
		Component comp;
		comp.mol_weight = 44.0;
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
			return reaction_const * surf_init * (1.0 - m) / (1 - m0) *
				exp(-activation_energy / Component::R / Component::T);
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
		};
	};

	static const int dolomite_components_num = 5;
	struct DolomiteReaction : Reaction<dolomite_components_num>
	{
	};

};

#endif /* ACID2D_REACTIONS_HPP_ */