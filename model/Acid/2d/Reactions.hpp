#ifndef ACID2D_REACTIONS_HPP_
#define ACID2D_REACTIONS_HPP_

#include <array>

#include "adolc/adouble.h"
#include "adolc/taping.h"

namespace acid2d
{
	struct Component
	{
		static double R;
		static double T;

		double mol_weight;
		double rho_stc;
	};

	// Laboratory
	Component getCaCO3()
	{
		Component comp;
		comp.mol_weight = 100.0;
		return comp;
	};
	Component getHCl()
	{
		Component comp;
		comp.mol_weight = 36.0;
		return comp;
	};
	Component getCaCl2()
	{
		Component comp;
		comp.mol_weight = 110;
		return comp;
	};
	Component getH2O()
	{
		Component comp;
		comp.mol_weight = 18.0;
		return comp;
	};
	Component getCO2()
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

		double surf_init;
		double m_init;
		double activation_energy;
		double reaction_const;
	};

	static const int calcite_components_num = 5;
	struct CalciteReaction : Reaction<calcite_components_num>
	{
		static const int calcite = 0;
		static const int acid = 1;
		static const int salt = 2;
		static const int water = 3;
		static const int co2 = 4;

		CalciteReaction()
		{
			comps[calcite	] = getCaCO3();		indices[calcite	] = -1.0;
			comps[acid		] = getHCl();		indices[acid	] = -2.0;
			comps[salt		] = getCaCl2();		indices[salt	] = 1.0;
			comps[water		] = getH2O();		indices[water	] = 1.0;
			comps[co2		] = getCO2();		indices[co2		] = 1.0;
		};
	};

	static const int dolomite_components_num = 5;
	struct DolomiteReaction : Reaction<dolomite_components_num>
	{
	};

};

#endif /* ACID2D_REACTIONS_HPP_ */