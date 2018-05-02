#ifndef WAXNIT1D_PROPERTIES_HPP_
#define WAXNIT1D_PROPERTIES_HPP_

#include <cassert>
#include "model/Basic1d/Properties.hpp"

namespace wax_nit1d
{
	// ADOLC stencil ids
	const int mid = basic1d::mid;
	const int left = basic1d::left;
	const int right = basic1d::right;

	struct Skeleton_Props : public basic1d::Skeleton_Props
	{
		// Initial values
		double p_init;
		double p_sat;
		double so_init;
		double m_init;
		double a_init;

		inline adouble getPermCoseni(adouble m) const
		{
			//return d_pore_r * d_pore_r * m * m * m / (1 - m) / (1 - m) / 150.0;
			return perm * (m * m * m / (1 - m) / (1 - m)) /
				(m_init * m_init * m_init / (1 - m_init) / (1 - m_init));
		};
	};
	struct Oil_Props : public basic1d::Liquid_Props
	{
		// Fluid volume factor
		Interpolate* b;
		inline adouble getB(adouble p) const
		{
			return b->Solve(p_ref) * exp((adouble)beta * (p_ref - p));
		};
		inline adouble getViscosity(const adouble p) const
		{
			return (adouble)(visc);
		};
		// Relative fluid permeability
		Interpolate* kr;
		inline adouble getKr(adouble s_o, const Skeleton_Props* props) const
		{
			adouble isAboveZero = (s_o - props->s_oc > 0.0) ? true : false;
			adouble isAboveCritical = (s_o > 1.0 - props->s_wc) ? true : false;
			adouble tmp;
			//condassign(tmp, isAboveZero, pow((s_o - (adouble)props->s_oc) / (adouble)(1.0 - props->s_wc - props->s_oc), 3.0), (adouble)0.0);
			//condassign(tmp, isAboveCritical, (adouble)1.0);
			return s_o;
			//return kr->Solve(s_w);
		};
		// Density of gas in STC [kg/m3]
		double dens_gas_stc;
		// Density of wax in STC [kg/m3]
		double dens_wax_stc;
		// Settling coeff
		double gamma;
		// Wax-oil ratio
		Interpolate* lp;
		inline adouble getRhoTilde(adouble p) const
		{
			return dens_stc / getB(p);
		};
		inline adouble getRho(adouble p) const
		{
			return getRhoTilde(p);
		};
	};
	struct Wax_Props
	{
		// Density of wax in STC [kg/m3]
		double dens_stc;
		inline adouble getRho() const
		{
			return dens_stc;
		};
	};

	struct Properties : public basic1d::Properties
	{
		Skeleton_Props props_sk;
		Oil_Props props_oil;
		Wax_Props props_wax;
		double L;

		std::vector< std::pair<double, double> > lp;
	};
};

#endif /* WAXNIT1D_PROPERTIES_HPP_ */
