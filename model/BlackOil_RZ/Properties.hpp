#ifndef BLACKOILRZ_PROPERTIES_HPP_
#define BLACKOILRZ_PROPERTIES_HPP_

#include <array>
#include "model/Basic2d/Properties.hpp"

namespace blackoil_rz
{
	const int mid = basic2d::mid;
	const int left = basic2d::left;
	const int right = basic2d::right;
	const int vertical = basic2d::vertical;

	struct Skeleton_Props : public basic2d::Skeleton_Props
	{
		// Initial values
		double p_init;
		double p_sat;
		double so_init;
		double sw_init;
	};
	struct Water_Props : public basic2d::Liquid_Props
	{
		// Fluid volume factor
		Interpolate* b;
		inline adouble getB(adouble p, adouble p_bub, adouble SATUR) const
		{
			//adouble tmp;
			//condassign(tmp, SATUR, b->Solve(p), b->Solve(p_bub));
			return exp((adouble)beta * (p - p_ref));
		};
		inline adouble getViscosity(adouble p) const
		{
			return (adouble)(visc);
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
	};
	struct Oil_Props : public basic2d::Liquid_Props
	{
		// Fluid volume factor
		Interpolate* b;
		inline adouble getB(adouble p, adouble p_bub, adouble SATUR) const
		{
			adouble tmp;
			condassign(tmp, SATUR, b->Solve(p), b->Solve(p_bub) * exp((adouble)beta * (p_bub - p)));
			return tmp;
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
	};
	struct Gas_Props : public basic2d::Gas_Props
	{
		inline adouble getViscosity(const adouble p) const
		{
			return (adouble)(visc);
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
	};

	struct Properties : public basic2d::Properties
	{
		std::vector<Skeleton_Props> props_sk;
		Oil_Props props_oil;
		Water_Props props_wat;
		Gas_Props props_gas;
	};
};

#endif /* BLACKOILRZ_PROPERTIES_HPP_ */