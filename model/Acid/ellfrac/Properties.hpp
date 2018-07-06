#ifndef ACIDELLFRAC_PROPERTIES_HPP_
#define ACIDELLFRAC_PROPERTIES_HPP_

#include <array>

#include "model/Basic1d/Properties.hpp"
#include "model/Acid/ellfrac/Reactions.hpp"

namespace acidellfrac
{
	// ADOLC stencil ids
	const int mid = basic1d::mid;
	const int left = basic1d::left;
	const int right = basic1d::right;

	struct FracProperties
	{
		double w2, l2, height;
		double p_init, c_init;
	};
	struct Skeleton_Props : public basic1d::Skeleton_Props
	{
		//SolidComponent cur_mineral;
		double xa_eqbm;

		// Initial values
		double p_init;
		double t_init;
		double sw_init;
		double so_init;
		double xa_init;
		double xw_init;

		inline adouble getPoro(adouble m0, adouble p) const
		{
			return m0;
			//return m0 * (1.0 + beta * (p - p_init));
		};
		inline adouble getPermCoseni(adouble m0, adouble p) const
		{
			//return d_pore_r * d_pore_r * m * m * m / (1 - m) / (1 - m) / 150.0;
			adouble m = getPoro(m0, p);
			return perm * (m * m * m / (1 - m) / (1 - m)) / 
							(m_init * m_init * m_init / (1 - m_init) / (1 - m_init));
		};
		inline double getInitDiam(double m_init, double k0)
		{
			return sqrt(150.0 * k0 * (1.0 - m_init) * (1.0 - m_init) / m_init / m_init / m_init);
		};

		inline adouble getDensity(adouble p) const
		{
			return dens_stc;
		};
	};
	struct Water_Props : public basic1d::Liquid_Props
	{
		//LiquidComponent acid;
		//SolidComponent salt;
		//LiquidComponent water;
		double D_e;

		Interpolate* kr;
		inline adouble getKr(adouble sw, adouble m, const Skeleton_Props* props) const
		{
			adouble isAboveZero = (sw - props->s_wc > 0.0) ? true : false;
			adouble isAboveCritical = (sw > 1.0 - props->s_oc) ? true : false;
			adouble tmp;
			condassign(tmp, isAboveZero, ((1.0 - m) * pow((sw - props->s_wc) / (1.0 - props->s_wc - props->s_oc), 3.0) + (m - props->m_init) * (sw - props->s_wc) / (1.0 - props->s_wc - props->s_oc)) / (1.0 - props->m_init), (adouble)0.0);
			//condassign(tmp, isAboveZero, pow((sw - props->s_wc) / (1.0 - props->s_wc - props->s_oc), 3.0), (adouble)0.0);
			condassign(tmp, isAboveCritical, (adouble)1.0);
			return tmp;
		};
		inline adouble getViscosity(adouble p, adouble xa, adouble xw, adouble xs) const
		{
			return visc;
		};
		inline adouble getDensity(adouble p, adouble xa, adouble xw, adouble xs) const
		{
			return dens_stc;
		};
	};
	struct Oil_Props : public basic1d::Liquid_Props
	{
		double gas_dens_stc;
		//LiquidComponent oil;
		Interpolate* b;
		inline adouble getB(adouble p) const
		{
			return exp((adouble)beta * (p_ref - p));
		};
		inline adouble getDensity(adouble p) const
		{
			return dens_stc / getB(p);
		};

		Interpolate* kr;
		inline adouble getKr(adouble sw, adouble m, const Skeleton_Props* props) const
		{
			adouble isAboveZero = (1 - sw - props->s_oc > 0.0) ? true : false;
			adouble isAboveCritical = (1 - sw > 1.0 - props->s_wc) ? true : false;
			adouble tmp;
			condassign(tmp, isAboveZero, ((1.0 - m) * pow((1 - sw - props->s_oc) / (1.0 - props->s_wc - props->s_oc), 3.0) + (m - props->m_init) * (1 - sw - props->s_oc) / (1.0 - props->s_wc - props->s_oc)) / (1.0 - props->m_init), (adouble)0.0);
			//condassign(tmp, isAboveZero, pow((1 - sw - props->s_oc) / (1.0 - props->s_wc - props->s_oc), 3.0), (adouble)0.0);
			condassign(tmp, isAboveCritical, (adouble)1.0);
			return tmp;
		};
		inline adouble getViscosity(adouble p) const
		{
			return visc;
		};
	};
	struct Gas_Props : public basic1d::Gas_Props
	{
		GasComponent co2;

		Interpolate* kr;
		inline adouble getKr(adouble sw, adouble so, const Skeleton_Props* props) const
		{
			adouble isAboveZero = (1.0 - so - sw - props->s_gc > 0.0) ? true : false;
			adouble isAboveCritical = (1.0 - so - sw > 1.0 - props->s_wc - props->s_oc) ? true : false;
			adouble tmp;
			condassign(tmp, isAboveZero, 0.5 * pow(((adouble)(1.0 - props->s_gc) - sw - so) / (adouble)(1.0 - props->s_wc - props->s_oc - props->s_gc), 3.0), (adouble)0.0);
			condassign(tmp, isAboveCritical, (adouble)0.5);
			return tmp;
		};
		inline adouble getViscosity(adouble p) const
		{
			return visc;
		};
		Interpolate* rho;
		inline adouble getDensity(adouble p) const
		{
			return rho->Solve(p);
			//return co2.getDensity(p);
		};
	};
	struct Properties : public basic1d::Properties
	{
		FracProperties props_frac;
		std::vector<Skeleton_Props> props_sk;
		Water_Props props_w;
		Oil_Props props_o;
		Gas_Props props_g;

		size_t cellsNum_mu_frac, cellsNum_z;
		size_t cellsNum_mu_poro;
		double re;

		std::vector<double> cs;
		std::vector<bool> LeftBoundIsRate;
		std::vector< std::pair<double, double> > rho_co2;
	};
};

#endif /* ACIDELLFRAC_PROPERTIES_HPP_ */