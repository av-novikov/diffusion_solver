#ifndef ACID2D_PROPERTIES_HPP_
#define ACID2D_PROPERTIES_HPP_

#include <array>

#include "model/Acid/2d/Reactions.hpp"

namespace acid2d
{
	// ADOLC stencil ids
	const int mid = basic2d::mid;
	const int left = basic2d::left;
	const int right = basic2d::right;
	const int vertical = basic2d::vertical;

	struct Skeleton_Props : public basic2d::Skeleton_Props
	{
		//SolidComponent cur_mineral;
		double xa_eqbm;

		// Initial values
		double p_init;
		double sw_init;
		double so_init;
		double xa_init;
		double xw_init;

		//double d_pore_r, d_pore_z;
		inline adouble getPermCoseni_r(adouble m) const
		{
			//return d_pore_r * d_pore_r * m * m * m / (1 - m) / (1 - m) / 150.0;
			return perm_r * (m * m * m / (1 - m) / (1 - m)) / 
							(m_init * m_init * m_init / (1 - m_init) / (1 - m_init));
		};
		inline adouble getPermCoseni_z(adouble m) const
		{
			//return d_pore_z * d_pore_z * m * m * m / (1 - m) / (1 - m) / 150.0;
			return perm_z * (m * m * m / (1 - m) / (1 - m)) /
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
	struct Water_Props : public basic2d::Liquid_Props
	{
		//LiquidComponent acid;
		//SolidComponent salt;
		//LiquidComponent water;

		Interpolate* kr;
		inline adouble getKr(adouble sw, adouble so, const Skeleton_Props* props) const
		{
			adouble isAboveZero = (sw - props->s_wc > 0.0) ? true : false;
			adouble isAboveCritical = (sw > 1.0 - props->s_oc - props->s_gc) ? true : false;
			adouble tmp;
			condassign(tmp, isAboveZero, pow((sw - (adouble)props->s_wc) / (adouble)(1.0 - props->s_wc - props->s_oc - props->s_gc), 3.0), (adouble)0.0);
			condassign(tmp, isAboveCritical, (adouble)1.0);
			return tmp;
		};
		inline adouble getViscosity(adouble p, adouble xa, adouble xw) const
		{
			return visc;
		};
		inline adouble getDensity(adouble p, adouble xa, adouble xw) const
		{
			return dens_stc;
		};
	};
	struct Oil_Props : public basic2d::Liquid_Props
	{
		//LiquidComponent oil;

		Interpolate* kr;
		inline adouble getKr(adouble sw, adouble so, const Skeleton_Props* props) const
		{
			adouble isAboveZero = (so - props->s_oc > 0.0) ? true : false;
			adouble isAboveCritical = (so > 1.0 - props->s_wc - props->s_gc) ? true : false;
			adouble tmp;
			condassign(tmp, isAboveZero, pow((so - (adouble)props->s_oc) / (adouble)(1.0 - props->s_wc - props->s_oc - props->s_gc), 3.0), (adouble)0.0);
			condassign(tmp, isAboveCritical, (adouble)1.0);
			return tmp;
		};
		inline adouble getViscosity(adouble p) const
		{
			return visc;
		};
		inline adouble getDensity(adouble p) const
		{
			return dens_stc;
		}
	};
	struct Gas_Props : public basic2d::Gas_Props
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
		inline adouble getDensity(adouble p) const
		{
			return co2.getDensity(p);
		};
	};
	struct Properties : public basic2d::Properties
	{
		std::vector<Skeleton_Props> props_sk;
		Water_Props props_w;
		Oil_Props props_o;
		Gas_Props props_g;

		std::vector<double> xa;
	};
};

#endif /* ACID2D_PROPERTIES_HPP_ */