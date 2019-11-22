#ifndef ACID2DREC_PROPERTIES_HPP_
#define ACID2DREC_PROPERTIES_HPP_

#include <array>
#include <valarray>

#include "model/Basic1d/Properties.hpp"
#include "model/Acid/2drec/Reactions.hpp"

namespace acid2drec
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
	struct Skeleton_Props
	{
		int id;
		//SolidComponent cur_mineral;
		double xa_eqbm;
		// Initial values
		double p_init;
		double t_init;
		double sw_init;
		double so_init;
		double xa_init;
		double xw_init;
		// Maximum porosity
		double m_max;
		// Coefficient in petrophysical relationsship
		double A;

		inline adouble getPoro(const adouble& m0, const adouble& p) const
		{
			return m0;// *exp(beta * (p - p_ref));
		};
		/*inline adouble getPermCoseni(const adouble& m0, const adouble& p) const
		{
			//adouble m = getPoro(m0, p);
			return perm * (m0 * m0 * m0 / (1 - m0) / (1 - m0)) / 
							(m_init * m_init * m_init / (1 - m_init) / (1 - m_init));
		};*/
		inline adouble getPermCoseni(const adouble& m0, const adouble& p) const
		{
			//adouble m = getPoro(m0, p);
			return perm * exp(A * (m0 - m_init));
		};
		inline adouble getDensity(const adouble& p) const
		{
			return dens_stc;
		};

		double beta;
		double dens_stc;
		double perm;
		double m_init;
		double p_ref;
		// Top and bottom depth of perforation
		double h1, h2;
		// Height of formation [m]
		double height;
		// Connate saturations
		double s_wc, s_oc, s_gc;

		double p_out;

		double hx, hy, hz;
		int cellsNum;
	};
	struct Water_Props : public basic1d::Liquid_Props
	{
		//LiquidComponent acid;
		//SolidComponent salt;
		//LiquidComponent water;
		double D_e;

		Interpolate* kr;
		inline adouble getKr(const adouble& sw, const adouble& m, const Skeleton_Props* props) const
		{
			adouble isAboveZero = (sw - props->s_wc > 0.0) ? true : false;
			adouble isAboveCritical = (sw > 1.0 - props->s_oc) ? true : false;
			adouble tmp;
			condassign(tmp, isAboveZero, ((props->m_max - m) * pow((sw - props->s_wc) / (1.0 - props->s_wc - props->s_oc), 3.0) + (m - props->m_init) * (sw - props->s_wc) / (1.0 - props->s_wc - props->s_oc)) / (props->m_max - props->m_init), (adouble)0.0);
			//condassign(tmp, isAboveZero, pow((sw - props->s_wc) / (1.0 - props->s_wc - props->s_oc), 3.0), (adouble)0.0);
			condassign(tmp, isAboveCritical, (adouble)1.0);
			return tmp;
		};
		inline adouble getViscosity(const adouble& p, const adouble& xa, const adouble& xw, const adouble& xs) const
		{
			return visc;
		};
		inline adouble getDensity(const adouble& p, const adouble& xa, const adouble& xw, const adouble& xs) const
		{
			return getRho(p);// dens_stc;
		};
		inline adouble getRho(const adouble& p) const
		{
			return dens_stc * exp(beta * (p - p_ref));
		};
	};
	struct Oil_Props : public basic1d::Liquid_Props
	{
		double gas_dens_stc;
		//LiquidComponent oil;
		Interpolate* b;
		inline adouble getB(const adouble& p) const
		{
			return exp((adouble)beta * (p_ref - p));
		};
		inline adouble getDensity(const adouble& p) const
		{
			return dens_stc / getB(p);
		};
		inline adouble getRho(const adouble& p) const
		{
			return dens_stc * exp(beta * (p - p_ref));
		};

		Interpolate* kr;
		inline adouble getKr(const adouble& sw, const adouble& m, const Skeleton_Props* props) const
		{
			adouble isAboveZero = (1 - sw - props->s_oc > 0.0) ? true : false;
			adouble isAboveCritical = (1 - sw > 1.0 - props->s_wc) ? true : false;
			adouble tmp;
			condassign(tmp, isAboveZero, ((props->m_max - m) * pow((1 - sw - props->s_oc) / (1.0 - props->s_wc - props->s_oc), 3.0) + (m - props->m_init) * (1 - sw - props->s_oc) / (1.0 - props->s_wc - props->s_oc)) / (props->m_max - props->m_init), (adouble)0.0);
			//condassign(tmp, isAboveZero, pow((1 - sw - props->s_oc) / (1.0 - props->s_wc - props->s_oc), 3.0), (adouble)0.0);
			condassign(tmp, isAboveCritical, (adouble)1.0);
			return tmp;
		};
		inline adouble getViscosity(const adouble& p) const
		{
			return visc;
		};
	};
	struct Gas_Props : public basic1d::Gas_Props
	{
		GasComponent co2;

		Interpolate* kr;
		inline adouble getKr(const adouble& sw, const adouble& so, const Skeleton_Props* props) const
		{
			adouble isAboveZero = (1.0 - so - sw - props->s_gc > 0.0) ? true : false;
			adouble isAboveCritical = (1.0 - so - sw > 1.0 - props->s_wc - props->s_oc) ? true : false;
			adouble tmp;
			condassign(tmp, isAboveZero, 0.5 * pow(((adouble)(1.0 - props->s_gc) - sw - so) / (adouble)(1.0 - props->s_wc - props->s_oc - props->s_gc), 3.0), (adouble)0.0);
			condassign(tmp, isAboveCritical, (adouble)0.5);
			return tmp;
		};
		inline adouble getViscosity(const adouble& p) const
		{
			return visc;
		};
		Interpolate* rho;
		inline adouble getDensity(const adouble& p) const
		{
			return rho->Solve(p);
			//return co2.getDensity(p);
		};
	};
    struct ProdProps
    {
        double x_size, y_size, z_size, R_dim;
        int nx, ny;
    };
	struct Properties : public basic1d::Properties
	{
		std::vector<Skeleton_Props> props_sk;
		Water_Props props_w;
		Oil_Props props_o;
		Gas_Props props_g;

		size_t cellsNum_x, cellsNum_y;
		double hx, hy, hz;

		double max_sol_volume, max_acid_volume;
		std::vector<double> cs;
		std::vector<bool> LeftBoundIsRate;
		std::vector< std::pair<double, double> > rho_co2;
        double R_dim;

        std::string prefix, permFile;
		bool fieldData;
		bool permFromFile;
	};
};

#endif /* ACID2DREC_PROPERTIES_HPP_ */