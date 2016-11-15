#ifndef ACID2D_HPP_
#define ACID2D_HPP_

#include <vector>

#include "model/cells/Variables.hpp"
#include "model/cells/CylCell2D.h"
#include "model/AbstractModel.hpp"
#include "util/Interpolate.h"

namespace acid2d
{
	struct Skeleton_Props
	{
		// Density of skeleton matter in STC [kg/m3]
		double dens_stc;
		// Compessibility [1/Pa]
		double beta;

		// Average pore diameter [m]
		double d_pore;

		// Permeability if needed
		Interpolate* K;
		double anisotropy;
		inline double getPerm_r(const double m) const
		{
			return d_pore * d_pore * m * m * m / (1 - m) / (1 - m) / 150.0;
		};
		inline double getPerm_r_dm(const double m) const
		{
			return d_pore * d_pore / 150.0 *
					(3.0 * m * m * (1 - m) + 2 * m * m * m) / (1 - m) / (1 - m) / (1 - m);
		};
		inline double getPerm_z(const double m) const
		{
			return getPerm_r(m) * anisotropy;
		};
		inline double getPerm_z_dm(const double m) const
		{
			return getPerm_r_dm(m) * anisotropy;
		};

		// Top and bottom depth of perforation
		double h1, h2;
		// Height of formation [m]
		double height;

		int cellsNum_z;

		double p_out;

		// Initial values
		double m_init;
		double p_init;
		double s_init;
		double Ya_init;
		double Ys_init;
	};
	struct Component
	{

	};
	struct Liquid_Props
	{
		Component acid;
		Component salt;
		Component water;

		// Viscosity [cP]
		double visc;
		Interpolate* visc_table;
		inline double getViscosity(const double p) const
		{
			return visc;
		};
		inline double getViscosity_dp(const double p) const
		{
			return 0.0;
		};

		// Fluid volume factor
		Interpolate* b;
		// Density of fluid in STC [kg/m3]
		double dens_stc;
		// Compessibility [1/Pa]
		double beta;
		// Reference pressure [Pa]
		double p_ref;
		inline double getDensity(const double p) const
		{
			return dens_stc * (1.0 + beta * (p - p_ref));
		};
		inline double getDensity_dp() const
		{
			return dens_stc * beta;
		};

		// Relative fluid permeability
		Interpolate* kr;
		inline double getKr_l(const double s_l) const
		{
			if (s_l > 1.0)
				return 1.0;
			else
				return kr->Solve(s_l);
		};
		inline double getKr_l_ds(const double s_l) const
		{
			if (s_l > 1.0)
				return kr->DSolve(1.0);
			else
				return kr->DSolve(s_l);
		};


	};
	struct Gas_Props
	{
		// Viscosity [cP]
		double visc;
		Interpolate* visc_table;
		inline double getViscosity(const double p) const
		{
			return visc;
		};
		inline double getViscosity_dp(const double p) const
		{
			return 0.0;
		};

		// Fluid volume factor
		Interpolate* b;
		// Density of fluid in STC [kg/m3]
		double dens_stc;
		// Reference pressure [Pa]
		double p_ref;
		inline double getDensity(const double p) const
		{
			return dens_stc * p / p_ref;
		};
		inline double getDensity_dp(const double p) const
		{
			return dens_stc / p_ref;
		};

		// Relative fluid permeability
		Interpolate* kr;
		inline double getKr_g(const double s_l) const
		{
			if (s_l > 1.0)
				return 0.0;
			else
				return kr->Solve(s_l);
		};
		inline double getKr_g_ds(const double s_l) const
		{
			if (s_l > 1.0)
				return kr->DSolve(1.0);
			else
				return kr->DSolve(s_l);
		};
	};
	struct Properties
	{
	};

	typedef VarSimpleAcid Variable;
	typedef NewCylCell2D<VarSimpleAcid, Skeleton_Props> Cell;
	template <typename TVariable> using TCell = NewCylCell2D<TVariable, Skeleton_Props>;

	class Acid2d : public AbstractModel<Variable, Properties, TCell, Acid2d>
	{
		template<typename> friend class Snapshotter;
		template<typename> friend class GRDECLSnapshotter;
		template<typename> friend class VTKSnapshotter;
		template<typename> friend class AbstractMethod;
		friend class Acid2dSolver;

	protected:
		// Continuum properties
		int skeletonsNum;
		std::vector<Skeleton_Props> props_sk;
		Liquid_Props props_l;
		Gas_Props props_g;

		// Number of cells in radial direction
		int cellsNum_r;
		// Number of cells in vertical direction
		int cellsNum_z;

		// BHP will be converted to the depth
		double depth_point;

		void setInitialState();
		// Set all properties
		void setProps(Properties& props);
		// Make all properties dimensionless
		void makeDimLess();
		// Build grid
		void buildGridLog();
		// Set perforated cells
		void setPerforated();
		// Set some deviation to rate distribution
		void setRateDeviation(int num, double ratio);
		// Check formations properties
		void checkSkeletons(const std::vector<Skeleton_Props>& props);

		// Service functions
		inline double upwindIsCur(int cur, int beta)
		{
			if (cells[cur].u_next.p < cells[beta].u_next.p)
				return 0.0;
			else
				return 1.0;
		};
		inline int getUpwindIdx(int cur, int beta)
		{
			if (cells[cur].u_next.p < cells[beta].u_next.p)
				return beta;
			else
				return cur;
		};
		inline void getNeighborIdx(int cur, int* const neighbor)
		{
			neighbor[0] = cur - cellsNum_z - 2;
			neighbor[1] = cur + cellsNum_z + 2;
			neighbor[2] = cur - 1;
			neighbor[3] = cur + 1;
		};

		inline double getTrans(Cell& cell, Cell& beta)
		{
			double k1, k2, S;

			if (abs(cell.num - beta.num) == 1) {
				k1 = cell.props->getPerm_z(cell.u_next.m);
				k2 = beta.props->getPerm_z(beta.u_next.m);
				if (k1 == 0.0 && k2 == 0.0)
					return 0.0;
				S = 2.0 * M_PI * cell.r * cell.hr;
				return 2.0 * k1 * k2 * S / (k1 * beta.hz + k2 * cell.hz);
			}
			else {
				k1 = cell.props->getPerm_r(cell.u_next.m);
				k2 = beta.props->getPerm_r(beta.u_next.m);
				S = 2.0 * M_PI * cell.hz * (cell.r + sign(beta.num - cell.num) * cell.hr / 2.0);
				return 2.0 * k1 * k2 * S / (k1 * beta.hr + k2 * cell.hr);
			}
		};
		// First eqn
		/*double solve_eq1(int cur);
		double solve_eq1_dm(int cur);
		double solve_eq1_dp(int cur);
		double solve_eq1_ds(int cur);
		double solve_eq1_dYs(int cur);
		double solve_eq1_dYa(int cur);
		double solve_eq1_dm_beta(int cur, int beta);
		double solve_eq1_dp_beta(int cur, int beta);
		double solve_eq1_ds_beta(int cur, int beta);
		double solve_eq1_dYs_beta(int cur, int beta);
		double solve_eq1_dYa_beta(int cur, int beta);*/

	public:
		Acid2d();
		~Acid2d();

		void setPeriod(int period);
		double getRate(int cur);
	};
};

#endif /* ACID2D_HPP_ */