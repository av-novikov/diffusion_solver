#ifndef BINGHAM2D_HPP_
#define BINGHAM2D_HPP_

#include "model/cells/Variables.hpp"
#include "model/AbstractModel.hpp"
#include "util/utils.h"
#include "model/cells/RadialCell.hpp"

#include "adolc/adouble.h"
#include "adolc/taping.h"

namespace bing1d
{
	static const int stencil = 3;
	static const int Lstencil = 3;
	static const int Rstencil = 2;
	static const int mid = 1;
	static const int left = 2;
	static const int right = 3;

	struct Skeleton_Props
	{
		// Porosity in STC
		double m;
		inline adouble getPoro(adouble p) const
		{
			return (adouble)(m)* ((adouble)(1.0) + (adouble)(beta)* (p /*- cell.props->p_init*/));
		};

		// Density of skeleton matter in STC [kg/m3]
		double dens_stc;
		// Compessibility [1/Pa]
		double beta;
		// Permeability along radial direction [mD]
		double perm_r;
		// Permeability along vertical direction [mD]
		double perm_z;

		inline double getPerm_r(const double r) const
		{
			return (r > radius_eff ? perm_r : perm_eff);
		};

		// Permeability of colmatage zone [mD]
		std::vector<double> perms_eff;
		// Radius of colmatage zone [m]
		std::vector<double> radiuses_eff;
		// Vector of skins
		std::vector<double> skins;
		double perm_eff;
		double radius_eff;
		double skin;

		// Height of formation [m]
		double height;

		int cellsNum_z;

		double p_out;

		// Initial values
		double p_init;
	};
	struct Oil_Props
	{
		// Viscosity [cP]
		double visc;
		Interpolate* visc_table;
		inline adouble getViscosity(const adouble p) const
		{
			return (adouble)(visc);
		};

		// Yield stress
		double tau0;

		// Density of fluid in STC [kg/m3]
		double dens_stc;
		// Compessibility [1/Pa]
		double beta;
		// Reference pressure [Pa]
		double p_ref;
		inline adouble getDensity(const adouble p) const
		{
			return (adouble)(dens_stc)* ((adouble)(1.0) + (adouble)(beta)* (p - (adouble)(p_ref)));
		};
		inline adouble getB(const adouble p) const
		{
			return (adouble)(1.0) - (adouble)(beta)* (p - (adouble)(p_ref));
		};
	};
	struct Properties
	{
		// Vector of start times of periods [sec]
		std::vector<double> timePeriods;
		// Vector of rates [m3/day]
		std::vector<double> rates;
		// Vector of BHPs [Pa]
		std::vector<double> pwf;

		// If left boundary condition would be 2nd type
		bool leftBoundIsRate;
		// If right boundary condition would be 1st type
		bool rightBoundIsPres;

		// Time step limits
		// Initial time step [sec]
		double ht;
		// Minimal time step [sec]
		double ht_min;
		// Maximum time step [sec]
		double ht_max;

		// Inner radius of well [m]
		double r_w;
		// Radius of formation [m]
		double r_e;

		// Number of cells in radial direction
		int cellsNum_r;

		Skeleton_Props props_sk;
		Oil_Props props_oil;

		double depth_point;

		// Data set (pressure, oil viscosity) ([Pa], [cP])
		std::vector< std::pair<double, double> > visc_oil;
	};

	typedef Var1phase Variable;
	typedef TapeVar1Phase TapeVariable;
	typedef RadialCell<Variable> Cell;

	class Bingham1d : public AbstractModel<Variable, Properties, RadialCell, Bingham1d>
	{
		template<typename> friend class Snapshotter;
		template<typename> friend class GRDECLSnapshotter;
		template<typename> friend class VTKSnapshotter;
		template<typename> friend class AbstractMethod;
		friend class Bing1dSolver;

		Skeleton_Props props_sk;
		Oil_Props props_oil;

		// Number of cells in radial direction
		int cellsNum_r;

		void setInitialState();
		void setProps(Properties& props);
		void makeDimLess();
		void buildGridLog();
		void setPerforated();

		inline double upwindIsCur(int cur, int beta)
		{
			if (cells[cur].u_next.p < cells[beta].u_next.p)
				return 0.0;
			else
				return 1.0;
		};
		inline int getUpwindIdx(int cur, int beta) const
		{
			if (cells[cur].u_next.p < cells[beta].u_next.p)
				return beta;
			else
				return cur;
		};
		inline void getNeighborIdx(int cur, int* const neighbor)
		{
			neighbor[0] = cur - 1;
			neighbor[1] = cur + 1;
		};
		inline double getTrans(const Cell& cell, const Cell& beta)
		{
			double k1, k2;
			k1 = (cell.r > props_sk.radius_eff ? props_sk.perm_r : props_sk.perm_eff);
			k2 = (beta.r > props_sk.radius_eff ? props_sk.perm_r : props_sk.perm_eff);
			return 4.0 * M_PI * k1 * k2 * (cell.r + sign(beta.num - cell.num) * cell.hr / 2.0) * props_sk.height / (k2 * cell.hr + k1 * beta.hr);
		};
		inline double getNablaP(const Cell& cell)
		{
			int neighbor[2];
			getNeighborIdx(cell.num, neighbor);

			Cell* nebr1 = &cells[neighbor[0]];
			Cell* nebr2 = &cells[neighbor[1]];

			double r_eff = props_sk.radius_eff;
			if ((nebr1->r < r_eff) && (nebr2->r > r_eff))
			{
				if (cell.r > r_eff)
					nebr1 = &const_cast<Cell&>(cell);
				else
					nebr2 = &const_cast<Cell&>(cell);
			}
			return (nebr2->u_next.p - nebr1->u_next.p) / (nebr2->r - nebr1->r);
		};
		inline double getNablaP(const Cell& cell1, const Cell& cell2)
		{
			return (cell2.u_next.p - cell1.u_next.p) / (cell2.r - cell1.r);
		};
		inline double getOilVelocity(const Cell& cell)
		{
			double gradp = getNablaP(cell);

			if (fabs(gradp) >= props_oil.tau0)
				return -props_sk.getPerm_r(cell.r) / props_oil.getViscosity(cell.u_next.p).value() *
							(1.0 - props_oil.tau0 / fabs(gradp)) * gradp;
			else
				return 0.0;
		};
		inline adouble linearInterp(adouble f1, const Cell& cell1, adouble f2, const Cell& cell2)
		{
			double r = cell1.r + sign(cell2.num - cell1.num) * cell1.hr / 2.0;
			return ((adouble)(fabs(r - cell1.r)) * f2 + (adouble)(fabs(r - cell2.r)) * f1) /
				(adouble)(cell2.r - cell1.r);
		};

		void solve_eqMiddle(int cur);
		void solve_eqLeft();
		void solve_eqRight();
		void setVariables(int cur);

	public:
		Bingham1d();
		~Bingham1d();

		double *x, *y;

		void setPeriod(int period);
		double getRate(int cur);
	};
};

#endif /* BINGHAM2D_HPP_ */
