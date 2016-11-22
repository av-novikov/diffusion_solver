#ifndef VPP2D_HPP_
#define VPP2D_HPP_

#include <vector>

#include "model/cells/Variables.hpp"
#include "model/cells/CylCell2D.h"
#include "model/AbstractModel.hpp"
#include "util/utils.h"

namespace vpp2d
{
	struct Skeleton_Props
	{
		// Porosity in STC
		double m;
		// Density of skeleton matter in STC [kg/m3]
		double dens_stc;
		// Compessibility [1/Pa]
		double beta;
		// Permeability along radial direction [mD]
		double perm_r;
		// Permeability along vertical direction [mD]
		double perm_z;

		// Permeability of colmatage zone [mD]
		std::vector<double> perms_eff;
		// Radius of colmatage zone [m]
		std::vector<double> radiuses_eff;
		// Vector of skins
		std::vector<double> skins;
		double perm_eff;
		double radius_eff;
		double skin;

		// Top and bottom depth of perforation
		double h1, h2;
		// Height of formation [m]
		double height;

		int cellsNum_z;

		double p_out;

		// Initial values
		double p_init;
		double s_init;
		double c_init;
	};
	struct Component
	{

	};
	struct Water_Props
	{
		Component active;

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
		double b_bore;
		// Density of fluid in STC [kg/m3]
		double dens_stc;
		// Compessibility [1/Pa]
		double beta;
		inline double getB(double p) const
		{
			return b->Solve(p);
		};
		inline double getBoreB(double p) const
		{
			return b_bore;
		};
		inline double getB_dp(double p) const
		{
			return b->DSolve(p);
		};
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
		inline double getKr(const double s_l) const
		{
			if (s_l > 1.0)
				return 1.0;
			else
				return kr->Solve(s_l);
		};
		inline double getKr_ds(const double s_l) const
		{
			if (s_l > 1.0)
				return kr->DSolve(1.0);
			else
				return kr->DSolve(s_l);
		};
	};
	struct Oil_Props
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
		double b_bore;
		inline double getB(double p) const
		{
			return b->Solve(p);
		};
		inline double getBoreB_oil(double p) const
		{
			return b_bore;
		};
		inline double getB_dp(double p) const
		{
			return b->DSolve(p);
		};
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
		inline double getKr(const double s_l) const
		{
			if (s_l > 1.0)
				return 0.0;
			else
				return kr->Solve(s_l);
		};
		inline double getKr_ds(const double s_l) const
		{
			if (s_l > 1.0)
				return kr->DSolve(1.0);
			else
				return kr->DSolve(s_l);
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

		// Perforated intervals
		std::vector<std::pair<int, int> > perfIntervals;
		// Time step limits
		// Initial time step [sec]
		double ht;
		// Minimal time step [sec]
		double ht_min;
		// Maximum time step [sec]
		double ht_max;
		// During the time flow rate decreases 'e' times in well test [sec] 
		double alpha;

		// Inner radius of well [m]
		double r_w;
		// Radius of formation [m]
		double r_e;

		// Number of cells in radial direction
		int cellsNum_r;
		// Number of cells in vertical direction
		int cellsNum_z;

		std::vector<Skeleton_Props> props_sk;
		Water_Props props_w;
		Oil_Props props_o;

		double depth_point;

		// Data set (saturation, relative oil permeability)
		std::vector< std::pair<double, double> > kr_o;
		// Data set (saturation, relative gas permeability)
		std::vector< std::pair<double, double> > kr_w;

		// Data set (pressure, oil volume factor) ([Pa], [m3/m3])
		std::vector< std::pair<double, double> > B_o;
		// Data set (pressure, gas volume factor) ([Pa], [m3/m3])
		std::vector< std::pair<double, double> > B_w;

		// Data set (pressure, oil volume factor) ([Pa], [cP])
		std::vector< std::pair<double, double> > visc_o;
		// Data set (pressure, gas volume factor) ([Pa], [cP])
		std::vector< std::pair<double, double> > visc_w;

		// Data set (pressure, gas volume factor) ([Pa], [1])
		std::vector< std::pair<double, double> > a;
	};

	typedef VarSimpleVPP Variable;
	typedef NewCylCell2D<VarSimpleVPP, Skeleton_Props> Cell;
	template <typename TVariable> using TCell = NewCylCell2D<TVariable, Skeleton_Props>;

	class VPP2d : public AbstractModel<Variable, Properties, TCell, VPP2d>
	{
		template<typename> friend class Snapshotter;
		template<typename> friend class GRDECLSnapshotter;
		template<typename> friend class VTKSnapshotter;
		template<typename> friend class AbstractMethod;
		friend class VPP2dSolver;

	protected:
		// Continuum properties
		int skeletonsNum;
		std::vector<Skeleton_Props> props_sk;
		Water_Props props_w;
		Oil_Props props_o;
		Interpolate* a;

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
		inline int getSkeletonIdx(const Cell& cell) const
		{
			int idx = 0;
			while (idx < props_sk.size())
			{
				if (cell.z <= props_sk[idx].h2 + EQUALITY_TOLERANCE)
					return idx;
				idx++;
			}
			exit(-1);
		};
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
				k1 = cell.props->perm_z;
				k2 = beta.props->perm_z;
				if (k1 == 0.0 && k2 == 0.0)
					return 0.0;
				S = 2.0 * M_PI * cell.r * cell.hr;
				return 2.0 * k1 * k2 * S / (k1 * beta.hz + k2 * cell.hz);
			}
			else {
				k1 = (cell.r > cell.props->radius_eff ? cell.props->perm_r : cell.props->perm_eff);
				k2 = (beta.r > beta.props->radius_eff ? beta.props->perm_r : beta.props->perm_eff);
				S = 2.0 * M_PI * cell.hz * (cell.r + sign(beta.num - cell.num) * cell.hr / 2.0);
				return 2.0 * k1 * k2 * S / (k1 * beta.hr + k2 * cell.hr);
			}
		};
		inline double getPerm_r(const Cell& cell) const
		{
			return (cell.r > cell.props->radius_eff ? cell.props->perm_r : cell.props->perm_eff);
		};
		inline double getNablaP(Cell& cell, int varNum, int axis)
		{
			Cell* nebr1;
			Cell* nebr2;
			double h, r_eff;

			switch (axis)
			{
			case R_AXIS:
				nebr1 = &cells[cell.num - cellsNum_z - 2];
				nebr2 = &cells[cell.num + cellsNum_z + 2];

				r_eff = cell.props->radius_eff;
				if ((nebr1->r < r_eff) && (nebr2->r > r_eff))
				{
					if (cell.r > r_eff)
						nebr1 = &cell;
					else
						nebr2 = &cell;
				}
				h = nebr2->r - nebr1->r;

				break;
			case Z_AXIS:
				nebr1 = &cells[cell.num - 1];
				nebr2 = &cells[cell.num + 1];
				h = nebr2->z - nebr1->z;
				break;
			}

			switch (varNum)
			{
			case PREV:
				return (nebr2->u_prev.p - nebr1->u_prev.p) / h;
			case ITER:
				return (nebr2->u_iter.p - nebr1->u_iter.p) / h;
			case NEXT:
				return (nebr2->u_next.p - nebr1->u_next.p) / h;
			}
		};
		inline double getWaterVelocity(Cell& cell, int varNum, int axis)
		{
			Variable* var;
			switch (varNum)
			{
			case PREV:
				var = &cell.u_prev;
				break;
			case ITER:
				var = &cell.u_iter;
				break;
			case NEXT:
				var = &cell.u_next;
				break;
			}

			switch (axis)
			{
			case R_AXIS:
				return -getPerm_r(cell) * props_w.getKr(var->s) / props_w.visc * getNablaP(cell, varNum, axis);
			case Z_AXIS:
				return -cell.props->perm_z * props_w.getKr(var->s) / props_w.visc * getNablaP(cell, varNum, axis);
			}
		};
		inline double getOilVelocity(Cell& cell, int varNum, int axis)
		{
			Variable* var;
			switch (varNum)
			{
			case PREV:
				var = &cell.u_prev;
				break;
			case ITER:
				var = &cell.u_iter;
				break;
			case NEXT:
				var = &cell.u_next;
				break;
			}

			switch (axis)
			{
			case R_AXIS:
				return -getPerm_r(cell) * props_o.getKr(var->s) / props_o.visc * getNablaP(cell, varNum, axis);
			case Z_AXIS:
				return -cell.props->perm_z * props_o.getKr(var->s) / props_o.visc * getNablaP(cell, varNum, axis);
			}
		};

		// Finds functional
		double solveH();

	public:
		VPP2d();
		~VPP2d();

		void setPeriod(int period);
		double getRate(int cur);
	};

};

#endif /* VPP2D_HPP_ */
