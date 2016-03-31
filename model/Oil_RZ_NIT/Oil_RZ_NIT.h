#ifndef OIL_RZ_NIT_H_
#define OIL_RZ_NIT_H_

#include <vector>
#include <map>
#include <cstdlib>

#include "util/utils.h"
#include "model/cells/Variables.hpp"
#include "model/cells/CylCell2D.h"
#include "model/AbstractModel.hpp"

namespace oil_rz_nit
{
	typedef CylCell2D<Var1phaseNIT> Cell;

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
		double p_init;
		double t_init;

		// Thermal properties

		// Thermal diffusivity coefficient [m2/sec]
		double kappa_eff;
		// Mass heat capacity [J/kg/K]
		double c;
		// Thermal conductivity coefficient [W/m/K]
		double lambda_r;
		// Thermal conductivity coefficient [W/m/K]
		double lambda_z;
	};

	struct Fluid_Props
	{
		// Viscosity [cP]
		double visc;
		// Density of fluid in STC [kg/m3]
		double dens_stc;
		// Volume factor for well bore
		double b_bore;
		// Compessibility [1/Pa]
		double beta;
		// Relative fluid permeability
		Interpolate* kr;
		// Fluid volume factor
		Interpolate* b;

		// Thermal properties

		// Mass heat capacity [J/kg/K]
		double c;
		// Thermal conductivity coefficient [W/m/K]
		double lambda;
		// Joule-thompson coefficient [K/Pa]
		double jt;
		// Adiabatic coefficient [K/Pa]
		double ad;
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
		// During the time flow rate decreases 'e' times in well test [sec] 
		double alpha;

		// Perforated intervals
		std::vector<std::pair<int, int> > perfIntervals;
		// Inner radius of well [m]
		double r_w;
		// Radius of formation [m]
		double r_e;

		// Number of cells in radial direction
		int cellsNum_r;
		// Number of cells in vertical direction
		int cellsNum_z;

		std::vector<Skeleton_Props> props_sk;
		Fluid_Props props_oil;

		// BHP will be converted to the depth
		double depth_point;
	};

	class Oil_RZ_NIT : public AbstractModel<Var1phaseNIT, Properties, CylCell2D, Oil_RZ_NIT>
	{
		template<typename> friend class Snapshotter;
		template<typename> friend class GRDECLSnapshotter;
		template<typename> friend class VTKSnapshotter;
		template<typename> friend class AbstractMethod;
		friend class OilRZNITSolver;	
	protected:

		// Continuum properties
		int skeletonsNum;
		std::vector<Skeleton_Props> props_sk;
		Fluid_Props props_oil;

		int cellsNum_r;
		int cellsNum_z;

		double alpha;

		double depth_point;

		void setInitialState();
		void setProps(Properties& props);
		void makeDimLess();
		void buildGridLog();
		void setPerforated();
		void setRateDeviation(int num, double ratio);
		void checkSkeletons(const std::vector<Skeleton_Props>& props);
		
		// Service functions
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
		inline int getSkeletonIdx(const Cell& cell) const
		{
			int idx = 0;
			while(idx < props_sk.size())
			{
				if(cell.z <= props_sk[idx].h2 + EQUALITY_TOLERANCE)
					return idx;
				idx++;
			}
			exit(-1);
		};
		// Solving coefficients
		inline double getPoro(double p, Cell& cell) const
		{
			const int idx = getSkeletonIdx(cell);
			return props_sk[idx].m * (1.0 + props_sk[idx].beta * (p /*- props_sk[idx].p_init*/) );
		};
		inline double getPoro_dp(Cell& cell) const
		{
			const int idx = getSkeletonIdx(cell);
			return props_sk[idx].m * props_sk[idx].beta;
		};
		inline double getRho(double p) const
		{
			return props_oil.dens_stc * (1.0 + props_oil.beta * p);
		};
		inline double getRho_dp() const
		{
			return props_oil.dens_stc * props_oil.beta;
		};
		inline double getTrans(Cell& cell, Cell& beta)
		{
			double k1, k2, S;
			const int idx1 = getSkeletonIdx(cell);
			const int idx2 = getSkeletonIdx(beta);

			if( abs(cell.num - beta.num) == 1) {
				k1 = props_sk[idx1].perm_z;
				k2 = props_sk[idx2].perm_z;
				if(k1 == 0.0 && k2 == 0.0)
					return 0.0;
				S = 2.0 * M_PI * cell.r * cell.hr;
				return 2.0 * k1 * k2 * S / (k1 * beta.hz + k2 * cell.hz);
			} else {
				k1 = (cell.r > props_sk[idx1].radius_eff ? props_sk[idx1].perm_r : props_sk[idx1].perm_eff);
				k2 = (beta.r > props_sk[idx2].radius_eff ? props_sk[idx2].perm_r : props_sk[idx2].perm_eff);
				if(k1 == 0.0 && k2 == 0.0)
					return 0.0;
				S = 2.0 * M_PI * cell.hz * (cell.r + sign(beta.num - cell.num) * cell.hr / 2.0);
				return 2.0 * k1 * k2 * S / (k1 * beta.hr + k2 * cell.hr);
			}
		};
		inline double getRho(Cell& cell, Cell& beta) const
		{
			if( abs(cell.num - beta.num) == 1 )
				return ( getRho(beta.u_next.p) * cell.hz + beta.hz * getRho(cell.u_next.p) ) / (beta.hz + cell.hz);
			else
				return ( getRho(beta.u_next.p) * cell.hr + beta.hr * getRho(cell.u_next.p) ) / (beta.hr + cell.hr);
		};
		inline double getRho_dp(Cell& cell, Cell& beta) const
		{
			if( abs(cell.num - beta.num) == 1 )
				return beta.hz * getRho_dp() / (beta.hz + cell.hz);
			else
				return beta.hr * getRho_dp() / (beta.hr + cell.hr);
		};
		inline double getRho_dp_beta(Cell& cell, Cell& beta) const
		{
			if( abs(cell.num - beta.num) == 1 )
				return cell.hz * getRho_dp() / (beta.hz + cell.hz);
			else
				return cell.hr * getRho_dp() / (beta.hr + cell.hr);
		};

		inline double getPerm_r(const Cell& cell) const
		{
			const int idx = getSkeletonIdx(cell);
			return (cell.r > props_sk[idx].radius_eff ? props_sk[idx].perm_r : props_sk[idx].perm_eff);
		};
		inline double getNablaP(Cell& cell, int varNum, int axis)
		{
			Cell* nebr1;
			Cell* nebr2;
			double h, r_eff;
			const int idx = getSkeletonIdx(cell);

			switch(axis)
			{
			case R_AXIS:
				nebr1 = &cells[cell.num - cellsNum_z - 2];
				nebr2 = &cells[cell.num + cellsNum_z + 2];
				
				r_eff = props_sk[idx].radius_eff;
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

			switch(varNum)
			{
			case PREV:
				return (nebr2->u_prev.p - nebr1->u_prev.p ) / h;
			case ITER:
				return (nebr2->u_iter.p - nebr1->u_iter.p ) / h;
			case NEXT:
				return (nebr2->u_next.p - nebr1->u_next.p ) / h;
			}
		};
		inline double getOilVelocity(Cell& cell, int varNum, int axis)
		{
			const int idx = getSkeletonIdx( cell );

			Var1phaseNIT* var;
			switch(varNum)
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

			switch(axis)
			{
			case R_AXIS:
				return -getPerm_r(cell) / props_oil.visc * getNablaP(cell, varNum, axis);
			case Z_AXIS:
				return -props_sk[idx].perm_z / props_oil.visc * getNablaP(cell, varNum, axis);
			}
		};

		inline double getCn(Cell& cell) const
		{
			const int idx = getSkeletonIdx(cell);
			return getPoro(cell.u_next.p, cell) * getRho(cell.u_next.p) * props_oil.c +
				(1.0 - getPoro(cell.u_next.p, cell)) * props_sk[idx].dens_stc * props_sk[idx].c;
		};
		inline double getAd(Cell& cell) const
		{
			return getPoro(cell.u_next.p, cell) * getRho(cell.u_next.p) * props_oil.c * props_oil.ad;
		};
		inline double getLambda(Cell& cell, int axis)
		{
			const int idx = getSkeletonIdx(cell);
			switch (axis)
			{
			case R_AXIS:
				return getPoro(cell.u_next.p, cell) * props_oil.lambda + (1.0 - getPoro(cell.u_next.p, cell)) * props_sk[idx].lambda_r;
			case Z_AXIS:
				return getPoro(cell.u_next.p, cell) * props_oil.lambda + (1.0 - getPoro(cell.u_next.p, cell)) * props_sk[idx].lambda_z;
			}
		};
		inline double getLambda(Cell& cell1, Cell& cell2)
		{
			double cm;
			if (abs(cell1.num - cell2.num) == 1)
				return (cell1.hz * getLambda(cell2, Z_AXIS) + cell2.hz * getLambda(cell1, Z_AXIS)) / (cell1.hz + cell2.hz);
			else 
				return (cell1.hr * getLambda(cell2, R_AXIS) + cell2.hr * getLambda(cell1, R_AXIS)) / (cell1.hr + cell2.hr);
		};
		inline double getJT(Cell& cell, int varNum, int axis)
		{
			Var1phaseNIT* var;
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

			return getRho(var->p) * props_oil.c * props_oil.jt * getOilVelocity(cell, varNum, axis);
		};
		inline double getA(Cell& cell, int varNum, int axis)
		{
			Var1phaseNIT* var;
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

			return getRho(var->p) * props_oil.c * getOilVelocity(cell, varNum, axis);
		};

		double solve_eq(int cur);
		double solve_eq_dp(int cur);
		double solve_eq_dp_beta(int cur, int beta);
		
		double solve_eqLeft(int cur);
		double solve_eqLeft_dp(int cur);
		double solve_eqLeft_dp_beta(int cur);
		
		double solve_eqRight(int cur);
		double solve_eqRight_dp(int cur);
		double solve_eqRight_dp_beta(int cur);

		double solve_eqTop(int cur);
		double solve_eqTop_dp(int cur);
		double solve_eqTop_dp_beta(int cur);

		double solve_eqBot(int cur);
		double solve_eqBot_dp(int cur);
		double solve_eqBot_dp_beta(int cur);

		double solveH();
		double getRate(int cur);

	public:
		Oil_RZ_NIT();
		~Oil_RZ_NIT();

		void setPeriod(int period);
	};

};

#endif /* OIL_RZ_NIT_H_ */