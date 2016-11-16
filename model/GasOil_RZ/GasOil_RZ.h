#ifndef GASOIL_RZ_H_
#define GASOIL_RZ_H_

#include <vector>
#include <map>
#include <string>

#include "model/cells/Variables.hpp"
#include "model/cells/CylCell2D.h"
#include "model/AbstractModel.hpp"
#include "util/Interpolate.h"
#include "util/utils.h"

namespace gasOil_rz
{
	typedef CylCell2D<Var2phase> Cell;

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
		double p_bub;
		double s_init;
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
		std::vector<std::pair<int,int> > perfIntervals;
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
		Fluid_Props props_oil;
		Fluid_Props props_gas;

		double depth_point;

		// Data set (saturation, relative oil permeability)
		std::vector< std::pair<double,double> > kr_oil;
		// Data set (saturation, relative gas permeability)
		std::vector< std::pair<double,double> > kr_gas;

		// Data set (pressure, oil volume factor) ([Pa], [m3/m3])
		std::vector< std::pair<double,double> > B_oil;
		// Data set (pressure, gas volume factor) ([Pa], [m3/m3])
		std::vector< std::pair<double,double> > B_gas;

		// Data set (pressure, gas content in oil) ([Pa], [m3/m3])
		std::vector< std::pair<double,double> > Rs;
	};

	class GasOil_RZ : public AbstractModel<Var2phase, Properties, CylCell2D, GasOil_RZ>
	{
		template<typename> friend class Snapshotter;
		template<typename> friend class GRDECLSnapshotter;
		template<typename> friend class VTKSnapshotter;
		template<typename> friend class AbstractMethod;
		friend class GasOil2DSolver;	

	protected:
		// Continuum properties
		int skeletonsNum;
		std::vector<Skeleton_Props> props_sk;
		Fluid_Props props_oil;
		Fluid_Props props_gas;

		// Number of cells in radial direction
		int cellsNum_r;
		// Number of cells in vertical direction
		int cellsNum_z;

		// Gas content in oil
		Interpolate* Rs;
		Interpolate* Prs;

		// BHP will be converted to the depth
		double depth_point;
		// During the time flow rate decreases 'e' times in well test [sec] 
		double alpha;

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
			if(cells[cur].u_next.p < cells[beta].u_next.p)
				return 0.0;
			else
				return 1.0;
		};
		inline int getUpwindIdx(int cur, int beta)
		{
			if(cells[cur].u_next.p < cells[beta].u_next.p)
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
		inline double getPerm_r(const Cell& cell) const
		{
			const int idx = getSkeletonIdx(cell);
			return (cell.r > props_sk[idx].radius_eff ? props_sk[idx].perm_r : props_sk[idx].perm_eff);
		};
		inline double getKr_oil(double sat_oil) const
		{
			if(sat_oil > 1.0)
				return 1.0;
			else
				return props_oil.kr->Solve(sat_oil);
		};
		inline double getKr_oil_ds(double sat_oil) const
		{
			if(sat_oil > 1.0)
				return props_oil.kr->DSolve(1.0);
			else
				return props_oil.kr->DSolve(sat_oil);
		};
		inline double getKr_gas(double sat_oil) const
		{
			if(sat_oil > 1.0)
				return 0.0;
			else
				return props_gas.kr->Solve(sat_oil);
		};
		inline double getKr_gas_ds(double sat_oil) const
		{
			if(sat_oil > 1.0)
				return props_gas.kr->DSolve(1.0);
			else
				return props_gas.kr->DSolve(sat_oil);
		};
		inline double getB_oil(double p, double p_bub, bool SATUR) const
		{
			if(SATUR)
				return props_oil.b->Solve(p);
			else
				return props_oil.b->Solve(p_bub) * (1.0 + props_oil.beta * (p_bub - p));
		};
		inline double getBoreB_oil(double p, double p_bub, bool SATUR) const
		{
			//return props_oil.b_bore;
			return getB_oil(p, p_bub, SATUR);
		};
		inline double getB_oil_dp(double p, double p_bub, bool SATUR) const
		{
			if(SATUR)
				return props_oil.b->DSolve(p);
			else
				return -props_oil.b->Solve(p_bub) * props_oil.beta;
		};
		inline double getB_oil_dp_bub(double p_bub, bool SATUR) const
		{
			if (SATUR)
				return 0.0;
			else
				return props_oil.b->Solve(p_bub) * props_oil.beta;
		};
		inline double getB_gas(double p) const
		{
			return props_gas.b->Solve(p);
		};
		inline double getB_gas_dp(double p) const
		{
			return props_gas.b->DSolve(p);
		};
		inline double getRs(double p, double p_bub, bool SATUR) const
		{
			if(SATUR)
				return Rs->Solve(p);
			else
				return Rs->Solve(p_bub);
		};
		inline double getRs_dp(double p, double p_bub, bool SATUR) const
		{
			if(SATUR)
				return Rs->DSolve(p);
			else
				return 0.0;
		};
		inline double getRs_dp_bub() const
		{
				return 0.0;
		};
		inline double getPresFromRs(double rs) const
		{
			return Prs->Solve(rs);
		};
		inline void solveP_bub()
		{
			/*int idx;
			double factRs, dissGas;

			for(int i = 0; i < cellsNum_r+2; i++)
				for(int j = 0; j < cellsNum_z+2; j++)
				{
					idx = i * (cellsNum_z + 2) + j;

					Var2phase& next = cells[idx].u_next;
					Var2phase& prev = cells[idx].u_prev;

					if(next.s > 1.0)
						next.s = 1.0;

					dissGas = (1.0 - next.s) * getB_oil(next.p, next.p_bub, next.SATUR) / ( (1.0 - next.s) * getB_oil(next.p, next.p_bub, next.SATUR) + next.s * getB_gas(next.p));
					factRs = getRs(prev.p, prev.p_bub, next.SATUR) + dissGas;

					if(getRs(next.p, next.p, next.SATUR) > factRs)
					{
						next.p_bub = getPresFromRs(factRs);
						next.SATUR = false;
					} else {
						next.p_bub = next.p;
						next.SATUR = true;
					}
				}*/
			int idx;

			for (int i = 0; i < cellsNum_r + 2; i++)
				for (int j = 0; j < cellsNum_z + 2; j++)
				{
					idx = i * (cellsNum_z + 2) + j;

					Var2phase& next = cells[idx].u_next;
					if (next.s > 1.0 + EQUALITY_TOLERANCE)
						next.s = 1.0;

					if (next.p > next.p_bub + EQUALITY_TOLERANCE)
					{
						next.SATUR = false;
						next.s = 1.0;
					}
					else
						next.SATUR = true;
				}
		};

		// Thermal functions
		inline double getRho_oil(double p, double p_bub, bool SATUR) const
		{
			return (props_oil.dens_stc + getRs(p, p_bub, SATUR) * props_gas.dens_stc) / getB_oil(p, p_bub, SATUR);
		};
		inline double getRho_gas(double p) const
		{
			return props_gas.dens_stc / getB_gas(p);
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

			Var2phase* var;
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
				return -getPerm_r(cell) * getKr_oil(var->s) / props_oil.visc * getNablaP(cell, varNum, axis);
			case Z_AXIS:
				return -props_sk[idx].perm_z * getKr_oil(var->s) / props_oil.visc * getNablaP(cell, varNum, axis);
			}
		};
		inline double getGasVelocity(Cell& cell, int varNum, int axis)
		{
			const int idx = getSkeletonIdx( cell );

			Var2phase* var;
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
				return -getPerm_r(cell) * getKr_gas(var->s) / props_gas.visc * getNablaP(cell, varNum, axis);
			case Z_AXIS:
				return -props_sk[idx].perm_z * getKr_gas(var->s) / props_gas.visc * getNablaP(cell, varNum, axis);
			}
		};

		// First eqn
		double solve_eq1(int cur);
		double solve_eq1_dp(int cur);
		double solve_eq1_ds(int cur);
		double solve_eq1_dp_beta(int cur, int beta);
		double solve_eq1_ds_beta(int cur, int beta);

		// Second eqn
		double solve_eq2(int cur);
		double solve_eq2_dp(int cur);
		double solve_eq2_ds(int cur);
		double solve_eq2_dp_beta(int cur, int beta);
		double solve_eq2_ds_beta(int cur, int beta);

		// Left boundary condition
		double solve_eqLeft(int cur);
		double solve_eqLeft_dp(int cur);
		double solve_eqLeft_ds(int cur);
		double solve_eqLeft_dp_beta(int cur);
		double solve_eqLeft_ds_beta(int cur);

		// Right boundary condition
		double solve_eqRight(int cur);
		double solve_eqRight_dp(int cur);
		double solve_eqRight_ds(int cur);
		double solve_eqRight_dp_beta(int cur);
		double solve_eqRight_ds_beta(int cur);

		// Finds functional
		double solveH();

	public:
		GasOil_RZ();
		~GasOil_RZ();
	
		void setPeriod(int period);
		double getRate(int cur);
	};
};

#endif /* GASOIL_RZ_H_ */