#ifndef GASOIL_RZ_H_
#define GASOIL_RZ_H_

#include <vector>
#include <map>
#include <string>

#include "model/GasOil_RZ/Properties.hpp"
#include "model/AbstractModel.hpp"

namespace gasOil_rz
{
	typedef Var2phase Variable;
	typedef NewCylCell2D<Variable, Skeleton_Props> Cell;
	template <typename TVariable> using TCell = NewCylCell2D<TVariable, Skeleton_Props>;

	class GasOil_RZ : public AbstractModel<Variable, Properties, TCell, GasOil_RZ>
	{
		template<typename> friend class Snapshotter;
		template<typename> friend class GRDECLSnapshotter;
		template<typename> friend class VTKSnapshotter;
		template<typename> friend class AbstractMethod;
		friend class GasOil2DSolver;	

	protected:
		static const int size;
		static const int schemeVarNum;
		static const int boundVarNum;

		// Continuum properties
		int skeletonsNum;
		std::vector<Skeleton_Props> props_sk;
		Oil_Props props_oil;
		Gas_Props props_gas;

		// Number of cells in radial direction
		int cellsNum_r;
		// Number of cells in vertical direction
		int cellsNum_z;

		// Gas content in oil
		//Interpolate* Rs;
		//Interpolate* Prs;

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
		inline int getUpwindIdx(int cur, int beta) const
		{
			if(cells[cur].u_next.p < cells[beta].u_next.p)
				return beta;
			else
				return cur;
		};
		inline int getUpwindIdxTape(int cur, int beta, int tapeIdx) const
		{
			if (cells[cur].u_next.p < cells[beta].u_next.p)
				return tapeIdx;
			else
				return 0;
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
		inline double getB_oil_dp_bub(double p, double p_bub, bool SATUR) const
		{
			if (SATUR)
				return 0.0;
			else
				return props_oil.b->Solve(p_bub) * props_oil.beta + 
						(1 - props_oil.beta * (p - p_bub)) * props_oil.b->DSolve(p_bub);
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
				return props_oil.Rs->Solve(p);
			else
				return props_oil.Rs->Solve(p_bub);
		};
		inline double getRs_dp(double p, double p_bub, bool SATUR) const
		{
			if(SATUR)
				return props_oil.Rs->DSolve(p);
			else
				return 0.0;
		};
		inline double getRs_dp_bub(double p_bub, bool SATUR) const
		{
			if (SATUR)
				return 0.0;
			else
				return props_oil.Rs->DSolve(p_bub);
		};
		inline void solveP_bub()
		{
			int idx;

			for (int i = 0; i < cellsNum_r + 2; i++)
				for (int j = 0; j < cellsNum_z + 2; j++)
				{
					idx = i * (cellsNum_z + 2) + j;

					Var2phase& next = cells[idx].u_next;

					if (next.SATUR)
					{
						if (next.s > 1.0 + EQUALITY_TOLERANCE)
						{
							next.SATUR = false;
//							next.s = 1.0;
							next.p_bub -= 0.01 * next.p_bub;
						}
						else
							next.p_bub = next.p;
					}
					else 
					{
						if (next.p_bub > next.p + EQUALITY_TOLERANCE)
						{
							next.SATUR = true;
//							next.p_bub = next.p;
							next.s -= 0.01;
						}
						else
							next.s = 1.0;
					}
				}
		};

		// Thermal functions
		inline double getRho_oil(double p, double p_bub, bool SATUR) const
		{
			return (props_oil.dens_stc + props_oil.getRs(p, p_bub, SATUR) * props_gas.dens_stc) / props_oil.getB(p, p_bub, SATUR);
		};
		inline double getRho_gas(double p) const
		{
			return props_gas.dens_stc / props_gas.getB(p);
		};
		inline double getNablaP(Cell& cell, int varNum, int axis)
		{
			Cell* nebr1;
			Cell* nebr2;
			double h, r_eff;

			switch(axis)
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
			Variable* var;
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
				return -cell.props->getPerm_r(cell.r) * props_oil.getKr(var->s) / props_oil.visc * getNablaP(cell, varNum, axis);
			case Z_AXIS:
				return -cell.props->perm_z * props_oil.getKr(var->s) / props_oil.visc * getNablaP(cell, varNum, axis);
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
				return -cell.props->getPerm_r(cell.r) * props_gas.getKr(var->s) / props_gas.visc * getNablaP(cell, varNum, axis);
			case Z_AXIS:
				return -props_sk[idx].perm_z * props_gas.getKr(var->s) / props_gas.visc * getNablaP(cell, varNum, axis);
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
		void setMiddleIndependent(double* x, int cur);

		// Left boundary condition
		double solve_eqLeft(int cur);
		double solve_eqLeft_dp(int cur);
		double solve_eqLeft_ds(int cur);
		double solve_eqLeft_dp_beta(int cur);
		double solve_eqLeft_ds_beta(int cur);
		void setLeftIndependent(double* x, int cur);

		// Right boundary condition
		double solve_eqRight(int cur);
		double solve_eqRight_dp(int cur);
		double solve_eqRight_ds(int cur);
		double solve_eqRight_dp_beta(int cur);
		double solve_eqRight_ds_beta(int cur);
		void setRightIndependent(double* x, int cur);

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