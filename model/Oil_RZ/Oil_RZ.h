#ifndef OIL_RZ_H_
#define OIL_RZ_H_

#include <vector>
#include <map>
#include <cstdlib>

#include "util/utils.h"
#include "model/cells/Variables.hpp"
#include "model/cells/CylindricalCell2D.h"
#include "model/AbstractModel.hpp"

namespace oil_rz
{
	typedef CylCell2D<Var1phase> Cell;

	struct Skeleton_Props
	{
		double m;
		double beta;
		double dens_stc;

		double perm_r;
		double perm_z;

		std::vector<double> perms_eff;
		std::vector<double> radiuses_eff;
		std::vector<double> skins;
		double perm_eff;
		double radius_eff;
		double skin;

		double h1, h2;
		double height;

		int cellsNum_z;

		double p_out;
		double p_init;
	};

	struct Oil_Props
	{
		double visc;
		double dens_stc;
		double beta;
		double b_bore;
	};

	struct Properties
	{
		std::vector<double> timePeriods;
		std::vector<double> rates;
		std::vector<double> pwf;

		bool leftBoundIsRate;
		bool rightBoundIsPres;

		std::vector<std::pair<int,int> > perfIntervals;

		double ht;
		double ht_min;
		double ht_max;

		double alpha;

		int cellsNum_r;
		int cellsNum_z;

		double r_w;
		double r_e;

		std::vector<Skeleton_Props> props_sk;
		double depth_point;
		
		Oil_Props props_oil;
	};
	class Oil_RZ : public AbstractModel<Var1phase, Properties, CylCell2D, Oil_RZ>
	{
		template<typename> friend class Snapshotter;
		template<typename> friend class GRDECLSnapshotter;
		template<typename> friend class VTKSnapshotter;
		template<typename> friend class AbstractMethod;
		friend class OilRZSolver;	
	protected:
		int skeletonsNum;
		std::vector<Skeleton_Props> props_sk;
		Oil_Props props_oil;

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
			double h;

			switch(axis)
			{
			case R_AXIS:
				nebr1 = &cells[cell.num - cellsNum_z - 2];
				nebr2 = &cells[cell.num + cellsNum_z + 2];
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

			Var1phase* var;
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

		double solve_eq(int cur);
		double solve_eq_dp(int cur);
		double solve_eq_dp_beta(int cur, int beta);
		
		double solve_eqLeft(int cur);
		double solve_eqLeft_dp(int cur);
		double solve_eqLeft_dp_beta(int cur);
		
		double solve_eqRight(int cur);
		double solve_eqRight_dp(int cur);
		double solve_eqRight_dp_beta(int cur);

		double solveH();

	public:
		Oil_RZ();
		~Oil_RZ();

		void setPeriod(int period);
	};

};

#endif /* OIL_RZ_H_ */