#ifndef GASOIL_RZ_H_
#define GASOIL_RZ_H_

#define ADOLC_ADVANCED_BRANCHING

#include <vector>
#include <map>
#include <string>

#include "model/cells/Variables.hpp"
#include "model/GasOil_RZ/Properties.hpp"
#include "model/AbstractModel.hpp"

namespace gasOil_rz
{
	static const int stencil = 5;
	static const int Lstencil = 3;
	static const int Rstencil = 2;
	static const int Vstencil = 2;

	typedef Var2phase Variable;
	typedef TapeVarGasOil TapeVariable;
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

		// Continuum properties
		int skeletonsNum;
		std::vector<Skeleton_Props> props_sk;
		Oil_Props props_oil;
		Gas_Props props_gas;

		// Number of cells in radial direction
		int cellsNum_r;
		// Number of cells in vertical direction
		int cellsNum_z;

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
		inline void getNeighborIdx(int cur, int* const neighbor)
		{
			neighbor[0] = cur - cellsNum_z - 2; 
			neighbor[1] = cur + cellsNum_z + 2;
			neighbor[2] = cur - 1;
			neighbor[3] = cur + 1;
		};
		inline int isNotSatur(int cur) const
		{
			return !cells[cur].u_next.SATUR;
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
		inline double getTrans(const Cell& cell, const Cell& beta) const 
		{
			double k1, k2, S;

			if( abs(cell.num - beta.num) == 1) {
				k1 = cell.props->perm_z;
				k2 = beta.props->perm_z;
				if(k1 == 0.0 && k2 == 0.0)
					return 0.0;
				S = 2.0 * M_PI * cell.r * cell.hr;
				return 2.0 * k1 * k2 * S / (k1 * beta.hz + k2 * cell.hz);
			} else {
				k1 = (cell.r > cell.props->radius_eff ? cell.props->perm_r : cell.props->perm_eff);
				k2 = (beta.r > beta.props->radius_eff ? beta.props->perm_r : beta.props->perm_eff);
				if(k1 == 0.0 && k2 == 0.0)
					return 0.0;
				S = 2.0 * M_PI * cell.hz * (cell.r + sign(beta.num - cell.num) * cell.hr / 2.0);
				return 2.0 * k1 * k2 * S / (k1 * beta.hr + k2 * cell.hr);
			}
		};

		inline void solveP_bub()
		{
			int idx;

			for (int i = 0; i < cellsNum_r + 2; i++)
				for (int j = 0; j < cellsNum_z + 2; j++)
				{
					idx = i * (cellsNum_z + 2) + j;

					Skeleton_Props* props = cells[idx].props;
					Variable& next = cells[idx].u_next;

					if (next.SATUR)
					{
						if ((next.s > 1.0 + EQUALITY_TOLERANCE) || (next.p > props_oil.p_sat))
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
							next.s -= 0.01;
//							next.p_bub = next.p;
						}
						else
							next.s = 1.0;
					}
				}
		};
		// Thermal functions
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
				return -cell.props->getPerm_r(cell.r) * props_oil.getKr(var->s).value() / props_oil.getViscosity(var->p).value() * getNablaP(cell, varNum, axis);
			case Z_AXIS:
				return -cell.props->perm_z * props_oil.getKr(var->s).value() / props_oil.getViscosity(var->p).value() * getNablaP(cell, varNum, axis);
			}
		};
		inline double getGasVelocity(Cell& cell, int varNum, int axis)
		{
			const int idx = getSkeletonIdx( cell );

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
				return -cell.props->getPerm_r(cell.r) * props_gas.getKr(var->s).value() / props_gas.getViscosity(var->p).value() * getNablaP(cell, varNum, axis);
			case Z_AXIS:
				return -props_sk[idx].perm_z * props_gas.getKr(var->s).value() / props_gas.getViscosity(var->p).value() * getNablaP(cell, varNum, axis);
			}
		};

		void setVariables(int cur);
		void solve_eqMiddle(int cur);
		void solve_eqLeft(int cur);
		void solve_eqRight(int cur);
		void solve_eqVertical(int cur);

		// Finds functional
		double solveH();

	public:
		GasOil_RZ();
		~GasOil_RZ();
	
		double* x;
		double* y;

		void setPeriod(int period);
		double getRate(int cur);
	};
};

#endif /* GASOIL_RZ_H_ */