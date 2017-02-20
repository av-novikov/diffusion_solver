#ifndef VPP2D_HPP_
#define VPP2D_HPP_

#include <vector>

#include "model/cells/Variables.hpp"
#include "model/cells/CylCell2D.h"
#include "model/AbstractModel.hpp"
#include "util/utils.h"
#include "model/VPP2d/Properties.hpp"

namespace vpp2d
{
	static const int R_AXIS = 0;
	static const int Z_AXIS = 1;

	typedef VarSimpleVPP<double> Variable;
	typedef TapeVarSimpleVPP TapeVariable;
	typedef NewCylCell2D<Variable, Skeleton_Props> Cell;
	template <typename TVariable> using TCell = NewCylCell2D<TVariable, Skeleton_Props>;

	class VPP2d : public AbstractModel<Variable, Properties, TCell, VPP2d>
	{
		template<typename> friend class Snapshotter;
		template<typename> friend class GRDECLSnapshotter;
		template<typename> friend class VTKSnapshotter;
		template<typename> friend class AbstractMethod;
		friend class VPPSolver;

	protected:
		// Concentration
		std::vector<double> c;
		double Conc;

		// Continuum properties
		int skeletonsNum;
		std::vector<Skeleton_Props> props_sk;
		Water_Props props_wat;
		Oil_Props props_oil;

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
		inline adouble upwindIsCur(int cur, int beta)
		{
			if (cells[cur].u_next.p < cells[beta].u_next.p)
				return adouble(0.0);
			else
				return adouble(1.0);
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
		inline double getTrans(const Cell& cell, const Cell& beta) const
		{
			double k1, k2, Square;

			if (abs(cell.num - beta.num) == 1) {
				k1 = cell.props->perm_z;
				k2 = beta.props->perm_z;
				if (k1 == 0.0 && k2 == 0.0)
					return 0.0;
				Square = 2.0 * M_PI * cell.r * cell.hr;
				return 2.0 * k1 * k2 * Square / (k1 * beta.hz + k2 * cell.hz);
			}
			else {
				k1 = (cell.r > cell.props->radius_eff ? cell.props->perm_r : cell.props->perm_eff);
				k2 = (beta.r > beta.props->radius_eff ? beta.props->perm_r : beta.props->perm_eff);
				if (k1 == 0.0 && k2 == 0.0)
					return 0.0;
				Square = 2.0 * M_PI * cell.hz * (cell.r + sign(beta.num - cell.num) * cell.hr / 2.0);
				return 2.0 * k1 * k2 * Square / (k1 * beta.hr + k2 * cell.hr);
			}
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
				return -cell.props->getPerm_r(cell.r) * props_wat.getKr(var->s).value() / props_wat.getViscosity(var->p).value() * getNablaP(cell, varNum, axis);
			case Z_AXIS:
				return -cell.props->perm_z * props_wat.getKr(var->s).value() / props_wat.getViscosity(var->p).value() * getNablaP(cell, varNum, axis);
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
				return -cell.props->getPerm_r(cell.r) * props_oil.getKr(var->s).value() / props_oil.getViscosity(var->p).value() * getNablaP(cell, varNum, axis);
			case Z_AXIS:
				return -cell.props->perm_z * props_oil.getKr(var->s).value() / props_oil.getViscosity(var->p).value() * getNablaP(cell, varNum, axis);
			}
		};

		// Finds functional
		double solveH();

		void setVariables(int cur);
		void solve_eqMiddle(int cur);
		void solve_eqLeft(int cur);
		void solve_eqRight(int cur);
		void solve_eqVertical(int cur);
		void solveFixVar();

	public:
		VPP2d();
		~VPP2d();

		double* x;
		double* y;

		void setPeriod(int period);
		double getRate(int cur);
	};

};

#endif /* VPP2D_HPP_ */
