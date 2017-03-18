#ifndef BLACKOILRZ_HPP_
#define BLACKOILRZ_HPP_

#include "model/BlackOil_RZ/Properties.hpp"
#include "model/cells/Variables.hpp"
#include "model/cells/CylCell2D.h"
#include "model/Basic2d/Basic2d.hpp"
#include "model/Basic2d/Basic2dSolver.hpp"
#include "util/Interpolate.h"

namespace blackoil_rz
{
	static const int R_AXIS = 0;
	static const int Z_AXIS = 1;

	static const int stencil = basic2d::stencil;
	static const int Lstencil = basic2d::Lstencil;
	static const int Rstencil = basic2d::Rstencil;
	static const int Vstencil = basic2d::Vstencil;

	typedef VarBlackOil Variable;
	typedef TapeVarBlackOil TapeVariable;
	typedef NewCylCell2D<Variable, Skeleton_Props> Cell;
	template <typename TVariable> using TCell = NewCylCell2D<TVariable, Skeleton_Props>;

	class BlackOil_RZ : public basic2d::Basic2d<Variable, Properties, Skeleton_Props, TCell, BlackOil_RZ>
	{
		template<typename> friend class Snapshotter;
		template<typename> friend class VTKSnapshotter;
		template<typename> friend class basic2d::Basic2dSolver;
		friend class BlackOil2dSolver;
	protected:
		std::vector<Skeleton_Props> props_sk;
		Liquid_Props props_wat, props_oil;
		Gas_Props props_gas;

		void setProps(Properties& props);
		void makeDimLess();

		// Service functions
		inline double getTrans(const Cell& cell, const Cell& beta) const
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
				if (k1 == 0.0 && k2 == 0.0)
					return 0.0;
				S = 2.0 * M_PI * cell.hz * (cell.r + sign(beta.num - cell.num) * cell.hr / 2.0);
				return 2.0 * k1 * k2 * S / (k1 * beta.hr + k2 * cell.hr);
			}
		};
		inline void solveP_bub()
		{

		};

		void solve_eqMiddle(const Cell& cell);
		void solve_eqLeft(const Cell& cell);
		void solve_eqRight(const Cell& cell);
		void solve_eqVertical(const Cell& cell);
		void setVariables(const Cell& cell);

	public:
		BlackOil_RZ();
		~BlackOil_RZ();

		double getRate(int cur) const;

		double* x;
		double* y;
		double** jac;
	};
};

#endif /* BLACKOILRZ_HPP_ */
