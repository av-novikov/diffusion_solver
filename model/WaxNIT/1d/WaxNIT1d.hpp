#ifndef WAXNIT1D_HPP_
#define WAXNIT1D_HPP_

#include "model/WaxNIT/1d/Properties.hpp"
#include "model/cells/Variables.hpp"
#include "model/cells/LinearCell.hpp"
#include "model/Basic1d/Basic1d.hpp"
#include "util/Interpolate.h"

namespace wax_nit1d
{
	static const int stencil = basic1d::stencil;
	static const int Lstencil = basic1d::Lstencil;
	static const int Rstencil = basic1d::Rstencil;

	typedef VarWaxNIT1d Variable;
	typedef TapeVarWaxNIT1d TapeVariable;
	typedef LinearCell<Variable, Skeleton_Props> Cell;
	template <typename TVariable> using TCell = LinearCell<TVariable, Skeleton_Props>;
	typedef Cell::Type Type;

	class WaxNIT1d : public basic1d::Basic1d<Variable, Properties, Skeleton_Props, TCell, WaxNIT1d>
	{
		template<typename> friend class Snapshotter;
		template<typename> friend class VTKSnapshotter;
		template<typename> friend class basic1d::Basic1dSolver;
		friend class WaxNIT1dSolver;
	protected:
		Water_Props props_wat;
		Oil_Props props_oil;
		Gas_Props props_gas;
		Wax_Props props_wax;
		double L;

		void setProps(Properties& props);
		void makeDimLess();
		void setInitialState();

		// Service functions
		inline adouble getOilVelocity(const Cell& cell) const
		{
			const TapeVariable& next = var[0];
			return -cell.props->getPermCoseni(next.m) * props_oil.getKr(next.s_w, next.s_o, cell.props) / props_oil.getViscosity(next.p) * 
				(var[2].p - var[1].p) / (cells[cell.num + 1].x - cells[cell.num - 1].x);
		};
		inline adouble getWatVelocity(const Cell& cell) const
		{
			const TapeVariable& next = var[0];
			return -cell.props->getPermCoseni(next.m) * props_wat.getKr(next.s_w, next.s_o, cell.props) / props_wat.getViscosity(next.p) *
				(var[2].p - var[1].p) / (cells[cell.num + 1].x - cells[cell.num - 1].x);
		};
		// Just for snapshotter
		inline double getOilVel(const Cell& cell) const
		{
			const auto& next = cell.u_next;
			const Cell& cell1 = cells[cell.num - 1];
			const Cell& cell2 = cells[cell.num + 1];
			return -cell.props->getPermCoseni(next.m).value() * props_oil.getKr(next.s_w, next.s_o, cell.props).value() /
					props_oil.getViscosity(next.p).value() * (cell2.u_next.p - cell1.u_next.p) / (cell2.x - cell1.x);
		};
		inline double getWatVel(const Cell& cell) const
		{
			const auto& next = cell.u_next;
			const Cell& cell1 = cells[cell.num - 1];
			const Cell& cell2 = cells[cell.num + 1];
			return -cell.props->getPermCoseni(next.m).value() * props_wat.getKr(next.s_w, next.s_o, cell.props).value() /
					props_wat.getViscosity(next.p).value() * (cell2.u_next.p - cell1.u_next.p) / (cell2.x - cell1.x);
		};

		inline adouble getOilVelocityAbs(const Cell& cell) const 
		{
			return fabs(getOilVelocity(cell));
		};
		inline adouble getTrans(const Cell& cell, adouble m_cell, const Cell& beta, adouble m_beta) const
		{
			adouble k1, k2, S;
			k1 = cell.props->getPermCoseni(m_cell);
			k2 = beta.props->getPermCoseni(m_beta);
			if (k1 == 0.0 && k2 == 0.0)
				return 0.0;
			S = props_sk.height;
			return 2.0 * k1 * k2 * S / (k1 * beta.hx + k2 * cell.hx);
		};
		inline adouble getAverage(adouble p1, const Cell& cell1, adouble p2, const Cell& cell2) const
		{
			return (p1 * (adouble)cell2.hx + p2 * (adouble)cell1.hx) / (adouble)(cell1.hx + cell2.hx);
		};

		TapeVariable* var;
		adouble* h;
		inline double getDistance(const Cell& cell1, const Cell& cell2) const
		{
			return cell2.x - cell1.x;
		};

		void solve_eqMiddle(const Cell& cell);
		void solve_eqLeft(const Cell& cell);
		void solve_eqRight(const Cell& cell);
		void setVariables(const Cell& cell);

	public:
		WaxNIT1d();
		~WaxNIT1d();

		double getRate(int cur) const;

		double* x;
		double* y;
		double** jac;

		static const int var_size = Variable::size;
	};
};


#endif /* WAXNIT1D_HPP_ */