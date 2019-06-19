#ifndef ACID1D_HPP_
#define ACID1D_HPP_

#include <vector>

#include "model/Acid/1d/Properties.hpp"
#include "model/cells/AcidVariables.hpp"
#include "model/cells/LinearCell.hpp"
#include "model/Basic1d/Basic1d.hpp"
#include "util/Interpolate.h"

namespace acid1d
{
	static const int stencil = basic1d::stencil;
	static const int Lstencil = basic1d::Lstencil;
	static const int Rstencil = basic1d::Rstencil;

	typedef JustAcid Variable;
	typedef TapeJustAcid TapeVariable;
	typedef LinearCell<Variable, Skeleton_Props> Cell;
	template <typename TVariable> using TCell = LinearCell<TVariable, Skeleton_Props>;
	typedef CalciteReaction CurrentReaction;
	typedef CurrentReaction::REACTS REACTS;
	typedef Cell::Type Type;

	class Acid1d : public basic1d::Basic1d<Variable, Properties, Skeleton_Props, TCell, Acid1d>
	{
		template<typename> friend class Snapshotter;
		template<typename> friend class VTKSnapshotter;
		friend class Acid1dSolver;

	protected:
		CurrentReaction reac;
		Water_Props props_w;
		Oil_Props props_o;
		Gas_Props props_g;

		double xa;
		std::vector<double> xas;

		void setProps(Properties& props);
		void makeDimLess();
		void setInitialState();

		// Service functions
		inline adouble getAverage(adouble p1, const Cell& cell1, adouble p2, const Cell& cell2) const
		{
			return (p1 * (adouble)cell2.hx + p2 * (adouble)cell1.hx) / (adouble)(cell1.hx + cell2.hx);
		};
		inline adouble getReactionRate(TapeVariable& var) const
		{
			return var.sw * props_w.getDensity(var.p, var.xa, var.xw, var.xs) *
					(var.xa - props_sk.xa_eqbm) * 
					reac.getReactionRate(props_sk.m_init, var.m) / reac.comps[REACTS::ACID].mol_weight;
		};
		double getNablaP(Cell& cell, int varNum)
		{
			Cell* nebr1;
			Cell* nebr2;
			double h, r_eff;

			nebr1 = &cells[cell.num - 1];
			nebr2 = &cells[cell.num + 1];

			r_eff = props_sk.radius_eff;
			if ((nebr1->x < r_eff) && (nebr2->x > r_eff))
			{
				if (cell.x > r_eff)
					nebr1 = &cell;
				else
					nebr2 = &cell;
			}
			h = nebr2->x - nebr1->x;
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

		inline adouble getTrans(const Cell& cell, adouble m_cell, const Cell& beta, adouble m_beta) const
		{
			adouble k1, k2, S;
			k1 = props_sk.getPermCoseni(m_cell);
			k2 = props_sk.getPermCoseni(m_beta);
			if (k1 == 0.0 && k2 == 0.0)
				return 0.0;
			S = props_sk.height;
			return 2.0 * k1 * k2 * S / (k1 * beta.hx + k2 * cell.hx);
		};
		inline double getWaterVelocity(Cell& cell)
		{
			const Variable next = cell.u_next;
			return -props_sk.getPermCoseni(next.m).value() * props_w.getKr(next.sw, &props_sk).value() / props_w.getViscosity(next.p, next.xa, next.xw, next.xs).value() * getNablaP(cell, NEXT);
		};
		inline double getOilVelocity(Cell& cell)
		{
			const Variable next = cell.u_next;
			return -props_sk.getPermCoseni(next.m).value() * props_o.getKr(next.sw, &props_sk).value() / props_o.getViscosity(next.p).value() * getNablaP(cell, NEXT);
		};

		void solve_eqMiddle(const Cell& cell);
		void solve_eqLeft(const Cell& cell);
		void solve_eqRight(const Cell& cell);
		void setVariables(const Cell& cell);

	public:
		Acid1d();
		~Acid1d();

		double* x;
		double* y;
		double** jac;

		void setPeriod(int period);
		double getRate(int cur);
		static const int var_size = Variable::size;
	};
};

#endif /* ACID1D_HPP_ */