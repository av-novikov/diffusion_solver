#ifndef ACID2D_HPP_
#define ACID2D_HPP_

#include <vector>

#include "model/Acid/2d/Properties.hpp"
#include "model/cells/AcidVariables.hpp"
#include "model/cells/CylCell2D.h"
#include "model/Basic2d/Basic2d.hpp"
#include "util/Interpolate.h"

namespace acid2d
{
	static const int R_AXIS = 0;
	static const int Z_AXIS = 1;

	static const int stencil = basic2d::stencil;
	static const int Lstencil = basic2d::Lstencil;
	static const int Rstencil = basic2d::Rstencil;
	static const int Vstencil = basic2d::Vstencil;

	typedef FirstAcid Variable;
	typedef TapeFirstAcid TapeVariable;
	typedef NewCylCell2D<Variable, Skeleton_Props> Cell;
	template <typename TVariable> using TCell = NewCylCell2D<TVariable, Skeleton_Props>;
	typedef CalciteReaction CurrentReaction;
	typedef CurrentReaction::REACTS REACTS;

	class Acid2d : public basic2d::Basic2d<Variable, Properties, Skeleton_Props, TCell, Acid2d>
	{
		template<typename> friend class Snapshotter;
		template<typename> friend class GRDECLSnapshotter;
		template<typename> friend class VTKSnapshotter;
		friend class Acid2dSolver;

	protected:
		CurrentReaction reac;
		std::vector<Skeleton_Props> props_sk;
		Water_Props props_w;
		Oil_Props props_o;
		Gas_Props props_g;

		double xa;
		std::vector<double> xas;

		void setProps(Properties& props);
		void makeDimLess();

		// Service functions
		inline adouble getAvarage(adouble p1, const Cell& cell1, adouble p2, const Cell& cell2) const
		{
			double r1, r2;
			if (abs(cell1.num - cell2.num) == 1) 
			{
				r1 = cell1.hz;		r2 = cell2.hz;
			}
			else
			{
				r1 = cell1.hr;		r2 = cell2.hr;
			}

			return (p1 * (adouble)r2 + p2 * (adouble)r1) / (adouble)(r1 + r2);
		};
		inline adouble getReactionRate(TapeVariable& var, const Skeleton_Props& props) const
		{
			return var.sw * props_w.getMolarDensity(var.p, var.xa, var.xw) * 
					pow((var.xa - (adouble)props.xa_eqbm), (adouble)reac.alpha) * 
					reac.getReactionRate(props.m_init, var.m);
		};
		inline adouble getTrans(const Cell& cell, adouble m_cell, const Cell& beta, adouble m_beta) const
		{
			adouble k1, k2, S;

			if (abs(cell.num - beta.num) == 1) {
				k1 = cell.props->getPermCoseni_z(m_cell);
				k2 = beta.props->getPermCoseni_z(m_beta);
				if (k1 == 0.0 && k2 == 0.0)
					return 0.0;
				S = 2.0 * M_PI * cell.r * cell.hr;
				return 2.0 * k1 * k2 * S / (k1 * beta.hz + k2 * cell.hz);
			}
			else {
				k1 = cell.props->getPermCoseni_r(m_cell);
				k2 = beta.props->getPermCoseni_r(m_beta);
				if (k1 == 0.0 && k2 == 0.0)
					return 0.0;
				S = 2.0 * M_PI * cell.hz * (cell.r + sign(beta.num - cell.num) * cell.hr / 2.0);
				return 2.0 * k1 * k2 * S / (k1 * beta.hr + k2 * cell.hr);
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
				return -cell.props->getPermCoseni_r(cell.u_next.m).value() * props_w.getKr(var->sw, var->so, cell.props).value() / props_w.getViscosity(var->p, var->xa, var->xw).value() * getNablaP(cell, varNum, axis);
			case Z_AXIS:
				return -cell.props->getPermCoseni_z(cell.u_next.m).value() * props_w.getKr(var->sw, var->so, cell.props).value() / props_w.getViscosity(var->p, var->xa, var->xw).value() * getNablaP(cell, varNum, axis);
			}
		};
		inline double getGasVelocity(Cell& cell, int varNum, int axis)
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
				return -cell.props->getPermCoseni_r(cell.u_next.m).value() * props_g.getKr(var->sw, var->so, cell.props).value() / props_g.getViscosity(var->p).value() * getNablaP(cell, varNum, axis);
			case Z_AXIS:
				return -cell.props->getPermCoseni_z(cell.u_next.m).value() * props_g.getKr(var->sw, var->so, cell.props).value() / props_g.getViscosity(var->p).value() * getNablaP(cell, varNum, axis);
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
				return -cell.props->getPermCoseni_r(cell.u_next.m).value() * props_o.getKr(var->sw, var->so, cell.props).value() / props_o.getViscosity(var->p).value() * getNablaP(cell, varNum, axis);
			case Z_AXIS:
				return -cell.props->getPermCoseni_z(cell.u_next.m).value() * props_o.getKr(var->sw, var->so, cell.props).value() / props_o.getViscosity(var->p).value() * getNablaP(cell, varNum, axis);
			}
		};

		void solve_eqMiddle(const Cell& cell);
		void solve_eqLeft(const Cell& cell);
		void solve_eqRight(const Cell& cell);
		void solve_eqVertical(const Cell& cell);
		void setVariables(const Cell& cell);

	public:
		Acid2d();
		~Acid2d();

		double* x;
		double* y;
		double** jac;

		//double getRate(int cur);
	};
};

#endif /* ACID2D_HPP_ */