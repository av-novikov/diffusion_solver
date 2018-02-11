#ifndef ACID2DNIT_HPP_
#define ACID2DNIT_HPP_

#include <vector>

#include "model/Acid/2dnit/Properties.hpp"
#include "model/cells/AcidVariables.hpp"
#include "model/cells/CylCell2D.h"
#include "model/Basic2d/Basic2d.hpp"
#include "util/Interpolate.h"

namespace acid2dnit
{
	static const int R_AXIS = 0;
	static const int Z_AXIS = 1;

	static const int stencil = basic2d::stencil;
	static const int Lstencil = basic2d::Lstencil;
	static const int Rstencil = basic2d::Rstencil;
	static const int Vstencil = basic2d::Vstencil;

	typedef JustAcidNIT Variable;
	typedef TapeJustAcidNIT TapeVariable;
	typedef NewCylCell2D<Variable, Skeleton_Props> Cell;
	template <typename TVariable> using TCell = NewCylCell2D<TVariable, Skeleton_Props>;
	typedef CalciteReaction CurrentReaction;
	typedef CurrentReaction::REACTS REACTS;
	typedef Cell::Type Type;

	class Acid2dNIT : public basic2d::Basic2d<Variable, Properties, Skeleton_Props, TCell, Acid2dNIT>
	{
		template<typename> friend class Snapshotter;
		template<typename> friend class GRDECLSnapshotter;
		template<typename> friend class VTKSnapshotter;
		friend class Acid2dNITSolver;

	protected:
		CurrentReaction reac;
		Water_Props props_w;
		Oil_Props props_o;
		Gas_Props props_g;

		double xa;
		std::vector<double> xas;
		double temp;
		std::vector<double> temps;

		void setProps(Properties& props);
		void makeDimLess();
		void setInitialState();

		// Service functions
		inline adouble getAverage(adouble p1, const Cell& cell1, adouble p2, const Cell& cell2) const
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
			return var.sw * props_w.getDensity(var.p, var.xa, var.xw, var.xs) *
					(var.xa - props.xa_eqbm) * 
					reac.getReactionRate(props.m_init, var.m, var.t) / reac.comps[REACTS::ACID].mol_weight;
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
		inline adouble getOilVelocity(const Cell& cell, const int axis) const
		{
			const TapeVariable& next = var[0];
			adouble tmp;
			adouble is_R_axis = (axis == R_AXIS) ? true : false;
			adouble is_Z_axis = (axis == Z_AXIS) ? true : false;
			condassign(tmp, is_R_axis,
				-cell.props->getPermCoseni_r(next.m) * props_o.getKr(next.sw, cell.props) / props_o.getViscosity(next.p) *
				(var[2].p - var[1].p) / (cells[cell.num + cellsNum_z + 2].r - cells[cell.num - cellsNum_z - 2].r));
			condassign(tmp, is_Z_axis,
				-cell.props->getPermCoseni_z(next.m) * props_o.getKr(next.sw, cell.props) / props_o.getViscosity(next.p) *
				((var[4].p - var[3].p) / (cells[cell.num + 1].z - cells[cell.num - 1].z) -
				grav * props_o.getDensity(next.p)));
			return tmp;
		};
		inline adouble getWatVelocity(const Cell& cell, const int axis) const
		{
			const TapeVariable& next = var[0];
			adouble tmp;
			adouble is_R_axis = (axis == R_AXIS) ? true : false;
			adouble is_Z_axis = (axis == Z_AXIS) ? true : false;
			condassign(tmp, is_R_axis,
				-cell.props->getPermCoseni_r(next.m) * props_w.getKr(next.sw, cell.props) / 
				props_w.getViscosity(next.p, next.xa, next.xw, next.xs) *
				(var[2].p - var[1].p) / (cells[cell.num + cellsNum_z + 2].r - cells[cell.num - cellsNum_z - 2].r));
			condassign(tmp, is_Z_axis,
				-cell.props->getPermCoseni_z(next.m) * props_w.getKr(next.sw, cell.props) / 
				props_w.getViscosity(next.p, next.xa, next.xw, next.xs) *
				((var[4].p - var[3].p) / (cells[cell.num + 1].z - cells[cell.num - 1].z) - 
				grav * props_w.getDensity(next.p, next.xa, next.xw, next.xs)));
			return tmp;
		};
		inline double getOilVel(const Cell& cell, const int axis) const
		{
			double tmp;
			const auto& next = cell.u_next;
			if (axis == R_AXIS)
			{
				const Cell& cell1 = cells[cell.num - cellsNum_z - 2];
				const Cell& cell2 = cells[cell.num + cellsNum_z + 2];
				tmp = -cell.props->getPermCoseni_r(next.m).value() * props_o.getKr(next.sw, cell.props).value() / 
					props_o.getViscosity(next.p).value() *
					(cell2.u_next.p - cell1.u_next.p) / (cell2.r - cell1.r);
			}
			else if (axis == Z_AXIS)
			{
				const Cell& cell1 = cells[cell.num - 1];
				const Cell& cell2 = cells[cell.num + 1];
				tmp = -cell.props->getPermCoseni_z(next.m).value() * props_o.getKr(next.sw, cell.props).value() / 
					props_o.getViscosity(next.p).value() *
					((cell2.u_next.p - cell1.u_next.p) / (cell2.z - cell1.z) - grav * props_o.getDensity(next.p).value());
			}
			return tmp;
		};
		inline double getWatVel(const Cell& cell, const int axis) const
		{
			double tmp;
			const auto& next = cell.u_next;
			if (axis == R_AXIS)
			{
				const Cell& cell1 = cells[cell.num - cellsNum_z - 2];
				const Cell& cell2 = cells[cell.num + cellsNum_z + 2];
				tmp = -cell.props->getPermCoseni_r(next.m).value() * props_w.getKr(next.sw, cell.props).value() /
					props_w.getViscosity(next.p, next.xa, next.xw, next.xs).value() *
					((cell2.u_next.p - cell1.u_next.p) / (cell2.r - cell1.r) - grav * props_w.getDensity(next.p, next.xa, next.xw, next.xs).value());
			}
			else if (axis == Z_AXIS)
			{
				const Cell& cell1 = cells[cell.num - 1];
				const Cell& cell2 = cells[cell.num + 1];
				tmp = -cell.props->getPermCoseni_z(next.m).value() * props_w.getKr(next.sw, cell.props).value() /
					props_w.getViscosity(next.p, next.xa, next.xw, next.xs).value() *
					((cell2.u_next.p - cell1.u_next.p) / (cell2.z - cell1.z) - grav * props_w.getDensity(next.p, next.xa, next.xw, next.xs).value());
			}
			return tmp;
		};

		TapeVariable* var;
		adouble* h;
		inline adouble getCn(const Cell& cell) const
		{
			const TapeVariable& next = var[0];
			return next.m *
				((1.0 - next.sw) * props_o.getDensity(next.p) * props_o.c +
					next.sw * props_w.getDensity(next.p, next.xa, next.xw, next.xs) * props_w.c) +
					(1.0 - next.m) * cell.props->dens_stc * cell.props->c;
		};
		inline adouble getAd(const Cell& cell) const
		{
			const TapeVariable& next = var[0];
			return next.m * ((1.0 - next.sw) * props_o.getDensity(next.p) * props_o.ad * props_o.c +
					next.sw * props_w.getDensity(next.p, next.xa, next.xw, next.xs) * props_w.ad * props_w.c);
		};
		inline adouble getLambda(const Cell& cell, const int axis) const
		{
			const TapeVariable& next = var[0];
			adouble tmp;
			adouble is_R_axis = (axis == R_AXIS) ? true : false;
			adouble is_Z_axis = (axis == Z_AXIS) ? true : false;
			condassign(tmp, is_R_axis, next.m *	((1.0 - next.sw) * props_o.lambda + next.sw * props_w.lambda) +
				(1.0 - next.m) * cell.props->lambda_r);
			condassign(tmp, is_Z_axis, next.m * ((1.0 - next.sw) * props_o.lambda + next.sw * props_w.lambda) +
				(1.0 - next.m) * cell.props->lambda_z);
			return tmp;
		};
		inline adouble getJT(const Cell& cell, const int axis) const
		{
			const TapeVariable& next = var[0];
			return	props_o.getDensity(next.p) * props_o.c * props_o.jt * getOilVelocity(cell, axis) +
				props_w.getDensity(next.p, next.xa, next.xw, next.xs) * props_w.c * props_w.jt * getWatVelocity(cell, axis);
		};
		inline adouble getA(const Cell& cell, const int axis) const
		{
			const TapeVariable& next = var[0];
			return	props_o.getDensity(next.p) * props_o.c * getOilVelocity(cell, axis) +
				props_w.getDensity(next.p, next.xa, next.xw, next.xs) * props_w.c * getWatVelocity(cell, axis);
		};
		struct DivIndices
		{
			adouble ther;
			adouble pres;
			DivIndices(adouble _ther, adouble _pres) : ther(_ther), pres(_pres) {};
		};
		inline DivIndices getDivCoeff(const Cell& cell, const Cell& beta) const
		{
			adouble r1, r2, lambda, A, a, jt, JT;
			if (abs(cell.num - beta.num) == 1)
			{
				A = getA(cell, Z_AXIS);
				a = fmax((adouble)0.0, signA(cell.z - beta.z) * A);
				JT = getJT(cell, Z_AXIS);
				jt = fmax(0.0, sign(cell.z - beta.z) * JT);
				r1 = cell.hz;		r2 = beta.hz;
				lambda = (r1 * getLambda(beta, Z_AXIS) + r2 * getLambda(cell, Z_AXIS)) / (r1 + r2) / r1;
			}
			else
			{
				A = getA(cell, R_AXIS);
				a = fmax((adouble)0.0, signA(cell.r - beta.r) * A);
				JT = getJT(cell, R_AXIS);
				jt = fmax(0.0, sign(cell.r - beta.r) * JT);
				r1 = cell.hr;		r2 = beta.hr;
				lambda = (r1 * getLambda(beta, R_AXIS) + r2 * getLambda(cell, R_AXIS)) / (r1 + r2) / r1;
			}
			DivIndices coeff(a + lambda, jt);
			return coeff;
		};
		inline double getDistance(const Cell& cell1, const Cell& cell2) const
		{
			double dist;
			if (fabs(cell2.z - cell1.z) > EQUALITY_TOLERANCE)
				dist = cell2.z - cell1.z;
			else
				dist = cell2.r - cell1.r;

			return dist;
		};

		void solve_eqMiddle(const Cell& cell);
		void solve_eqLeft(const Cell& cell);
		void solve_eqRight(const Cell& cell);
		void solve_eqVertical(const Cell& cell);
		void setVariables(const Cell& cell);

	public:
		Acid2dNIT();
		~Acid2dNIT();

		double* x;
		double* y;
		double** jac;

		void setPeriod(int period);
		double getRate(int cur);
		static const int var_size = Variable::size;
	};
};

#endif /* ACID2DNIT_HPP_ */